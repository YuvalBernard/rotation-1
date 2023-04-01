data {
int N; // num. of measurements
  real<lower=0> R1a;
  real<lower=0> R2a;
  real<lower=0> w1;
  real dw;
  real<lower=0> tp;
  vector[N] xZ; // offsets measured in spectrum.
  vector[N] Z; // Z-spectrum values.
}
transformed data {
  real sigma = sd(Z);
  vector[N] theta = atan(w1./xZ);
  vector[N] R_eff = R1a * cos(theta).^2 + R2a * sin(theta).^2;
}
parameters {
  real<lower=0> R2b;
  real<lower=0> k;
  real<lower=0,upper=1> f;
}
transformed parameters {
  real gamma_squared = w1^2 * (k + R2b)/k + (k + R2b)^2;
  vector[N] R_ex_max = f*k*sin(theta).^2 .* ((xZ - dw).^2 + R2b*(w1^2 + xZ.^2)/k + R2b*(k + R2b))/gamma_squared; 
  vector[N] R_ex = R_ex_max * gamma_squared ./ (gamma_squared + (xZ - dw).^2);
  vector[N] R1p = R_ex + R_eff;
  vector[N] Z_ss = cos(theta).^2 * R1a ./ R1p;
  vector[N] Z_hat = (cos(theta).^2 - Z_ss) .* exp(-R1p * tp) + Z_ss;
}
model {
  R2b ~ normal(27000,2000);
  f ~ normal(0.02,0.01);
  k ~ normal(285,10);
  Z ~ normal(Z_hat,sigma);
}
generated quantities {
  real Z_pred[N];
  Z_pred = normal_rng(Z_hat, sigma);
}
