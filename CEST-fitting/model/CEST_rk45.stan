functions {
  vector cest(real t,vector M,real xZ,real R1a,real R2a,real R1b,real R2b,real k,real f,real w1,real dw) {
    real ka = k*f;     // pool a exchange rate        
    vector[6] dMdt;
    dMdt[1] = -(R2a + ka)*M[1] + -xZ*2*pi()*M[2] + k*M[4];
    dMdt[2] = xZ*2*pi()*M[1] + -(R2a + ka)*M[2] + w1*2*pi()*M[3] + k*M[5];
    dMdt[3] = -w1*2*pi()*M[2] + -(R1a + ka)*M[3] + k*M[6] + R1a;
    dMdt[4] = ka*M[1] + -(R2b + k)*M[4] + (dw - xZ)*2*pi()*M[5];
    dMdt[5] = ka*M[2] + -(dw - xZ)*2*pi()*M[4] + -(R2b + k)*M[5] + w1*2*pi()*M[6];
    dMdt[6] = k*M[3] + -w1*2*pi()*M[5] + -(R1b + k)*M[6] + f*R1b;
    return dMdt;
  }
}
data {
  int N; // num. of measurements
  real<lower=0> R1a;
  real<lower=0> R2a;
  real<lower=0> w1;
  real dw;
  array[1] real<lower=0> tp;
  real xZ[N]; // offsets measured in spectrum.
  real Z[N]; // Z-spectrum values.
}
parameters {
  real<lower=0> R1b;
  real<lower=0> R2b;
  real<lower=0> k;
  real<lower=0,upper=1> f;
  real<lower=0> sigma;
}
transformed parameters {
  vector[6] M0 = to_vector([0,0,1,0,0,f]);
  vector[N] Z_hat;
  array[1] vector[6] magnetization;
  for (i in 1:N) {
    magnetization = ode_rk45(cest,M0,0,tp,xZ[i],R1a,R2a,R1b,R2b,k,f,w1,dw);
    Z_hat[i] = magnetization[1,3];
  }
}
model {
  R1b ~ normal(0,50);
  R2b ~ normal(27000,2000);
  f ~ normal(0.02,0.01);
  k ~ normal(285,10);
  sigma ~ cauchy(0,0.01);
  Z ~ normal(Z_hat,sigma);
}
generated quantities {
  real Z_pred[N];
  Z_pred = normal_rng(Z_hat, sigma);
}