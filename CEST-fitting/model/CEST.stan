// Fit a 2 pool Z-spectrum

// Given parameters are R1a, R2a, dw, w1.
// Parameters to fit are R1b, R2b, f, k.

data {
  int N; // num. of measurements
  real<lower=0> R1a;
  real<lower=0> R2a;
  real<lower=0> w1;
  real dw;
  real<lower=0> tp;
  real xZ[N]; // offsets measured in spectrum.
  real Z[N]; // Z-spectrum values.
}

transformed data{
  real<lower=0> sigma = sd(Z); // Standard deviation of Z.
}

parameters {
  real<lower=0> R1b;
  real<lower=0> R2b;
  real<lower=0> k;
  real<lower=0,upper=1> f;
}

transformed parameters {
  matrix[6,6] A = rep_matrix(0,6,6);
  vector[6] b = to_vector([0,0,R1a,0,0,f*R1b]);
  vector[6] M; // Auxiliary vector
  real Z_hat[N];
  
  A[1,1] = -(R2a + f*k);
  A[1,4] = k;
  A[2,2] = -(R2a + f*k);
  A[2,3] = w1*2*pi();
  A[2,5] = k;
  A[3,2] = -w1*2*pi();
  A[3,3] = -(R1a + f*k);
  A[3,6] = k;
  A[4,1] = f*k;
  A[4,4] = -(R2b + k);
  A[5,2] = f*k;
  A[5,5] = -(R2b + k);
  A[5,6] = w1*2*pi();
  A[6,3] = f*k;
  A[6,5] = -w1*2*pi();
  A[6,6] = -(R1b + k);
  
  for(i in 1:N){
    A[1,2] = -xZ[i]*2*pi();
    A[2,1] = xZ[i]*2*pi();
    A[4,5] = (-xZ[i] + dw)*2*pi();
    A[5,4] = (xZ[i] - dw)*2*pi();
    
    M = matrix_exp(A * tp) * (to_vector([0,0,1,0,0,f]) + A\b) - A\b;
    
    Z_hat[i] = M[3];
  }
}

model {
   R1b ~ normal(0,50);
   R2b ~ normal(27000,10000);
   f ~ cauchy(0.02,0.01);
   k ~ normal(285,10);
   
   Z ~ normal(Z_hat, sigma);
}

generated quantities {
  real Z_pred[N];
  
  Z_pred = normal_rng(Z_hat, sigma);
}
