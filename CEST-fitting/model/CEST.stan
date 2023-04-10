data {
  int N; // num. of measurements
  real<lower=0> R1a;
  real<lower=0> R2a;
  real<lower=0> w1;
  real dw;
  real<lower=0> tp;
  array[N] real xZ; // offsets measured in spectrum.
  array[N] real Z; // Z-spectrum values.
}

parameters {
  real<lower=0> R1b;
  real<lower=0> R2b;
  real<lower=0> k;
  real<lower=0, upper=1> f;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[6,6] A = rep_matrix(0,6,6);
  vector[6] b = to_vector([0,0,R1a,0,0,f*R1b]);
  vector[6] M0 = to_vector([0,0,1,0,0,f]);
  vector[6] M; // auxiliary vector
  array[N] real Z_tilde;
  
  A[1,1] = -(R2a + f * k);
  A[1,4] = k;
  A[2,2] = -(R2a + f * k);
  A[2,3] = w1 * 2 * pi();
  A[2,5] = k;
  A[3,2] = -w1 * 2 * pi();
  A[3,3] = -(R1a + f * k);
  A[3,6] = k;
  A[4,1] = f * k;
  A[4,4] = -(R2b + k);
  A[5,2] = f * k;
  A[5,5] = -(R2b + k);
  A[5,6] = w1 * 2 * pi();
  A[6,3] = f * k;
  A[6,5] = -w1 * 2 * pi();
  A[6,6] = -(R1b + k);
  
  for(i in 1:N){
    A[1,2] = -xZ[i] * 2 * pi();
    A[2,1] = xZ[i] * 2 * pi();
    A[4,5] = (-xZ[i] + dw) * 2 * pi();
    A[5,4] = (xZ[i] - dw) * 2 * pi();
    
    M = matrix_exp(A * tp) * (M0 + A\b) - A\b;
    
    Z_tilde[i] = M[3];
  }
}

model {
   R1b ~ cauchy(0,50);
   R2b ~ normal(27000,3000);
   f ~ cauchy(0.02,0.005);
   k ~ cauchy(285,10);
   sigma ~ cauchy(0,0.05);
   
   Z ~ normal(Z_tilde, sigma);
}

generated quantities {
  array[N] real Z_rep;
  
  Z_rep = normal_rng(Z_tilde, sigma);
}
