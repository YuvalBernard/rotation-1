data {
  int N; // num. of measurements
  real<lower=0> R1a;
  real<lower=0> R2a;
  real<lower=0> w1;
  real dw;
  real<lower=0> tp;
  vector[N] xZ; // offsets measured in spectrum.
  vector<lower=0, upper=1>[N] Z; // Z-spectrum values.
}

transformed data {
  real w1_rad = w1 * 2 * pi();
  real dw_rad = dw * 2 * pi();
  vector[N] xZ_rad = xZ * 2 * pi();
}

parameters {
  real<lower=0, upper = pi() / 2> R1b_unif;
  real<lower=0> R2b_std;
  real<lower=0, upper = pi() / 2> k_unif;
  real<lower=0, upper = pi() / 2> f_unif;
  real<lower=0, upper = pi() / 2> sigma_unif;
}

transformed parameters {
  real R1b = 50 * tan(R1b_unif); // R1b ~ cauchy(0,50)
  real R2b = 27000 + 3000 * R2b_std; // R2b ~ normal(27000,3000)
  real k = 285 + 10 * tan(k_unif); // k ~ cauchy(285,10)
  real f = 0.02 + 0.005 * tan(f_unif); // f ~ cauchy(0.02,0.005) 
  real sigma = 0.05 * tan(sigma_unif); // sigma ~ cauchy(0,0.05)
  
  vector<lower=0, upper=1.05>[N] Z_tilde;
  
  // Calculation process of Z_tilde need not be exposed to other blocks
  {
    matrix[6,6] A = rep_matrix(0,6,6);
  
    A[1,1] = -(R2a + f * k);
    A[1,4] = k;
    A[2,2] = -(R2a + f * k);
    A[2,3] = w1_rad;
    A[2,5] = k;
    A[3,2] = -w1_rad;
    A[3,3] = -(R1a + f * k);
    A[3,6] = k;
    A[4,1] = f * k;
    A[4,4] = -(R2b + k);
    A[5,2] = f * k;
    A[5,5] = -(R2b + k);
    A[5,6] = w1_rad;
    A[6,3] = f * k;
    A[6,5] = -w1_rad;
    A[6,6] = -(R1b + k);
  
    for(i in 1:N){
      vector[6] b = to_vector([0,0,R1a,0,0,f*R1b]);
      vector[6] M0 = to_vector([0,0,1,0,0,f]);
      vector[6] M; // auxiliary vector
      vector[6] A_inv_b = A\b;
    
      A[1,2] = -xZ_rad[i];
      A[2,1] = xZ_rad[i];
      A[4,5] = (-xZ_rad[i] + dw_rad);
      A[5,4] = (xZ_rad[i] - dw_rad);
    
      M = matrix_exp(A * tp) * (M0 + A_inv_b) - A_inv_b;
    
      Z_tilde[i] = M[3];
    }
  }

}

model {
  R2b_std ~ std_normal(); 
  
  Z ~ normal(Z_tilde, sigma) T[0, 1];
}

generated quantities {
  vector[N] Z_rep;
  // generate random numbers from normal(Z_tilde, sigma) truncated between 0 and 1
  {
    vector[N] p_lb;
    vector[N] p_ub;
  
    for(i in 1:N){
      p_lb[i] = normal_cdf(0 | Z_tilde[i], sigma);
      p_ub[i] = normal_cdf(1 | Z_tilde[i], sigma);
    }
    vector[N] u = to_vector(uniform_rng(p_lb, p_ub));
    Z_rep = Z_tilde + sigma * std_normal_qf(u);
  }
}
