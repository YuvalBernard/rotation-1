functions {
  vector cest(vector pars, vector w1,
              data array[] real x_r, data array[] int x_i) {
    matrix[6,6] A = rep_matrix(0,6,6);
    A[1,1] = -(R2a + f * k);
    A[1,4] = k;
    A[2,2] = -(R2a + f * k);
    A[2,5] = k;
    A[3,3] = -(R1a + f * k);
    A[3,6] = k;
    A[4,1] = f * k;
    A[4,4] = -(R2b + k);
    A[5,2] = f * k;
    A[5,5] = -(R2b + k);
    A[6,3] = f * k;
    A[6,6] = -(R1b + k);
    
  }
}

pars = to_vector([R2b,R1b,f,k,sigma,dw_rad]); // or append_row

data {
  int K; // number of experiments.
  array[K] int N; // Specifies number of measurements in K different experiments.
  real<lower=0> R1a;
  real<lower=0> R2a;
  vector[K] w1; // Different RF field intensities for different experiments.
  real dw;
  real<lower=0> tp;
  vector[sum(N)] xZ; // offsets measured in spectrum.
  vector[sum(N)] Z; // Z-spectrum values.
}

transformed data {
  vector[K] w1_rad = w1 * 2 * pi();
  real dw_rad = dw * 2 * pi();
  int sum_N = sum(N);
  vector[sum_N] xZ_rad = xZ * 2 * pi();
  real max_Z = max(Z);
}

parameters {
  real<lower=0> R2b_std;
  real<lower=0> R1b_std;
  real<lower=0> k_std;
  // real<lower=0> f_std;
  // real<lower=0> sigma_std;
  real<lower=0> sigma;
  real<lower=0> f;
}

transformed parameters {
  real R2b = 27000 + 5000 * R2b_std;
  real R1b = 5 * R1b_std;
  real k = 150 + 100 * k_std;
  // real f = 0.015 + 0.005 * f_std;
  // real sigma = 0.05 * sigma_std;
  
  vector<lower=0, upper=max_Z>[sum_N] Z_tilde;
  
  // Calculation process of Z_tilde need not be exposed to other blocks
  {
    matrix[6,6] A = rep_matrix(0,6,6);
    int ii = 0; // auxiliary count variable
    
    A[1,1] = -(R2a + f * k);
    A[1,4] = k;
    A[2,2] = -(R2a + f * k);
    A[2,5] = k;
    A[3,3] = -(R1a + f * k);
    A[3,6] = k;
    A[4,1] = f * k;
    A[4,4] = -(R2b + k);
    A[5,2] = f * k;
    A[5,5] = -(R2b + k);
    A[6,3] = f * k;
    A[6,6] = -(R1b + k);
    
    for(j in 1:K){
      A[2,3] = w1_rad[j];
      A[3,2] = -w1_rad[j];
      A[5,6] = w1_rad[j];
      A[6,5] = -w1_rad[j];
      
      for(i in 1:N[j]){
        vector[6] b = to_vector([0,0,R1a,0,0,f*R1b]);
        vector[6] M0 = to_vector([0,0,1,0,0,f]);
        vector[6] M; // auxiliary vector
    
        A[1,2] = -xZ_rad[i+ii];
        A[2,1] = xZ_rad[i+ii];
        A[4,5] = (-xZ_rad[i+ii] + dw_rad);
        A[5,4] = (xZ_rad[i+ii] - dw_rad);
        
        vector[6] A_inv_b = A\b;
    
        M = matrix_exp(A * tp) * (M0 + A_inv_b) - A_inv_b;
    
        Z_tilde[i+ii] = M[3];
      }
      ii += N[j];
    }
  }

}

model {
  R2b_std ~ std_normal(); 
  R1b_std ~ std_normal();
  // f_std ~ std_normal();
  f ~ std_normal();
  k_std ~ std_normal();
  // sigma_std ~ std_normal();
  sigma ~ inv_gamma(3, 0.5);
  Z ~ normal(Z_tilde, sigma) T[0, max_Z];
}

generated quantities {
  vector[sum_N] Z_rep;
  // generate random numbers from normal(Z_tilde, sigma) truncated between 0 and max_Z
  {
    vector[sum_N] p_lb;
    vector[sum_N] p_ub;
  
    for(i in 1:sum_N){
      p_lb[i] = normal_cdf(0 | Z_tilde[i], sigma);
      p_ub[i] = normal_cdf(max_Z | Z_tilde[i], sigma);
    }
    vector[sum_N] u = to_vector(uniform_rng(p_lb, p_ub));
    Z_rep = Z_tilde + sigma * std_normal_qf(u);
  }
}
