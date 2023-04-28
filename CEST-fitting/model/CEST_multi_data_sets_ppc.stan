// Strategy: Same as CEST for one data set, however...
// all data sets are appended to each other,
// and inference is done over multiple measurements together.

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
}

generated quantities {
  real R2b = 27000 + 5000 * inv_Phi(uniform_rng(normal_cdf(0 | 27000, 5000),1));
  real R1b = uniform_rng(0.1,10);
  real k = 200 + 100 * inv_Phi(uniform_rng(normal_cdf(0 | 200, 100),1));
  real sigma = inv_gamma_rng(3,0.5);
  real f_mu = uniform_rng(0,0.1);
  real f_kappa = uniform_rng(2,100);
  real f = beta_proportion_rng(f_mu,f_kappa);

  vector[sum_N] Z_tilde;
  
  // Calculation process of Z_tilde
  {
    matrix[6,6] A = rep_matrix(0,6,6);
    int ii = 0; // auxiliary count variable
    vector[6] b = to_vector([0,0,R1a,0,0,f*R1b]);
    vector[6] M0 = to_vector([0,0,1,0,0,f]);
    vector[6] M; // auxiliary vector
    
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
  
  array[sum_N] real Z_sim = normal_rng(Z_tilde,sigma);
}
