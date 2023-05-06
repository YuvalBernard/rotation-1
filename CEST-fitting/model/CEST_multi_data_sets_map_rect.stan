/**
  * Fit each data set (with different w1) to Bloch-McConnell
  * equations. Use Map-rect function to parallelize calculation.
  * Each experiment is contained in a shard.
  * Warning: Program assumes each experiment has
  *          the same number of measurement points.
*/

functions {
  vector cest(vector pars, vector w1_rad,
              data array[] real x_r, data array[] int x_i) {
    // Specify parameters shared accros shards.
    real R2b = pars[1];
    real R1b = pars[2];
    real f = pars[3];
    real k = pars[4];
    real sigma = pars[5];
    real dw_rad = pars[6];
    real R1a = pars[7];
    real R2a = pars[8];
    real tp = pars[9];
    // Extract variates and covariates.
    int N = x_i[1];
    vector[N] xZ_rad = to_vector(x_r[1:N]);
    vector[N] Z = to_vector(x_r[N+1:2*N]);
    // Model data to Bloch-McConnell equations.
    // Define relevant matrices and vectors.
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
    A[2,3] = w1_rad[1];
    A[3,2] = -w1_rad[1];
    A[5,6] = w1_rad[1];
    A[6,5] = -w1_rad[1];
    vector[6] b = to_vector([0,0,R1a,0,0,f*R1b]);
    vector[6] M0 = to_vector([0,0,1,0,0,f]);
    vector[6] M; // auxiliary vector.
    vector[N] Z_tilde;
    // Solve BM equations using matrix exponential.
    for (i in 1:N) {
      A[1,2] = -xZ_rad[i];
      A[2,1] = xZ_rad[i];
      A[4,5] = (-xZ_rad[i] + dw_rad);
      A[5,4] = (xZ_rad[i] - dw_rad);
      vector[6] A_inv_b = A\b;
      M = matrix_exp(A * tp) * (M0 + A_inv_b) - A_inv_b;
      Z_tilde[i] = M[3];
    }
    // Calculate log-posterior
    // vector[N] lp = Z_tilde + sigma * inv_Phi(uniform_lpdf(Z | normal_cdf(0 | Z_tilde, sigma),
    //                                             normal_cdf(max(Z) | Z_tilde, sigma)));
    real lp = normal_lpdf(Z | Z_tilde, sigma);
    return [lp]';
  }
}

data {
  int K; // number of experiments.
  int N; // Specifies number of measurements in each experiment.
  array[K,1] int x_i; // created for bullshit semantics of map_rect.
  real<lower=0> R1a;
  real<lower=0> R2a;
  vector[K] w1; // Different RF field intensities for different experiments.
  real dw;
  real<lower=0> tp;
  vector[N] xZ; // offsets measured in spectrum.
  vector[N * K] Z; // Z-spectrum values.
}

transformed data {
  // Convert data to rad/s and define useful constants.
  array[K] vector[1] w1_rad;
  for (i in 1:K) {
    w1_rad[i][1] = w1[i] * 2 * pi();
  }
  real dw_rad = dw * 2 * pi();
  vector[N] xZ_rad = xZ * 2 * pi();
  // Separate data into shards.
  array[K, 2 * N] real x_r; // holds xZ and Z
  {
    int pos = 1;
    for (k in 1:K) {
      int end = pos + N - 1;
      x_r[k] = to_array_1d(append_row(xZ, Z[pos:end]));
      pos += N;
    }
  }
}

parameters {
  real R2b_std;
  real<lower=0> R1b_std;
  real k_std;
  real<lower=0> sigma;
  real<lower=0> f;
}

transformed parameters {
  real R2b = 27000 + 5000 * R2b_std;
  real R1b = 5 * R1b_std;
  real k = 200 + 100 * k_std;
}

model {
  R2b_std ~ std_normal(); 
  R1b_std ~ std_normal();
  f ~ std_normal();
  k_std ~ std_normal();
  sigma ~ std_normal();
  target += sum(map_rect(cest, 
                         to_vector([R2b,R1b,f,k,sigma,dw_rad,R1a,R2a,tp]),
                         w1_rad, x_r, x_i));
}

// generated quantities {
//   vector[sum_N] Z_rep;
//   // generate random numbers from normal(Z_tilde, sigma) truncated between 0 and max_Z
//   {
//     vector[sum_N] p_lb;
//     vector[sum_N] p_ub;
//   
//     for(i in 1:sum_N){
//       p_lb[i] = normal_cdf(0 | Z_tilde[i], sigma);
//       p_ub[i] = normal_cdf(max_Z | Z_tilde[i], sigma);
//     }
//     vector[sum_N] u = to_vector(uniform_rng(p_lb, p_ub));
//     Z_rep = Z_tilde + sigma * std_normal_qf(u);
//   }
// }