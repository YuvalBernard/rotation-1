functions { 
   // function declarations and definitions
   real evaluate_Z(real offset, real R1a, real R2a, real R1b, real R2b, real f, real k, real dw, real tp, real w1) {
       real ka;
       vector[6] M;  
       matrix[6,6] A;  
       vector[6] b;
       ka = f*k;
       b = [0, 0, R1a, 0, 0, f*R1b]'; 
       A = [ [-(R2a + ka), offset, 0, k, 0, 0],
             [-offset, -(R2a + ka), w1, 0, k, 0],
             [0, -w1, -(R1a + ka), 0, 0, k],              
             [ka, 0, 0, -(R2b + k), offset + dw, 0],   
             [0, ka, 0, -(offset + dw), -(R2b + k), w1],
             [0, 0, ka, 0, -w1, -(R1b + k)] ];
       M = matrix_exp(A*tp) * (b + A\b) - A\b;
       return M[3];
   }
}
data {
   int<lower=0> N; // number of measurements   
   real<lower=0> w1; // RF field intensity (in Hz)   
   real<lower=0> tp; // saturation time (in s)
   vector[N] xZ; // offsets measured in spectrum   
   vector[N] Z; // experimental Z spectrum data    
}
parameters {
   real<lower=0> R1a;
   real<lower=0> R2a;
   real<lower=0> R1b;
   real<lower=0> R2b;
   real<lower=0, upper=1> f;
   real<lower=0> k;
   real<upper=0> dw;
   real<lower=0> sigma;
   vector[N] ZPred;
}

model { 
  //priors and likelihood go here.
   vector[N] Zmean;
   R1a ~ normal(8,4.0/3);    
   R2a ~ normal(393,10.0/3);    
   R1b ~ uniform(0.01,1);   
   R2b ~ normal(28000,4000.0/3); 
   f ~ uniform(0,1);   
   k ~ normal(285,5);  
   dw ~ normal(-260,2);
   sigma ~ cauchy(0,1);
   for(i in 1:N){    
       Zmean[i] = evaluate_Z(xZ[i]*2*pi(),R1a,R2a,R1b,R2b,f,k,dw*2*pi(),tp,w1*2*pi());
       Z[i] ~ normal(Zmean[i],sigma);
       ZPred[i] ~ normal(Zmean[i],sigma);
   }
}
