functions{ // squared exponential kernel
  matrix cov_exp(vector x,
                 real alpha,
                 real rho){
    int n = rows(x);
    matrix[n, n] K;
    // KP
    for (i in 1:(n)){
      K[i,i] = pow(alpha, 0.2e1);
      for (j in (i+1):n){
        K[i,j] = exp(-pow(x[i] - x[j], 0.2e1) * pow(rho, -0.2e1));
        K[i,j] = pow(alpha, 0.2e1) * K[i,j];
        K[j,i] = K[i,j];
      }
      K[n,n] = pow(alpha, 0.2e1);
    }
    return(K);
  }
}

data {
  int<lower=0> N; 
  vector[N] x;
  vector[N] y;
  
  int<lower=0> N_pred;
  vector[N_pred] x_pred;
}
parameters {
  real<lower=0,upper=5> u; // physical parameter
  real<lower=0, upper=2> sigma; // noise parameter
  real<lower=0, upper=4> alpha; // marginal sd (delta)
  real<lower=0, upper=6> rho; // length scale (delta)
}

model {
  matrix[N, N] cov = cov_exp(x, alpha, rho)+diag_matrix(rep_vector(sigma^2, N));
  matrix[N, N] L_cov = cholesky_decompose(cov);
  // priors
  u ~ normal(1,2); 
  y ~ multi_normal_cholesky(5*exp(-u*x), L_cov);
}
