functions{
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
  int<lower=1> N; //number of total obervations
  int<lower=1> Ns; //number of individuals
  // int<lower=1> n_ns; //number of observations per individual 
  int<lower=1,upper=Ns> id[Ns]; //individual ids
  int<lower=1> id_ind[Ns]; //indicator for the segment function
  int<lower=1> id_size[Ns]; //individual data sizes
  // vector[N] x; //same across all individuals (e.g. time)
  vector[N] x; //input data
  row_vector[N] y; //output data
}
parameters {
  // non-centered parameterization parameters
  real<lower=0,upper=3.0> u_tilde[Ns];
  real rho_tilde[Ns];   // non-centered sd of rho (delta process)
  real alpha_tilde[Ns]; // non-centered sd of alpha (delta precess)
  real<lower=0> sigma;  // same noise across individuals
  
  // Global-level parameters for delta 
  real<lower=0> rho_m;   // median of individual log-normal
  real<lower=0> rho_s;   //sd of of individual log-normal
  real<lower=0> alpha_m; // median of alpha log-normal
  real<lower=0> alpha_s; //sd of alpha log-normal
  
  // Global-level parameters for u
  real<lower=0.5, upper=1.8> mu; 
  real<lower=1, upper=2> tau;
}
transformed parameters {
  real<lower=0> u[Ns];   // physical parameters
  real<lower=0> rho[Ns];   // length scale
  real<lower=0> alpha[Ns];  // marginal standard deviation
  // Non-centered parameterization of individual parameters
  for (s in 1:Ns) {
    rho[s] = exp(log(rho_m) + rho_s * rho_tilde[s]);
    alpha[s] = exp(log(alpha_m) + alpha_s * alpha_tilde[s]);
    u[s] = mu + tau * u_tilde[s];
  }
}
model {
  // matrix[N, N] cov[Ns];
  // matrix[N, N] L_cov[Ns];
  // for (s in 1:Ns) { // individual covariance 
  //   cov[s] = cov_exp(x[s], alpha[s], rho[s])+diag_matrix(rep_vector(sigma^2, N));
  //   L_cov[s] = cholesky_decompose(cov[s]);
  // }
  
  // priors
  // Global parameters
  rho_m ~ inv_gamma(2, 0.5);
  alpha_m ~ normal(0,2);
  rho_s ~ normal(0, 0.5);
  alpha_s ~ normal(0, 0.5);

  // non-centered parameterization of individual parameters
  rho_tilde ~ normal(0, 1); 
  alpha_tilde ~ normal(0, 1); 
  u_tilde ~ normal(0, 1);  

  // likelihood
  // for (i in 1:Ns){
  //   y[i] ~ multi_normal_cholesky(5*exp(-u[id[i]]*x[i]), L_cov[id[i]]);
  // }
  //likelihood
  for (i in 1:Ns){
    vector[id_size[i]] x_s;
    row_vector[id_size[i]] y_s;
    matrix[id_size[i], id_size[i]] cov;
    matrix[id_size[i], id_size[i]] L_cov;
    
    x_s = segment(x, id_ind[i], id_size[i]);
    y_s = segment(y, id_ind[i], id_size[i]);
    cov = cov_exp(x_s, alpha[i], rho[i])+diag_matrix(rep_vector(sigma^2, id_size[i]));
    L_cov = cholesky_decompose(cov);
    y_s ~ multi_normal_cholesky(5*exp(-u[id[i]]*x_s), L_cov);
  }
}
