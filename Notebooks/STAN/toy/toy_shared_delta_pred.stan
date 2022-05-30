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
  int<lower=1> N; // number of observations per individual
  int<lower=1> Ns; // number of individuals 
  int<lower=1,upper=Ns> id[Ns]; // individual id
  matrix[Ns, N] x; // individual input vector
  matrix[Ns, N] y; // matrix of all individual outputs
  // predictions
  int<lower=0> N_pred;
  matrix[Ns, N_pred] x_pred;
}
transformed data{
  int<lower=1> N_tot;
  N_tot=N+N_pred;
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
  // predictions
  matrix[Ns,N_pred] y_pred;
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
  matrix[Ns, N_tot] x_tot;
  matrix[N_tot, N_tot] cov[Ns];
  matrix[N_tot, N_tot] L_cov[Ns];
  matrix[Ns, N_tot] z;
  
  x_tot=append_col(x,x_pred);
  z= append_col(y,y_pred);
  for (s in 1:Ns) {
    cov[s] = cov_exp(to_vector(x_tot[s]), alpha[s], rho[s])+diag_matrix(rep_vector(sigma^2, N_tot));
    L_cov[s] = cholesky_decompose(cov[s]);
  }
  
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
  for (i in 1:Ns){
    z[i] ~ multi_normal_cholesky(5*exp(-u[id[i]]*x_tot[i]), L_cov[id[i]]);
  }
}
