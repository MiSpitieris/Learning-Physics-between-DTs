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
  int<lower=1> N; //number of observations per individual
  int<lower=1> Ns;
  int<lower=1,upper=Ns> id[Ns]; //number of individuals
  matrix[Ns, N] x; //different across individuals (e.g. time)
  matrix[Ns, N] y; //input data
  // predictions
  int<lower=0> N_pred;
  matrix[Ns, N_pred] x_pred;
}
transformed data{
  int<lower=1> N_tot;
  N_tot=N+N_pred;
}
parameters {
  real<lower=0, upper=3.0> u_tilde[Ns]; // non-centered parameterization 

  real<lower=0, upper=4> alpha; // marginal sd (delta)
  real<lower=0, upper=6> rho;  // length scale (delta)
  
  real<lower=0> sigma; // noise sd
  real<lower=0.5, upper=1.8> mu; // Global mean for u
  real<lower=1, upper=2> tau; // Global sd for u
  // predictions
  matrix[Ns,N_pred] y_pred;
}
transformed parameters {
  real<lower=0> u[Ns];   // individual u
  // Non-centered parameterization 
  for (s in 1:Ns) {
    u[s] = mu + tau * u_tilde[s];
  }
}
model {
  matrix[Ns,N_tot] x_tot;
  matrix[N_tot, N_tot] cov[Ns];
  matrix[N_tot, N_tot] L_cov[Ns];
  matrix[Ns,N_tot] z;
  
  x_tot=append_col(x,x_pred);
  z= append_col(y,y_pred);
  for (s in 1:Ns) {
    cov[s] = cov_exp(to_vector(x_tot[s]), alpha, rho)+diag_matrix(rep_vector(sigma^2, N_tot));
    L_cov[s] = cholesky_decompose(cov[s]);
  }

  u_tilde ~ normal(0, 1);  // non-centered
  // likelihood
  for (i in 1:Ns){
        z[i] ~ multi_normal_cholesky(5*exp(-u[id[i]]*x_tot[i]), L_cov[id[i]]);

  }
}
