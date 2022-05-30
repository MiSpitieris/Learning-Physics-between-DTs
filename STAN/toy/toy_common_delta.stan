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
  // vector[N] x; //same across all individuals (e.g. time)
  vector[N] x[Ns]; //different across individuals (e.g. time)
  row_vector[N] y[Ns];
}
parameters {
  real<lower=0, upper=3.0> u_tilde[Ns]; // non-centered parameterization 

  real<lower=0, upper=4> alpha; // marginal sd (delta)
  real<lower=0, upper=6> rho;  // length scale (delta)
  
  real<lower=0> sigma; // noise sd
  real<lower=0.5, upper=1.8> mu; // Global mean for u
  real<lower=1, upper=2> tau; // Global sd for u
}
transformed parameters {
  real<lower=0> u[Ns];   // individual u
  // Non-centered parameterization 
  for (s in 1:Ns) {
    u[s] = mu + tau * u_tilde[s];
  }
}
model {
  matrix[N, N] cov[Ns];
  matrix[N, N] L_cov[Ns];
  for (s in 1:Ns) {
    cov[s] = cov_exp(x[s], alpha, rho)+diag_matrix(rep_vector(sigma^2, N));
    L_cov[s] = cholesky_decompose(cov[s]);
  }

  u_tilde ~ normal(0, 1);  // non-centered
  // likelihood
  for (i in 1:Ns){
    y[i] ~ multi_normal_cholesky(5*exp(-u[id[i]]*x[i]), L_cov[i]);
  }
}
