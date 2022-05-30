data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
  
  int<lower=0>N_pred;
  vector[N_pred] x_pred;
}
parameters {
  real<lower=0,upper=5> u;
  real<lower=0, upper=2> sigma;
}

model {
  // priors
  u ~ normal(1,2);
  // likelihood
  y ~ normal(5*exp(-u*x), sigma);
}
generated quantities {
  vector[N_pred] y_pred;
  for(n in 1:N_pred){
    y_pred[n] = normal_rng(5*exp(-u*x_pred[n]), sigma);
  }
}

