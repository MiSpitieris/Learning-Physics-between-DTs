data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
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
