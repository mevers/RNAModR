// Stan model

data {
  // Number of observations
  int<lower=1> N;
  // Number of parameters
  int<lower=1> K;
  // Data
  matrix[N,K] X;
  // Response
  vector[N] y;
}

parameters {
  // Define parameters to estimate
  vector[K] beta;
  real<lower=0> sigma;
}

model {
  y ~ lognormal(X * beta, sigma);
}
