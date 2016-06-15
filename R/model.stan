// Stan model

data {
  // Number of observations
  int<lower=1> N;
  // Data
  int X[N];
  // Response
  int y[N];
}

parameters {
  // Define parameters to estimate
  vector[K] beta;
}

model {
  // Priors
  beta[0] <- gamma(1, 1);
  beta[1] <- gamma(1, 1);

  y ~ poisson(beta[0] + beta[1] * X);
}
