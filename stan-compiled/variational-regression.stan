data {
  int<lower=0> n;
  vector[n] y;
  vector[n] x;
  real tau;
}

parameters {
  real b;
  real<lower=0> sigma;
}

model {
  target += -log(sigma);
  target += normal_lpdf(b | 0, sigma*tau);
  target += normal_lpdf(y | b*x, sigma);
}
