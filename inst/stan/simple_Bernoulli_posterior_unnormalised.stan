data{
  int<lower=1> N0;
  int<lower=0> y0;
  real<lower=0> c;
  real<lower=0> d;
  real<lower=0> eta;
  real<lower=0> nu;
  int<lower=1> N;
  int<lower=0> y;
}
parameters{
  real<lower=0, upper=1> theta;
  real<lower=0, upper=1> a_0;
}
transformed parameters{
  real logL = (y0 -1) * log(theta) + ((N0 - y0) - 1) * log1m(theta);
}
model{
  /*Power prior*/
  target += a_0 * logL;
  target += beta_lpdf(theta | c, d);
  target += beta_lpdf(a_0 | eta, nu);
  /* Likelihood */
  target += y * log(theta) + (N - y)  * log1m(theta);
}
