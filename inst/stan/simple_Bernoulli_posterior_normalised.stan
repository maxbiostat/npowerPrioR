data{
  int<lower=1> N0;
  int<lower=0> y0;
  real<lower=0> c;
  real<lower=0> d;
  real<lower=0> nu;
  real<lower=0> eta;
  int<lower=1> N;
  int<lower=0> y;
}
parameters{
  real<lower=0,upper=1> theta;
  real<lower=0,upper=1> a_0;
}
model{
  /*Power prior*/
  target += beta_lpdf(theta | a_0*y0 + c, a_0 *(N0-y0) + d);
  target += beta_lpdf(a_0 | eta, nu);
  /* Likelihood */
  target += y * log(theta) + (N - y) * log1m(theta);
}
