functions{
  int which_min(real [] y ){
    int ans = sort_indices_asc(y)[1];
    return(ans);
  }
  real approximate_ca0(real x, real[] x_pred, real[] y_pred){
    int K = size(x_pred);
    real deltas [K];
    real ans;
    int i;
    if(size(y_pred) != K) reject("x_pred and y_pred aren't of the same size");
    for(k in 1:K) deltas[k] = fabs(x_pred[k] - x);
    i = which_min(deltas);
    if(i != 1){
      real x1 = x_pred[i];
      real x2 = x_pred[i + 1];
      real y1 = y_pred[i];
      real y2 = y_pred[i + 1];
      ans = y1 + (y2-y1) * (x-x1)/(x2-x1);
    }else{
      ans = y_pred[i];
    }
    return(ans);
  }
}
data{
  int<lower=1> N0;
  int<lower=0> y0;
  real<lower=0> c;
  real<lower=0> d;
  real<lower=0> eta;
  real<lower=0> nu;
  int<lower=1> N;
  int<lower=0> y;
  /* Approximation stuff*/
  int<lower=0> K;
  real pred_grid_x[K];
  real pred_grid_y[K];
}
parameters{
  real<lower=0, upper=1> theta;
  real<lower=0, upper=1> a_0;
}
model{
  /*Power prior*/
  target += -approximate_ca0(a_0, pred_grid_x, pred_grid_y);
  target += (a_0 * y0 + c -1) * log(theta) + (a_0 * (N0 - y0) + d - 1) * log1m(theta);
  target += beta_lpdf(theta | c, d);
  target += beta_lpdf(a_0 | eta, nu);
  /* Likelihood */
  target += y * log(theta) + (N - y)  * log1m(theta);
}
