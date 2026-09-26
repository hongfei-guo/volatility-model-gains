// Restricted boundary AR(1)-GJR-GARCH(1,1), labelled GJR-R.
// lambda = 2 * (1 - alpha - beta), hence alpha + beta + lambda/2 = 1.
// This specification is treated as a separate boundary benchmark.
// Forecasting must still use the observed sign of the latest innovation.
data {
  int<lower=2> T;
  vector[T] y;
  int<lower=1> K;
  real<lower=1e-8> y_scale;
  real<lower=1e-8> prior_mean_sd;
  real<lower=1e-8> prior_ar_sd;
  real<lower=1e-8> prior_log_scale_sd;
  real<lower=1e-8> prior_persistence_concentration;
  real<lower=1e-8> prior_nu_rate;
}
parameters {
  real mu_y;
  real<lower=-0.999, upper=0.999> phi_y;
  real log_omega;
  simplex[3] persistence;
  real<lower=0> nu_minus_two;
}
transformed parameters {
  real<lower=0> alpha = persistence[1];
  real<lower=0> beta = persistence[2];
  real<lower=0> lambda = 2 * persistence[3];
  real<lower=0> omega = exp(log_omega);
  real<lower=2> nu = 2 + nu_minus_two;
  vector[T] mean_y;
  vector[T] resid;
  vector<lower=0>[T] sigma2;

  mean_y[1] = mu_y;
  resid[1] = y[1] - mean_y[1];
  sigma2[1] = square(y_scale);
  for (t in 2:T) {
    real neg = resid[t - 1] < 0 ? 1.0 : 0.0;
    mean_y[t] = mu_y + phi_y * (y[t - 1] - mu_y);
    sigma2[t] = omega
                + alpha * square(resid[t - 1])
                + lambda * neg * square(resid[t - 1])
                + beta * sigma2[t - 1];
    resid[t] = y[t] - mean_y[t];
  }
}
model {
  mu_y ~ normal(0, prior_mean_sd);
  phi_y ~ normal(0, prior_ar_sd);
  log_omega ~ normal(log(0.01 * square(y_scale) + 1e-8), prior_log_scale_sd);
  persistence ~ dirichlet(rep_vector(prior_persistence_concentration, 3));
  nu_minus_two ~ exponential(prior_nu_rate);

  for (t in 1:T) {
    target += K * student_t_lpdf(
      y[t] | nu, mean_y[t], sqrt(sigma2[t]) * sqrt((nu - 2) / nu)
    );
  }
}
