// Standard stationary AR(1)-GJR-GARCH(1,1), labelled GJR-U.
// The positive- and negative-shock ARCH coefficients are separately nonnegative,
// while lambda is free to be positive or negative. The simplex guarantees
// beta + (alpha_pos + alpha_neg)/2 < 1, equivalently
// alpha + beta + lambda/2 < 1, and omega > 0.
// The observed sign of the latest innovation is always used.
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
  real<lower=0> sigma2_bar;
  simplex[4] persistence;
  real<lower=0> nu_minus_two;
}
transformed parameters {
  // alpha_pos is the coefficient after a nonnegative shock; alpha_neg is the
  // total coefficient after a negative shock. Their difference is lambda.
  real<lower=0> alpha_pos = 2 * persistence[1];
  real<lower=0> alpha_neg = 2 * persistence[2];
  real<lower=0> beta = persistence[3];
  real alpha = alpha_pos;
  real lambda = alpha_neg - alpha_pos;
  real<lower=0> omega = sigma2_bar * persistence[4];
  real<lower=2> nu = 2 + nu_minus_two;
  vector[T] mean_y;
  vector[T] resid;
  vector<lower=0>[T] sigma2;

  mean_y[1] = mu_y;
  resid[1] = y[1] - mean_y[1];
  sigma2[1] = sigma2_bar;
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
  sigma2_bar ~ lognormal(log(square(y_scale)), prior_log_scale_sd);
  persistence ~ dirichlet(rep_vector(prior_persistence_concentration, 4));
  nu_minus_two ~ exponential(prior_nu_rate);

  for (t in 1:T) {
    target += K * student_t_lpdf(
      y[t] | nu, mean_y[t], sqrt(sigma2[t]) * sqrt((nu - 2) / nu)
    );
  }
}
