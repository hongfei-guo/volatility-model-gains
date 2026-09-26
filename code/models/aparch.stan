functions {
  // E|Z|^r for a standardized Student-t variable Z with variance one.
  real std_t_abs_moment(real nu, real r) {
    return exp(
      0.5 * r * log(nu - 2)
      + lgamma(0.5 * (r + 1))
      + lgamma(0.5 * (nu - r))
      - 0.5 * log(pi())
      - lgamma(0.5 * nu)
    );
  }
}
// AR(1)-APARCH(1,1). The simplex is defined on effective persistence:
// alpha * E[(|Z|-gamma Z)^delta] + beta < 1.
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
  real<lower=1e-8> prior_asym_sd;
  real<lower=0.2> prior_delta_mean;
  real<lower=1e-8> prior_delta_sd;
  real<lower=0.21> delta_max;
}
parameters {
  real mu_y;
  real<lower=-0.999, upper=0.999> phi_y;
  real log_sd_bar;
  real<lower=-1, upper=1> gamma;
  real<lower=0.2, upper=delta_max> delta;
  simplex[3] effective_persistence;
  real<lower=0> nu_gap;
}
transformed parameters {
  // Smooth approximation to max(2, delta).  The softplus term is strictly
  // above max(0, delta-2), so nu is strictly above both 2 and delta without
  // imposing the unnecessarily strong restriction nu > delta + 2.
  real<lower=2> nu_floor = 2 + log1p_exp(20 * (delta - 2)) / 20;
  real<lower=2> nu = nu_floor + nu_gap;
  real<lower=0> abs_moment = std_t_abs_moment(nu, delta);
  real<lower=0> asym_moment = abs_moment
    * (pow(1 - gamma, delta) + pow(1 + gamma, delta)) / 2;
  real<lower=0> alpha = effective_persistence[1] / asym_moment;
  real<lower=0> beta = effective_persistence[2];
  real<lower=0> sd_delta_bar = exp(delta * log_sd_bar);
  real<lower=0> omega = sd_delta_bar * effective_persistence[3];
  vector[T] mean_y;
  vector[T] resid;
  vector<lower=0>[T] sd_delta;
  vector<lower=0>[T] sigma;

  mean_y[1] = mu_y;
  resid[1] = y[1] - mean_y[1];
  sd_delta[1] = sd_delta_bar;
  sigma[1] = exp(log_sd_bar);
  for (t in 2:T) {
    mean_y[t] = mu_y + phi_y * (y[t - 1] - mu_y);
    sd_delta[t] = omega
      + alpha * pow(abs(resid[t - 1]) - gamma * resid[t - 1], delta)
      + beta * sd_delta[t - 1];
    sigma[t] = pow(sd_delta[t], 1 / delta);
    resid[t] = y[t] - mean_y[t];
  }
}
model {
  mu_y ~ normal(0, prior_mean_sd);
  phi_y ~ normal(0, prior_ar_sd);
  log_sd_bar ~ normal(log(y_scale), prior_log_scale_sd);
  gamma ~ normal(0, prior_asym_sd);
  delta ~ normal(prior_delta_mean, prior_delta_sd);
  effective_persistence ~ dirichlet(rep_vector(prior_persistence_concentration, 3));
  nu_gap ~ exponential(prior_nu_rate);

  for (t in 1:T) {
    target += K * student_t_lpdf(
      y[t] | nu, mean_y[t], sigma[t] * sqrt((nu - 2) / nu)
    );
  }
}
