// AR(1)-STAR-SV extension used in the empirical comparison.
//
// The shock regime is classified from raw lagged returns outside Stan:
//   0 = non-negative, 1 = mild negative, 2 = severe negative.
// The volatility transition, however, uses the standardized innovation from
// the common AR(1) return equation. This extends the published zero-mean
// specification with the conditional mean used by all eight models.
//
// Prediction is not performed in generated quantities. Monthly structural
// point estimates feed the independent daily particle filter, which preserves
// latent-state uncertainty without mixing the artificial data-clone posterior.
data {
  int<lower=2> T;
  vector[T] y;
  array[T] int<lower=0, upper=2> trigger;
  int<lower=1> K;
  real<lower=1e-8> y_scale;
  real<lower=1e-8> prior_mean_sd;
  real<lower=1e-8> prior_ar_sd;
  real<lower=1e-8> prior_log_scale_sd;
  real<lower=1e-8> prior_phi_h_sd;
  real<lower=1e-8> prior_sigma_eta_sd;
  real<lower=1e-8> prior_nu_rate;
  real<lower=1e-8> prior_asym_sd;
}
parameters {
  real mu_y;
  real<lower=-0.999, upper=0.999> phi_y;
  real mu_h;
  real<lower=-0.999, upper=0.999> phi_h;
  real<lower=0> sigma_eta;
  real<lower=0, upper=1> gamma_severe_magnitude;
  real<lower=0, upper=1> gamma_mild_fraction;
  real<lower=0> nu_minus_two;
  array[K] vector[T] z_h;
}
transformed parameters {
  real<lower=-1, upper=0> gamma_severe = -gamma_severe_magnitude;
  real<lower=-1, upper=0> gamma_mild = gamma_severe
    + gamma_severe_magnitude * gamma_mild_fraction;
  real<lower=2> nu = 2 + nu_minus_two;
}
model {
  vector[T] mean_y;
  vector[T] resid;

  mu_y ~ normal(0, prior_mean_sd);
  phi_y ~ normal(0, prior_ar_sd);
  mu_h ~ normal(2 * log(y_scale), prior_log_scale_sd);
  phi_h ~ normal(0.95, prior_phi_h_sd);
  sigma_eta ~ normal(0, prior_sigma_eta_sd);
  gamma_severe_magnitude ~ normal(0, prior_asym_sd);
  gamma_mild_fraction ~ beta(1, 1);
  nu_minus_two ~ exponential(prior_nu_rate);
  for (k in 1:K) z_h[k] ~ std_normal();

  mean_y[1] = mu_y;
  resid[1] = y[1] - mean_y[1];
  for (t in 2:T) {
    mean_y[t] = mu_y + phi_y * (y[t - 1] - mu_y);
    resid[t] = y[t] - mean_y[t];
  }

  for (k in 1:K) {
    vector[T] h;
    // Marin-Veiga's initialization is used because the exact stationary law
    // of the regime-dependent STAR transition has no closed Gaussian form.
    h[1] = mu_h + sigma_eta * z_h[k][1];
    for (t in 2:T) {
      real standardized_residual_previous = resid[t - 1] / exp(0.5 * h[t - 1]);
      real gamma_selected = trigger[t - 1] == 2
        ? gamma_severe
        : (trigger[t - 1] == 1 ? gamma_mild : 0.0);
      h[t] = mu_h + phi_h * (h[t - 1] - mu_h)
        + gamma_selected * standardized_residual_previous
        + sigma_eta * z_h[k][t];
    }
    y ~ student_t(
      nu, mean_y, exp(0.5 * h) * sqrt((nu - 2) / nu)
    );
  }
}
