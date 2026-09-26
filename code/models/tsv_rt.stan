// Return-triggered threshold stochastic volatility under the common AR(1)
// conditional-mean specification used by the model-comparison design.
// The regime is data-defined: y[t-1] < 0. Exact zero is non-negative.
data {
  int<lower=2> T;
  vector[T] y;
  array[T] int<lower=0, upper=1> negative_return;
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
  real<lower=0> delta_tsv;
  real<lower=0> nu_minus_two;
  array[K] vector[T] z_h;
}
transformed parameters {
  real<lower=2> nu = 2 + nu_minus_two;
}
model {
  vector[T] mean_y;
  array[K] vector[T] h;

  mu_y ~ normal(0, prior_mean_sd);
  phi_y ~ normal(0, prior_ar_sd);
  mu_h ~ normal(2 * log(y_scale), prior_log_scale_sd);
  phi_h ~ normal(0.95, prior_phi_h_sd);
  sigma_eta ~ normal(0, prior_sigma_eta_sd);
  delta_tsv ~ normal(0, prior_asym_sd);
  nu_minus_two ~ exponential(prior_nu_rate);

  mean_y[1] = mu_y;
  for (t in 2:T) {
    mean_y[t] = mu_y + phi_y * (y[t - 1] - mu_y);
  }

  for (k in 1:K) {
    z_h[k] ~ std_normal();
    // The centered sign term has variance 1/4 under a symmetric return law.
    // This Gaussian initialization matches its contribution to unconditional
    // state variance while leaving the transition equation exact.
    h[k][1] = mu_h
      + sqrt(square(sigma_eta) + square(delta_tsv) / 4)
        / sqrt(1 - square(phi_h) + 1e-8) * z_h[k][1];
    for (t in 2:T) {
      real threshold = negative_return[t - 1] - 0.5;
      h[k][t] = mu_h + phi_h * (h[k][t - 1] - mu_h)
        + delta_tsv * threshold
        + sigma_eta * z_h[k][t];
    }
    y ~ student_t(
      nu, mean_y, exp(0.5 * h[k]) * sqrt((nu - 2) / nu)
    );
  }
}
