# DAX parameter-uncertainty sensitivity

This directory contains the machine-readable results supporting Online
Appendix A9.1. The supplementary analysis compares the 2018--2023 DAX
predictive-density forecasts from STAR-SV and ARSV after allowing for local
structural-parameter uncertainty. It uses only the main-OOS sample and does
not enter the holdout inference.

## Files

- `sensitivity_summary.csv`: paper-facing raw scores, STAR-SV advantages,
  relative changes, sampling and Monte Carlo uncertainty components, combined
  descriptive bands, and direction indicators.
- `daily_density_scores.csv`: date-level LPDS and CRPS values for the main,
  plug-in, and parameter-integrated forecast variants.
- `parameter_draws.csv`: inspection material recording the parameter nodes,
  block assignments, seeds, and model coordinates used in the integration.
  The default reproduction starts from `daily_density_scores.csv` and does
  not consume this file.
## Reported descriptive band

For score (s), let \(\widehat\rho_s\) denote the relative change in the
STAR-SV advantage. The reported combined descriptive band is

\[
\widehat\rho_s \pm 1.96
\sqrt{\widehat{\operatorname{se}}_{\mathrm{HAC},s}^{2}
      +\widehat{\operatorname{se}}_{\mathrm{MC},s}^{2}}.
\]

The HAC component uses the joint covariance of the integrated-minus-plug-in
effect and the plug-in advantage. The Monte Carlo component uses the
12-block delete-one-block jackknife. The resulting band is a descriptive
summary; it is neither a frequentist confidence interval nor a decision
threshold.

The project-level run instructions are in `../../../README.md`. The summary is
rebuilt from `daily_density_scores.csv` by
`../../../code/reproduce_results.R`, using the software recorded under
`../../../environment/`.

The prediction-generating code and fitted inputs are described in
`../../../design/parameter_uncertainty.md`.
