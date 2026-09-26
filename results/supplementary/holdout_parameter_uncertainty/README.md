# Holdout parameter-mixture comparisons

These files support Appendix A9.2 and its summary in Section 6.1.

- `daily_predictions.csv`: VaR and ES by market, date, model, prediction
  construction and node scope. Return observations are not included.
- `sensitivity_summary.csv`: paired means, HAC intervals, block Monte Carlo
  standard errors and precision assessments for P, I, I−P and P−C.
- `nested_estimates.csv`: comparisons across the five nested node counts.

Run `Rscript code/reproduce_results.R` from the archive root after preparing
the returns. The original C forecasts come from `forecasts/holdout.csv`.
Scientific definitions, fitted inputs, seeds and optional prediction-generation
commands are documented in `design/parameter_uncertainty.md` at the archive root.
The NIKKEI sensitivity increment retains its limited-precision qualification.
