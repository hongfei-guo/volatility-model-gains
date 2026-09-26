# Holdout design

The scientific rules below were fixed on 6 August 2026, before holdout
forecast generation.

## Purpose

The 2024–2025 holdout assesses whether forecast-performance differences
observed during the 2018–2023 main out-of-sample period persist in a later,
non-overlapping sample. The model set, forecast construction, evaluation
criteria, and treatment of unavailable forecasts were fixed before holdout
forecast generation began.

## Models and forecasts

The analysis covers DAX, FTSE 100, Nikkei 225, the SSE Composite, and the
S&P 500. It compares GARCH, GJR-R, GJR-U, APARCH, ARSV, AARSV, TSV-RT, and
STAR-SV. Each model uses the most recent 2,000 observations, monthly
structural-parameter refits, daily state updates, and one-day-ahead forecasts.
The tail levels are 95%, 97.5%, and 99%.

## Primary inference

The two primary holdout tests compare STAR-SV with ARSV at the 97.5% tail
level in DAX and NIKKEI. For each comparison, the reported quantity is the
mean daily FZ0 difference, STAR-SV minus ARSV. Negative values favour STAR-SV.
Paired HAC inference yields a two-sided p-value and a 95% confidence interval.
The two p-values form one family and receive Holm adjustment at the 5% level.

Only a negative mean difference with a Holm-adjusted p-value below 0.05
corroborates the corresponding main-period advantage. A positive mean
difference with a Holm-adjusted p-value below 0.05 is a statistical reversal.
All other outcomes are non-corroborating; retention or change of the sign is
reported descriptively.

## Additional evidence

The remaining model comparisons are secondary. Their paired HAC p-values are
nominal and unadjusted and do not alter the primary holdout conclusion. The
95% and 99% holdout results are supplementary.

Three forecast properties are evaluated separately:

- relative VaR–ES accuracy using FZ0, rankings, and Model Confidence Sets;
- model-specific calibration using coverage, independence, shortfall, and
  conditional-calibration diagnostics;
- predictive-density accuracy using LPDS and CRPS.

No composite winner is formed across these criteria.

## Evaluation samples and unavailable forecasts

Cross-model comparisons, rankings, and Model Confidence Sets use the
within-market intersection of valid forecast dates across all eight models.
Model-specific calibration uses every valid forecast date for the model being
assessed. A model-month without a valid forecast remains unavailable. Missing
forecasts are not interpolated, backfilled, or assigned an artificial loss.

The archive records the exact evaluation dates, unavailable forecast dates,
random seeds, bootstrap settings, and simulated Z2 reference values used in
the reported analysis.
