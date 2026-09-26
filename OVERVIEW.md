# Testing Whether Volatility Model Gains Persist

Research by Hongfei Guo, J. Miguel Marín and Helena Veiga.

This study compares eight volatility models across five equity indices and tests whether forecast advantages persist in a prespecified later holdout. Evaluation covers joint VaR–ES accuracy, calibration, predictive density and forecast availability.

[UC3M Working Paper 2026-13](https://hdl.handle.net/10016/50798)

## Explore the code

- `code/models/`: eight Stan model specifications.
- `code/evaluation_functions.R`: forecast evaluation calculations.
- `code/reproduce_results.R`: result reproduction entry point.
- `code/render_exhibits.R`: manuscript exhibits.

See [README.md](README.md) for dependencies, run order and data requirements. The default reproduction starts from supplied forecasts and lawfully obtained source index observations; it does not rerun monthly estimation of all models. Third-party index observations are not redistributed.
