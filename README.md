# Testing Whether Volatility Model Gains Persist

*Testing Whether Volatility Model Gains Persist: A Prespecified Holdout in Tail Risk Forecasting*

Hongfei Guo, J. Miguel Marín, and Helena Veiga

Do forecast advantages survive in a later sample? This study compares eight
volatility models across five equity indices, separating the discovery of
forecast gains from their evaluation in a prespecified subsequent holdout.
It assesses joint VaR–ES accuracy, calibration, predictive density and forecast
availability. This repository contains the R evaluation code, Stan model
specifications and saved forecasts used to reproduce the reported comparisons.

[Research overview](OVERVIEW.md) · [Working paper](https://hdl.handle.net/10016/50798) · [Archived replication materials](https://doi.org/10.5281/zenodo.22657288)

## Scope

The code reproduces the forecast-evaluation results, main-paper tables,
Figures 1 and 2, and evaluation-based online-appendix results from the supplied daily
forecast panels. It reconstructs the common and model-specific evaluation
samples; recomputes FZ0, paired HAC inference, Holm adjustment, calibration
diagnostics, density summaries, rankings, and Model Confidence Sets; and
renders the main exhibits.

The executable reproduction begins with the model forecasts used in the paper
and source index data obtained independently by the replicator. Monthly
estimation of all eight models is computationally intensive and is not part of
the default run. The Stan files under `code/models/` document the estimated
specifications. The supplied monthly parameter estimates, diagnostics,
forecast seeds, and forecast panels provide the numerical inputs needed to
inspect the reported analysis.

## Repository contents

- `data/`: source-data documentation and the required user-supplied file
  structure. Third-party index observations are not redistributed.
- `forecasts/`: model forecasts for the 2018–2023 main out-of-sample period
  and the 2024–2025 holdout. The `crps` column contains the CRPS values used in
  the paper; row-level realized returns are not included.
- `design/`: the holdout design, MCS settings, and simulated Z2 reference
  values.
- `code/`: one return-construction check, one result-reproduction entry point,
  supporting functions, exhibit code, the eight Stan specifications, and the
  parameter-mixture prediction calculations.
- `results/`: machine-readable values reported in the paper and online
  appendix.
- `inputs/sensitivity/`: fitted point estimates, joint covariance matrices,
  and training dates for the two parameter-mixture analyses.
- `environment/`: the dependency lockfile for the executable reproduction.

## Data availability

The daily index observations were obtained from Yahoo Finance and are not
redistributed because they remain subject to the provider's terms. The ticker
map, sample periods, required input schema, and return transformation are
documented in `data/README.md`. Replicators must obtain the source observations
directly from Yahoo Finance or another licensed provider.

The released LPDS, CRPS, and PIT columns are model-evaluation outputs computed
using the observations used in the paper. Return-dependent FZ0 and coverage
statistics are reconstructed from the replicator's source data. Exact numerical
agreement therefore requires the same historical closing levels; later provider
revisions can prevent bitwise reproduction.

## Software

The R dependency environment is recorded in `environment/renv.lock`. To
restore it, install `renv` and run the following command from the repository root:

```text
Rscript -e 'renv::restore(lockfile = "environment/renv.lock")'
```

Monthly model estimation used CmdStan 2.36.0 through cmdstanr 0.9.0. These
dependencies are not required for the default forecast-evaluation reproduction.

## Reproduction

Obtain the five source files described in `data/README.md` and place them under
`data/source/`. Then run the following commands from the repository root:

```text
Rscript code/build_returns.R
Rscript code/reproduce_results.R
```

The first command constructs the return series from the user-supplied closing
index levels.
The second command writes all regenerated outputs to `reproduced/`. It also
checks selected headline numerical results and stops if any of them differ.

The result reproduction performs the following steps:

1. merges the reconstructed returns with the main-period and holdout forecast
   panels;
2. constructs the within-market common samples and the model-specific samples;
3. recomputes joint VaR–ES scores, paired HAC comparisons, Holm-adjusted
   primary tests, calibration diagnostics, density summaries, and rankings;
4. recomputes the Model Confidence Sets using the recorded seeds and 5,000
   moving-block bootstrap replications;
5. rebuilds Tables 1–4 and Figures 1 and 2;
6. validates the released pre-main-OOS clone-count and prior-sensitivity
   summary for GJR-U and AARSV in SP500 and SHCOMP;
7. rebuilds the SHCOMP common-sample, DAX density and holdout parameter-mixture,
   and GJR-R versus GJR-U summaries;
8. aggregates the monthly asymmetry estimates and calculates the power and
   sample-size quantities discussed in Section 6.1.

## Generating the parameter-mixture forecasts

The default run starts from saved predictions and scores. The calculations
that generated the two parameter-mixture analyses are also supplied. After
constructing the returns, a market-month can be generated with:

```text
Rscript code/generate_sensitivity.R density DAX 2018-01
Rscript code/generate_sensitivity.R holdout DAX 2024-01
Rscript code/generate_sensitivity.R holdout NIKKEI 2024-01
```

Each command processes both ARSV and STAR-SV and writes to `recomputed/`.
Available market-months are listed in `inputs/sensitivity/meta.csv`.
These calculations use 768 nodes and 5,000 particles per node and are
computationally intensive. They do not run MCMC. The compact fitted inputs
are sufficient to reconstruct the nodes and forecasts conditional on the
reported estimates; they do not reproduce estimation of those estimates.
Details and comparison instructions are in `design/parameter_uncertainty.md`.

## Other estimation-based evidence

The eight Stan files specify the models. The default reproduction does not
repeat monthly estimation, initialization searches, or pre-2018 model fitting.
The STAR-SV window file records the selected windows and validation periods;
the candidate-level validation score panels are not supplied. The clone-count
and prior-sensitivity file records the reported changes, diagnostic values,
and assessment criteria. Its validation step checks those stored quantities;
it does not reconstruct them from alternative fits or forecast panels. The
K-specific joint draws and complete underlying comparison inputs are not
included.

## Random-number settings

The MCS settings file records the bootstrap size, realized block length, and
seed for each market, tail level, and confidence level. The Z2 reference file
records the simulation size, seed, and simulated cutoffs for each sample size
and tail level. The DAX date-level density scores retain the block structure
used to reproduce the reported descriptive uncertainty bands; the accompanying
parameter-draw file records the stochastic nodes and seeds for inspection.

## Paper correspondence

| Paper item | Machine-readable location |
|---|---|
| Table 1 | `results/exhibits/table1_*.csv` |
| Table 2 | `results/exhibits/table2.csv` |
| Table 3 | `results/exhibits/table3.csv` |
| Table 4 | `results/exhibits/table4_*.csv` |
| Figure 1 | `results/exhibits/table2.csv`, `results/exhibits/figure1.pdf`, `results/exhibits/figure1.png` |
| Figure 2 | `results/exhibits/figure2_data.csv`, `results/exhibits/figure2.pdf`, `results/exhibits/figure2.png` |
| Online Appendix A1 | `results/main_oos/market_characteristics.csv`, `results/holdout/market_characteristics.csv` |
| Online Appendix A2 | period-specific evaluation-sample and unavailable-forecast files |
| Online Appendix A3 | `code/models/`, `results/supplementary/star_window_selection.csv` |
| Online Appendix A4 | `results/supplementary/prior_sensitivity.csv` (GJR-U and AARSV fits in SP500 and SHCOMP) |
| Online Appendix A5 | period-specific comparison, ranking, and MCS files |
| Online Appendix A6 | period-specific calibration files and `design/z2_reference_values.csv` |
| Online Appendix A7 | period-specific predictive-density and ranking files |
| Online Appendix A8 | `results/supplementary/shcomp_common_sample/` |
| Online Appendix A9.1 | `results/supplementary/dax_parameter_uncertainty/` |
| Online Appendix A9.2 | `results/supplementary/holdout_parameter_uncertainty/` |
| Online Appendix A10 | `results/supplementary/gjr_comparison/fz0_table.csv` |
| Online Appendix A11 | `results/supplementary/asymmetry_parameter_summary.csv`, aggregated from `monthly_parameter_estimates.csv` |
| Section 6.1 power calculations | `results/exhibits/power_calculations.csv` |

## Key numerical checks

- DAX holdout STAR-SV minus ARSV mean FZ0 difference: `0.0126942`;
  Holm-adjusted p-value: `0.8223942`.
- NIKKEI holdout STAR-SV minus ARSV mean FZ0 difference: `−0.1145950`;
  Holm-adjusted p-value: `0.5749785`.
- Neither primary holdout comparison is corroborated.
- Figure 1 contains 60 main-OOS contrast cells: six contrasts for five markets
  at the 95% and 97.5% levels.
- The pre-main-OOS clone-count and prior-sensitivity summary contains 12 rows
  covering GJR-U and AARSV fits in SP500 and SHCOMP.
- DAX parameter-uncertainty relative changes: `−0.7137%` for LPDS and
  `−0.0599%` for CRPS.

- Holdout parameter-mixture mean FZ0 differences: `0.0127398421` for DAX
  and `−0.1142469154` for NIKKEI.
- NIKKEI's sensitivity increment has MC/HAC `0.1224359` and retains its
  limited-precision qualification.

## Notes on interpretation

Cross-model comparisons and rankings use each market’s eight-model common
evaluation sample. Calibration uses every valid date available for the model
being assessed. Unavailable forecasts are reported and are never imputed.
The 95% and 99% holdout results and all comparisons outside the two-test
primary family are supplementary or secondary evidence as described in
`design/holdout_protocol.md`.
