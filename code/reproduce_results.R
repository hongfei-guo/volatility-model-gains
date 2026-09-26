#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve script path.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
archive_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(archive_root, "code", "evaluation_functions.R"))
source(file.path(archive_root, "code", "render_exhibits.R"))
source(file.path(archive_root, "code", "sensitivity", "holdout_calculation.R"))
source(file.path(archive_root, "code", "appendix_outputs.R"))

font_cache <- file.path(tempdir(), "fontconfig-cache")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache)

output_root <- file.path(archive_root, "reproduced")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
for (directory in c("main_oos", "holdout", "exhibits", "supplementary")) {
  path <- file.path(output_root, directory)
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
}

z2_reference <- utils::read.csv(
  file.path(archive_root, "design", "z2_reference_values.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
mcs_settings <- utils::read.csv(
  file.path(archive_root, "design", "mcs_settings.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)

period_spec <- list(
  main_oos = list(
    forecast_file = "main_oos.csv",
    start = as.Date("2018-01-01"),
    end = as.Date("2023-12-31")
  ),
  holdout = list(
    forecast_file = "holdout.csv",
    start = as.Date("2024-01-01"),
    end = as.Date("2025-12-31")
  )
)

period_results <- list()
returns_dir <- file.path(output_root, "returns")
if (!dir.exists(returns_dir)) {
  stop(
    "Reconstructed returns are missing. Run Rscript code/build_returns.R first.",
    call. = FALSE
  )
}
for (period in names(period_spec)) {
  message("Reproducing ", period, " results...")
  spec <- period_spec[[period]]
  model_specific_panel <- read_forecast_panel(
    file.path(archive_root, "forecasts", spec$forecast_file)
  )
  model_specific_panel <- attach_realized_returns(model_specific_panel, returns_dir)
  model_specific_panel <- add_loss_columns(model_specific_panel)
  common_panel <- common_evaluation_panel(model_specific_panel)
  samples <- evaluation_samples(
    model_specific_panel, common_panel,
    returns_dir, spec$start, spec$end
  )
  summary <- evaluation_summary(common_panel, period, z2_reference)
  rankings <- ranking_table(summary)
  calibration <- calibration_table(
    evaluation_summary(model_specific_panel, period, z2_reference)
  )
  density <- unique(summary[c("market", "model", "n", "mean_lpds", "mean_crps")])
  comparisons <- within_family_comparisons(common_panel, period)
  mcs <- run_mcs_battery(common_panel, period, mcs_settings)
  cumulative <- cumulative_loss_differences(common_panel)
  characteristics <- market_characteristics(model_specific_panel)

  period_output <- file.path(output_root, period)
  write_csv(samples$support, file.path(period_output, "evaluation_samples.csv"))
  write_csv(samples$unavailable, file.path(period_output, "unavailable_forecasts.csv"))
  write_csv(samples$common_dates, file.path(period_output, "common_evaluation_dates.csv"))
  write_csv(summary, file.path(period_output, "evaluation_summary.csv"))
  write_csv(rankings, file.path(period_output, "ranking_summary.csv"))
  write_csv(calibration, file.path(period_output, "model_specific_calibration_summary.csv"))
  write_csv(density, file.path(period_output, "predictive_density_summary.csv"))
  write_csv(comparisons, file.path(period_output, "within_family_comparisons.csv"))
  write_csv(mcs, file.path(period_output, "mcs_results.csv"))
  write_csv(cumulative, file.path(period_output, "cumulative_fz0_differences.csv"))
  write_csv(characteristics, file.path(period_output, "market_characteristics.csv"))

  period_results[[period]] <- list(
    model_specific_panel = model_specific_panel,
    common_panel = common_panel,
    samples = samples,
    summary = summary,
    rankings = rankings,
    calibration = calibration,
    density = density,
    comparisons = comparisons,
    mcs = mcs,
    cumulative = cumulative,
    characteristics = characteristics
  )
}

message("Rebuilding main-paper exhibits...")
exhibit_output <- file.path(output_root, "exhibits")
build_table1_outputs(
  period_results$main_oos$samples$support,
  period_results$holdout$samples$support,
  period_results$main_oos$samples$unavailable,
  period_results$holdout$samples$unavailable,
  exhibit_output
)
build_table2_output(period_results$main_oos$comparisons, exhibit_output)
build_table3_output(
  period_results$main_oos$comparisons,
  period_results$holdout$comparisons,
  exhibit_output
)
build_table4_outputs(period_results, exhibit_output)
build_main_oos_contrast_figure(period_results$main_oos$comparisons, exhibit_output)
build_cumulative_primary_figure(period_results, exhibit_output)

message("Rebuilding supplementary evaluation outputs...")
supplementary_output <- file.path(output_root, "supplementary")

# Pre-main-OOS clone-count and prior sensitivity summary. Monthly estimation is
# outside the default reproduction; this step validates the released summary
# against the scope and scenarios reported in Online Appendix A4.
prior_sensitivity <- utils::read.csv(
  file.path(archive_root, "results", "supplementary", "prior_sensitivity.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
required_sensitivity_columns <- c(
  "market", "model", "training_start", "training_end", "validation_start",
  "validation_end", "K", "prior_name", "max_rhat", "min_ess_bulk",
  "min_ess_tail", "max_standardized_parameter_difference_from_K5_baseline",
  "max_relative_VaR_ES_difference_from_K5_baseline",
  "mean_lpds_difference_from_K5_baseline",
  "relative_crps_difference_from_K5_baseline", "parameter_tolerance_se",
  "risk_relative_tolerance", "lpds_absolute_tolerance",
  "crps_relative_tolerance", "assessment"
)
assert_true(
  all(required_sensitivity_columns %in% names(prior_sensitivity)),
  "The released sensitivity summary is missing required columns."
)
assert_true(
  nrow(prior_sensitivity) == 12L &&
    setequal(unique(prior_sensitivity$market), c("SP500", "SHCOMP")) &&
    setequal(unique(prior_sensitivity$model), c("GJR_U", "AARSV")),
  "The released sensitivity summary must contain 12 GJR_U/AARSV rows for SP500 and SHCOMP."
)
sensitivity_key <- paste(
  prior_sensitivity$market, prior_sensitivity$model,
  prior_sensitivity$K, prior_sensitivity$prior_name, sep = "|"
)
expected_sensitivity_key <- as.vector(outer(
  paste(rep(c("SP500", "SHCOMP"), each = 2L), c("GJR_U", "AARSV"), sep = "|"),
  c("2|baseline", "10|baseline", "5|diffuse"), paste, sep = "|"
))
assert_true(
  !anyDuplicated(sensitivity_key) && setequal(sensitivity_key, expected_sensitivity_key),
  "The released sensitivity summary does not contain the expected scenarios."
)
assert_true(
  all(prior_sensitivity$assessment[prior_sensitivity$K == 2L] ==
        "low_clone_diagnostic") &&
    all(prior_sensitivity$assessment[prior_sensitivity$K != 2L] ==
          "within_tolerance"),
  "The released sensitivity assessments do not match Online Appendix A4."
)
assert_true(
  all(is.finite(prior_sensitivity$max_rhat)) &&
    all(prior_sensitivity$max_rhat <= 1.01) &&
    all(prior_sensitivity$min_ess_bulk >= 400) &&
    all(prior_sensitivity$min_ess_tail >= 400),
  "A released sensitivity fit does not satisfy the reported diagnostics."
)
evaluated_sensitivity <- prior_sensitivity[prior_sensitivity$K != 2L, , drop = FALSE]
assert_true(
  nrow(evaluated_sensitivity) == 8L &&
    all(abs(evaluated_sensitivity[[
      "max_standardized_parameter_difference_from_K5_baseline"
    ]]) <= evaluated_sensitivity$parameter_tolerance_se) &&
    all(abs(evaluated_sensitivity[[
      "max_relative_VaR_ES_difference_from_K5_baseline"
    ]]) <= evaluated_sensitivity$risk_relative_tolerance) &&
    all(abs(evaluated_sensitivity[[
      "mean_lpds_difference_from_K5_baseline"
    ]]) <= evaluated_sensitivity$lpds_absolute_tolerance) &&
    all(abs(evaluated_sensitivity[[
      "relative_crps_difference_from_K5_baseline"
    ]]) <= evaluated_sensitivity$crps_relative_tolerance),
  "A released K=10 or diffuse-prior sensitivity result exceeds its reported tolerance."
)
write_csv(
  prior_sensitivity,
  file.path(supplementary_output, "prior_sensitivity.csv")
)

# SHCOMP seven-model common-sample sensitivity.
main_model_specific <- period_results$main_oos$model_specific_panel
main_common_dates <- unique(
  period_results$main_oos$common_panel$date[
    period_results$main_oos$common_panel$market == "SHCOMP"
  ]
)
seven_models <- setdiff(model_order, "STAR_SV")
shcomp_rows <- list()
shcomp_summary_rows <- list()
for (confidence in tail_levels) {
  suffix <- risk_suffix(confidence)
  cells <- lapply(seven_models, function(model) {
    full <- main_model_specific[
      main_model_specific$market == "SHCOMP" & main_model_specific$model == model,
      , drop = FALSE
    ]
    restricted <- full[full$date %in% main_common_dates, , drop = FALSE]
    excluded <- full[!full$date %in% main_common_dates, , drop = FALSE]
    data.frame(
      market = "SHCOMP",
      confidence = confidence,
      model = model,
      n_full = nrow(full),
      mean_fz0_full = mean(full[[paste0("fz0_", suffix)]]),
      n_restricted = nrow(restricted),
      mean_fz0_restricted = mean(restricted[[paste0("fz0_", suffix)]]),
      n_excluded = nrow(excluded),
      mean_fz0_excluded = mean(excluded[[paste0("fz0_", suffix)]]),
      stringsAsFactors = FALSE
    )
  })
  cells <- do.call(rbind, cells)
  cells$rank_full <- rank(cells$mean_fz0_full, ties.method = "min")
  cells$rank_restricted <- rank(cells$mean_fz0_restricted, ties.method = "min")
  cells$mean_shift_full_minus_restricted <-
    cells$mean_fz0_full - cells$mean_fz0_restricted
  cells$rank_change_full_minus_restricted <- cells$rank_full - cells$rank_restricted
  cells <- cells[, c(
    "market", "confidence", "model", "n_full", "mean_fz0_full", "rank_full",
    "n_restricted", "mean_fz0_restricted", "rank_restricted", "n_excluded",
    "mean_fz0_excluded", "mean_shift_full_minus_restricted",
    "rank_change_full_minus_restricted"
  )]
  shcomp_rows[[length(shcomp_rows) + 1L]] <- cells
  full_order <- paste(cells$model[order(cells$rank_full)], collapse = ">")
  restricted_order <- paste(cells$model[order(cells$rank_restricted)], collapse = ">")
  shcomp_summary_rows[[length(shcomp_summary_rows) + 1L]] <- data.frame(
    market = "SHCOMP",
    confidence = confidence,
    n_full = unique(cells$n_full),
    n_restricted = unique(cells$n_restricted),
    n_excluded = unique(cells$n_excluded),
    winner_full = cells$model[which.min(cells$rank_full)],
    winner_restricted = cells$model[which.min(cells$rank_restricted)],
    winner_unchanged = cells$model[which.min(cells$rank_full)] ==
      cells$model[which.min(cells$rank_restricted)],
    complete_rank_order_full = full_order,
    complete_rank_order_restricted = restricted_order,
    complete_rank_order_unchanged = identical(full_order, restricted_order),
    models_with_rank_change = sum(cells$rank_full != cells$rank_restricted),
    max_absolute_mean_shift = max(abs(cells$mean_shift_full_minus_restricted)),
    stringsAsFactors = FALSE
  )
}
shcomp_dir <- file.path(supplementary_output, "shcomp_common_sample")
write_csv(do.call(rbind, shcomp_rows), file.path(shcomp_dir, "fz0_rank_sensitivity.csv"))
write_csv(
  do.call(rbind, shcomp_summary_rows),
  file.path(shcomp_dir, "fz0_rank_stability_summary.csv")
)

# Main-OOS GJR-R versus GJR-U FZ0 and calibration comparison.
gjr_rows <- list()
gjr_calibration_rows <- list()
main_common <- period_results$main_oos$common_panel
main_summary <- period_results$main_oos$summary
for (market in market_order) {
  market_panel <- main_common[main_common$market == market, , drop = FALSE]
  first <- market_panel[market_panel$model == "GJR_R", , drop = FALSE]
  second <- market_panel[market_panel$model == "GJR_U", , drop = FALSE]
  merged <- merge(first, second, by = "date", suffixes = c("_r", "_u"))
  for (confidence in tail_levels) {
    suffix <- risk_suffix(confidence)
    difference <- merged[[paste0("fz0_", suffix, "_r")]] -
      merged[[paste0("fz0_", suffix, "_u")]]
    test <- hac_mean_test(difference)
    row_r <- main_summary[
      main_summary$market == market & main_summary$model == "GJR_R" &
        abs(main_summary$confidence - confidence) <= 1e-12,
      , drop = FALSE
    ]
    row_u <- main_summary[
      main_summary$market == market & main_summary$model == "GJR_U" &
        abs(main_summary$confidence - confidence) <= 1e-12,
      , drop = FALSE
    ]
    gjr_rows[[length(gjr_rows) + 1L]] <- data.frame(
      market = market,
      confidence = confidence,
      n = nrow(merged),
      hac_lag = as.integer(test[["lag"]]),
      mean_fz0_gjr_r = row_r$mean_fz0,
      mean_fz0_gjr_u = row_u$mean_fz0,
      mean_difference_gjr_r_minus_gjr_u = unname(test[["mean"]]),
      hac_se = unname(test[["se"]]),
      hac_statistic = unname(test[["statistic"]]),
      hac_p_value = unname(test[["p_value"]]),
      fraction_days_gjr_r_lower_fz0 = mean(difference < 0),
      lower_mean_model = if (row_r$mean_fz0 < row_u$mean_fz0) "GJR_R" else "GJR_U",
      two_sided_hac_evidence_at_5pct = unname(test[["p_value"]]) < 0.05,
      stringsAsFactors = FALSE
    )
    gjr_calibration_rows[[length(gjr_calibration_rows) + 1L]] <- data.frame(
      market = market,
      confidence = confidence,
      n = nrow(merged),
      gjr_r_z2 = row_r$z2,
      gjr_r_zone = row_r$z2_zone,
      gjr_u_z2 = row_u$z2,
      gjr_u_zone = row_u$z2_zone,
      z2_green_cutoff = row_r$z2_green_cutoff,
      z2_red_cutoff = row_r$z2_red_cutoff,
      zone_comparison = if (
        row_r$z2_zone == row_u$z2_zone
      ) "same_zone" else if (
        match(row_r$z2_zone, c("red", "yellow", "green")) >
          match(row_u$z2_zone, c("red", "yellow", "green"))
      ) "gjr_r_better_zone" else "gjr_u_better_zone",
      stringsAsFactors = FALSE
    )
  }
}
gjr_dir <- file.path(supplementary_output, "gjr_comparison")
write_csv(do.call(rbind, gjr_rows), file.path(gjr_dir, "fz0_comparison.csv"))
write_csv(
  do.call(rbind, gjr_calibration_rows),
  file.path(gjr_dir, "calibration_comparison.csv")
)

# DAX density sensitivity: rebuild the published summary from its released
# date-level scores. The parameter nodes are supplied for deeper inspection.
dax_summary <- rebuild_dax_sensitivity(file.path(
  archive_root, "results", "supplementary", "dax_parameter_uncertainty",
  "daily_density_scores.csv"
))
write_csv(
  dax_summary,
  file.path(supplementary_output, "dax_parameter_uncertainty", "sensitivity_summary.csv")
)


# Holdout parameter-mixture comparison and appendix summaries.
holdout_sensitivity <- rebuild_holdout_sensitivity(archive_root,
  period_results$holdout$common_panel,period_results$holdout$samples$common_dates)
holdout_directory <- file.path(supplementary_output,'holdout_parameter_uncertainty')
for(name in c('summary','nested')) {
  filename <- if(name=='summary') 'sensitivity_summary.csv' else 'nested_estimates.csv'
  key <- if(name=='summary') c('market','metric') else c('market','B')
  check_reference(holdout_sensitivity[[name]],file.path(archive_root,'results/supplementary/holdout_parameter_uncertainty',filename),key)
  write_csv(holdout_sensitivity[[name]],file.path(holdout_directory,filename))
}
write_csv(parameter_summary(file.path(archive_root,'results/supplementary/monthly_parameter_estimates.csv')),
  file.path(supplementary_output,'asymmetry_parameter_summary.csv'))
write_csv(power_summary(read.csv(file.path(exhibit_output,'table3.csv'))),file.path(exhibit_output,'power_calculations.csv'))
write_csv(gjr_table(period_results),file.path(supplementary_output,'gjr_comparison','fz0_table.csv'))

message("Checking the primary reported results...")
table3 <- utils::read.csv(
  file.path(exhibit_output, "table3.csv"), stringsAsFactors = FALSE
)
assert_true(
  abs(table3$holdout_mean_fz0_difference[table3$market == "DAX"] - 0.0126942137098656) <= 1e-12,
  "DAX primary holdout mean differs from the paper."
)
assert_true(
  abs(table3$holdout_holm_adjusted_p_value[table3$market == "DAX"] - 0.822394226045449) <= 1e-12,
  "DAX primary Holm p-value differs from the paper."
)
assert_true(
  abs(table3$holdout_mean_fz0_difference[table3$market == "NIKKEI"] + 0.114595015179796) <= 1e-12,
  "NIKKEI primary holdout mean differs from the paper."
)
assert_true(
  abs(table3$holdout_holm_adjusted_p_value[table3$market == "NIKKEI"] - 0.574978489530676) <= 1e-12,
  "NIKKEI primary Holm p-value differs from the paper."
)
assert_true(all(table3$outcome == "not_corroborated"),
            "A primary holdout outcome is no longer non-corroborating.")
assert_true(
  max(abs(dax_summary$relative_change_pct - c(-0.713681587602221, -0.0598577001375825))) <= 1e-10,
  "DAX parameter-uncertainty relative changes differ from the paper."
)

message("Reproduction completed successfully. Outputs: ", output_root)
