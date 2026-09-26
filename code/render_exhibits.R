build_table1_outputs <- function(main_support, holdout_support,
                                 main_unavailable, holdout_unavailable,
                                 output_dir) {
  design <- data.frame(
    item = c(
      "Markets and returns", "Evaluation periods", "Forecast timing",
      "Common estimation target", "Evaluation and samples"
    ),
    design = c(
      "DAX, FTSE, NIKKEI, SHCOMP, and SP500; 100 times daily log close-to-close returns on local trading calendars",
      "Main OOS: 2018-01-01 to 2023-12-31; holdout: 2024-01-01 to 2025-12-31",
      "Most recent 2,000 observations; monthly structural refits; daily state updates and one-day-ahead forecasts",
      "Stationary AR(1) conditional mean; unit-variance standardized Student-t innovations; data cloning K=5 with plug-in structural estimates",
      "Main OOS: 95% and 97.5% comparison levels; holdout: two primary 97.5% tests; common and model-specific evaluation samples; no imputation"
    ),
    stringsAsFactors = FALSE
  )
  models <- data.frame(
    model = c("GARCH", "GJR-R", "GJR-U", "APARCH", "ARSV", "AARSV", "TSV-RT", "STAR-SV"),
    family = c(rep("GARCH", 4L), rep("SV", 4L)),
    comparison_role = c(
      "Symmetric benchmark", "Restricted boundary benchmark", "Benchmark extension", "Benchmark extension",
      "Symmetric benchmark", "Benchmark extension", "Benchmark extension", "Benchmark extension"
    ),
    mechanism = c(
      "Symmetric squared-shock response", "Observed-sign squared-shock response",
      "Observed-sign squared-shock response", "Power and leverage response",
      "Symmetric latent log-volatility", "Continuous standardized-return leverage response",
      "Raw-negative-return threshold response",
      "Mild/severe raw-return regime with standardized-residual response"
    ),
    qualification = c(
      "Stationary baseline", "Unit-persistence boundary",
      "Free asymmetry with stationary effective persistence",
      "Powered-moment and effective-persistence restrictions", "SV benchmark",
      "Equity-leverage sign restriction", "Return-triggered threshold specification",
      "AR(1)-mean extension; trailing window 5 for DAX and 10 for the other markets"
    ),
    stringsAsFactors = FALSE
  )
  sample_rows <- lapply(market_order, function(market) {
    main <- main_support[main_support$market == market, , drop = FALSE]
    holdout <- holdout_support[holdout_support$market == market, , drop = FALSE]
    missing_text <- function(frame, label) {
      missing <- frame[frame$market == market, , drop = FALSE]
      if (!nrow(missing)) return(character())
      groups <- split(missing, factor(missing$model, levels = model_order), drop = TRUE)
      vapply(groups, function(group) {
        months <- paste(sort(unique(group$refit_month)), collapse = ", ")
        paste0(
          label, ": ", display_model(group$model[[1]]), " ", months,
          " (", nrow(group), " days)"
        )
      }, character(1))
    }
    missing <- c(
      missing_text(main_unavailable, "Main"),
      missing_text(holdout_unavailable, "Holdout")
    )
    data.frame(
      market = market,
      main_calendar_n = unique(main$calendar_n),
      main_common_n = unique(main$common_n),
      holdout_calendar_n = unique(holdout$calendar_n),
      holdout_common_n = unique(holdout$common_n),
      unavailable_forecasts = if (length(missing)) paste(missing, collapse = "; ") else "None",
      stringsAsFactors = FALSE
    )
  })
  samples <- do.call(rbind, sample_rows)
  write_csv(design, file.path(output_dir, "table1_design.csv"))
  write_csv(models, file.path(output_dir, "table1_models.csv"))
  write_csv(samples, file.path(output_dir, "table1_samples.csv"))
}

build_table2_output <- function(main_comparisons, output_dir) {
  contrast_order <- c(
    "GJR_R|GARCH", "GJR_U|GARCH", "APARCH|GARCH",
    "AARSV|ARSV", "TSV_RT|ARSV", "STAR_SV|ARSV"
  )
  out <- main_comparisons[main_comparisons$confidence %in% c(0.95, 0.975), , drop = FALSE]
  out$contrast <- paste(out$model, out$benchmark, sep = "|")
  out$contrast_order <- match(out$contrast, contrast_order)
  out$market_order <- match(out$market, market_order)
  out <- out[order(out$confidence, out$contrast_order, out$market_order), , drop = FALSE]
  out$contrast_id <- paste0(out$model, "_MINUS_", out$benchmark)
  out$display_text <- paste(display_model(out$model), "minus", display_model(out$benchmark))
  out$direction <- ifelse(
    out$mean_fz0_difference < 0, "first_model_favoured",
    ifelse(out$mean_fz0_difference > 0, "benchmark_favoured", "tie")
  )
  out$nominal_hac_evidence_5pct <- out$raw_two_sided_p_value < 0.05
  out <- out[, c(
    "confidence", "contrast_order", "contrast_id", "display_text", "market",
    "market_order", "n", "mean_fz0_difference", "hac_se", "ci_lower_95",
    "ci_upper_95", "raw_two_sided_p_value", "fraction_days_better", "direction",
    "nominal_hac_evidence_5pct"
  )]
  names(out)[names(out) == "raw_two_sided_p_value"] <- "hac_p_value"
  write_csv(out, file.path(output_dir, "table2.csv"))
}

build_main_oos_contrast_figure <- function(main_comparisons, output_dir) {
  contrast_order <- c(
    "GJR_R|GARCH", "GJR_U|GARCH", "APARCH|GARCH",
    "AARSV|ARSV", "TSV_RT|ARSV", "STAR_SV|ARSV"
  )
  contrast_labels <- c(
    "GJR-R - GARCH", "GJR-U - GARCH", "APARCH - GARCH",
    "AARSV - ARSV", "TSV-RT - ARSV", "STAR-SV - ARSV"
  )
  out <- main_comparisons[
    main_comparisons$confidence %in% c(0.95, 0.975),
    , drop = FALSE
  ]
  out$contrast <- paste(out$model, out$benchmark, sep = "|")
  out <- out[out$contrast %in% contrast_order, , drop = FALSE]
  out$contrast_order <- match(out$contrast, contrast_order)
  out$market_order <- match(out$market, market_order)
  out$contrast_id <- paste0(out$model, "_MINUS_", out$benchmark)
  out$display_text <- contrast_labels[out$contrast_order]
  out$tail_level <- ifelse(out$confidence == 0.95, "95%", "97.5%")
  out$interval_lower_95 <- out$mean_fz0_difference - 1.96 * out$hac_se
  out$interval_upper_95 <- out$mean_fz0_difference + 1.96 * out$hac_se
  out <- out[order(out$confidence, out$market_order, out$contrast_order), , drop = FALSE]

  key <- paste(out$market, out$contrast_id, out$confidence, sep = "|")
  assert_true(nrow(out) == 60L, "The main-OOS contrast figure requires 60 cells.")
  assert_true(!anyDuplicated(key), "The main-OOS contrast figure contains duplicate cells.")
  assert_true(!anyNA(out[c("market_order", "contrast_order")]),
              "The main-OOS contrast figure contains an unknown market or contrast.")
  assert_true(all(is.finite(out$mean_fz0_difference)) &&
                all(is.finite(out$hac_se)) && all(out$hac_se >= 0),
              "The main-OOS contrast figure contains invalid estimates.")

  figure_data <- out[, c(
    "market", "market_order", "confidence", "tail_level", "model", "benchmark",
    "contrast_order", "contrast_id", "display_text", "n",
    "mean_fz0_difference", "hac_se", "interval_lower_95", "interval_upper_95"
  )]
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is unavailable; Figure 1 image files were skipped.")
    return(invisible(figure_data))
  }

  plot_data <- figure_data
  plot_data$market <- factor(plot_data$market, levels = market_order)
  plot_data$display_text <- factor(plot_data$display_text, levels = contrast_labels)
  plot_data$tail_level <- factor(
    plot_data$tail_level, levels = c("95%", "97.5%")
  )
  palette <- c("#000000", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
  shapes <- c(16, 17, 15, 18, 3, 4)
  dodge <- ggplot2::position_dodge(width = 0.72)
  figure <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = market, y = mean_fz0_difference,
      colour = display_text, shape = display_text, group = display_text
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0, linewidth = 0.35, linetype = "dashed", colour = "grey35"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = interval_lower_95, ymax = interval_upper_95),
      width = 0, linewidth = 0.45, position = dodge
    ) +
    ggplot2::geom_point(size = 2.35, stroke = 0.9, position = dodge) +
    ggplot2::facet_wrap(~tail_level, ncol = 1L) +
    ggplot2::scale_colour_manual(values = stats::setNames(palette, contrast_labels)) +
    ggplot2::scale_shape_manual(values = stats::setNames(shapes, contrast_labels)) +
    ggplot2::scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
    ggplot2::labs(
      x = "Market", y = "Mean FZ0 difference",
      colour = "Contrast", shape = "Contrast"
    ) +
    ggplot2::theme_minimal(base_size = 10.5) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.25, colour = "grey88"),
      strip.text = ggplot2::element_text(face = "bold", size = 10.5),
      axis.title = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(size = 9),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(face = "bold", size = 9),
      legend.text = ggplot2::element_text(size = 8.2),
      legend.box.spacing = grid::unit(0.1, "cm"),
      plot.margin = ggplot2::margin(5.5, 7, 4, 7)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2L, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 2L, byrow = TRUE)
    )
  ggplot2::ggsave(
    file.path(output_dir, "figure1.pdf"), figure,
    width = 7.4, height = 6.0, units = "in", device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    file.path(output_dir, "figure1.png"), figure,
    width = 7.4, height = 6.0, units = "in", dpi = 300, bg = "white"
  )
  invisible(figure_data)
}

build_table3_output <- function(main_comparisons, holdout_comparisons, output_dir) {
  select_primary <- function(frame) {
    frame[
      frame$market %in% c("DAX", "NIKKEI") & frame$model == "STAR_SV" &
        frame$benchmark == "ARSV" & abs(frame$confidence - 0.975) <= 1e-12,
      , drop = FALSE
    ]
  }
  main <- select_primary(main_comparisons)
  holdout <- select_primary(holdout_comparisons)
  main <- main[match(c("DAX", "NIKKEI"), main$market), , drop = FALSE]
  holdout <- holdout[match(c("DAX", "NIKKEI"), holdout$market), , drop = FALSE]
  direction_status <- ifelse(
    sign(main$mean_fz0_difference) != sign(holdout$mean_fz0_difference),
    "sign_changed",
    ifelse(holdout$mean_fz0_difference < 0, "favourable_direction_retained",
           "unfavourable_direction_retained")
  )
  direction_display <- ifelse(
    direction_status == "sign_changed",
    "Sign changed; not significant",
    ifelse(
      direction_status == "favourable_direction_retained",
      "Negative direction retained", "Positive direction retained"
    )
  )
  out <- data.frame(
    market = main$market,
    contrast = "STAR_SV_MINUS_ARSV",
    confidence = 0.975,
    development_n = main$n,
    development_mean_fz0_difference = main$mean_fz0_difference,
    development_hac_se = main$hac_se,
    development_ci_lower_95 = main$ci_lower_95,
    development_ci_upper_95 = main$ci_upper_95,
    development_raw_two_sided_p_value = main$raw_two_sided_p_value,
    holdout_n = holdout$n,
    holdout_mean_fz0_difference = holdout$mean_fz0_difference,
    holdout_hac_se = holdout$hac_se,
    holdout_ci_lower_95 = holdout$ci_lower_95,
    holdout_ci_upper_95 = holdout$ci_upper_95,
    holdout_raw_two_sided_p_value = holdout$raw_two_sided_p_value,
    holdout_holm_adjusted_p_value = holdout$holm_adjusted_p_value,
    direction_status = direction_status,
    direction_display = direction_display,
    outcome = holdout$outcome,
    stringsAsFactors = FALSE
  )
  write_csv(out, file.path(output_dir, "table3.csv"))
}

unique_leader <- function(frame, value_column, direction) {
  value <- frame[[value_column]]
  best <- if (direction == "min") min(value) else max(value)
  rows <- frame[abs(value - best) <= 1e-14, , drop = FALSE]
  assert_true(nrow(rows) == 1L, "An assessment cell does not have a unique leader.")
  display_model(rows$model[[1]])
}

build_table4_outputs <- function(period_results, output_dir) {
  summary_rows <- list()
  gjr_rows <- list()
  for (period in names(period_results)) {
    result <- period_results[[period]]
    rankings <- result$rankings
    mcs <- result$mcs
    calibration <- result$calibration
    density <- result$density
    for (market in market_order) {
      ranking_cell <- rankings[
        rankings$market == market & abs(rankings$confidence - 0.975) <= 1e-12,
        , drop = FALSE
      ]
      mcs_cell <- mcs[
        mcs$market == market & abs(mcs$tail_confidence - 0.975) <= 1e-12 &
          abs(mcs$mcs_confidence - 0.95) <= 1e-12,
        , drop = FALSE
      ]
      calibration_cell <- calibration[
        calibration$market == market & abs(calibration$confidence - 0.975) <= 1e-12,
        , drop = FALSE
      ]
      density_cell <- density[density$market == market, , drop = FALSE]
      zones <- tolower(calibration_cell$z2_zone)
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        period = period,
        period_label = if (period == "main_oos") "Main-OOS, 2018--2023" else "Holdout, 2024--2025",
        market = market,
        common_support_n = unique(ranking_cell$n),
        fz0_leader_975 = unique_leader(ranking_cell, "mean_fz0", "min"),
        mcs_95_set_size_975 = sum(mcs_cell$included),
        christoffersen_cc_nonrejections_5pct = sum(calibration_cell$christoffersen_cc_p_value >= 0.05),
        de_conditional_nonrejections_5pct = sum(calibration_cell$de_conditional_p_value >= 0.05),
        z2_green = sum(zones == "green"),
        z2_yellow = sum(zones == "yellow"),
        z2_red = sum(zones == "red"),
        lpds_leader = unique_leader(density_cell, "mean_lpds", "max"),
        crps_leader = unique_leader(density_cell, "mean_crps", "min"),
        stringsAsFactors = FALSE
      )
    }

    fz0_counts <- c(GJR_R = 0L, GJR_U = 0L, ties = 0L)
    for (market in market_order) for (confidence in tail_levels) {
      cell <- rankings[
        rankings$market == market & rankings$model %in% c("GJR_R", "GJR_U") &
          abs(rankings$confidence - confidence) <= 1e-12,
        , drop = FALSE
      ]
      values <- stats::setNames(cell$mean_fz0, cell$model)
      if (abs(values[["GJR_R"]] - values[["GJR_U"]]) <= 1e-14) {
        fz0_counts[["ties"]] <- fz0_counts[["ties"]] + 1L
      } else if (values[["GJR_R"]] < values[["GJR_U"]]) {
        fz0_counts[["GJR_R"]] <- fz0_counts[["GJR_R"]] + 1L
      } else {
        fz0_counts[["GJR_U"]] <- fz0_counts[["GJR_U"]] + 1L
      }
    }
    density_count <- function(column, direction) {
      counts <- c(GJR_R = 0L, GJR_U = 0L, ties = 0L)
      for (market in market_order) {
        cell <- density[density$market == market & density$model %in% c("GJR_R", "GJR_U"), ]
        values <- stats::setNames(cell[[column]], cell$model)
        difference <- values[["GJR_R"]] - values[["GJR_U"]]
        if (abs(difference) <= 1e-14) counts[["ties"]] <- counts[["ties"]] + 1L
        else if ((direction == "max" && difference > 0) ||
                 (direction == "min" && difference < 0)) {
          counts[["GJR_R"]] <- counts[["GJR_R"]] + 1L
        } else counts[["GJR_U"]] <- counts[["GJR_U"]] + 1L
      }
      counts
    }
    lpds_counts <- density_count("mean_lpds", "max")
    crps_counts <- density_count("mean_crps", "min")
    gjr_rows[[length(gjr_rows) + 1L]] <- data.frame(
      period = period,
      period_label = if (period == "main_oos") "Main-OOS, 2018--2023" else "Holdout, 2024--2025",
      fz0_cells = 15L,
      fz0_gjr_r_better = fz0_counts[["GJR_R"]],
      fz0_gjr_u_better = fz0_counts[["GJR_U"]],
      fz0_ties = fz0_counts[["ties"]],
      density_markets = 5L,
      lpds_gjr_r_better = lpds_counts[["GJR_R"]],
      lpds_gjr_u_better = lpds_counts[["GJR_U"]],
      lpds_ties = lpds_counts[["ties"]],
      crps_gjr_r_better = crps_counts[["GJR_R"]],
      crps_gjr_u_better = crps_counts[["GJR_U"]],
      crps_ties = crps_counts[["ties"]],
      stringsAsFactors = FALSE
    )
  }
  write_csv(do.call(rbind, summary_rows), file.path(output_dir, "table4_assessments.csv"))
  write_csv(do.call(rbind, gjr_rows), file.path(output_dir, "table4_gjr_counts.csv"))
}

build_cumulative_primary_figure <- function(period_results, output_dir) {
  rows <- list()
  period_labels <- c(main_oos = "Main OOS 2018--2023", holdout = "Holdout 2024--2025")
  for (period in names(period_results)) {
    cumulative <- period_results[[period]]$cumulative
    comparisons <- period_results[[period]]$comparisons
    for (market in c("DAX", "NIKKEI")) {
      path <- cumulative[
        cumulative$market == market & cumulative$model == "STAR_SV" &
          cumulative$benchmark == "ARSV" & abs(cumulative$confidence - 0.975) <= 1e-12,
        , drop = FALSE
      ]
      path$date <- as.Date(path$date)
      path <- path[order(path$date), , drop = FALSE]
      comparison <- comparisons[
        comparisons$market == market & comparisons$model == "STAR_SV" &
          comparisons$benchmark == "ARSV" & abs(comparisons$confidence - 0.975) <= 1e-12,
        , drop = FALSE
      ]
      assert_true(nrow(path) == comparison$n,
                  "The cumulative FZ0 path and comparison sample sizes differ.")
      assert_true(
        abs(tail(path$cumulative_difference, 1L) / nrow(path) -
              comparison$mean_fz0_difference) <= 1e-12,
        "The cumulative FZ0 endpoint does not reproduce the mean comparison."
      )
      baseline <- data.frame(
        period = period,
        segment_label = period_labels[[period]],
        market = market,
        date = min(path$date) - 1,
        evaluation_index = 0L,
        is_baseline = TRUE,
        daily_difference = NA_real_,
        cumulative_difference = 0,
        stringsAsFactors = FALSE
      )
      data_rows <- data.frame(
        period = period,
        segment_label = period_labels[[period]],
        market = market,
        date = path$date,
        evaluation_index = seq_len(nrow(path)),
        is_baseline = FALSE,
        daily_difference = path$daily_difference,
        cumulative_difference = path$cumulative_difference,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1L]] <- rbind(baseline, data_rows)
    }
  }
  plot_data <- do.call(rbind, rows)
  write_csv(plot_data, file.path(output_dir, "figure2_data.csv"))

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning(paste(
      "Package 'ggplot2' is unavailable; cumulative FZ0 data were written but",
      "image files were skipped."
    ))
    return(invisible(plot_data))
  }
  plot_data$market <- factor(plot_data$market, levels = c("DAX", "NIKKEI"))
  plot_data$segment_label <- factor(
    plot_data$segment_label,
    levels = c("Main OOS 2018--2023", "Holdout 2024--2025")
  )
  endpoints <- do.call(rbind, lapply(split(
    plot_data[!plot_data$is_baseline, , drop = FALSE],
    interaction(plot_data$market[!plot_data$is_baseline],
                plot_data$segment_label[!plot_data$is_baseline], drop = TRUE)
  ), function(x) x[which.max(x$date), , drop = FALSE]))
  endpoints$end_label <- sprintf("%+.1f", endpoints$cumulative_difference)
  spans <- vapply(as.character(endpoints$market), function(market) {
    x <- plot_data[plot_data$market == market, ]
    diff(range(x$cumulative_difference))
  }, numeric(1))
  negative <- endpoints$cumulative_difference < 0
  endpoints$label_y <- endpoints$cumulative_difference +
    pmax(ifelse(negative, 8, 4), ifelse(negative, 0.09, 0.045) * spans)
  boundary <- do.call(rbind, lapply(c("DAX", "NIKKEI"), function(market) {
    x <- plot_data[plot_data$market == market, ]
    data.frame(
      market = factor(market, levels = c("DAX", "NIKKEI")),
      date = as.Date("2024-01-01"),
      y = max(x$cumulative_difference) + 0.04 * diff(range(x$cumulative_difference)),
      label = "Holdout reset"
    )
  }))
  palette <- c("Main OOS 2018--2023" = "#0072B2", "Holdout 2024--2025" = "#D55E00")
  figure <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = date, y = cumulative_difference, color = segment_label,
      group = interaction(market, segment_label)
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, color = "grey45") +
    ggplot2::geom_vline(
      xintercept = as.Date("2024-01-01"), linewidth = 0.35,
      linetype = "dotted", color = "grey50"
    ) +
    ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
    ggplot2::geom_text(
      data = endpoints, ggplot2::aes(y = label_y, label = end_label),
      hjust = 1.08, vjust = 0.5, size = 3.0, show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = boundary, ggplot2::aes(x = date, y = y, label = label),
      inherit.aes = FALSE, hjust = -0.08, vjust = 1, size = 2.9, color = "grey35"
    ) +
    ggplot2::facet_wrap(~market, ncol = 1, scales = "fixed") +
    ggplot2::scale_color_manual(values = palette, drop = FALSE) +
    ggplot2::scale_x_date(
      date_breaks = "1 year", date_labels = "%Y",
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::labs(
      x = NULL, y = "Cumulative FZ0 difference (STAR-SV minus ARSV)", color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position = "top", legend.justification = "left",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.25, color = "grey85"),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      strip.background = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      plot.margin = ggplot2::margin(6, 12, 6, 6)
    )
  ggplot2::ggsave(
    file.path(output_dir, "figure2.pdf"), figure,
    width = 7.2, height = 6.0, units = "in", device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    file.path(output_dir, "figure2.png"), figure,
    width = 7.2, height = 6.0, units = "in", dpi = 300, bg = "white"
  )
  invisible(plot_data)
}
