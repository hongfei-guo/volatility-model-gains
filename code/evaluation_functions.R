assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

assert_columns <- function(frame, required, label) {
  missing <- setdiff(required, names(frame))
  assert_true(
    length(missing) == 0L,
    paste0(label, " is missing columns: ", paste(missing, collapse = ", "))
  )
}

write_csv <- function(frame, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  options(digits = 17, scipen = 999)
  utils::write.csv(frame, path, row.names = FALSE, quote = TRUE, na = "")
}

model_order <- c(
  "GARCH", "GJR_R", "GJR_U", "APARCH",
  "ARSV", "AARSV", "TSV_RT", "STAR_SV"
)
market_order <- c("DAX", "FTSE", "NIKKEI", "SHCOMP", "SP500")
tail_levels <- c(0.95, 0.975, 0.99)

risk_suffix <- function(confidence) sprintf("%04d", round(1000 * confidence))

display_model <- function(model) {
  out <- as.character(model)
  out[out == "GJR_R"] <- "GJR-R"
  out[out == "GJR_U"] <- "GJR-U"
  out[out == "STAR_SV"] <- "STAR-SV"
  out[out == "TSV_RT"] <- "TSV-RT"
  out
}

read_forecast_panel <- function(path) {
  panel <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "date", "market", "model", "family", "refit_month",
    "training_start", "training_end", "training_n",
    "VaR_0950", "ES_0950", "VaR_0975", "ES_0975",
    "VaR_0990", "ES_0990", "lpds", "crps", "pit"
  )
  assert_columns(panel, required, basename(path))
  panel$date <- as.Date(panel$date)
  assert_true(!anyNA(panel$date), paste0(basename(path), " has invalid dates."))
  assert_true(!anyDuplicated(panel[c("market", "model", "date")]),
              paste0(basename(path), " has duplicate scientific keys."))
  assert_true(setequal(unique(panel$market), market_order),
              paste0(basename(path), " has an unexpected market set."))
  assert_true(setequal(unique(panel$model), model_order),
              paste0(basename(path), " has an unexpected model set."))
  panel
}

attach_realized_returns <- function(panel, returns_dir) {
  return_rows <- lapply(market_order, function(market) {
    path <- file.path(returns_dir, paste0(market, ".csv"))
    if (!file.exists(path)) {
      stop("Missing reconstructed return series: ", path, call. = FALSE)
    }
    frame <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    assert_columns(frame, c("date", "return"), basename(path))
    frame$date <- as.Date(frame$date)
    frame$return <- suppressWarnings(as.numeric(frame$return))
    if (anyNA(frame$date) || any(!is.finite(frame$return)) || anyDuplicated(frame$date)) {
      stop(market, " return series must have valid, unique dates and finite values.")
    }
    data.frame(
      market = market,
      date = frame$date,
      realized_return = frame$return,
      stringsAsFactors = FALSE
    )
  })
  returns <- do.call(rbind, return_rows)
  panel_key <- paste(panel$market, panel$date, sep = "\r")
  return_key <- paste(returns$market, returns$date, sep = "\r")
  index <- match(panel_key, return_key)
  if (anyNA(index)) {
    stop("The supplied returns do not cover every released forecast date.", call. = FALSE)
  }
  panel$realized_return <- returns$realized_return[index]
  panel$realized_loss <- -panel$realized_return
  panel
}

fz0_loss <- function(loss, var, es, confidence) {
  alpha <- 1 - confidence
  finite <- is.finite(loss) & is.finite(var) & is.finite(es)
  if (any(!finite)) {
    stop("FZ0 inputs must be finite on every released forecast date.", call. = FALSE)
  }
  if (any(finite & (es <= 0 | es < var))) {
    stop("FZ0 requires positive ES forecasts with ES not below VaR.", call. = FALSE)
  }
  out <- rep(NA_real_, length(loss))
  ok <- finite
  out[ok] <- as.numeric(loss[ok] > var[ok]) * (loss[ok] - var[ok]) /
    (alpha * es[ok]) + var[ok] / es[ok] + log(es[ok]) - 1
  out
}

add_loss_columns <- function(panel) {
  for (confidence in tail_levels) {
    suffix <- risk_suffix(confidence)
    var <- panel[[paste0("VaR_", suffix)]]
    es <- panel[[paste0("ES_", suffix)]]
    panel[[paste0("fz0_", suffix)]] <- fz0_loss(
      panel$realized_loss, var, es, confidence
    )
    panel[[paste0("hit_", suffix)]] <- as.integer(panel$realized_loss > var)
  }
  panel
}

common_evaluation_panel <- function(panel) {
  rows <- lapply(market_order, function(market) {
    market_panel <- panel[panel$market == market, , drop = FALSE]
    dates <- lapply(model_order, function(model) {
      unique(market_panel$date[market_panel$model == model])
    })
    common_dates <- Reduce(intersect, dates)
    out <- market_panel[market_panel$date %in% common_dates, , drop = FALSE]
    out <- out[order(match(out$model, model_order), out$date), , drop = FALSE]
    assert_true(
      nrow(out) == 8L * length(common_dates),
      paste0("The common evaluation panel is incomplete for ", market, ".")
    )
    out
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

evaluation_samples <- function(panel, common_panel, returns_dir, period_start, period_end) {
  support_rows <- list()
  unavailable_rows <- list()
  common_rows <- list()
  for (market in market_order) {
    returns <- utils::read.csv(
      file.path(returns_dir, paste0(market, ".csv")),
      stringsAsFactors = FALSE
    )
    dates <- as.Date(returns$date)
    calendar <- sort(dates[dates >= period_start & dates <= period_end])
    market_panel <- panel[panel$market == market, , drop = FALSE]
    common_dates <- sort(unique(common_panel$date[common_panel$market == market]))
    common_rows[[length(common_rows) + 1L]] <- data.frame(
      market = market, date = as.character(common_dates), stringsAsFactors = FALSE
    )
    for (model in model_order) {
      valid <- sort(unique(market_panel$date[market_panel$model == model]))
      unavailable <- setdiff(calendar, valid)
      support_rows[[length(support_rows) + 1L]] <- data.frame(
        market = market,
        model = model,
        calendar_n = length(calendar),
        valid_n = length(valid),
        common_n = length(common_dates),
        unavailable_n = length(unavailable),
        first_valid_date = as.character(min(valid)),
        last_valid_date = as.character(max(valid)),
        stringsAsFactors = FALSE
      )
      if (length(unavailable)) {
        unavailable_rows[[length(unavailable_rows) + 1L]] <- data.frame(
          market = market,
          model = model,
          date = as.character(unavailable),
          refit_month = format(unavailable, "%Y-%m"),
          reason = "forecast_unavailable",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  list(
    support = do.call(rbind, support_rows),
    unavailable = if (length(unavailable_rows)) do.call(rbind, unavailable_rows) else
      data.frame(market = character(), model = character(), date = character(),
                 refit_month = character(), reason = character()),
    common_dates = do.call(rbind, common_rows)
  )
}

kupiec_test <- function(loss, var, confidence) {
  ok <- is.finite(loss) & is.finite(var)
  hit <- as.integer(loss[ok] > var[ok])
  n <- length(hit)
  assert_true(n > 0L, "Kupiec test has no finite observations.")
  x <- sum(hit)
  alpha <- 1 - confidence
  phat <- x / n
  ll_null <- stats::dbinom(x, size = n, prob = alpha, log = TRUE)
  ll_alt <- stats::dbinom(x, size = n, prob = phat, log = TRUE)
  lr <- max(0, -2 * (ll_null - ll_alt))
  c(hits = x, n = n, hit_rate = phat, statistic = lr,
    p_value = stats::pchisq(lr, 1, lower.tail = FALSE))
}

christoffersen_tests <- function(loss, var, confidence) {
  ok <- is.finite(loss) & is.finite(var)
  hit <- as.integer(loss[ok] > var[ok])
  assert_true(length(hit) >= 2L, "Christoffersen test has too few observations.")
  prev <- head(hit, -1L)
  next_ <- tail(hit, -1L)
  n00 <- sum(prev == 0 & next_ == 0)
  n01 <- sum(prev == 0 & next_ == 1)
  n10 <- sum(prev == 1 & next_ == 0)
  n11 <- sum(prev == 1 & next_ == 1)
  xlogp <- function(count, probability) {
    if (count == 0L) return(0)
    if (probability <= 0) return(-Inf)
    count * log(probability)
  }
  p01 <- if (n00 + n01) n01 / (n00 + n01) else 0
  p11 <- if (n10 + n11) n11 / (n10 + n11) else 0
  p1 <- (n01 + n11) / max(n00 + n01 + n10 + n11, 1)
  ll_iid <- xlogp(n00 + n10, 1 - p1) + xlogp(n01 + n11, p1)
  ll_markov <- xlogp(n00, 1 - p01) + xlogp(n01, p01) +
    xlogp(n10, 1 - p11) + xlogp(n11, p11)
  lr_ind <- max(0, -2 * (ll_iid - ll_markov))
  uc <- kupiec_test(loss[ok], var[ok], confidence)
  lr_cc <- unname(uc[["statistic"]]) + lr_ind
  c(
    lr_ind = lr_ind,
    p_ind = stats::pchisq(lr_ind, 1, lower.tail = FALSE),
    lr_cc = lr_cc,
    p_cc = stats::pchisq(lr_cc, 2, lower.tail = FALSE)
  )
}

acerbi_z2 <- function(loss, var, es, confidence) {
  alpha <- 1 - confidence
  ok <- is.finite(loss) & is.finite(var) & is.finite(es) & es > 0
  assert_true(any(ok), "Z2 statistic has no finite observations.")
  1 - mean(loss[ok] * as.numeric(loss[ok] > var[ok]) / (alpha * es[ok]))
}

simulate_z2_cutoffs <- function(confidence, n, simulations, seed,
                                probabilities = c(0.05, 0.0001),
                                chunk_size = 10000L) {
  alpha <- 1 - confidence
  q <- stats::qnorm(confidence)
  es <- stats::dnorm(q) / alpha
  z2 <- numeric(simulations)
  set.seed(seed)
  first <- 1L
  while (first <= simulations) {
    last <- min(simulations, first + chunk_size - 1L)
    size <- last - first + 1L
    exceedances <- stats::rbinom(size, n, alpha)
    total <- sum(exceedances)
    x <- if (total) stats::qnorm(stats::runif(total, confidence, 1)) else numeric()
    cumulative <- c(0, cumsum(x / (alpha * es)))
    end_index <- cumsum(exceedances)
    start_index <- end_index - exceedances
    sums <- cumulative[end_index + 1L] - cumulative[start_index + 1L]
    z2[first:last] <- 1 - sums / n
    first <- last + 1L
  }
  stats::quantile(z2, probabilities, names = FALSE, type = 8)
}

du_escanciano_tests <- function(pit, confidence, lags = 4L) {
  alpha <- 1 - confidence
  u <- pmin(pmax(as.numeric(pit), 0), 1)
  u <- u[is.finite(u)]
  n <- length(u)
  assert_true(n > lags + 1L, "Du-Escanciano test has too few observations.")
  H <- (alpha - u) / alpha * as.numeric(u <= alpha)
  expected <- alpha / 2
  variance <- alpha * (1 / 3 - alpha / 4)
  z <- sqrt(n) * (mean(H) - expected) / sqrt(variance)
  centered <- H - expected
  gamma0 <- mean(centered^2)
  rho <- vapply(seq_len(lags), function(k) {
    if (gamma0 <= 0) return(0)
    mean(centered[(k + 1L):n] * centered[1L:(n - k)]) / gamma0
  }, numeric(1))
  conditional <- n * sum(rho^2)
  c(
    unconditional_statistic = z,
    unconditional_p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
    conditional_statistic = conditional,
    conditional_p_value = stats::pchisq(conditional, lags, lower.tail = FALSE)
  )
}

lookup_z2_reference <- function(reference, period, confidence, n) {
  row <- reference[
    reference$period == period & abs(reference$confidence - confidence) <= 1e-12 &
      reference$n == n,
    , drop = FALSE
  ]
  assert_true(nrow(row) == 1L, "Z2 reference-value lookup did not return one row.")
  row
}

evaluation_summary <- function(panel, period, z2_reference) {
  groups <- split(panel, interaction(panel$market, panel$model, drop = TRUE))
  rows <- list()
  for (group in groups) {
    for (confidence in tail_levels) {
      suffix <- risk_suffix(confidence)
      var <- group[[paste0("VaR_", suffix)]]
      es <- group[[paste0("ES_", suffix)]]
      uc <- kupiec_test(group$realized_loss, var, confidence)
      cc <- christoffersen_tests(group$realized_loss, var, confidence)
      z2 <- acerbi_z2(group$realized_loss, var, es, confidence)
      reference <- lookup_z2_reference(z2_reference, period, confidence, nrow(group))
      zone <- if (z2 >= reference$green_cutoff) {
        "green"
      } else if (z2 >= reference$red_cutoff) {
        "yellow"
      } else {
        "red"
      }
      de <- du_escanciano_tests(group$pit, confidence)
      rows[[length(rows) + 1L]] <- data.frame(
        market = as.character(group$market[[1]]),
        model = as.character(group$model[[1]]),
        confidence = confidence,
        n = nrow(group),
        average_var = mean(var),
        average_es = mean(es),
        hits = unname(uc[["hits"]]),
        hit_rate = unname(uc[["hit_rate"]]),
        kupiec_statistic = unname(uc[["statistic"]]),
        kupiec_p_value = unname(uc[["p_value"]]),
        christoffersen_ind_statistic = unname(cc[["lr_ind"]]),
        christoffersen_ind_p_value = unname(cc[["p_ind"]]),
        christoffersen_cc_statistic = unname(cc[["lr_cc"]]),
        christoffersen_cc_p_value = unname(cc[["p_cc"]]),
        z2 = z2,
        z2_green_cutoff = reference$green_cutoff,
        z2_red_cutoff = reference$red_cutoff,
        z2_zone = zone,
        de_unconditional_statistic = unname(de[["unconditional_statistic"]]),
        de_unconditional_p_value = unname(de[["unconditional_p_value"]]),
        de_conditional_statistic = unname(de[["conditional_statistic"]]),
        de_conditional_p_value = unname(de[["conditional_p_value"]]),
        mean_fz0 = mean(group[[paste0("fz0_", suffix)]]),
        mean_lpds = mean(group$lpds),
        mean_crps = mean(group$crps),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out <- out[order(match(out$market, market_order), match(out$model, model_order), out$confidence), ]
  row.names(out) <- NULL
  out
}

calibration_table <- function(summary) {
  columns <- c(
    "market", "model", "confidence", "n", "average_var", "average_es",
    "hits", "hit_rate", "kupiec_statistic", "kupiec_p_value",
    "christoffersen_ind_statistic", "christoffersen_ind_p_value",
    "christoffersen_cc_statistic", "christoffersen_cc_p_value",
    "z2", "z2_green_cutoff", "z2_red_cutoff", "z2_zone",
    "de_unconditional_statistic", "de_unconditional_p_value",
    "de_conditional_statistic", "de_conditional_p_value"
  )
  summary[columns]
}

hac_mean_test <- function(x, lag = NULL) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  assert_true(n >= 3L, "HAC mean test has too few observations.")
  centered <- x - mean(x)
  if (is.null(lag)) lag <- max(1L, floor(4 * (n / 100)^(2 / 9)))
  lag <- min(as.integer(lag), n - 1L)
  lrv <- mean(centered^2)
  for (k in seq_len(lag)) {
    gamma <- mean(centered[(k + 1L):n] * centered[1L:(n - k)])
    lrv <- lrv + 2 * (1 - k / (lag + 1)) * gamma
  }
  se <- sqrt(max(lrv, 0) / n)
  statistic <- if (se > 0) mean(x) / se else NA_real_
  c(
    n = n, lag = lag, mean = mean(x), se = se, statistic = statistic,
    p_value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  )
}

within_family_comparisons <- function(panel, period) {
  pairs <- list(
    c("GJR_R", "GARCH"), c("GJR_U", "GARCH"), c("APARCH", "GARCH"),
    c("AARSV", "ARSV"), c("TSV_RT", "ARSV"), c("STAR_SV", "ARSV")
  )
  if (period == "holdout") pairs <- c(pairs, list(c("STAR_SV", "AARSV")))
  rows <- list()
  for (market in market_order) {
    market_panel <- panel[panel$market == market, , drop = FALSE]
    for (pair in pairs) {
      first <- market_panel[market_panel$model == pair[[1]], , drop = FALSE]
      second <- market_panel[market_panel$model == pair[[2]], , drop = FALSE]
      merged <- merge(first, second, by = "date", suffixes = c("_model", "_benchmark"))
      for (confidence in tail_levels) {
        suffix <- risk_suffix(confidence)
        difference <- merged[[paste0("fz0_", suffix, "_model")]] -
          merged[[paste0("fz0_", suffix, "_benchmark")]]
        test <- hac_mean_test(difference)
        rows[[length(rows) + 1L]] <- data.frame(
          market = market,
          model = pair[[1]],
          benchmark = pair[[2]],
          contrast_id = paste0(pair[[1]], "_MINUS_", pair[[2]]),
          confidence = confidence,
          n = length(difference),
          mean_fz0_difference = unname(test[["mean"]]),
          hac_lag = as.integer(test[["lag"]]),
          hac_se = unname(test[["se"]]),
          ci_lower_95 = unname(test[["mean"]] - 1.96 * test[["se"]]),
          ci_upper_95 = unname(test[["mean"]] + 1.96 * test[["se"]]),
          hac_statistic = unname(test[["statistic"]]),
          raw_two_sided_p_value = unname(test[["p_value"]]),
          fraction_days_better = mean(difference < 0),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$primary_test <- FALSE
  out$holm_adjusted_p_value <- NA_real_
  out$outcome <- ""
  if (period == "holdout") {
    primary <- out$market %in% c("DAX", "NIKKEI") & out$model == "STAR_SV" &
      out$benchmark == "ARSV" & abs(out$confidence - 0.975) <= 1e-12
    assert_true(sum(primary) == 2L, "The primary holdout family is incomplete.")
    out$primary_test[primary] <- TRUE
    out$holm_adjusted_p_value[primary] <- stats::p.adjust(
      out$raw_two_sided_p_value[primary], method = "holm"
    )
    out$outcome[primary] <- ifelse(
      out$mean_fz0_difference[primary] < 0 & out$holm_adjusted_p_value[primary] < 0.05,
      "corroborated",
      ifelse(
        out$mean_fz0_difference[primary] > 0 & out$holm_adjusted_p_value[primary] < 0.05,
        "statistical_reversal", "not_corroborated"
      )
    )
  }
  out[c(
    "market", "model", "benchmark", "contrast_id", "confidence", "n",
    "mean_fz0_difference", "hac_lag", "hac_se", "ci_lower_95",
    "ci_upper_95", "hac_statistic", "raw_two_sided_p_value",
    "holm_adjusted_p_value", "fraction_days_better", "primary_test", "outcome"
  )]
}

ranking_table <- function(summary) {
  summary$rank_fz0 <- ave(
    summary$mean_fz0,
    interaction(summary$market, summary$confidence, drop = TRUE),
    FUN = function(x) rank(x, ties.method = "min")
  )
  density <- unique(summary[c("market", "model", "mean_crps", "mean_lpds")])
  density$rank_crps <- ave(
    density$mean_crps, density$market,
    FUN = function(x) rank(x, ties.method = "min")
  )
  density$rank_lpds <- ave(
    -density$mean_lpds, density$market,
    FUN = function(x) rank(x, ties.method = "min")
  )
  out <- merge(
    summary, density[c("market", "model", "rank_crps", "rank_lpds")],
    by = c("market", "model"), all.x = TRUE, sort = FALSE
  )
  out[order(match(out$market, market_order), match(out$model, model_order), out$confidence), ]
}

run_mcs_battery <- function(panel, period, settings) {
  assert_true(requireNamespace("MCS", quietly = TRUE),
              "Package 'MCS' is required for Model Confidence Sets.")
  rows <- list()
  for (market in market_order) {
    for (confidence in tail_levels) {
      suffix <- risk_suffix(confidence)
      market_data <- panel[panel$market == market, c("date", "model", paste0("fz0_", suffix))]
      names(market_data)[[3]] <- "loss"
      wide <- reshape(market_data, idvar = "date", timevar = "model", direction = "wide")
      wide <- wide[order(wide$date), , drop = FALSE]
      wide <- wide[stats::complete.cases(wide), , drop = FALSE]
      loss_matrix <- as.matrix(wide[paste0("loss.", model_order)])
      colnames(loss_matrix) <- model_order
      for (mcs_confidence in c(0.95, 0.80)) {
        setting <- settings[
          settings$period == period & settings$market == market &
            abs(settings$tail_confidence - confidence) <= 1e-12 &
            abs(settings$mcs_confidence - mcs_confidence) <= 1e-12,
          , drop = FALSE
        ]
        assert_true(nrow(setting) == 1L, "MCS setting lookup did not return one row.")
        set.seed(setting$seed)
        result <- withCallingHandlers(
          MCS::MCSprocedure(
            Loss = loss_matrix,
            alpha = 1 - mcs_confidence,
            B = setting$bootstrap_replications,
            statistic = "Tmax",
            k = NULL,
            min.k = 5L,
            verbose = FALSE
          ),
          warning = function(condition) {
            if (grepl(
              "no non-missing arguments to max; returning -Inf",
              conditionMessage(condition), fixed = TRUE
            )) invokeRestart("muffleWarning")
          }
        )
        assert_true(as.integer(result@Info$k) == as.integer(setting$block_length),
                    "The MCS automatic block length differs from the recorded value.")
        shown <- as.data.frame(result@show, stringsAsFactors = FALSE)
        shown$model <- rownames(shown)
        rownames(shown) <- NULL
        table <- data.frame(
          model = model_order,
          included = model_order %in% shown$model,
          stringsAsFactors = FALSE
        )
        if (nrow(shown)) table <- merge(table, shown, by = "model", all.x = TRUE, sort = FALSE)
        table$market <- market
        table$tail_confidence <- confidence
        table$mcs_confidence <- mcs_confidence
        table$bootstrap_scheme <- setting$bootstrap_scheme
        table$bootstrap_replications <- setting$bootstrap_replications
        table$block_length <- as.integer(result@Info$k)
        table$block_selection_rule <- setting$block_selection_rule
        table$mcs_p_value <- result@Info$mcs_pvalue
        table$seed <- setting$seed
        rows[[length(rows) + 1L]] <- table
      }
    }
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

cumulative_loss_differences <- function(panel) {
  benchmarks <- c(
    GJR_R = "GARCH", GJR_U = "GARCH", APARCH = "GARCH",
    AARSV = "ARSV", TSV_RT = "ARSV", STAR_SV = "ARSV"
  )
  rows <- list()
  for (market in market_order) {
    market_panel <- panel[panel$market == market, , drop = FALSE]
    for (model in names(benchmarks)) {
      benchmark <- benchmarks[[model]]
      merged <- merge(
        market_panel[market_panel$model == model, ],
        market_panel[market_panel$model == benchmark, ],
        by = "date", suffixes = c("_model", "_benchmark")
      )
      merged <- merged[order(merged$date), , drop = FALSE]
      for (confidence in tail_levels) {
        suffix <- risk_suffix(confidence)
        difference <- merged[[paste0("fz0_", suffix, "_model")]] -
          merged[[paste0("fz0_", suffix, "_benchmark")]]
        rows[[length(rows) + 1L]] <- data.frame(
          date = as.character(merged$date),
          market = market,
          model = model,
          benchmark = benchmark,
          confidence = confidence,
          daily_difference = difference,
          cumulative_difference = cumsum(difference),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

market_characteristics <- function(panel) {
  rows <- lapply(market_order, function(market) {
    anchor <- panel[panel$market == market & panel$model == "GARCH", , drop = FALSE]
    anchor <- anchor[order(anchor$date), , drop = FALSE]
    x <- anchor$realized_return
    centered <- x - mean(x)
    sd_x <- stats::sd(x)
    skew <- mean(centered^3) / sd_x^3
    kurt <- mean(centered^4) / sd_x^4
    acf_sq <- as.numeric(stats::acf(x^2, plot = FALSE, lag.max = 20)$acf[-1])
    jb <- length(x) / 6 * (skew^2 + (kurt - 3)^2 / 4)
    data.frame(
      market = market,
      first_date = as.character(min(anchor$date)),
      last_date = as.character(max(anchor$date)),
      n = length(x),
      mean = mean(x),
      sd = sd_x,
      skewness = skew,
      kurtosis = kurt,
      excess_kurtosis = kurt - 3,
      minimum = min(x),
      maximum = max(x),
      negative_fraction = mean(x < 0),
      absolute_return_acf_1 = as.numeric(stats::acf(abs(x), plot = FALSE, lag.max = 1)$acf[[2]]),
      squared_return_acf_1 = acf_sq[[1]],
      squared_return_acf_5 = acf_sq[[5]],
      squared_return_acf_20 = acf_sq[[20]],
      jarque_bera = jb,
      jarque_bera_p_value = stats::pchisq(jb, 2, lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

dax_joint_hac <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  lag <- min(n - 1L, max(1L, floor(4 * (n / 100)^(2 / 9))))
  center <- colMeans(x)
  centered <- sweep(x, 2L, center, "-")
  lrv <- crossprod(centered) / n
  for (k in seq_len(lag)) {
    gamma <- crossprod(
      centered[(k + 1L):n, , drop = FALSE],
      centered[1L:(n - k), , drop = FALSE]
    ) / (n - k)
    lrv <- lrv + (1 - k / (lag + 1)) * (gamma + t(gamma))
  }
  list(n = n, lag = lag, mean = center, covariance = (lrv + t(lrv)) / (2 * n))
}

dax_contrast_series <- function(daily, score, scope, variant) {
  subset <- daily[daily$scope == scope & daily$variant == variant, , drop = FALSE]
  star <- subset[subset$model == "STAR_SV", c("date", score)]
  arsv <- subset[subset$model == "ARSV", c("date", score)]
  names(star)[[2]] <- "star"
  names(arsv)[[2]] <- "arsv"
  joined <- merge(star, arsv, by = "date", sort = TRUE)
  assert_true(nrow(joined) == 1523L, "DAX sensitivity support changed.")
  if (score == "lpds") joined$star - joined$arsv else joined$arsv - joined$star
}

jackknife_se <- function(values) {
  blocks <- length(values)
  sqrt((blocks - 1) / blocks * sum((values - mean(values))^2))
}

rebuild_dax_sensitivity <- function(path) {
  daily <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  rows <- lapply(c("lpds", "crps"), function(score) {
    plugin <- dax_contrast_series(daily, score, "all_parameter_blocks", "matched_plugin")
    integrated <- dax_contrast_series(daily, score, "all_parameter_blocks", "parameter_integrated")
    delta <- integrated - plugin
    joint <- dax_joint_hac(cbind(delta = delta, advantage_plugin = plugin))
    advantage_plugin <- mean(plugin)
    relative_change <- mean(delta) / advantage_plugin
    gradient <- c(1 / advantage_plugin, -mean(delta) / advantage_plugin^2)
    relative_hac_se <- sqrt(sum(outer(gradient, gradient) * joint$covariance))
    relative_loo <- vapply(seq_len(12L), function(block) {
      scope <- paste0("leave_block_out_", block)
      p <- dax_contrast_series(daily, score, scope, "matched_plugin")
      i <- dax_contrast_series(daily, score, scope, "parameter_integrated")
      mean(i - p) / mean(p)
    }, numeric(1))
    relative_mc_se <- jackknife_se(relative_loo)
    relative_combined_se <- sqrt(relative_hac_se^2 + relative_mc_se^2)
    mean_score <- function(model, variant, scope) {
      x <- daily[[score]][daily$model == model & daily$variant == variant & daily$scope == scope]
      assert_true(length(x) == 1523L, "DAX sensitivity mean has incomplete support.")
      mean(x)
    }
    arsv_main <- mean_score("ARSV", "reported_main", "reported_main")
    arsv_plugin <- mean_score("ARSV", "matched_plugin", "all_parameter_blocks")
    arsv_integrated <- mean_score("ARSV", "parameter_integrated", "all_parameter_blocks")
    star_main <- mean_score("STAR_SV", "reported_main", "reported_main")
    star_plugin <- mean_score("STAR_SV", "matched_plugin", "all_parameter_blocks")
    star_integrated <- mean_score("STAR_SV", "parameter_integrated", "all_parameter_blocks")
    advantage <- function(star, arsv) if (score == "lpds") star - arsv else arsv - star
    data.frame(
      score = toupper(score),
      n = 1523L,
      hac_lag = joint$lag,
      arsv_main = arsv_main,
      arsv_plugin = arsv_plugin,
      arsv_integrated = arsv_integrated,
      star_sv_main = star_main,
      star_sv_plugin = star_plugin,
      star_sv_integrated = star_integrated,
      advantage_main = advantage(star_main, arsv_main),
      advantage_plugin = advantage(star_plugin, arsv_plugin),
      advantage_integrated = advantage(star_integrated, arsv_integrated),
      relative_change_pct = 100 * relative_change,
      relative_hac_se_pct = 100 * relative_hac_se,
      relative_mc_se_pct = 100 * relative_mc_se,
      relative_combined_se_pct = 100 * relative_combined_se,
      combined_band_lower_pct = 100 * (relative_change - 1.96 * relative_combined_se),
      combined_band_upper_pct = 100 * (relative_change + 1.96 * relative_combined_se),
      direction_stable = sign(advantage_plugin) == sign(advantage(star_integrated, arsv_integrated)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
