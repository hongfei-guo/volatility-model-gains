hpu_spec <- function () 
list(markets = c("DAX", "NIKKEI"), models = c("ARSV", "STAR_SV"), months = format(seq(as.Date("2024-01-01"), 
    as.Date("2025-12-01"), by = "month"), "%Y-%m"), excluded_nikkei_months = c("2024-04", "2024-10", 
    "2025-08"), n_dates = c(DAX = 507L, NIKKEI = 425L), star_m = c(DAX = 5L, NIKKEI = 10L), confidence = 0.975, 
    K = 5L, training_n = 2000L, B = 768L, P = 5000L, blocks = 12L, block_size = 64L, nested_B = c(128L, 
        256L, 384L, 512L, 768L), base_seed = 20260904L, covariance_tolerance = 1e-10, resample_ess_fraction = 0.5, 
    dead_fraction_limit = 0.05)

hpu_assert <- function (ok, message) 
{
    if (length(ok) != 1L || is.na(ok) || !ok) 
        stop(message, call. = FALSE)
}

hpu_tasks <- function () 
{
    s <- hpu_spec()
    x <- expand.grid(month = s$months, market = s$markets, stringsAsFactors = FALSE)
    x$task_index <- seq_len(nrow(x))
    x <- x[!(x$market == "NIKKEI" & x$month %in% s$excluded_nikkei_months), ]
    x$task_id <- paste(x$market, x$month, sep = "_")
    rownames(x) <- NULL
    x[, c("task_id", "task_index", "market", "month")]
}

hpu_outer_shocks <- function () 
{
    s <- hpu_spec()
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    out <- matrix(NA_real_, s$B, 8L)
    for (pair in seq_len(s$B/2L)) {
        set.seed(s$base_seed + pair)
        z <- rnorm(8L)
        out[2L * pair - 1L, ] <- z
        out[2L * pair, ] <- -z
    }
    out
}

hpu_inner_seed <- function (task_index, node, step) 
{
    s <- hpu_spec()
    hpu_assert(task_index %in% hpu_tasks()$task_index && node %in% seq_len(s$B) && step %in% 0:31, "Invalid seed coordinate.")
    as.integer(s$base_seed + 10000L + ((task_index - 1L) * s$B + node - 1L) * 32L + step)
}

hpu_scope_blocks <- function () 
{
    s <- hpu_spec()
    out <- list(pooled = seq_len(s$blocks))
    for (g in seq_len(s$blocks)) out[[paste0("loo_", g)]] <- setdiff(seq_len(s$blocks), g)
    for (B in s$nested_B[s$nested_B < s$B]) out[[paste0("nested_", B)]] <- seq_len(B/s$block_size)
    out
}

hpu_score_mixture <- function (means, sds, nus, realized, scopes = hpu_scope_blocks(), block_size = hpu_spec()$block_size) 
{
    hpu_assert(is.matrix(sds) && nrow(sds) > 0L && ncol(sds) == length(means) && length(nus) == length(means) && 
        all(is.finite(means)) && all(is.finite(sds)) && all(sds > 0) && all(is.finite(nus)) && all(nus > 
        2) && length(realized) == 1L && is.finite(realized), "Invalid predictive mixture.")
    hpu_assert(length(scopes) > 0L && !is.null(names(scopes)) && !anyDuplicated(names(scopes)), "Invalid scope labels.")
    rows <- lapply(names(scopes), function(scope) {
        blocks <- scopes[[scope]]
        hpu_assert(length(blocks) > 0L && !anyDuplicated(blocks) && all(blocks %in% seq_len(ncol(sds)/block_size)), 
            "Invalid block selection.")
        cols <- unlist(lapply(blocks, function(g) ((g - 1L) * block_size + 1L):(g * block_size)))
        risk <- mixture_loss_var_es(hpu_spec()$confidence, rep(means[cols], each = nrow(sds)), as.vector(sds[, 
            cols, drop = FALSE]), rep(nus[cols], each = nrow(sds)))
        fz <- fz0_loss(-realized, risk[["VaR"]], risk[["ES"]], hpu_spec()$confidence)
        hpu_assert(all(is.finite(c(risk, fz))) && risk[["ES"]] > 0, "Nonfinite risk/score or nonpositive ES.")
        data.frame(scope, B = length(cols), VaR = unname(risk[["VaR"]]), ES = unname(risk[["ES"]]), FZ0 = unname(fz))
    })
    do.call(rbind, rows)
}

hpu_validate_series <- function (data) 
{
    hpu_assert(is.data.frame(data) && identical(names(data), c("date", "return")), "Return input must contain exactly date, return.")
    hpu_assert(inherits(data$date, "Date") && !anyNA(data$date) && !anyDuplicated(data$date) && all(diff(as.numeric(data$date)) > 
        0) && is.numeric(data$return) && all(is.finite(data$return)), "Invalid return calendar.")
    invisible(TRUE)
}

hpu_training <- function (data, fit, month) 
{
    hpu_validate_series(data)
    first <- which(format(data$date, "%Y-%m") == month)[1L]
    hpu_assert(!is.na(first) && first > hpu_spec()$training_n, "Missing training window.")
    train <- data[seq.int(first - hpu_spec()$training_n, first - 1L), , drop = FALSE]
    hpu_assert(as.character(min(train$date)) == as.character(fit$training_start) && as.character(max(train$date)) == 
        as.character(fit$training_end), "Training dates differ from accepted fit.")
    train
}

hpu_check_task_result <- function (x, task, dates) 
{
    s <- hpu_spec()
    keys <- c("date", "model", "variant", "scope")
    expected <- expand.grid(date = as.character(dates), model = s$models, variant = c("P", "I"), scope = names(hpu_scope_blocks()), 
        stringsAsFactors = FALSE)
    key <- function(z) do.call(paste, c(z[keys], sep = "|"))
    hpu_assert(is.data.frame(x) && all(c(keys, "task_id", "market", "VaR", "ES", "FZ0", "B") %in% names(x)), 
        "Incomplete task result schema.")
    hpu_assert(!anyDuplicated(key(x)) && setequal(key(x), key(expected)) && all(x$task_id == task$task_id) && 
        all(x$market == task$market), "Task output membership mismatch.")
    hpu_assert(all(is.finite(as.matrix(x[c("VaR", "ES", "FZ0", "B")]))) && all(x$ES > 0) && all(x$ES >= 
        x$VaR - 1e-10), "Invalid task results.")
    expected_B <- vapply(hpu_scope_blocks(), length, integer(1)) * s$block_size
    hpu_assert(all(x$B == expected_B[x$scope]), "Task scope has wrong outer-node count.")
    invisible(TRUE)
}

hpu_run_task <- function (task_id, data, fits, support_dates) 
{
    s <- hpu_spec()
    tasks <- hpu_tasks()
    task <- tasks[tasks$task_id == task_id, ]
    hpu_assert(nrow(task) == 1L && identical(names(fits), s$models), "Invalid task or fit pair.")
    hpu_validate_series(data)
    dates <- data$date[format(data$date, "%Y-%m") == task$month]
    hpu_assert(length(dates) > 0L && length(dates) <= 31L && identical(as.character(dates), as.character(support_dates[format(support_dates, 
        "%Y-%m") == task$month])), "Common support does not equal this retained month's trading calendar.")
    shocks <- hpu_outer_shocks()
    out <- list()
    for (model in s$models) {
        fit <- fits[[model]]
        params <- hpu_fit_nodes(fit, task$market, model, task$month, shocks)
        training <- hpu_training(data, fit, task$month)
        states <- list(P = vector("list", s$B), I = vector("list", s$B))
        for (b in seq_len(s$B)) for (variant in c("P", "I")) {
            theta <- if (variant == "P") 
                params$point
            else params$nodes[[b]]
            states[[variant]][[b]] <- initialize_sv_filter(model, theta, training$return, s$P, hpu_inner_seed(task$task_index, 
                b, 0L), s$resample_ess_fraction, star_m = s$star_m[[task$market]], training_dates = training$date, 
                dead_fraction_limit = s$dead_fraction_limit)
        }
        for (j in seq_along(dates)) {
            y <- data$return[match(dates[[j]], data$date)]
            for (variant in c("P", "I")) {
                means <- nus <- numeric(s$B)
                sds <- matrix(NA_real_, s$P, s$B)
                for (b in seq_len(s$B)) {
                  theta <- if (variant == "P") 
                    params$point
                  else params$nodes[[b]]
                  step <- sv_predict_components_update(states[[variant]][[b]], theta, y, seed = hpu_inner_seed(task$task_index, 
                    b, j), resample_ess_fraction = s$resample_ess_fraction, observation_id = as.character(dates[[j]]), 
                    dead_fraction_limit = s$dead_fraction_limit)
                  states[[variant]][[b]] <- step$state
                  means[[b]] <- step$components$mean[[1L]]
                  nus[[b]] <- step$components$nu[[1L]]
                  hpu_assert(length(step$components$sd) == s$P, "Incorrect particle supply.")
                  sds[, b] <- step$components$sd
                }
                scores <- hpu_score_mixture(means, sds, nus, y)
                scores$date <- as.character(dates[[j]])
                scores$model <- model
                scores$variant <- variant
                scores$task_id <- task_id
                scores$market <- task$market
                out[[length(out) + 1L]] <- scores
            }
        }
    }
    result <- do.call(rbind, out)
    hpu_check_task_result(result, task, dates)
    result
}

hpu_summary <- function (daily, original, support) 
{
    s <- hpu_spec()
    tasks <- hpu_tasks()
    hpu_assert(identical(names(support), c("market", "date")) && !anyDuplicated(paste(support$market, 
        support$date)) && setequal(unique(support$market), s$markets), "Invalid support table.")
    for (market in s$markets) {
        dates <- as.Date(support$date[support$market == market])
        hpu_assert(length(dates) == s$n_dates[[market]] && !anyNA(dates) && all(diff(as.numeric(dates)) > 
            0), "Wrong support calendar.")
    }
    hpu_assert(setequal(unique(daily$task_id), tasks$task_id), "Missing or extraneous tasks.")
    for (i in seq_len(nrow(tasks))) {
        t <- tasks[i, ]
        ds <- as.Date(support$date[support$market == t$market])
        ds <- ds[format(ds, "%Y-%m") == t$month]
        hpu_check_task_result(daily[daily$task_id == t$task_id, ], t, ds)
    }
    ckey <- paste(original$market, original$date, original$model)
    expected <- unlist(lapply(s$models, function(m) paste(support$market, support$date, m)))
    hpu_assert(!anyDuplicated(ckey) && setequal(ckey, expected) && all(is.finite(original$FZ0)), "Original anchor membership or values invalid.")
    rows <- list()
    nested <- list()
    for (market in s$markets) {
        dates <- support$date[support$market == market]
        contrast <- function(frame) {
            a <- frame[frame$market == market & frame$model == "ARSV", ]
            b <- frame[frame$market == market & frame$model == "STAR_SV", ]
            hpu_assert(!anyDuplicated(a$date) && !anyDuplicated(b$date) && setequal(a$date, dates) && 
                setequal(b$date, dates), "Contrast date mismatch.")
            b$FZ0[match(dates, b$date)] - a$FZ0[match(dates, a$date)]
        }
        anchor <- contrast(original)
        scopes <- names(hpu_scope_blocks())
        series <- lapply(scopes, function(sc) {
            P <- contrast(daily[daily$scope == sc & daily$variant == "P", ])
            I <- contrast(daily[daily$scope == sc & daily$variant == "I", ])
            cbind(P = P, I = I, delta = I - P, P_minus_C = P - anchor)
        })
        names(series) <- scopes
        for (B in s$nested_B) {
            scope <- if (B == s$B) 
                "pooled"
            else paste0("nested_", B)
            values <- colMeans(series[[scope]])
            nested[[length(nested) + 1L]] <- data.frame(market, B, mean_d_P = unname(values[["P"]]), 
                mean_d_I = unname(values[["I"]]), delta = unname(values[["delta"]]), mean_P_minus_C = unname(values[["P_minus_C"]]))
        }
        loo <- t(vapply(series[paste0("loo_", seq_len(s$blocks))], colMeans, numeric(4)))
        mc <- sqrt((s$blocks - 1)/s$blocks * colSums(sweep(loo, 2L, colMeans(loo), "-")^2))
        for (metric in colnames(series$pooled)) {
            h <- hac_mean_test(series$pooled[, metric])
            ratio <- if (h[["se"]] > 0) 
                mc[[metric]]/h[["se"]]
            else if (mc[[metric]] == 0) 
                0
            else Inf
            rows[[length(rows) + 1L]] <- data.frame(market, metric, n = length(dates), estimate = unname(h[["mean"]]), 
                hac_se = unname(h[["se"]]), hac_lower = unname(h[["mean"]] - 1.96 * h[["se"]]), hac_upper = unname(h[["mean"]] + 1.96 * h[["se"]]), 
                mc_se = unname(mc[[metric]]), mc_to_hac = unname(ratio), block_stability = if (ratio <= 
                  0.1) 
                  "small_relative_to_HAC"
                else "limited", original_difference = mean(anchor))
        }
    }
    list(summary = do.call(rbind, rows), nested = do.call(rbind, nested))
}

hpu_original_scores <- function (risks, data_by_market, support) 
{
    s <- hpu_spec()
    hpu_assert(all(c("market", "model", "date", "VaR_0975", "ES_0975") %in% names(risks)), "Missing original risk columns.")
    key <- paste(risks$market, risks$date, risks$model)
    expected <- unlist(lapply(s$models, function(m) paste(support$market, support$date, m)))
    hpu_assert(!anyDuplicated(key) && setequal(key, expected), "Original risks do not match common support.")
    out <- risks[, c("market", "model", "date")]
    out$FZ0 <- NA_real_
    for (market in s$markets) {
        data <- data_by_market[[market]]
        hpu_validate_series(data)
        ii <- which(risks$market == market)
        pos <- match(as.Date(risks$date[ii]), data$date)
        hpu_assert(!anyNA(pos) && all(is.finite(risks$VaR_0975[ii])) && all(is.finite(risks$ES_0975[ii])) && 
            all(risks$ES_0975[ii] > 0) && all(risks$ES_0975[ii] >= risks$VaR_0975[ii] - 1e-10), "Invalid original risk inputs.")
        out$FZ0[ii] <- fz0_loss(-data$return[pos], risks$VaR_0975[ii], risks$ES_0975[ii], s$confidence)
    }
    hpu_assert(all(is.finite(out$FZ0)), "Original scores are nonfinite.")
    out
}

