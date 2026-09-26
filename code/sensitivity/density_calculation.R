dax_sens_assert <- function (condition, message) 
{
    if (!isTRUE(condition)) 
        stop(message, call. = FALSE)
    invisible(TRUE)
}

dax_sens_model_index <- function (model) 
{
    match(toupper(model), dax_sens_spec()$models) - 1L
}

dax_sens_month_index <- function (refit_month) 
{
    match(as.character(refit_month), dax_sens_spec()$months) - 1L
}

dax_sens_outer_seed <- function (block_index, pair_index) 
{
    spec <- dax_sens_spec()
    g <- as.integer(block_index) - 1L
    q <- as.integer(pair_index) - 1L
    dax_sens_assert(g %in% 0:(spec$blocks - 1L) && q %in% 0:(spec$pairs_per_block - 1L), "Invalid outer seed index.")
    as.integer(spec$base_seed + 1L + spec$pairs_per_block * g + q)
}

dax_sens_inner_seed <- function (model, refit_month, block_index, pair_index, sign, operation_index) 
{
    spec <- dax_sens_spec()
    m <- dax_sens_model_index(model)
    r <- dax_sens_month_index(refit_month)
    g <- as.integer(block_index) - 1L
    q <- as.integer(pair_index) - 1L
    s <- match(as.character(sign), c("positive", "negative")) - 1L
    o <- as.integer(operation_index)
    dax_sens_assert(!anyNA(c(m, r, g, q, s, o)), "Invalid inner seed key.")
    dax_sens_assert(m %in% 0:1 && r %in% 0:71 && g %in% 0:(spec$blocks - 1L) && q %in% 0:(spec$pairs_per_block - 
        1L) && s %in% 0:1 && o %in% 0:31, "Inner seed key is outside the specified radix.")
    inner_id <- 1L + (((((m * 72L + r) * spec$blocks + g) * spec$pairs_per_block + q) * 2L + s) * 32L + 
        o)
    as.integer(spec$base_seed + 100000L + inner_id)
}

dax_sens_crps_seed <- function (date_index, block_index) 
{
    spec <- dax_sens_spec()
    d <- as.integer(date_index)
    g <- as.integer(block_index)
    dax_sens_assert(length(d) == 1L && d >= 1L && d <= spec$n_dates, "Invalid CRPS date index.")
    dax_sens_assert(length(g) == 1L && g >= 1L && g <= spec$blocks, "Invalid CRPS block index.")
    as.integer(spec$base_seed + 1000L + (d - 1L) * spec$blocks + g)
}

dax_sens_outer_shocks <- function () 
{
    spec <- dax_sens_spec()
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    rows <- vector("list", spec$B)
    cursor <- 0L
    for (g in 1:spec$blocks) for (q in 1:spec$pairs_per_block) {
        seed <- dax_sens_outer_seed(g, q)
        set.seed(seed)
        positive <- stats::rnorm(8L)
        for (sign in c("positive", "negative")) {
            cursor <- cursor + 1L
            shock <- if (sign == "positive") 
                positive
            else -positive
            row <- data.frame(outer_index = cursor, block_index = g, pair_index = q, sign = sign, outer_seed = seed, 
                stringsAsFactors = FALSE)
            for (j in seq_len(8L)) row[[paste0("shock_", j)]] <- shock[[j]]
            rows[[cursor]] <- row
        }
    }
    out <- bind_rows_base(rows)
    dax_sens_assert(nrow(out) == spec$B, "Outer shock table has the wrong size.")
    out
}

dax_sens_theta_order <- function (model) 
{
    common <- c("mu_y", "phi_y", "mu_h", "phi_h", "sigma_eta", "nu")
    if (toupper(model) == "STAR_SV") 
        c(common, "gamma_severe", "gamma_mild")
    else common
}

dax_sens_working_names <- function (model) 
{
    common <- c("mu_y", "z_phi_y", "mu_h", "z_phi_h", "log_sigma_eta", "log_nu_minus_2")
    if (toupper(model) == "STAR_SV") 
        c(common, "logit_severity", "logit_mild_fraction")
    else common
}

dax_sens_working_map <- function (model, theta) 
{
    model <- toupper(model)
    theta <- as.numeric(theta[dax_sens_theta_order(model)])
    names(theta) <- dax_sens_theta_order(model)
    dax_sens_assert(all(is.finite(theta)), "Natural-scale parameter vector is nonfinite.")
    z <- c(mu_y = theta[["mu_y"]], z_phi_y = atanh(theta[["phi_y"]]/0.999), mu_h = theta[["mu_h"]], z_phi_h = atanh(theta[["phi_h"]]/0.999), 
        log_sigma_eta = log(theta[["sigma_eta"]]), log_nu_minus_2 = log(theta[["nu"]] - 2))
    if (model == "STAR_SV") {
        a <- -theta[["gamma_severe"]]
        f <- (theta[["gamma_mild"]] - theta[["gamma_severe"]])/a
        z <- c(z, logit_severity = stats::qlogis(a), logit_mild_fraction = stats::qlogis(f))
    }
    dax_sens_assert(all(is.finite(z)), "Working-scale parameter vector is nonfinite.")
    z
}

dax_sens_working_inverse <- function (model, z) 
{
    model <- toupper(model)
    z <- as.numeric(z[dax_sens_working_names(model)])
    names(z) <- dax_sens_working_names(model)
    theta <- c(mu_y = z[["mu_y"]], phi_y = 0.999 * tanh(z[["z_phi_y"]]), mu_h = z[["mu_h"]], phi_h = 0.999 * 
        tanh(z[["z_phi_h"]]), sigma_eta = exp(z[["log_sigma_eta"]]), nu = 2 + exp(z[["log_nu_minus_2"]]))
    if (model == "STAR_SV") {
        a <- stats::plogis(z[["logit_severity"]])
        f <- stats::plogis(z[["logit_mild_fraction"]])
        theta <- c(theta, gamma_severe = -a, gamma_mild = -a * (1 - f))
    }
    theta <- theta[dax_sens_theta_order(model)]
    dax_sens_assert(all(is.finite(theta)), "Inverse-mapped parameter vector is nonfinite.")
    validate_sv_parameters(model, as.list(theta))
    if (model == "STAR_SV") {
        dax_sens_assert(theta[["gamma_severe"]] > -1 && theta[["gamma_severe"]] < theta[["gamma_mild"]] && 
            theta[["gamma_mild"]] < 0, "STAR-SV inverse map reached a parameter boundary.")
    }
    theta
}

dax_sens_working_jacobian <- function (model, theta) 
{
    model <- toupper(model)
    theta_order <- dax_sens_theta_order(model)
    theta <- as.numeric(theta[theta_order])
    names(theta) <- theta_order
    working <- dax_sens_working_names(model)
    J <- matrix(0, nrow = length(working), ncol = length(theta_order), dimnames = list(working, theta_order))
    J["mu_y", "mu_y"] <- 1
    J["z_phi_y", "phi_y"] <- 0.999/(0.999^2 - theta[["phi_y"]]^2)
    J["mu_h", "mu_h"] <- 1
    J["z_phi_h", "phi_h"] <- 0.999/(0.999^2 - theta[["phi_h"]]^2)
    J["log_sigma_eta", "sigma_eta"] <- 1/theta[["sigma_eta"]]
    J["log_nu_minus_2", "nu"] <- 1/(theta[["nu"]] - 2)
    if (model == "STAR_SV") {
        gs <- theta[["gamma_severe"]]
        gm <- theta[["gamma_mild"]]
        a <- -gs
        f <- (gm - gs)/a
        J["logit_severity", "gamma_severe"] <- -1/(a * (1 - a))
        J["logit_mild_fraction", "gamma_severe"] <- (gm/gs^2)/(f * (1 - f))
        J["logit_mild_fraction", "gamma_mild"] <- (-1/gs)/(f * (1 - f))
    }
    dax_sens_assert(all(is.finite(J)), "Working-scale Jacobian is nonfinite.")
    J
}

dax_sens_parameter_nodes <- function (fit, model, refit_month, outer_shocks) 
{
    checked <- dax_sens_validate_fit(fit, model, refit_month)
    dimension <- length(checked$z_hat)
    rows <- vector("list", nrow(outer_shocks))
    parameters <- vector("list", nrow(outer_shocks))
    for (b in seq_len(nrow(outer_shocks))) {
        epsilon <- as.numeric(unlist(outer_shocks[b, paste0("shock_", seq_len(dimension)), drop = FALSE], 
            use.names = FALSE))
        z <- checked$z_hat + as.numeric(checked$L %*% epsilon)
        names(z) <- names(checked$z_hat)
        theta <- dax_sens_working_inverse(model, z)
        parameters[[b]] <- as.list(theta)
        row <- cbind(data.frame(model = model, refit_month = refit_month, stringsAsFactors = FALSE), 
            outer_shocks[b, , drop = FALSE])
        for (name in c("mu_y", "phi_y", "mu_h", "phi_h", "sigma_eta", "nu", "gamma_severe", "gamma_mild")) {
            row[[name]] <- if (name %in% names(theta)) 
                theta[[name]]
            else NA_real_
        }
        rows[[b]] <- row
    }
    list(frame = bind_rows_base(rows), parameters = parameters, checked = checked)
}

dax_sens_sensitivity_scopes <- function () 
{
    spec <- dax_sens_spec()
    c("pooled", paste0("loo_", seq_len(spec$blocks)), paste0("nested_", spec$nested_B[spec$nested_B < 
        spec$B]))
}

dax_sens_scope_supply <- function (scope) 
{
    spec <- dax_sens_spec()
    block_size <- 2L * spec$pairs_per_block
    if (identical(scope, "pooled")) {
        blocks <- spec$blocks
    }
    else if (startsWith(scope, "loo_")) {
        omitted <- suppressWarnings(as.integer(sub("loo_", "", scope, fixed = TRUE)))
        dax_sens_assert(length(omitted) == 1L && omitted %in% seq_len(spec$blocks), "Invalid LOO scope.")
        blocks <- spec$blocks - 1L
    }
    else if (startsWith(scope, "nested_")) {
        nested_B <- suppressWarnings(as.integer(sub("nested_", "", scope, fixed = TRUE)))
        dax_sens_assert(length(nested_B) == 1L && nested_B %in% spec$nested_B && nested_B < spec$B, "Invalid nested scope.")
        blocks <- nested_B%/%block_size
    }
    else {
        stop("Unknown sensitivity scope.", call. = FALSE)
    }
    c(n_components = as.integer(blocks * block_size * spec$P), n_crps_draws = as.integer(blocks * spec$crps_draws_per_block))
}

dax_sens_sv_predict_components_update <- function (state, params, realized_return, seed, resample_ess_fraction = 0.5, observation_id = NA_character_, 
    dead_fraction_limit = 0.05) 
{
    core <- sv_predict_components_update(state = state, params = params, realized_return = realized_return, 
        seed = seed, resample_ess_fraction = resample_ess_fraction, observation_id = observation_id, 
        dead_fraction_limit = dead_fraction_limit)
    core$components$mean <- core$components$mean[[1L]]
    core$components$nu <- core$components$nu[[1L]]
    core
}

dax_sens_crps_plans <- function (date_index, n_draws = dax_sens_spec()$crps_draws_per_block) 
{
    spec <- dax_sens_spec()
    lapply(seq_len(spec$blocks), function(g) {
        set.seed(dax_sens_crps_seed(date_index, g))
        list(component_uniforms = stratified_uniforms(n_draws), t_uniforms = independently_permuted_stratified_uniforms(n_draws))
    })
}

dax_sens_stratified_component_indices <- function (n_components, uniforms) 
{
    uniforms <- as.numeric(uniforms)
    dax_sens_assert(length(uniforms) >= 1L && all(is.finite(uniforms)) && all(uniforms > 0) && all(uniforms <= 
        1), "CRPS component uniforms are invalid.")
    pmin(floor(uniforms * as.integer(n_components)) + 1L, as.integer(n_components))
}

dax_sens_score_components <- function (means, sds, nus, realized, crps_plans) 
{
    spec <- dax_sens_spec()
    block_size <- 2L * spec$pairs_per_block
    dax_sens_assert(is.matrix(sds) && nrow(sds) == spec$P && ncol(sds) == spec$B, "Sensitivity component matrix has the wrong dimensions.")
    dax_sens_assert(length(means) == spec$B && length(nus) == spec$B, "Sensitivity component metadata have the wrong dimensions.")
    dax_sens_assert(all(is.finite(means)) && all(is.finite(sds)) && all(sds > 0) && all(is.finite(nus)) && 
        all(nus > 2), "Sensitivity mixture contains invalid components.")
    dax_sens_assert(is.list(crps_plans) && length(crps_plans) == spec$blocks, "CRPS block plans have the wrong dimensions.")
    block_log_sums <- numeric(spec$blocks)
    block_draws <- matrix(NA_real_, nrow = spec$crps_draws_per_block, ncol = spec$blocks)
    for (g in seq_len(spec$blocks)) {
        columns <- ((g - 1L) * block_size + 1L):(g * block_size)
        log_density <- dstd_t(realized, mean = rep(means[columns], each = spec$P), sd = as.vector(sds[, 
            columns, drop = FALSE]), nu = rep(nus[columns], each = spec$P), log = TRUE)
        block_log_sums[[g]] <- log_sum_exp(log_density)
        plan <- crps_plans[[g]]
        dax_sens_assert(length(plan$component_uniforms) == spec$crps_draws_per_block && length(plan$t_uniforms) == 
            spec$crps_draws_per_block, "A CRPS block plan has the wrong draw supply.")
        n_block_components <- length(columns) * spec$P
        component_index <- dax_sens_stratified_component_indices(n_block_components, plan$component_uniforms)
        local_outer <- (component_index - 1L)%/%spec$P + 1L
        particle <- (component_index - 1L)%%spec$P + 1L
        outer <- columns[local_outer]
        selected_sd <- sds[cbind(particle, outer)]
        block_draws[, g] <- means[outer] + std_t_scale(selected_sd, nus[outer]) * stats::qt(plan$t_uniforms, 
            df = nus[outer])
    }
    scopes <- c("pooled", paste0("loo_", seq_len(spec$blocks)), paste0("nested_", spec$nested_B[spec$nested_B < 
        spec$B]))
    rows <- vector("list", length(scopes))
    for (i in seq_along(scopes)) {
        scope <- scopes[[i]]
        if (scope == "pooled") {
            included_blocks <- seq_len(spec$blocks)
        }
        else if (startsWith(scope, "loo_")) {
            omitted <- as.integer(sub("loo_", "", scope, fixed = TRUE))
            included_blocks <- setdiff(seq_len(spec$blocks), omitted)
        }
        else {
            nested_B <- as.integer(sub("nested_", "", scope, fixed = TRUE))
            included_blocks <- seq_len(nested_B%/%block_size)
        }
        columns <- unlist(lapply(included_blocks, function(g) {
            ((g - 1L) * block_size + 1L):(g * block_size)
        }), use.names = FALSE)
        n_components <- length(columns) * spec$P
        lpds <- log_sum_exp(block_log_sums[included_blocks]) - log(n_components)
        draws <- as.vector(block_draws[, included_blocks, drop = FALSE])
        crps <- empirical_crps(draws, realized)
        dax_sens_assert(all(is.finite(c(lpds, crps))), "Sensitivity score is nonfinite.")
        rows[[i]] <- data.frame(scope = scope, n_components = n_components, n_crps_draws = length(draws), 
            lpds = lpds, crps = crps, stringsAsFactors = FALSE)
    }
    bind_rows_base(rows)
}

dax_sens_training_slice <- function (inputs, refit_month) 
{
    spec <- dax_sens_spec()
    month_dates <- inputs$support_dates[format(inputs$support_dates, "%Y-%m") == refit_month]
    first <- match(month_dates[[1]], inputs$data$date)
    training <- inputs$data[seq.int(first - spec$estimation_window, first - 1L), ]
    list(training = training, dates = month_dates, realized = inputs$data$return[match(month_dates, inputs$data$date)])
}

dax_sens_month_task <- function (inputs, model, refit_month, outer_shocks) 
{
    spec <- dax_sens_spec()
    fit <- inputs$fit
    nodes <- dax_sens_parameter_nodes(fit, model, refit_month, outer_shocks)
    slice <- dax_sens_training_slice(inputs, refit_month)
    theta_hat <- as.list(nodes$checked$theta)
    state_p <- vector("list", spec$B)
    state_i <- vector("list", spec$B)
    for (b in seq_len(spec$B)) {
        key <- outer_shocks[b, ]
        seed <- dax_sens_inner_seed(model, refit_month, key$block_index, key$pair_index, key$sign, 0L)
        state_p[[b]] <- initialize_sv_filter(model, theta_hat, slice$training$return, spec$P, seed, spec$resample_ess_fraction, 
            star_m = if (model == "STAR_SV") 
                5L
            else 5L, training_dates = slice$training$date, dead_fraction_limit = spec$particle_dead_fraction_limit)
        state_i[[b]] <- initialize_sv_filter(model, nodes$parameters[[b]], slice$training$return, spec$P, 
            seed, spec$resample_ess_fraction, star_m = if (model == "STAR_SV") 
                5L
            else 5L, training_dates = slice$training$date, dead_fraction_limit = spec$particle_dead_fraction_limit)
    }
    daily <- vector("list", length(slice$dates) * 2L)
    max_dead <- c(P = 0, I = 0)
    min_ess <- c(P = Inf, I = Inf)
    date_index <- match(slice$dates, inputs$support_dates)
    dax_sens_assert(!anyNA(date_index), "A monthly forecast date is absent from specified support.")
    for (j in seq_along(slice$dates)) {
        means_p <- means_i <- nus_p <- nus_i <- numeric(spec$B)
        sds_p <- matrix(NA_real_, nrow = spec$P, ncol = spec$B)
        sds_i <- matrix(NA_real_, nrow = spec$P, ncol = spec$B)
        for (b in seq_len(spec$B)) {
            key <- outer_shocks[b, ]
            seed <- dax_sens_inner_seed(model, refit_month, key$block_index, key$pair_index, key$sign, 
                j)
            step_p <- dax_sens_sv_predict_components_update(state_p[[b]], theta_hat, slice$realized[[j]], 
                seed, spec$resample_ess_fraction, as.character(slice$dates[[j]]), spec$particle_dead_fraction_limit)
            step_i <- dax_sens_sv_predict_components_update(state_i[[b]], nodes$parameters[[b]], slice$realized[[j]], 
                seed, spec$resample_ess_fraction, as.character(slice$dates[[j]]), spec$particle_dead_fraction_limit)
            state_p[[b]] <- step_p$state
            state_i[[b]] <- step_i$state
            means_p[[b]] <- step_p$components$mean
            means_i[[b]] <- step_i$components$mean
            nus_p[[b]] <- step_p$components$nu
            nus_i[[b]] <- step_i$components$nu
            sds_p[, b] <- step_p$components$sd
            sds_i[, b] <- step_i$components$sd
            max_dead[["P"]] <- max(max_dead[["P"]], step_p$particle_diagnostics$max_dead_fraction)
            max_dead[["I"]] <- max(max_dead[["I"]], step_i$particle_diagnostics$max_dead_fraction)
            min_ess[["P"]] <- min(min_ess[["P"]], step_p$filter_ess)
            min_ess[["I"]] <- min(min_ess[["I"]], step_i$filter_ess)
        }
        plans <- dax_sens_crps_plans(date_index[[j]])
        score_p <- dax_sens_score_components(means_p, sds_p, nus_p, slice$realized[[j]], plans)
        score_i <- dax_sens_score_components(means_i, sds_i, nus_i, slice$realized[[j]], plans)
        score_p$variant <- "P"
        score_i$variant <- "I"
        base <- data.frame(date = as.character(slice$dates[[j]]), model = model, refit_month = refit_month, 
            stringsAsFactors = FALSE)
        daily[[2L * j - 1L]] <- cbind(base, score_p[, c("variant", "scope", "n_components", "n_crps_draws", 
            "lpds", "crps")])
        daily[[2L * j]] <- cbind(base, score_i[, c("variant", "scope", "n_components", "n_crps_draws", 
            "lpds", "crps")])
    }
    diagnostics <- data.frame(model = model, refit_month = refit_month, variant = c("P", "I"), max_particle_dead_fraction = unname(max_dead), 
        minimum_forecast_ess = unname(min_ess), stringsAsFactors = FALSE)
    value <- list(parameter_draws = nodes$frame, daily = bind_rows_base(daily), diagnostics = diagnostics)
    value
}

