effective_sample_size <- function (weights) 
{
    weights <- as.numeric(weights)
    if (!length(weights) || any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) 
        stop("Invalid particle weights.")
    weights <- weights/sum(weights)
    1/sum(weights^2)
}

systematic_resample <- function (weights, n = length(weights)) 
{
    weights <- as.numeric(weights)
    if (any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
        stop("Invalid particle weights.")
    }
    n <- as.integer(n)
    if (n < 1L) 
        stop("At least one particle is required.")
    weights <- weights/sum(weights)
    positions <- ((0:(n - 1L)) + stats::runif(1))/n
    pmin(findInterval(positions, cumsum(weights)) + 1L, length(weights))
}

normalize_log_weights <- function (log_likelihood, previous_weights = NULL) 
{
    log_likelihood <- as.numeric(log_likelihood)
    n <- length(log_likelihood)
    if (!n) 
        stop("No particle likelihoods were supplied.")
    if (is.null(previous_weights)) 
        previous_weights <- rep(1/n, n)
    previous_weights <- as.numeric(previous_weights)
    if (length(previous_weights) != n || any(!is.finite(previous_weights)) || any(previous_weights < 
        0) || sum(previous_weights) <= 0) {
        stop("Invalid previous particle weights.")
    }
    previous_weights <- previous_weights/sum(previous_weights)
    if (any(is.na(log_likelihood)) || any(log_likelihood == Inf)) {
        stop("Particle log likelihoods must be finite or negative infinity.")
    }
    log_previous_weight <- rep(-Inf, n)
    positive <- previous_weights > 0
    log_previous_weight[positive] <- log(previous_weights[positive])
    log_weight <- log_previous_weight + log_likelihood
    normalizer <- log_sum_exp(log_weight)
    if (!is.finite(normalizer)) {
        stop("All particle likelihoods underflowed or were non-finite.")
    }
    weights <- exp(log_weight - normalizer)
    list(weights = weights, log_normalizer = normalizer, ess = effective_sample_size(weights))
}

particle_log_density <- function (x, mean, h, nu) 
{
    h <- as.numeric(h)
    n <- length(h)
    if (!n) 
        stop("Particle log-density evaluation received no states.")
    x <- rep_len(as.numeric(x), n)
    mean <- rep_len(as.numeric(mean), n)
    nu <- rep_len(as.numeric(nu), n)
    out <- rep(-Inf, n)
    valid <- is.finite(x) & is.finite(mean) & is.finite(h) & is.finite(nu) & nu > 2
    if (!any(valid)) 
        return(out)
    residual <- x[valid] - mean[valid]
    nu_valid <- nu[valid]
    log_scale <- 0.5 * h[valid] + 0.5 * (log(nu_valid - 2) - log(nu_valid))
    log_kernel <- numeric(length(residual))
    nonzero <- residual != 0
    if (any(nonzero)) {
        log_abs_z <- log(abs(residual[nonzero])) - log_scale[nonzero]
        ratio_log <- 2 * log_abs_z - log(nu_valid[nonzero])
        stable_kernel <- ratio_log
        moderate <- ratio_log <= 50
        stable_kernel[moderate] <- log1p(exp(ratio_log[moderate]))
        log_kernel[nonzero] <- stable_kernel
    }
    constant <- lgamma((nu_valid + 1)/2) - lgamma(nu_valid/2) - 0.5 * (log(nu_valid) + log(pi))
    value <- constant - 0.5 * (nu_valid + 1) * log_kernel - log_scale
    value[is.na(value) | value == Inf] <- -Inf
    out[valid] <- value
    out
}

particle_diagnostics_update <- function (previous, h, weights, prior_weights, observation_id, dead_fraction_limit) 
{
    h <- as.numeric(h)
    weights <- as.numeric(weights)
    prior_weights <- as.numeric(prior_weights)
    if (length(h) != length(weights) || length(weights) != length(prior_weights) || any(!is.finite(h)) || 
        any(!is.finite(weights)) || any(weights < 0)) {
        stop("Particle diagnostic state is invalid.")
    }
    dead <- weights == 0
    newly_dead <- prior_weights > 0 & dead
    dead_count <- sum(dead)
    dead_fraction <- dead_count/length(weights)
    if (dead_fraction > dead_fraction_limit) {
        stop(sprintf(paste0("Particle-death gate failed at %s: %d/%d particles have zero mass ", "(limit %.3f)."), 
            as.character(observation_id %||% "unknown"), dead_count, length(weights), dead_fraction_limit))
    }
    live <- weights > 0
    if (!any(live)) 
        stop("Particle filter has no live particles.")
    previous <- previous %||% list(cumulative_dead = 0L, first_dead_observation = NA_character_, max_dead_fraction = 0)
    first_dead <- previous$first_dead_observation
    if (sum(newly_dead) > 0L && (length(first_dead) != 1L || is.na(first_dead) || !nzchar(first_dead))) {
        first_dead <- as.character(observation_id %||% "unknown")
    }
    list(last_dead_count = as.integer(dead_count), last_dead_fraction = as.numeric(dead_fraction), last_newly_dead_count = as.integer(sum(newly_dead)), 
        cumulative_dead = as.integer(previous$cumulative_dead + sum(newly_dead)), first_dead_observation = first_dead, 
        max_dead_fraction = max(previous$max_dead_fraction, dead_fraction), live_h_min = min(h[live]), 
        live_h_max = max(h[live]))
}

star_trigger_one <- function (previous_return, recent_returns, trailing = 5L) 
{
    recent_returns <- tail(as.numeric(recent_returns), as.integer(trailing))
    if (!is.finite(previous_return) || previous_return >= 0) 
        return(0L)
    reference <- mean(recent_returns, na.rm = TRUE)
    if (is.finite(reference) && previous_return <= reference) 
        2L
    else 1L
}

sv_transition_mean <- function (model, params, h_prev, resid_prev, previous_return, recent_returns, star_m = 5L) 
{
    model <- toupper(model)
    h_prev <- as.numeric(h_prev)
    base <- params[["mu_h"]] + params[["phi_h"]] * (h_prev - params[["mu_h"]])
    if (model == "ARSV") 
        return(base)
    if (model == "AARSV") {
        standardized_residual <- resid_prev/exp(0.5 * h_prev)
        return(base + params[["gamma_aarsv"]] * standardized_residual)
    }
    if (model == "STAR_SV") {
        trigger <- star_trigger_one(previous_return, recent_returns, trailing = star_m)
        gamma <- if (trigger == 2L) {
            params[["gamma_severe"]]
        }
        else if (trigger == 1L) {
            params[["gamma_mild"]]
        }
        else 0
        return(base + gamma * resid_prev/exp(0.5 * h_prev))
    }
    stop(sprintf("Unsupported SV model '%s'.", model))
}

validate_sv_parameters <- function (model, params) 
{
    required <- c("mu_y", "phi_y", "mu_h", "phi_h", "sigma_eta", "nu")
    extra <- switch(toupper(model), ARSV = character(), AARSV = "gamma_aarsv", STAR_SV = c("gamma_severe", 
        "gamma_mild"), stop("Unsupported SV model."))
    missing <- setdiff(c(required, extra), names(params))
    if (length(missing)) {
        stop(sprintf("%s parameters are missing: %s", model, paste(missing, collapse = ", ")))
    }
    values <- unlist(params[c(required, extra)], use.names = FALSE)
    if (any(!is.finite(values))) 
        stop("SV parameters must be finite.")
    if (abs(params[["phi_y"]]) >= 1 || abs(params[["phi_h"]]) >= 1 || params[["sigma_eta"]] < 0 || params[["nu"]] <= 
        2) {
        stop(sprintf("%s common parameter constraints failed.", model))
    }
    if (toupper(model) == "STAR_SV" && !(params[["gamma_severe"]] <= params[["gamma_mild"]] && params[["gamma_mild"]] <= 
        0)) {
        stop("STAR-SV coefficient ordering failed.")
    }
    invisible(TRUE)
}

sv_initial_state_sd <- function (model, params) 
{
    model <- toupper(model)
    innovation_variance <- params[["sigma_eta"]]^2
    if (model == "STAR_SV") {
        return(params[["sigma_eta"]])
    }
    if (model == "AARSV") {
        innovation_variance <- innovation_variance + params[["gamma_aarsv"]]^2
    }
    else if (model != "ARSV") {
        stop(sprintf("Unsupported SV model '%s'.", model))
    }
    denominator <- max(1 - params[["phi_h"]]^2, 1e-08)
    value <- sqrt(innovation_variance/denominator)
    if (!is.finite(value) || value < 0) 
        stop("Invalid SV initial-state scale.")
    value
}

sv_filter_observe <- function (state, params, realized_return, seed = NULL, resample_ess_fraction = 0.5, observation_id = NA_character_, 
    dead_fraction_limit = 0.05) 
{
    if (!is.null(seed)) 
        set.seed(as.integer(seed))
    n_particles <- length(state$h)
    threshold <- resample_ess_fraction * n_particles
    if (effective_sample_size(state$weights) < threshold) {
        ancestors <- systematic_resample(state$weights)
        h_prev <- state$h[ancestors]
        prior_weights <- rep(1/n_particles, n_particles)
    }
    else {
        h_prev <- state$h
        prior_weights <- state$weights/sum(state$weights)
    }
    star_m <- as.integer(state$star_m %||% 5L)
    if (any(prior_weights > 0 & !is.finite(h_prev))) {
        stop("A live particle has a non-finite latent state before propagation.")
    }
    live_before <- prior_weights > 0
    h_new <- h_prev
    valid_proposal <- rep(FALSE, n_particles)
    if (any(live_before)) {
        live_index <- which(live_before)
        transition_mean <- sv_transition_mean(state$model, params, h_prev[live_before], state$last_resid, 
            state$last_y, state$recent_returns, star_m)
        proposed <- transition_mean + params[["sigma_eta"]] * stats::rnorm(sum(live_before))
        finite_proposal <- is.finite(proposed)
        valid_index <- live_index[finite_proposal]
        h_new[valid_index] <- proposed[finite_proposal]
        valid_proposal[valid_index] <- TRUE
    }
    mean_new <- params[["mu_y"]] + params[["phi_y"]] * (state$last_y - params[["mu_y"]])
    log_likelihood <- rep(-Inf, n_particles)
    log_likelihood[valid_proposal] <- particle_log_density(realized_return, mean_new, h_new[valid_proposal], 
        params[["nu"]])
    normalized <- normalize_log_weights(log_likelihood, prior_weights)
    particle_diagnostics <- particle_diagnostics_update(state$particle_diagnostics, h_new, normalized$weights, 
        prior_weights, observation_id, dead_fraction_limit)
    filtered_mean <- sum(normalized$weights * h_new)
    filtered_sd <- sqrt(sum(normalized$weights * (h_new - filtered_mean)^2))
    new_state <- list(model = state$model, h = h_new, weights = normalized$weights, last_y = realized_return, 
        last_mean = mean_new, last_resid = realized_return - mean_new, recent_returns = tail(c(state$recent_returns, 
            realized_return), star_m), star_m = star_m, ess = normalized$ess, log_likelihood_increment = normalized$log_normalizer, 
        particle_diagnostics = particle_diagnostics)
    list(state = new_state, filtered_h_mean = filtered_mean, filtered_h_sd = filtered_sd, filter_ess = normalized$ess)
}

initialize_sv_filter <- function (model, params, training_returns, n_particles, seed, resample_ess_fraction = 0.5, star_m = 5L, 
    training_dates = NULL, dead_fraction_limit = 0.05) 
{
    model <- toupper(model)
    y <- safe_numeric(training_returns, "training_returns")
    validate_sv_parameters(model, params)
    if (is.null(training_dates)) {
        training_dates <- rep(NA_character_, length(y))
    }
    else {
        training_dates <- as.character(as.Date(training_dates))
        if (length(training_dates) != length(y) || anyNA(training_dates)) {
            stop("SV filter training dates must match the return series.")
        }
    }
    set.seed(as.integer(seed))
    stationary_sd <- sv_initial_state_sd(model, params)
    h <- stats::rnorm(as.integer(n_particles), params[["mu_h"]], stationary_sd)
    mean_t <- params[["mu_y"]]
    initial_prior_weights <- rep(1/length(h), length(h))
    first_weights <- normalize_log_weights(particle_log_density(y[1], mean_t, h, params[["nu"]]), initial_prior_weights)
    particle_diagnostics <- particle_diagnostics_update(NULL, h, first_weights$weights, initial_prior_weights, 
        training_dates[[1]], dead_fraction_limit)
    state <- list(model = model, h = h, weights = first_weights$weights, last_y = y[1], last_mean = mean_t, 
        last_resid = y[1] - mean_t, recent_returns = y[1], star_m = as.integer(star_m), ess = first_weights$ess, 
        log_likelihood_increment = first_weights$log_normalizer, particle_diagnostics = particle_diagnostics)
    if (length(y) >= 2L) {
        for (t in 2:length(y)) {
            observed <- sv_filter_observe(state, params, y[t], seed = NULL, resample_ess_fraction = resample_ess_fraction, 
                observation_id = training_dates[[t]], dead_fraction_limit = dead_fraction_limit)
            state <- observed$state
        }
    }
    state$recent_returns <- tail(y, state$star_m)
    state
}

sv_predict_components_update <- function (state, params, realized_return, seed, resample_ess_fraction = 0.5, observation_id = NA_character_, 
    dead_fraction_limit = 0.05, component_consumer = NULL) 
{
    n_particles <- length(state$h)
    set.seed(as.integer(seed))
    ancestors <- systematic_resample(state$weights, n_particles)
    h_prev <- state$h[ancestors]
    if (any(!is.finite(h_prev)) || any(state$weights[ancestors] <= 0)) {
        stop("Predictive resampling selected a dead or non-finite particle.")
    }
    star_m <- as.integer(state$star_m %||% 5L)
    transition_mean <- sv_transition_mean(state$model, params, h_prev, state$last_resid, state$last_y, 
        state$recent_returns, star_m)
    h_pred <- transition_mean + params[["sigma_eta"]] * stats::rnorm(n_particles)
    mean_next <- params[["mu_y"]] + params[["phi_y"]] * (state$last_y - params[["mu_y"]])
    sigma_pred <- exp(0.5 * h_pred)
    if (any(!is.finite(h_pred)) || any(!is.finite(sigma_pred)) || any(sigma_pred <= 0)) {
        stop(sprintf("A live predictive particle produced an invalid scale at %s.", as.character(observation_id %||% 
            "unknown")))
    }
    component_weights <- rep(1/n_particles, n_particles)
    components <- list(mean = rep(mean_next, n_particles), sd = sigma_pred, weights = component_weights, 
        nu = rep(params[["nu"]], n_particles))
    component_result <- if (is.null(component_consumer)) {
        NULL
    }
    else {
        component_consumer(components)
    }
    log_likelihood <- particle_log_density(realized_return, mean_next, h_pred, params[["nu"]])
    normalized <- normalize_log_weights(log_likelihood, component_weights)
    particle_diagnostics <- particle_diagnostics_update(state$particle_diagnostics, h_pred, normalized$weights, 
        component_weights, observation_id, dead_fraction_limit)
    filtered_mean <- sum(normalized$weights * h_pred)
    filtered_sd <- sqrt(sum(normalized$weights * (h_pred - filtered_mean)^2))
    h_filtered <- h_pred
    weights_filtered <- normalized$weights
    if (normalized$ess < resample_ess_fraction * n_particles) {
        resampled <- systematic_resample(normalized$weights)
        h_filtered <- h_pred[resampled]
        weights_filtered <- rep(1/n_particles, n_particles)
    }
    new_state <- list(model = state$model, h = h_filtered, weights = weights_filtered, last_y = realized_return, 
        last_mean = mean_next, last_resid = realized_return - mean_next, recent_returns = tail(c(state$recent_returns, 
            realized_return), star_m), star_m = star_m, ess = effective_sample_size(weights_filtered), 
        log_likelihood_increment = normalized$log_normalizer, particle_diagnostics = particle_diagnostics)
    list(state = new_state, components = components, component_result = component_result, filter_metrics = list(filter_ess = normalized$ess, 
        filtered_h_mean = filtered_mean, filtered_h_sd = filtered_sd, particle_dead_count = particle_diagnostics$last_dead_count, 
        particle_dead_fraction = particle_diagnostics$last_dead_fraction, particle_dead_cumulative = particle_diagnostics$cumulative_dead, 
        particle_first_dead_observation = particle_diagnostics$first_dead_observation, particle_max_dead_fraction = particle_diagnostics$max_dead_fraction, 
        particle_live_h_min = particle_diagnostics$live_h_min, particle_live_h_max = particle_diagnostics$live_h_max), 
        filter_ess = normalized$ess, particle_diagnostics = particle_diagnostics)
}

