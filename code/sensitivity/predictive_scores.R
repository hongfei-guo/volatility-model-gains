std_t_scale <- function (sd, nu) 
{
    sd <- as.numeric(sd)
    nu <- as.numeric(nu)
    if (any(!is.finite(sd)) || any(sd <= 0)) {
        stop("Student-t standard deviations must be positive and finite.")
    }
    if (any(!is.finite(nu)) || any(nu <= 2)) {
        stop("Student-t degrees of freedom must exceed 2.")
    }
    sd * sqrt((nu - 2)/nu)
}

dstd_t <- function (x, mean = 0, sd = 1, nu, log = FALSE) 
{
    scale <- std_t_scale(sd, nu)
    z <- (x - mean)/scale
    value <- stats::dt(z, df = nu, log = TRUE) - log(scale)
    if (log) 
        value
    else exp(value)
}

pstd_t <- function (q, mean = 0, sd = 1, nu) 
{
    scale <- std_t_scale(sd, nu)
    stats::pt((q - mean)/scale, df = nu)
}

qstd_t <- function (p, mean = 0, sd = 1, nu) 
{
    p <- as.numeric(p)
    if (any(!is.finite(p)) || any(p < 0 | p > 1)) {
        stop("Probabilities must lie in [0, 1].")
    }
    mean + std_t_scale(sd, nu) * stats::qt(p, df = nu)
}

normalize_mixture_inputs <- function (component_mean, component_sd, nu, weights = NULL) 
{
    component_sd <- as.numeric(component_sd)
    if (!length(component_sd) || any(!is.finite(component_sd)) || any(component_sd <= 0)) {
        stop("Mixture component standard deviations must be positive and finite.")
    }
    component_mean <- rep_len(as.numeric(component_mean), length(component_sd))
    if (any(!is.finite(component_mean))) 
        stop("Mixture component means must be finite.")
    nu <- rep_len(as.numeric(nu), length(component_sd))
    if (any(!is.finite(nu)) || any(nu <= 2)) 
        stop("Every mixture component requires nu > 2.")
    if (is.null(weights)) 
        weights <- rep(1/length(component_sd), length(component_sd))
    weights <- as.numeric(weights)
    if (length(weights) != length(component_sd) || any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 
        0) {
        stop("Mixture weights must be finite, nonnegative, and conformable.")
    }
    weights <- weights/sum(weights)
    list(mean = component_mean, sd = component_sd, nu = nu, weights = weights)
}

mixture_log_density <- function (x, component_mean, component_sd, nu, weights = NULL) 
{
    mix <- normalize_mixture_inputs(component_mean, component_sd, nu, weights)
    weighted_log_sum_exp(dstd_t(x, mix$mean, mix$sd, mix$nu, log = TRUE), mix$weights)
}

mixture_cdf <- function (x, component_mean, component_sd, nu, weights = NULL) 
{
    mix <- normalize_mixture_inputs(component_mean, component_sd, nu, weights)
    sum(mix$weights * pstd_t(x, mix$mean, mix$sd, mix$nu))
}

student_t_lower_partial_mean <- function (q, mean, sd, nu) 
{
    scale <- std_t_scale(sd, nu)
    z <- (q - mean)/scale
    mean * stats::pt(z, df = nu) - scale * ((nu + z^2)/(nu - 1)) * stats::dt(z, df = nu)
}

mixture_quantile <- function (probability, component_mean, component_sd, nu, weights = NULL, tolerance = 1e-10, max_iterations = 200L) 
{
    if (length(probability) != 1L || !is.finite(probability) || probability <= 0 || probability >= 1) {
        stop("Mixture quantile probability must lie strictly in (0, 1).")
    }
    mix <- normalize_mixture_inputs(component_mean, component_sd, nu, weights)
    lower_probability <- max(min(probability * 1e-04, probability/10), 1e-12)
    upper_probability <- min(max(1 - (1 - probability) * 1e-04, probability + 0.9 * (1 - probability)), 
        1 - 1e-12)
    lower <- min(qstd_t(lower_probability, mix$mean, mix$sd, mix$nu))
    upper <- max(qstd_t(upper_probability, mix$mean, mix$sd, mix$nu))
    objective <- function(x) mixture_cdf(x, mix$mean, mix$sd, mix$nu, mix$weights) - probability
    f_lower <- objective(lower)
    f_upper <- objective(upper)
    expansion <- max(diff(range(c(lower, upper))), max(mix$sd), 1)
    counter <- 0L
    while (f_lower > 0 && counter < 50L) {
        lower <- lower - expansion
        expansion <- expansion * 2
        f_lower <- objective(lower)
        counter <- counter + 1L
    }
    expansion <- max(diff(range(c(lower, upper))), max(mix$sd), 1)
    counter <- 0L
    while (f_upper < 0 && counter < 50L) {
        upper <- upper + expansion
        expansion <- expansion * 2
        f_upper <- objective(upper)
        counter <- counter + 1L
    }
    if (!is.finite(f_lower) || !is.finite(f_upper) || f_lower > 0 || f_upper < 0) {
        stop("Could not bracket the predictive-mixture quantile.")
    }
    stats::uniroot(objective, interval = c(lower, upper), tol = tolerance, maxiter = as.integer(max_iterations), 
        check.conv = TRUE)$root
}

mixture_loss_var_es <- function (confidence, component_mean, component_sd, nu, weights = NULL) 
{
    if (length(confidence) != 1L || confidence <= 0.5 || confidence >= 1) {
        stop("Loss confidence must lie in (0.5, 1).")
    }
    mix <- normalize_mixture_inputs(component_mean, component_sd, nu, weights)
    alpha <- 1 - confidence
    return_quantile <- mixture_quantile(alpha, mix$mean, mix$sd, mix$nu, mix$weights)
    lower_moment <- sum(mix$weights * student_t_lower_partial_mean(return_quantile, mix$mean, mix$sd, 
        mix$nu))
    var_loss <- -return_quantile
    es_loss <- -lower_moment/alpha
    validate_var_es(var_loss, es_loss, confidence)
    c(VaR = as.numeric(var_loss), ES = as.numeric(es_loss))
}

stratified_uniforms <- function (n) 
{
    n <- as.integer(n)
    if (n < 1L) 
        return(numeric())
    (seq_len(n) - stats::runif(n))/n
}

independently_permuted_stratified_uniforms <- function (n) 
{
    n <- as.integer(n)
    u <- stratified_uniforms(n)
    if (n <= 1L) 
        return(u)
    u[sample.int(n, size = n, replace = FALSE)]
}

empirical_crps <- function (draws, realized) 
{
    draws <- sort(as.numeric(draws))
    n <- length(draws)
    if (!n || any(!is.finite(draws)) || !is.finite(realized)) 
        return(NA_real_)
    if (n == 1L) 
        return(abs(draws - realized))
    first <- mean(abs(draws - realized))
    weights <- 2 * seq_len(n) - n - 1
    half_pair_expectation <- sum(weights * draws)/n^2
    first - half_pair_expectation
}
