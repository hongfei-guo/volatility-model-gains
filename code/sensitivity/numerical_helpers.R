`%||%` <- function (x, y) 
{
    if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) 
        y
    else x
}

safe_numeric <- function (x, name = deparse(substitute(x)), allow_na = FALSE) 
{
    x <- as.numeric(x)
    if (!length(x)) 
        stop(sprintf("%s is empty.", name))
    if (!allow_na && any(!is.finite(x))) 
        stop(sprintf("%s contains non-finite values.", name))
    x
}

log_sum_exp <- function (log_values) 
{
    log_values <- as.numeric(log_values)
    if (!length(log_values)) 
        return(-Inf)
    m <- max(log_values)
    if (!is.finite(m)) 
        return(m)
    m + log(sum(exp(log_values - m)))
}

weighted_log_sum_exp <- function (log_values, weights) 
{
    log_values <- as.numeric(log_values)
    weights <- as.numeric(weights)
    if (length(log_values) != length(weights)) 
        stop("log_values and weights must have equal length.")
    if (any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
        stop("Weights must be finite, nonnegative, and have positive sum.")
    }
    weights <- weights/sum(weights)
    positive <- weights > 0
    log_sum_exp(log_values[positive] + log(weights[positive]))
}

bind_rows_base <- function (objects) 
{
    objects <- Filter(function(x) !is.null(x) && NROW(x) > 0L, objects)
    if (!length(objects)) 
        return(data.frame())
    columns <- unique(unlist(lapply(objects, names), use.names = FALSE))
    filled <- lapply(objects, function(x) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
        missing <- setdiff(columns, names(x))
        for (name in missing) x[[name]] <- NA
        x[columns]
    })
    out <- do.call(rbind, filled)
    rownames(out) <- NULL
    out
}

validate_var_es <- function (var, es, confidence, tolerance = 1e-10) 
{
    if (any(!is.finite(c(var, es)))) 
        stop("VaR/ES forecasts must be finite.")
    if (es + tolerance < var) {
        stop(sprintf("ES must not be below VaR at confidence %.3f.", 
            confidence))
    }
    invisible(TRUE)
}
