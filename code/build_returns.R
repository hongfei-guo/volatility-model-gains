#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve script path.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
archive_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

input_dir <- file.path(archive_root, "data", "source")
output_dir <- file.path(archive_root, "reproduced", "returns")

registry <- data.frame(
  market = c("DAX", "FTSE", "NIKKEI", "SHCOMP", "SP500"),
  ticker = c("^GDAXI", "^FTSE", "^N225", "000001.SS", "^GSPC"),
  first_date = as.Date(c(
    "2000-01-04", "2000-01-05", "2000-01-05", "2000-01-05", "2000-01-04"
  )),
  last_date = as.Date(c(
    "2025-12-30", "2025-12-31", "2025-12-30", "2025-12-31", "2025-12-31"
  )),
  observations = c(6600L, 6566L, 6368L, 6294L, 6538L),
  stringsAsFactors = FALSE
)

normalize_names <- function(names) {
  gsub("[^a-z0-9]+", "", tolower(names))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(registry))) {
  market <- registry$market[[i]]
  input_path <- file.path(input_dir, paste0(market, ".csv"))
  if (!file.exists(input_path)) {
    stop(
      "Missing source data for ", market, ": ", input_path,
      ". See data/README.md.", call. = FALSE
    )
  }

  prices <- utils::read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
  normalized <- normalize_names(names(prices))
  date_column <- which(normalized == "date")
  close_column <- which(normalized == "close")
  if (length(date_column) != 1L || length(close_column) != 1L) {
    stop(market, " source file must contain exactly one Date column and one Close column.")
  }

  prices <- data.frame(
    date = as.Date(prices[[date_column]]),
    close = suppressWarnings(as.numeric(prices[[close_column]])),
    stringsAsFactors = FALSE
  )
  if (anyNA(prices$date) || any(!is.finite(prices$close)) || any(prices$close <= 0)) {
    stop(market, " source dates and closing levels must be complete and valid.")
  }
  if (anyDuplicated(prices$date)) stop(market, " source file contains duplicate dates.")
  prices <- prices[order(prices$date), , drop = FALSE]

  returns <- data.frame(
    date = prices$date[-1L],
    return = 100 * diff(log(prices$close)),
    stringsAsFactors = FALSE
  )
  returns <- returns[
    returns$date >= registry$first_date[[i]] &
      returns$date <= registry$last_date[[i]],
    , drop = FALSE
  ]

  if (nrow(returns) != registry$observations[[i]] ||
      min(returns$date) != registry$first_date[[i]] ||
      max(returns$date) != registry$last_date[[i]]) {
    stop(
      market, " reconstructed series does not match the required calendar. ",
      "The source file must include at least one trading day before ",
      registry$first_date[[i]], " and all observations through ",
      registry$last_date[[i]], ".", call. = FALSE
    )
  }

  utils::write.csv(
    data.frame(date = as.character(returns$date), return = returns$return),
    file.path(output_dir, paste0(market, ".csv")),
    row.names = FALSE, quote = TRUE
  )
}

message("Return construction completed successfully. Outputs: ", output_dir)
