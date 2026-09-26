# Computations supporting appendix tables and the power discussion.
rebuild_holdout_sensitivity <- function(root, holdout_panel, support) {
  daily <- read.csv(file.path(root,'results/supplementary/holdout_parameter_uncertainty/daily_predictions.csv'),stringsAsFactors=FALSE)
  assert_columns(daily,c('scope','B','VaR','ES','date','model','variant','market'),'Holdout sensitivity predictions')
  returns <- unique(holdout_panel[c('market','date','realized_return')])
  key <- function(x) paste(x$market,x$date)
  assert_true(!anyDuplicated(key(returns)), 'Ambiguous return matching.')
  i <- match(key(daily),key(returns))
  assert_true(!anyNA(i),'Sensitivity dates are absent from the supplied return series.')
  daily$FZ0 <- fz0_loss(-returns$realized_return[i],daily$VaR,daily$ES,0.975)
  daily$task_id <- paste(daily$market,format(as.Date(daily$date),'%Y-%m'),sep='_')
  support <- support[support$market %in% c('DAX','NIKKEI'),c('market','date')]
  support$date <- as.character(support$date)
  support <- support[order(match(support$market,c('DAX','NIKKEI')),support$date),]
  original <- holdout_panel[holdout_panel$market %in% c('DAX','NIKKEI') & holdout_panel$model %in% c('ARSV','STAR_SV'),c('market','date','model','fz0_0975')]
  original$date <- as.character(original$date)
  names(original)[4] <- 'FZ0'
  hpu_summary(daily,original,support)
}

parameter_summary <- function(path) {
  x <- read.csv(path,stringsAsFactors=FALSE)
  mapping <- c(AARSV='gamma_aarsv',APARCH='gamma',GJR_R='lambda',GJR_U='lambda',TSV_RT='delta_tsv')
  keep <- (x$model %in% names(mapping) & x$parameter==unname(mapping[x$model])) |
    (x$model=='STAR_SV' & x$parameter %in% c('gamma_mild','gamma_severe'))
  keep[is.na(keep)] <- FALSE
  x <- x[keep,]
  assert_true(!anyDuplicated(x[c('period','market','model','refit_month','parameter')]) && all(is.finite(x$estimate)),
              'Monthly parameter records are duplicated or nonfinite.')
  groups <- split(x,interaction(x$period,x$market,x$model,x$parameter,drop=TRUE))
  out <- do.call(rbind,lapply(groups,function(z) data.frame(z[1,c('period','market','model','parameter')],
    months=nrow(z),median=median(z$estimate),minimum=min(z$estimate),maximum=max(z$estimate))))
  rownames(out) <- NULL
  assert_true(nrow(out)==70L,'Asymmetry parameter summary has incomplete coverage.')
  out[order(out$market,out$model,out$parameter,out$period),]
}

power_summary <- function(table3) {
  z <- qnorm(0.975); zpower <- qnorm(0.80)
  mu <- table3$development_mean_fz0_difference; se <- table3$holdout_hac_se
  power <- function(se) pnorm(-z-mu/se)+pnorm(mu/se-z)
  data.frame(market=table3$market,main_oos_difference=mu,holdout_n=table3$holdout_n,
    holdout_hac_se=se,mde_80pct=(z+zpower)*se,
    rejection_probability=power(se),
    rejection_probability_displayed_se=power(round(se,3)),
    approximate_required_days=table3$holdout_n*((z+zpower)*se/abs(mu))^2)
}

gjr_table <- function(period_results) {
  out <- list()
  for(period in names(period_results)) {
    x <- period_results[[period]]$common_panel
    for(market in market_order) for(confidence in tail_levels) {
      a <- x[x$market==market & x$model=='GJR_R',]
      b <- x[x$market==market & x$model=='GJR_U',]
      assert_true(setequal(a$date,b$date),'GJR comparison dates differ.')
      loss <- paste0('fz0_',risk_suffix(confidence))
      d <- a[[loss]]-b[[loss]][match(a$date,b$date)]
      h <- hac_mean_test(d)
      out[[length(out)+1L]] <- data.frame(period,market,confidence,n=length(d),mean_difference=mean(d),
        hac_se=if(period=='main_oos') h[['se']] else NA_real_,
        p_value=if(period=='main_oos') h[['p_value']] else NA_real_)
    }
  }
  do.call(rbind,out)
}

check_reference <- function(actual,path,key,tolerance=1e-10) {
  ref <- read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
  make_key <- function(x) do.call(paste,c(x[key],sep='|'))
  assert_true(!anyDuplicated(make_key(actual)) && !anyDuplicated(make_key(ref)) &&
    setequal(make_key(actual),make_key(ref)) && setequal(names(actual),names(ref)),paste('Output keys differ:',basename(path)))
  actual <- actual[match(make_key(ref),make_key(actual)),names(ref),drop=FALSE]
  for(n in names(ref)) {
    if(is.numeric(ref[[n]])) assert_true(isTRUE(all.equal(actual[[n]],ref[[n]],tolerance=tolerance,check.attributes=FALSE)),paste('Output differs:',basename(path),n))
    else assert_true(identical(as.character(actual[[n]]),as.character(ref[[n]])),paste('Output differs:',basename(path),n))
  }
  invisible(TRUE)
}
