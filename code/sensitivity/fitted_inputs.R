# Fitted quantities used by the parameter-mixture calculations.
read_sensitivity_fit <- function(root, analysis, market, model, month) {
  directory <- file.path(root, 'inputs', 'sensitivity')
  select <- function(name) {
    x <- read.csv(file.path(directory, paste0(name, '.csv')), stringsAsFactors=FALSE)
    x[x$analysis == analysis & x$market == market & x$model == model & x$month == month, , drop=FALSE]
  }
  meta <- select('meta'); p <- select('parameters'); v <- select('covariance')
  assert_true(nrow(meta)==1L && meta$training_n==2000L && meta$K==5L,
              'Missing or inconsistent fitted-input metadata.')
  nm <- dax_sens_theta_order(model)
  assert_true(!anyDuplicated(p$parameter) && setequal(p$parameter,nm) && all(is.finite(p$estimate)),
              'Invalid parameter vector.')
  keys <- paste(v$row,v$column)
  grid <- expand.grid(row=nm,column=nm,stringsAsFactors=FALSE)
  assert_true(!anyDuplicated(keys) && setequal(keys,paste(grid$row,grid$column)) && all(is.finite(v$value)),
              'Invalid covariance entries.')
  covariance <- matrix(v$value[match(paste(grid$row,grid$column),keys)],length(nm),dimnames=list(nm,nm))
  assert_true(max(abs(covariance-t(covariance)))<1e-10,'Covariance is asymmetric.')
  list(market=market,model=model,period=if(analysis=='density') 'main' else 'holdout',
       refit_month=month,training_start=meta$training_start,training_end=meta$training_end,
       training_n=meta$training_n,K=meta$K,star_m=meta$star_m,
       point_estimate=setNames(p$estimate[match(nm,p$parameter)],nm),mle_covariance=covariance)
}

dax_sens_spec <- function() list(
  models=c('ARSV','STAR_SV'),
  months=format(seq(as.Date('2018-01-01'),as.Date('2023-12-01'),by='month'),'%Y-%m'),
  n_dates=1523L,estimation_window=2000L,K=5L,B=768L,blocks=12L,pairs_per_block=32L,
  nested_B=c(128L,256L,384L,512L,768L),P=5000L,crps_draws_per_block=5000L,
  base_seed=20260812L,resample_ess_fraction=0.5,particle_dead_fraction_limit=0.05)

dax_sens_validate_fit <- function(fit, model, refit_month) {
  assert_true(identical(fit$model,model) && identical(fit$refit_month,refit_month),
              'Fitted-input identity differs from the requested calculation.')
  theta <- fit$point_estimate[dax_sens_theta_order(model)]
  validate_sv_parameters(model,as.list(theta))
  J <- dax_sens_working_jacobian(model,theta)
  sigma <- J %*% fit$mle_covariance[names(theta),names(theta)] %*% t(J)
  list(theta=theta,z_hat=dax_sens_working_map(model,theta),L=t(chol(sigma,pivot=FALSE)))
}

hpu_fit_nodes <- function(fit, market, model, month, shocks=hpu_outer_shocks()) {
  assert_true(fit$market==market && fit$period=='holdout' && fit$K==5L &&
                (model!='STAR_SV' || fit$star_m==hpu_spec()$star_m[[market]]),
              'Fitted inputs do not match the holdout specification.')
  checked <- dax_sens_validate_fit(fit,model,month)
  assert_true(identical(dim(shocks),c(768L,8L)) && all(is.finite(shocks)), 'Invalid normal draw supply.')
  nodes <- lapply(seq_len(768L),function(b) {
    z <- checked$z_hat + as.numeric(checked$L %*% shocks[b,seq_along(checked$z_hat)])
    names(z) <- names(checked$z_hat)
    as.list(dax_sens_working_inverse(model,z))
  })
  list(point=as.list(checked$theta),nodes=nodes)
}
