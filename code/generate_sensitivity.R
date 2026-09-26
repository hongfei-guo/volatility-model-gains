#!/usr/bin/env Rscript
# Generate one market-month of parameter-mixture forecasts from fitted quantities.
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
root <- normalizePath(file.path(dirname(script),'..'),mustWork=TRUE)
source(file.path(root,'code/evaluation_functions.R'))
for(f in c('numerical_helpers','predictive_scores','particle_filter','density_calculation','holdout_calculation','fitted_inputs')) {
  source(file.path(root,'code/sensitivity',paste0(f,'.R')))
}
args <- commandArgs(TRUE)
if(length(args)!=3L) stop('Usage: Rscript code/generate_sensitivity.R density|holdout MARKET YYYY-MM',call.=FALSE)
analysis <- args[1]; market <- args[2]; month <- args[3]
assert_true(analysis %in% c('density','holdout'),'Unknown analysis.')
meta <- read.csv(file.path(root,'inputs/sensitivity/meta.csv'),stringsAsFactors=FALSE)
selected <- meta[meta$analysis==analysis & meta$market==market & meta$month==month,]
assert_true(nrow(selected)==2L && setequal(selected$model,c('ARSV','STAR_SV')),'Market-month is outside the reported sample.')
data <- read.csv(file.path(root,'reproduced/returns',paste0(market,'.csv')),stringsAsFactors=FALSE)
data$date <- as.Date(data$date); hpu_validate_series(data)
period <- if(analysis=='density') 'main_oos' else 'holdout'
support <- read.csv(file.path(root,'results',period,'common_evaluation_dates.csv'),stringsAsFactors=FALSE)
support_dates <- sort(as.Date(support$date[support$market==market]))
fits <- setNames(lapply(c('ARSV','STAR_SV'),function(model) read_sensitivity_fit(root,analysis,market,model,month)),c('ARSV','STAR_SV'))
for(fit in fits) hpu_training(data,fit,month)
output <- file.path(root,'recomputed',analysis,paste0(market,'_',month,'.csv'))
assert_true(!file.exists(output),'Output already exists; preserve it before repeating the calculation.')
RNGkind('Mersenne-Twister','Inversion','Rejection')
if(analysis=='holdout') {
  result <- hpu_run_task(paste(market,month,sep='_'),data,fits,support_dates)
  result$FZ0 <- NULL; result$task_id <- NULL
} else {
  assert_true(market=='DAX' && month %in% dax_sens_spec()$months,'Invalid density comparison month.')
  result <- do.call(rbind,lapply(names(fits),function(model) {
    dax_sens_month_task(list(fit=fits[[model]],data=data,support_dates=support_dates),model,month,dax_sens_outer_shocks())$daily
  }))
  result$variant <- c(P='matched_plugin',I='parameter_integrated')[result$variant]
  result$scope <- ifelse(result$scope=='pooled','all_parameter_blocks',sub('^loo_','leave_block_out_',result$scope))
}
write_csv(result,output)
message('Predictions written to ',output)
