# GetModelData.R
# script to add model to observed data

GetModelData <- function(pars) {
  
  source("SampleRun.R")

  all.out <- foreach(i = unique(input.all$treatment), .combine = 'rbind', 
                     .export = c(ls(envir = .GlobalEnv), "pars"),
                     .packages = c("deSolve")) %dopar% {
                       SampleRun(pars, input.all[input.all$treatment==i, ])
                     }
  
  return(all.out)
}
