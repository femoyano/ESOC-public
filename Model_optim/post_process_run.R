### ===================================== ###
### Post-process at end of run ###
### ===================================== ###

require(deSolve)
require(FME)
library(doParallel)
cores = detectCores()
# cores = 1
cat("Cores detected:", cores, "\n")
registerDoParallel(cores = cores)

# load('Run_Optim_0120-2357_decMM-upt1st-diffcubic_.RData')
## ------------------------------------ ##
##  Check parameter sensitivities     ----
## ------------------------------------ ##

# par_corr_plot <- pairs(Sfun, which = c("C_R"), col = c("blue", "green"))
# ident <- collin(Sfun)
# ident_plot <- plot(ident, ylim=c(0,20))
# ident[ident$N==9 & ident$collinearity<15,]


## ------------------------------------ ##
### -------- Get model data ----------- ##           
## ------------------------------------ ##

# pars_replace <- mcmcMod$bestpar
pars_replace <- fitMod$par
pars <- ParsReplace(pars_replace, pars_default)

source("GetModelData.R")

# Get model output with optimized parameters
system.time(mod.out <- GetModelData(pars))

source('post_process_T-fits_T-M-plots.R')
