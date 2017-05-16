#### Documentations ===========================================================
# Script used to prepare settings and run parameter optimization
# author(s):
# Fernando Moyano (fmoyano #at# uni-goettingen.de)
#### ==========================================================================

### ----------------------------------- ###
###      Non User Setup                 ###
### ----------------------------------- ###

### Libraries =================================================================

require(deSolve)
require(FME)
require(plyr)
require(reshape2)

starttime  <- format(Sys.time(), "%m%d-%H%M")

### Define time variables =====================================================
year     <- 31104000 # seconds in a year
hour     <- 3600     # seconds in an hour
sec      <- 1        # seconds in a second!

# ----- fixed model setup ----
ode.method <- "lsoda" # see ode function
flag_des   <- 1       # Must be 1: model unstable when doing stepwise.
tstep      <- hour    # Model time step
tsave      <- hour    # save time step (only relevant for stepwise model)

### Sourced required files ====================================================
source("flux_functions.R")
source("Model.R")
source("initial_state.R")
source(cost_fun)
source("AccumCalc.R")
source("ParsReplace.R")
source("SampleRun.R")

# Parameter setup =============================================================
pars_default <- read.csv(pars.default.file, row.names = 1)
pars_default <- setNames(pars_default[[1]], row.names(pars_default))

pars_optim       <- read.csv(pars.optim.file, row.names = 1)
pars_optim_init  <- setNames(pars_optim[[1]], row.names(pars_optim))
pars_optim_lower <- setNames(pars_optim[[2]], row.names(pars_optim))
pars_optim_upper <- setNames(pars_optim[[3]], row.names(pars_optim))

# Input Setup =================================================================
input_path    <- file.path(".")
input.all     <- read.csv(file.path(input_path, "mtdata_model_input.csv"))
obs.accum     <- read.csv(file.path(input_path, "mtdata_co2.csv"))
site.data.mz  <- read.csv(file.path(input_path, "site_Closeaux.csv"))
site.data.bf  <- read.csv(file.path(input_path, "site_BareFallow42p.csv"))

# Save text
savetxt2 <- paste0('_dec', dec_fun, '-upt', upt_fun, '-diff', diff_fun, '_')

### ----------------------------------- ###
###      Optimization/Calibration       ###
### ----------------------------------- ###

# Test model cost and computation time --------------
if(run.test) system.time(cost <- ModCost(pars_optim_init))

### Check sensitivity of parameters ---------------
if(run.sens) Sfun <- sensFun(ModCost, pars_optim_init)

## Optimize parameters
if (run.mfit) {
  fitMod <- modFit(f = ModCost, p = pars_optim_init, method = mf.method,
                   upper = pars_optim_upper, lower = pars_optim_lower)
}

## Saving work space
save.image(file = paste("Run_Optim_", starttime, savetxt, savetxt2, ".RData", sep = ""))
