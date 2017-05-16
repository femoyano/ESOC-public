#### run_optim_smp.R

#### Documentations ===========================================================
# Script used to run optimization as shared memory job (SMP)
# author(s):
# Fernando Moyano (fmoyano #at# uni-goettingen.de)
#### ==========================================================================

t0 <- Sys.time()

### ----------------------------------- ###
###   Settings for parallel processing  ###
### ----------------------------------- ###
library(doParallel)
cores = detectCores()
# cores = 1
cat("Cores detected:", cores, "\n")
registerDoParallel(cores = cores)

### ----------------------------------- ###
###       User Setup                    ###
### ----------------------------------- ###
source('setup.R')

### ----------------------------------- ###
###         Run optimization            ###
### ----------------------------------- ###
list2env(setup, envir = .GlobalEnv)

source("main_optim.R")

print(Sys.time() - t0)
