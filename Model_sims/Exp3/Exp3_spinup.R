### Exp3: spinup
### Running only the base cases to get the stating point

### Model Options --------------------------------------------------------------
library(deSolve)
library(foreach)

a <- commandArgs(trailingOnly = TRUE)
D_0 <- as.numeric(a)
savetext <- paste0("_D", D_0)

mode <- c("fixed", "dynamic")

parsname  <- "pars_bestfit_man_decMM_upt1st_power"
diff_fun  <- "power" # Options: 'hama', 'power'
flag_lea  <- 0  # C_D leakage flow
sim_years <- 50
t_save    <- "month" # "hour", "day", "month", "year"

D_amp   <- 1
Y_amp   <- 0
T_slope <- 0

eq_moist <- c(0.08, 0.10, 0.12, 0.15, 0.20, 0.25)
eq_litt  <- 0.05 # g m-2 h-1
eq_temp  <- c(5,10,15,20,25,30)

pars_file  <- paste0("parsets/", parsname, ".csv")
site_file  <- "input/site_Closeaux_toc012.csv"

flag_dif  = 1  # diffusion limitation on or off
dec_fun   = "MM"   # One of: 'MM', '2nd', '1st'
upt_fun   = "1st"  # One of: 'MM', '2nd', '1st'

# Time constants
year       <- 31104000 # seconds in a year
month      <- 2592000  # seconds in a month
day        <- 86400    # seconds in a day
hour       <- 3600     # seconds in an hour
sec        <- 1        # seconds in a second!

# Read files
pars <- read.csv(pars_file, row.names = 1)
pars <- setNames(pars[[1]], row.names(pars))
pars[["D_0"]] <- D_0
source("SitePars.R")
pars <- SitePars(site_file)

source("Create_T_Input.R")


### Run code -------------------------------------------------------------------
RunFun <- function(m,i,j) {

  source("flux_functions.R", local = TRUE)

  if(m == "dynamic") litt <- eq_litt * j^2.75 / 0.25^2.75 else if(m == "fixed") litt = eq_litt

  source("GetEquil.R", local = TRUE)
  initial_state <- GetEquil(eq_temp = i, eq_moist = j, eq_litt = litt)

  source("Create_T_Input.R", local = TRUE)
  input <- data.frame(
    temp = Create_T_Input(years = sim_years, T_mean = i, D_amp = D_amp,
                         Y_amp = Y_amp, slope = T_slope),
    moist = j,
    litt = litt)


  temp <- input$temp  + 273.15  # [K] soil temperature
  moist <- input$moist
  I_sl <- input$litt / 1000
  I_ml <- input$litt / 10 / 1000 # metabolic intput 1/10 of structural

  times_input <- seq(1:length(temp))

  end   <- tail(times_input, 1)
  times <- seq(0, end, by = get(t_save)/hour)

  # ----- Define input interpolation functions
  Approx_I_sl  <<- approxfun(times_input, I_sl , method = "linear", rule = 2)
  Approx_I_ml  <<- approxfun(times_input, I_ml , method = "linear", rule = 2)
  Approx_temp  <<- approxfun(times_input, temp  , method = "linear", rule = 2)
  Approx_moist <<- approxfun(times_input, moist , method = "linear", rule = 2)

  source("Model.R", local = TRUE)
  ode.method <- "lsoda"  # see ode function
  out <- ode(initial_state, times, Model_desolve, pars, method = ode.method)
  out <- as.data.frame(out)
  out$mode <- m; out$tempC <- i; out$VWC <- j
  return(out)
}

### Setings for parallel processing  -------------------------------------------
library(doParallel)
cores = detectCores()
registerDoParallel(cores = cores)

all_out <- foreach(m = mode, .combine='rbind',
                   .export = ls(envir = .GlobalEnv),
                   .packages = c('deSolve')) %:%
  foreach(i = eq_temp, .combine='rbind') %:%
  foreach(j = eq_moist, .combine='rbind') %dopar% {
    RunFun(m,i,j)
  }

save(all_out, file = paste0("spin3_", t_save, "_", parsname, "_years", sim_years,
                            savetext, ".RData"))

