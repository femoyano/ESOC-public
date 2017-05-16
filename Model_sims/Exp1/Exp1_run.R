### Exp1
### Running the different cases from the starting point

### Model Options --------------------------------------------------------------
library(deSolve)
library(foreach)

mode <- c("fixed", "dynamic")

parsname  <- "pars_bestfit_man_decMM_upt1st_power"
diff_fun  <- "power" # Options: 'hama', 'power'
flag_lea  <- 0  # C_D leakage flow
sim_years <- 50
t_save    <- "month" # "hour", "day", "month", "year"

D_amp   <- 1
Y_amp   <- 0
T_slope <- 0
nypr    <- c(3, 6, 12, 24, 45, 90, 180)
ypr     <- c(1000, 2000)

eq_litt  <- 0.05 # g m-2 h-1
eq_moist <- 0.2
eq_temp  <- c(5,15,25)

pars_file  <- paste0("parsets/", parsname, ".csv")
site_file  <- "input/site_Closeaux_toc012.csv"
start_file <- paste0("spin1_", t_save, "_", parsname, "_years", sim_years, ".RData")

flag_dif  <- 1  # diffusion limitation on or off
dec_fun   <- "MM"   # One of: 'MM', '2nd', '1st'
upt_fun   <- "1st"  # One of: 'MM', '2nd', '1st'

# Time constants
year       <- 31104000 # seconds in a year
month      <- 2592000  # seconds in a month
day        <- 86400    # seconds in a day
hour       <- 3600     # seconds in an hour
sec        <- 1        # seconds in a second!

# Read files
pars <- read.csv(pars_file, row.names = 1)
pars <- setNames(pars[[1]], row.names(pars))
source("SitePars.R")
pars <- SitePars(site_file)

load(start_file)
start <- all_out
rm(all_out)

source("Create_T_Input.R")
temp15 <- Create_T_Input(years = 1, T_mean = 15, D_amp = D_amp,
                         Y_amp = Y_amp, slope = T_slope)


### Run code -------------------------------------------------------------------
RunFun <- function(m,i,j,k) {

  source("flux_functions.R", local = TRUE)

  startvals <- tail(start[start$mode == m & start$tempC == i &
                            start$ypr == j, 2:9], 1)
  initial_state <- setNames(as.numeric(startvals), names(startvals))

  source("Create_T_Input.R", local = TRUE)
  temp <- Create_T_Input(years = 1, T_mean = i, D_amp = D_amp,
                         Y_amp = Y_amp, slope = T_slope)

  source("Create_ML_Input.R", local = TRUE)
  input <- Create_ML_Input(years = 1, nypr = k, ypr = j, clay = pars[['clay']],
                           sand = pars[['sand']], ps = pars[['ps']],
                           temp = temp)

  moist <- input$moist
  litt <- input$litt
  leak <- input$leak

  if(i != 15 & m == "fixed") {
    input <- Create_ML_Input(years = 1, nypr = k, ypr = j,
                             clay = pars[['clay']], sand = pars[['sand']],
                             ps = pars[['ps']], temp = temp15)
    moist <- input$moist
  }

  if(m == "fixed") litt <- rep(0.05, length(temp))

  temp <- rep(temp, sim_years)
  moist <- rep(moist, sim_years)
  litt <- rep(litt, sim_years)
  leak <- rep(leak, sim_years)

  temp <- temp  + 273.15  # [K] soil temperature
  I_sl <- litt / 1000
  I_ml <- litt/10 / 1000 # metabolic intput 1/10 of structural

  times_input <- seq(1:length(temp))

  end   <- tail(times_input, 1)
  times <- seq(0, end, by = get(t_save)/hour)

  # ----- Define input interpolation functions
  Approx_I_sl  <<- approxfun(times_input, I_sl , method = "linear", rule = 2)
  Approx_I_ml  <<- approxfun(times_input, I_ml , method = "linear", rule = 2)
  Approx_temp  <<- approxfun(times_input, temp  , method = "linear", rule = 2)
  Approx_moist <<- approxfun(times_input, moist , method = "linear", rule = 2)
  Approx_leak <<- approxfun(times_input, leak , method = "linear", rule = 2)

  source("Model.R", local = TRUE)
  ode.method <- "lsoda"  # see ode function
  out <- ode(initial_state, times, Model_desolve, pars, method = ode.method)
  out <- as.data.frame(out)
  out$mode <- m; out$tempC <- i; out$ypr <- j; out$nypr <- k
  return(out)
}

### Setings for parallel processing  -------------------------------------------
library(doParallel)registerDoMPI(cl)
cores = detectCores()
registerDoParallel(cores = cores)

all_out <- foreach(m = mode, .combine='rbind',
                   .export = ls(envir = .GlobalEnv),
                   .packages = c('deSolve')) %:%
  foreach(i = eq_temp, .combine='rbind') %:%
  foreach(j = ypr, .combine='rbind') %:%
  foreach(k = nypr, .combine='rbind') %dopar% {
    RunFun(m,i,j,k)
  }

save(all_out, file = paste0("exp1_", t_save, "_", parsname, "_years", sim_years, ".RData"))



