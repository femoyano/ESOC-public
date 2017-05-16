# Plot moisture effect at different precipitation regimes

### ----------------------------------- ###
###                Setup                ###
### ----------------------------------- ###

run_name = "sim"
t_save  = "month" # How often to save model results: either "hour", "day", "month", or "year"

# Get initial_state for variables: either defualt values (initial_state.R) or from a previous run
source("initial_state.R")
# load("initial_state.RData")

# Input file (can use make_input.R to create artificial input)
input_file <- "test1_input_T20_D1_Y20_P2000_N48.RData"

# Load the parameters used in the model
pars_file <- "parsets/pars_bestfit_man_decMM_upt1st_power.csv"

site_file <- "input/site_Closeaux_toc012.csv"

### === Model Options =====================

flag.dif  = 1  # diffusion limitation on or off
diff_fun  = "power" # Options: 'hama', 'power'
dec.fun   = "MM"   # One of: 'MM', '2nd', '1st'
upt.fun   = "1st"  # One of: 'MM', '2nd', '1st'

### ----------------------------------- ###
###       Non_user Setup                ###
### ----------------------------------- ###

# Libraries
require(deSolve)
require(plyr)
require(reshape2)

# Time constants
year       <- 31104000 # seconds in a year
month      <- 2592000  # seconds in a month
day        <- 86400    # seconds in a day
hour       <- 3600     # seconds in an hour
sec        <- 1        # seconds in a second!

# Sourced flux and Model functions
source("flux_functions.R")
source("Model.R")

# Loading pars
pars <- read.csv(pars_file, row.names = 1)
pars <- setNames(pars[[1]], row.names(pars))

# Loading input
load(input_file)

# Get data for making input ----------
# input_file <- "input/modelinput_1year_cycles.csv"
site.data  <- read.csv(site_file)
source("prepare_site_data.R")

# Prepare input -----------
temp <- input$temp
moist <- input$moist
litt <- input$litt

times_input <- seq(1:length(temp))

# Set equilibrium values --------
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

### ----------------------------------- ###
###              Run                    ###
### ----------------------------------- ###

ode.method <- "lsoda"  # see ode function
out <- ode(initial_state, times, Model_desolve, pars, method = ode.method)
out <- as.data.frame(out)

# ### ----------------------------------- ###
# ###        Analysis and Plotting        ###
# ### ----------------------------------- ###
# source('analysis_A.R')

# ### ----------------------------------- ###
# ###        Saving                       ###
# ### ----------------------------------- ###

# Save final state to use as inistial for another run
initial_state <- setNames(as.numeric(tail(out[2:7],1)), names(out[2:7]))
initial_state[5:6] <- 0
save(file = paste0(run_name, "_initial_state.RData"), initial_state)

# Save model output
save(file = paste0(run_name, "_out.RData"), out)

