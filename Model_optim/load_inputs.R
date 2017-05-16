### Climate forcing and litter input =========================================

### Documentation ============================================================
# Forcing data will be interpolated to model time step
# Soil T must be in kelvin and soil moisture in fraction of volume
# Litter input must be in gC m-2 h-1

### Spatial soil data =========================================================

site.data  <- read.csv(site.file)

sand   <- site.data$sand  # [g g^-1] clay fraction values
clay   <- site.data$clay  # [g g^-1] sand fraction values
silt   <- site.data$silt  # [g g^-1] silt fraction values
ps     <- site.data$ps    # [m^3 m^-3] soil pore space
depth  <- site.data$depth # [m] soil depth

### Climate and litter data ===================================================

input.data <- read.csv(input.file) # input data file

# Adjust time units and extract data
times.input <- input.data[,1]        # time vector of input data
input.tstep <- get(names(input.data)[1])
times_input <- times.input * input.tstep / tstep # convert input data to model step unit
temp        <- input.data$temp       # [K] soil temperature
moist       <- input.data$moist      # [m3 m-3] specific soil volumetric moisture

# Note! Model units changed to kg, but input still in g m-2, so converted here to kg.
I_ml  <- input.data$litter_met / hour * tstep  / 1000 # [kgC m^-2 tstep^-1] convert litter input rates to the model rate
I_sl  <- input.data$litter_str / hour * tstep  / 1000 # [kgC m^-2 tstep^-1] convert litter input rates to the model rate

rm(input.data, site.data, times.input, input.tstep)
