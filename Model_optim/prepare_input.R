#### ---- Prepares input required for each run ----

# Choose site
if (input$site[1] == "bare_fallow") {
  site.data <- site.data.bf
  V_D_ref <- pars[["V_D_bf"]]
  f_CD    <- pars[["f_CD_bf"]]
  f_CE    <- pars[["f_CE_bf"]]
  f_CM    <- pars[["f_CM_bf"]]
} else if (input$site[1] == "maize") {
  site.data <- site.data.mz
  V_D_ref <- pars[["V_D_mz"]]
  f_CD    <- pars[["f_CD_mz"]]
  f_CE    <- pars[["f_CE_mz"]]
  f_CM    <- pars[["f_CM_mz"]]
} else stop("no site name match in prepare_input.R")

sand   <- site.data$sand  # [kg kg^-1] clay fraction values
clay   <- site.data$clay  # [kg kg^-1] sand fraction values
silt   <- site.data$silt  # [kg kg^-1] silt fraction values
ps     <- site.data$ps    # [m^3 m^-3] soil pore space
depth  <- site.data$depth # [m] soil depth
toc    <- site.data$toc   # [kgC kgSoil-1]

### ----- Calculate variables and add to parameter list
b       <- 2.91 + 15.9 * clay                         # [] b parameter (Campbell 1974) as in Cosby  et al. 1984 - Alternatively: obtain from land model.
psi_sat <- exp(6.5 - 1.3 * sand) / 1000               # [kPa] saturation water potential (Cosby et al. 1984 after converting their data from cm H2O to Pa) - Alternatively: obtain from land model.
Rth     <- ps * (psi_sat / pars[["psi_Rth"]])^(1 / b) # [m3 m-3] Threshold relative water content for mic. respiration (water retention formula from Campbell 1984)
fc      <- ps * (psi_sat / pars[["psi_fc"]])^(1 / b)  # [m3 m-3] Field capacity relative water content (water retention formula from Campbell 1984) - Alternatively: obtain from land model.
D_e0    <- pars[["D_d0"]] #/ 10

# Add new parameters to pars
parameters <- c(pars, V_D_ref = V_D_ref, sand = sand, silt = silt, clay = clay,
                ps = ps, depth = depth, Rth = Rth, fc = fc, b = b, D_e0 = D_e0)

### ----- Calculate initial C pool sizes

# Assign the pool sizes
TOC <- toc * parameters[["pd"]] * (1 - parameters[["ps"]]) * parameters[["depth"]]
initial_state[["C_P"]]  <- TOC * (1 - f_CD - f_CE * 2 - f_CM)
initial_state[["C_D"]]  <- TOC * f_CD
initial_state[["C_E"]]  <- TOC * f_CE
initial_state[["C_Em"]] <- TOC * f_CE
initial_state[["C_M"]]  <- TOC * f_CM
initial_state[["C_Rg"]] <- 0
initial_state[["C_Rm"]] <- 0


### ----- Extract data
times_input <- input$hour       # time of input data
temp        <- input$temp       # [K] soil temperature
moist       <- input$moist      # [m3 m-3] specific soil volumetric moisture
I_ml        <- input$litter_met / hour * tstep # [kgC m^-2 tstep^-1] convert litter input rates to the model rate
I_sl        <- input$litter_str / hour * tstep # [kgC m^-2 tstep^-1] convert litter input rates to the model rate

### ----- Obtain data times: start and end
start <- times_input[1]
end   <- tail(times_input, 1)
times <- seq(start, end)

### ----- Define input interpolation functions
Approx_I_sl  <<- approxfun(times_input, I_sl , method = "linear", rule = 2)
Approx_I_ml  <<- approxfun(times_input, I_ml , method = "linear", rule = 2)
Approx_temp  <<- approxfun(times_input, temp  , method = "linear", rule = 2)
Approx_moist <<- approxfun(times_input, moist , method = "linear", rule = 2)
