# Get site data

SitePars <- function(site_file) {

  site.data  <- read.csv(site_file)

  sand   <- site.data$sand  # [kgC kgSoil-1] clay fraction values
  clay   <- site.data$clay  # [kgC kgSoil-1] sand fraction values
  silt   <- site.data$silt  # [kgC kgSoil-1] silt fraction values
  ps     <- site.data$ps    # [m^3 m^-3] soil pore space
  depth  <- site.data$depth # [m] soil depth
  toc    <- site.data$toc   # [kgC kgSoil-1]

  ### ----- Calculate variables and add to parameter list
  psi_fc  <- 33
  b       <- 2.91 + 15.9 * clay                         # [] b parameter (Campbell 1974) as in Cosby  et al. 1984 - Alternatively: obtain from land model.
  psi_sat <- exp(6.5 - 1.3 * sand) / 1000               # [kPa] saturation water potential (Cosby et al. 1984 after converting their data from cm H2O to Pa) - Alternatively: obtain from land model.
  Rth     <- ps * (psi_sat / pars[["psi_Rth"]])^(1 / b) # [m3 m-3] Threshold relative water content for mic. respiration (water retention formula from Campbell 1984)
  fc      <- ps * (psi_sat / pars[["psi_fc"]])^(1 / b)  # [m3 m-3] Field capacity relative water content (water retention formula from Campbell 1984) - Alternatively: obtain from land model.
  psi_ad  <- 100000 # kPa air dry water potential
  Mad     <- ps * (psi_sat / psi_ad)^(1 / b)  # [m3 m-3] Water content at air dryness

  D_d0    <- pars[["D_0"]]        # [h-1] Diffusion conductance for dissolved C
  D_e0    <- pars[["D_0"]] / 10   # [h-1] Diffusion conductance for enzymes
  # mc      <- pars[['mc_0']] * pars[["pd"]] * (1 - ps) * depth # [kgC m-2] basal microbial carbon

  # Add new parameters to pars
  pars <- c(pars, sand = sand, silt = silt, clay = clay, ps = ps, depth = depth,
            z = depth, Rth = Rth, fc = fc, D_d0 = D_d0, D_e0 = D_e0, b = b)
}
