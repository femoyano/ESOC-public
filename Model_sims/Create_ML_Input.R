# Create moisture and litter input

Create_ML_Input <- function(years, nypr, ypr, clay, sand, ps, temp) {

  # Potential evapotranspiration
  Gs <- 0 # ground heat flux
  cp <- 1.013 * 10^(-3) # (MJ kg-1) air specific heat capacity
  rho <- 1.225 # density of air (kg m-3) (average value)
  ra <- 100 # aerodynamic resistance (between crops and forest)
  gamma <- 0.067 # psychrometric constant
  lambda <- 2.3  # latent heat of vaporization MJ kg-1

  # Water potential calculations
  psi_ad  <- 100000 # kPa
  psi_fc  <- 33
  b       <- 2.91 + 15.9 * clay               # [] b parameter (Campbell 1974) as in Cosby  et al. 1984 - Alternatively: obtain from land model.
  # b <- 13.6 * clay + 3.5 # From Loll and Moldrup, Pore-size distribution and soil-water retention
  psi_sat <- exp(6.5 - 1.3 * sand) / 1000     # [kPa] saturation water potential (Cosby et al. 1984 after converting their data from cm H2O to Pa) - Alternatively: obtain from land model.
  Mad     <- ps * (psi_sat / psi_ad)^(1 / b)  # [m3 m-3] Water content at air dryness
  fc      <- ps * (psi_sat / psi_fc)^(1 / b)  # [m3 m-3] Field capacity relative water content (water retention formula from Campbell 1984) - Alternatively: obtain from land model.

  hour <- 3600
  volsoil <- 0.5 # soil volume m3
  Mdens   <- 1000 # water density kg m-3
  hours <- years * 360*24

  pr <- ypr * years
  npr <- nypr * years
  pr_freq <- round(hours/npr)
  pr_event <- pr / (volsoil * Mdens) / npr

  litt <- NA
  wp <- psi_fc
  leak <- 0    # soil water leakage
  M <- fc+pr_event

  makedata <- function(M) {
    for(i in 2:hours) {
      delta <- 4098 * (0.6108 * exp( 17.27 * temp[i] / (temp[i] + 237.3))) / (temp[i] + 237.3)^2  # (kPa ºC-1) calculated assuming Tair = Tsoil
      Rn <- 8.5 * temp[i] / 1000000 # (MJ m-2 s-1) very rough correlation using ST002 from Hainich 2016
      es <- 0.6108 * exp(17.27 * temp[i] / (temp[i] + 237.3)) # (kPa) assume Tair = Tsoil
      ea <- es*0.5 # (kPa) assume a constant RH of 50%
      Ep = (delta * (Rn - Gs) + rho * cp * (es - ea) / ra) / (lambda * (delta + gamma)) # (kg m−2 s−1)
      psi <- psi_sat * (M[i-1]/ps)^(-b)
      Ea <- Ep *  (log(psi) - log(psi_ad)) / (log(psi_fc) - log(psi_ad)) / (volsoil * Mdens) * hour
      M[i] <- M[i-1] - Ea
      if(i%%pr_freq == 0)  M[i] <- M[i] + pr_event
      leak[i] <- 0
      if(M[i] > fc) {
        leak[i] <- M[i] - fc
        M[i] <- fc
      }
      litt[i] <- Ea * 500 # litter in g m-2 h-1 scalled by Ea
      wp[i] <- psi
    }
    litt[1] <- litt[2]
    data.frame(moist = M, litt = litt, wp = wp, leak = leak)
  }

  input <- makedata(M)
  M <- min(min(input$moist) + pr_event, fc)
  input <- makedata(M)
}
