sumfun <- function(df) {
  v  <- df$C_R_r
  m  <- mean(v)
  n  <- length(v)
  sd <- sd(v)
  cv <- m/sd
  se <- sd/sqrt(n)
  data.frame(C_R_r = m, C_R_sd = sd, C_R_cv = cv, C_R_se = se, 
             site = df$site[1], moist_grav = mean(df$moist_grav),
             temp = mean(df$temp), C13 = mean(df$C13), C13sd = sd(df$C13),
             moist_vol = mean(df$moist_vol))
}

out <- ddply(obs.accum, .(meas_rep), sumfun)
