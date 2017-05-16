# Create Input

Create_T_Input <- function(years, T_mean, D_amp, Y_amp, slope) {
  hour = seq(1:(years * 360*24))
  D_pi = hour / length(hour) * pi * 2 * years * 360
  Y_pi = hour / length(hour) * pi * 2 * years
  D_out <- sin(D_pi) * (D_amp/2) + T_mean
  Y_out <- sin(Y_pi) * (Y_amp/2)
  out <- D_out + Y_out
  change <- approx(c(0,slope), n = length(out))[[2]]
  out <- out + change
}
