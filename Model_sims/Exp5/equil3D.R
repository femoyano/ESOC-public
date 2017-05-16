# See the equilibrium with T and M and I

source("exp5/GetEquil2.R")
temp = seq(from = 5, by = 1, length.out = 30)
t_matrix <- matrix(temp, nrow = length(temp))
moist = seq(from = 0.1, by =  0.01, length.out = 20)
m_matrix <- matrix(moist, ncol = length(moist))
d <- expand.grid(temp = temp, moist = moist)
p <- 2.75
d$litt <- 0.05 * moist^p / 0.25^p
d$C <- GetEquil(eq_temp = d$temp, eq_moist = d$moist, eq_litt = d$litt,
                var = "C_P")

data <- list(
  temp = temp,
  moist = moist,
  C = matrix(d$C, nrow = length(temp), ncol = length(moist)),
  type = "surface")



plot_ly(x = data$temp, y = data$moist, z = data$C) %>% add_surface()
