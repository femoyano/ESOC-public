# ParsReplace.R
# Adds or replaces parameters from the list of optimized parameters ----------------------

ParsReplace <- function(pars_optim, pars) {
  for(n in names(pars_optim)) pars[[n]] <- pars_optim[[n]]
  return(pars)
}
