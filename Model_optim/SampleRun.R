### Defines the function that runs the model for each sample -------------------
SampleRun <- function(pars, input) {

  require(deSolve)
  source("prepare_input.R", local = TRUE)

  out <- ode(initial_state, times, Model_desolve, parameters,
             method = ode.method)

  out <- with(as.list(parameters), {
    # converting to gC respired per kg soil
    out[, "C_Rg"] <- out[, "C_Rg"]/(depth * (1 - ps) * pd) * 1000
    out[, "C_Rm"] <- out[, "C_Rm"]/(depth * (1 - ps) * pd) * 1000
    out[, "C_dec_r"] <- out[, "C_dec_r"]/(depth * (1 - ps) * pd) * 1000
    return(out)
  })

  out <- cbind(out, C_R = out[, "C_Rg"] + out[, "C_Rm"],
               C_dec = cumsum(out[, "C_dec_r"]),
               treatment = rep(input$treatment[1], nrow(out)))
  return(out)
}
