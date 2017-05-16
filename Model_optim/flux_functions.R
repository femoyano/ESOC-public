# flux_functions.r

#### Documentation
#### ============================================================
#### Functions called by the model
#### Fernando Moyano (fmoyano #at# uni-goettingen.de)
#### ==========================================================================

## Functions to calculate diffusion depending on options -----
if (diff_fun == "hama") {
  get_D_sm <- function(moist, ps, Rth, b, p1, p2) {
    if (moist <= Rth) {D_sm <- 0} else
     D_sm <- (ps - Rth)^p1 * ((moist - Rth) / (ps - Rth))^p2 # p1=1.5, p2=2.5
  }
}
if (diff_fun == "archie") {
  get_D_sm <- function(moist, ps, Rth, b, p1, p2) D_sm <- ps^p1 * (moist/ps)^p2  # p1=1.3-2.25, p2=2-2.5
}
if (diff_fun == "power") {
  get_D_sm <- function(moist, ps, Rth, b, p1, p2) D_sm <- moist^p1  # p1=3
}
if (diff_fun == "MQ") {
  get_D_sm <- function(moist, ps, Rth, b, p1, p2) D_sm <- moist^p1 / ps^p2 # p1=3.33, p2=2
}
if (diff_fun == "PC") {
  get_D_sm <- function(moist, ps, Rth, b, p1, p2) D_sm <- p1*moist^p2  # p1=2.8, p2=3
}
if (diff_fun == "sade") {
  get_D_sm <- function(moist, ps, Rth, b, p1, p2) D_sm <- p1*(moist/ps)^p2  # p1=0.73, p2=1.98
}
if (diff_fun == "oles1") {
  get_D_sm <- function(moist, ps, Rth, b, p1, p2) {
    if (moist < 0.022*b) {D_sm <- 0} else
    D_sm <- p1 * moist * (moist - 0.022 * b) / (ps - 0.022 * b) # p1=0.45
  }
}
if (diff_fun == "oles2") {
  get_D_sm <- function(moist, ps, Rth, b, p1, p2) {
    D_sm <- p1 * moist * (moist / ps)^(p2*b) # p1=0.45, p2=0.3
  }
}

## Reaction kinetics ---------
ReactionMM <- function(S, E, V, K, depth) {
  S <- S/depth
  E <- E/depth
  Flux <- (V * E * S)/(K + S) * depth
}

Reaction2nd <- function(S, E, V, depth) {
  S <- S/depth
  E <- E/depth
  Flux <- (V * S * E) * depth
}

Reaction1st <- function(S, V, fc.mod) {
  Flux <- V * S
}

# ==============================================================================
# Temperature responses

# Temperature response for equilibrium reactions = Arrhenius
# (for K values)
TempRespEq <- function(k_ref, T, T_ref, E, R) {
  k_ref * exp(-E/R * (1/T - 1/T_ref))
}
