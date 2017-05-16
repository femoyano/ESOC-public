#### initial_state.r

#### Documentation ============================================================
# Sets the initial values for the state variables
# Needs to be a named vector for using ode
#### ==========================================================================

initial_state <- c(
  C_P  = 100  , # particulate carbon
  C_D  = 0.01 , # soluble carbon
  C_E  = 0.01 , # enzymes carbon
  C_Em = 0.01 , # enzyme carbon at microbes
  C_M  = 1    , # microbial carbon
  C_Rg  = 0   , # growth respiration
  C_Rm  = 0    # maintenance respiration
)
