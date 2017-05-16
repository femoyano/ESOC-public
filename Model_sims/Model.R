#### Model.R ================================================================

#### Documentation ==========================================================
# Main function running the model.
# This version is for use with deSolve ode function. It calculates and returns
# the state variable rates of change.
# author(s):
# Fernando Moyano (fmoyano #at# uni-goettingen.de)
#### ========================================================================

Model_desolve <- function(t, initial_state, pars) {
  # must be defined as: func <- function(t, y, parms,...) for use with ode

  with(as.list(c(initial_state, pars)), {

    # set time used for interpolating input data.
    t_i <- t
    # if (spinup) t_i <- t%%end  # this causes spinups to repeat the input data

    # Calculate the input and forcing at time t
    I_sl <- Approx_I_sl(t_i)
    I_ml <- Approx_I_ml(t_i)
    temp <- Approx_temp(t_i)
    moist <- Approx_moist(t_i)
    if(flag_lea) leak <- Approx_leak(t_i)

    # Calculate temporally changing variables
    K_D   <- TempRespEq(K_D_ref, temp, T_ref, E_K, R)
    if(upt_fun == "MM") K_U   <- TempRespEq(K_U_ref, temp, T_ref, E_K, R)
    V_D   <- TempRespEq(V_D_ref, temp, T_ref, E_V, R)
    V_U <- TempRespEq(V_U_ref, temp, T_ref, E_V, R)
    r_md  <- TempRespEq(r_md_ref, temp, T_ref, E_m , R)
    r_ed  <- TempRespEq(r_ed_ref, temp, T_ref, E_e , R)
    r_mr  <- TempRespEq(r_mr_ref, temp, T_ref, E_r , R)

    ## Diffusion calculations  --------------------------------------
    # Note: for diffusion fluxes, no need to divide by moist and depth to get specific
    # concentrations and multiply again for total since they cancel out.
    D_sm <- get_D_sm(moist, ps, Rth, b, p1, p2)
    D_d <- D_d0 * D_sm
    D_e <- D_e0 * D_sm

    # Diffusion calculations
    Diff_cd <- D_d * C_D

    ### Calculate all fluxes ------
    # Input rate
    F_slcp <- I_sl
    F_mlcd <- I_ml

    # Decomposition rate
    if(dec_fun == "MM") {
      F_cpcd <- ReactionMM(C_P, C_E, V_D, K_D, depth)
    }
    if(dec_fun == "2nd") {
      F_cpcd <- Reaction2nd(C_P, C_E, V_D, depth)
    }
    if(dec_fun == "1st") {
      F_cpcd <- Reaction1st(C_P, V_D)
    }

    # Calculate the uptake flux
    if(upt_fun == "MM") {
      Ucd <- ReactionMM(Diff_cd, C_M, V_U, K_U, depth)
    }
    if(upt_fun == "2nd") {
      Ucd <- Reaction2nd(Diff_cd, C_M, V_U, depth)
    }
    if(upt_fun == "1st") {
      Ucd <- Reaction1st(Diff_cd, V_U)
    }

    # Microbial growth, mortality, respiration and enzyme production
    F_cdcm <- Ucd * f_gr * (1 - f_ue)
    F_cdcr <- Ucd * (1 - f_gr)
    F_cdem <- Ucd * f_gr * f_ue

    F_cmcp <- C_M * r_md * f_mp
    F_cmcd <- C_M * r_md * (1 - f_mp)
    F_cmcr <- C_M * r_mr

    F_emce <- D_e * (C_Em - C_E)
    F_emcd <- C_Em * r_ed
    F_cecd <- C_E * r_ed

    # Leakage flow
    if(flag_lea) F_cdcl <- C_D * (leak/moist) else F_cdcl <- 0

    ## Rate of change calculation for state variables ---------------
    dC_P <- F_slcp + F_cmcp - F_cpcd
    dC_D <- F_mlcd + F_cpcd + F_cecd + F_emcd + F_cmcd - F_cdcm - F_cdcr - F_cdem - F_cdcl
    dC_E <- F_emce - F_cecd
    dC_Em <- F_cdem - F_emce - F_emcd
    dC_M <- F_cdcm - F_cmcp - F_cmcr - F_cmcd
    dC_Rg <- F_cdcr
    dC_Rm <- F_cmcr
    dC_L <- F_cdcl

    return(list(c(dC_P, dC_D, dC_E, dC_Em, dC_M, dC_Rg, dC_Rm, dC_L),
      C_dec_r = F_cpcd, temp = temp, moist = moist, D_d = D_d))

  })  # end of with(...

}  # end of Model
