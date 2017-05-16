### Equilibrium solutions
GetEquil <- function(eq_temp, eq_moist, eq_litt, var) {
  initial_state <- with(as.list(pars), {

    source("flux_functions.R")
    # Site data
    moist = eq_moist
    temp = eq_temp + 273.15
    I_sl = eq_litt / 1000
    I_ml = eq_litt/10 / 1000

    # Calculate intermediate variables
    b = 2.91 + 15.9 * clay
    # k_ads = Temp.Resp.Eq(k_ads_ref, temp, T_ref, E_V, R) k_des
    # = Temp.Resp.Eq(k_des_ref, temp, T_ref, E_V, R)
    psi_sat = exp(6.5 - 1.3 * sand)/1000
    Rth = ps * (psi_sat/psi_Rth)^(1/b)
    fc = ps * (psi_sat/psi_fc)^(1/b)
    if(diff_fun == "hama")  D_sm = (ps - Rth)^p1 * ((moist - Rth)/(ps - Rth))^p2
    if(diff_fun == "cubic") D_sm = moist^3
    if(diff_fun == "power") D_sm = moist^p1

    # Calculate end variables
    K_D = TempRespEq(K_D_ref, temp, T_ref, E_K, R)
    V_D = TempRespEq(V_D_ref, temp, T_ref, E_V, R)
    V_U = TempRespEq(V_U_ref, temp, T_ref, E_V, R)
    r_md = TempRespEq(r_md_ref, temp, T_ref, E_m, R)
    r_ed = TempRespEq(r_ed_ref, temp, T_ref, E_e, R)
    r_mr = TempRespEq(r_mr_ref, temp, T_ref, E_r, R)
    D_d = D_0*D_sm
    D_e = D_0 / 10*D_sm
    MD = 200*(100*clay)^0.6*pd*(1 - ps) / 1000000 #from mg kg-1 to kg m-3
    M_fc = 1 #sy.Min(1, M / fc)
    # Ka = k_ads/k_des


    # --------------------- Options: dec_fun = 'MM' upt_fun = '1st'

    C_P <- -K_D*r_ed*z*(2*D_e + r_ed)*(I_ml*f_gr*f_mp*f_ue*r_md -
                                         I_ml*f_gr*f_mp*r_md + I_sl*f_gr*f_mp*f_ue*r_md - I_sl*f_gr*f_mp*r_md +
                                         I_sl*f_gr*f_ue*r_mr + I_sl*f_gr*r_md - I_sl*r_md - I_sl*r_mr)/
      (D_e*I_ml*V_D*f_gr*f_ue*r_md + D_e*I_ml*V_D*f_gr*f_ue*r_mr +
         2*D_e*I_ml*f_gr*f_mp*f_ue*r_ed*r_md - 2*D_e*I_ml*f_gr*f_mp*r_ed*r_md +
         D_e*I_sl*V_D*f_gr*f_ue*r_md + D_e*I_sl*V_D*f_gr*f_ue*r_mr +
         2*D_e*I_sl*f_gr*f_mp*f_ue*r_ed*r_md - 2*D_e*I_sl*f_gr*f_mp*r_ed*r_md +
         2*D_e*I_sl*f_gr*f_ue*r_ed*r_mr + 2*D_e*I_sl*f_gr*r_ed*r_md -
         2*D_e*I_sl*r_ed*r_md - 2*D_e*I_sl*r_ed*r_mr +
         I_ml*f_gr*f_mp*f_ue*r_ed**2*r_md - I_ml*f_gr*f_mp*r_ed**2*r_md +
         I_sl*f_gr*f_mp*f_ue*r_ed**2*r_md - I_sl*f_gr*f_mp*r_ed**2*r_md +
         I_sl*f_gr*f_ue*r_ed**2*r_mr + I_sl*f_gr*r_ed**2*r_md -
         I_sl*r_ed**2*r_md - I_sl*r_ed**2*r_mr)

    C_D <- -(I_ml + I_sl)*(r_md + r_mr)/
      (D_d*V_U*(f_gr*f_ue*r_mr + f_gr*r_md - r_md - r_mr))

    C_M <- f_gr*(I_ml + I_sl)*(f_ue - 1)/
      (f_gr*f_ue*r_mr + f_gr*r_md - r_md - r_mr)

    C_E <- -D_e*f_gr*f_ue*(I_ml + I_sl)*(r_md + r_mr)/
      (r_ed*(2*D_e + r_ed)*(f_gr*f_ue*r_mr + f_gr*r_md - r_md - r_mr))

    C_Em <- -f_gr*f_ue*(D_e + r_ed)*(I_ml + I_sl)*(r_md + r_mr)/
      (r_ed*(2*D_e + r_ed)*(f_gr*f_ue*r_mr + f_gr*r_md - r_md - r_mr))

    C_T <- C_P + C_D + C_M + C_E

    return(get(var))
  })
}


