# -*- coding: utf-8 -*-
from __future__ import division
import time
import sympy as sy
from math import exp

"""
Created on Tue Dec 29 11:56:58 2015
Last modified: Dec 2016
@author: Fernando Moyano
Script for solving the steady state equations
"""

# Model options used in this script:
diff_fun  = "power" # Options: 'hama', 'cubic'
dec_fun   = "MM" # One of: 'MM', '2nd', '1st'
upt_fun   = "1st" # One of: 'MM', '2nd', '1st'

# Define varialbes
year = 31104000       # seconds in a year
month = 2592000       # seconds in a month
day = 86400           # seconds in a day
hour = 3600           # seconds in an hour
sec = 1               # seconds in a second!
tstep = hour

# Define functions

def T_resp_eq(k_ref, T, T_ref, E, R):
    return k_ref * sy.exp(-E/R * (1/T-1/T_ref))

# Define symbols
C_P, C_D, C_A, C_E, C_Em, C_M = \
    sy.symbols('C_P C_D C_A C_E C_Em C_M')


f_gr, f_mp, f_ue, r_ed, r_md, r_mr = \
(sy.symbols('f_gr f_mp f_ue r_ed r_md r_mr'))
V_D, K_D, V_U, K_U = (sy.symbols('V_D  K_D V_U K_U '))
D_d, D_e = (sy.symbols('D_d D_e'))
M, I_sl, I_ml, z = sy.symbols('M I_sl I_ml z')
Ka, k_ads, k_des, MD = sy.symbols('Ka k_ads k_des MD')

# Define fluxes

F_slcp = I_sl
F_mlcd = I_ml

D_diff = D_d * (C_D - 0)

if dec_fun == "MM":
      F_cpcd = z * V_D * C_P/z * C_E/z / (K_D + C_P/z)    
if dec_fun == "2nd":
      F_cpcd = z * V_D * C_P/z * C_E/z
if dec_fun == "1st":
      F_cpcd = V_D * C_P
      
if upt_fun == "MM":
      Ucd = z * V_U * D_diff/z * C_M/z / (K_U + C_D/z)
if upt_fun == "2nd":
      Ucd = z * V_U * D_diff/z * C_M/z
if upt_fun == "1st":
      Ucd = V_U * D_diff
      
# Microbial growth, mortality, respiration and enzyme production
F_cdcm = Ucd * f_gr * (1 - f_ue)
F_cdcr = Ucd * (1 - f_gr)
F_cdem = Ucd * f_gr * f_ue

F_cmcp = C_M * r_md * f_mp
F_cmcd = C_M * r_md * (1 - f_mp)
F_cmcr = C_M * r_mr

F_emce = D_e * (C_Em - C_E)
F_emcd = C_Em * r_ed
F_cecd = C_E * r_ed

dC_P = F_slcp + F_cmcp - F_cpcd
dC_D = F_mlcd + F_cpcd + F_cecd + F_emcd + F_cmcd - F_cdcm - F_cdcr - F_cdem
dC_E = F_emce - F_cecd
dC_Em = F_cdem - F_emce - F_emcd
dC_M = F_cdcm - F_cmcp - F_cmcr - F_cmcd

sol = sy.solve([dC_P, dC_D, dC_M, dC_E, dC_Em],
                        [C_P, C_D, C_M, C_E, C_Em], dict=True)

sol = sol[0]
E_m = 10 = sol[C_P]
sol_C_D = sol[C_D]
sol_C_M = sol[C_M]
sol_C_Em = sol[C_Em]
sol_C_E = sol[C_E]


#%%

# Site data
clay = 0.15
sand = 0.28
silt = 0.57
ps = 0.45
I_sl_v = 0.00005
I_ml_v = 0.000005
z_v = 0.3

# Intermediate parameter values
D_d0 = 1.37 / hour * tstep
D_e0 = 0.137 / hour * tstep
E_e = 10
E_K = 90
E_m = 10
E_r = 90
E_V = 90
K_D_ref = 60
K_U_ref = 1
k_ads_ref = 1.08e-6 / sec * tstep
k_des_ref = 1.19e-10 / sec * tstep
p1 = 2.75
p2 = 1.26
pd = 2700
psi_fc = 33
psi_Rth = 15000
R = 0.008314
r_ed_ref = 0.00017 / hour * tstep
r_md_ref = 0.002 / hour * tstep
r_mr_ref = 0.000036
T = 293.15
T_ref = 293.15
V_D_ref = 0.3 / hour * tstep
V_U_ref = 0.09 / hour * tstep

# End parameter values
f_gr_v = 0.7
f_ue_v = 0.025 / hour * tstep
M_v = 0.2

# Calculate intermediate variables
b = 2.91 + 15.9 * clay
k_ads = T_resp_eq(k_ads_ref, T, T_ref, E_V, R)
k_des = T_resp_eq(k_des_ref, T, T_ref, E_V, R)
psi_sat = exp(6.5 - 1.3 * sand) / 1000
Rth = ps * (psi_sat / psi_Rth)**(1 / b)
fc = ps * (psi_sat / psi_fc)**(1 / b)
if diff_fun == "hama":
    D_sm = (ps - Rth)**p1 * ((M_v - Rth)/(ps - Rth))**p2
if diff_fun == "power":
    D_sm = M_v**p1

# Calculate end variables
K_D_v = T_resp_eq(K_D_ref, T, T_ref, E_K, R)
V_D_v = T_resp_eq(V_D_ref, T, T_ref, E_V, R)
K_U_v = T_resp_eq(K_U_ref, T, T_ref, E_K, R)
V_U_v = T_resp_eq(V_U_ref, T, T_ref, E_V, R)
r_ed_v = T_resp_eq(r_ed_ref, T, T_ref, E_e, R)
r_mr_v = T_resp_eq(r_mr_ref, T, T_ref, E_r, R)
r_md_v = T_resp_eq(r_md_ref, T, T_ref, E_m, R)
D_d_v = D_d0 * D_sm
D_e_v = D_e0 * D_sm
MD_v = 200 * (100 * clay)**0.6 * pd * (1 - ps) / 1000000 #from mg kg-1 to kg m-3
M_fc_v = 1 #sy.Min(1, M / fc)
Ka_v = k_ads/k_des

# Substitute variables (parameters) with values

eq_C_P = sol_C_P.subs([
    (D_d, D_d_v), (D_e, D_e_v), (f_gr, f_gr_v), (f_ue, f_ue_v),
    (r_mr, r_mr_v),(I_ml, I_ml_v), (I_sl, I_sl_v), (K_D, K_D_v), (K_U, K_U_v), 
    (M, M_v), (r_ed, r_ed_v), (r_md, r_md_v),
    (V_D, V_D_v), (V_U, V_U_v), (z, z_v)
    ])

eq_C_D = sol_C_D.subs([
    (D_d, D_d_v), (D_e, D_e_v), (f_gr, f_gr_v), (f_ue, f_ue_v),
    (r_mr, r_mr_v),(I_ml, I_ml_v), (I_sl, I_sl_v), (K_D, K_D_v), (K_U, K_U_v), 
    (M, M_v), (r_ed, r_ed_v), (r_md, r_md_v),
    (V_D, V_D_v), (V_U, V_U_v), (z, z_v)
    ])

eq_C_E = sol_C_E.subs([
    (D_d, D_d_v), (D_e, D_e_v), (f_gr, f_gr_v), (f_ue, f_ue_v),
    (r_mr, r_mr_v),(I_ml, I_ml_v), (I_sl, I_sl_v), (K_D, K_D_v), (K_U, K_U_v), 
    (M, M_v), (r_ed, r_ed_v), (r_md, r_md_v),
    (V_D, V_D_v), (V_U, V_U_v), (z, z_v)
    ])

eq_C_Em = sol_C_Em.subs([
    (D_d, D_d_v), (D_e, D_e_v), (f_gr, f_gr_v), (f_ue, f_ue_v),
    (r_mr, r_mr_v),(I_ml, I_ml_v), (I_sl, I_sl_v), (K_D, K_D_v), (K_U, K_U_v), 
    (M, M_v), (r_ed, r_ed_v), (r_md, r_md_v),
    (V_D, V_D_v), (V_U, V_U_v), (z, z_v)
    ])

eq_C_M = sol_C_M.subs([
    (D_d, D_d_v), (D_e, D_e_v), (f_gr, f_gr_v), (f_ue, f_ue_v),
    (r_mr, r_mr_v),(I_ml, I_ml_v), (I_sl, I_sl_v), (K_D, K_D_v), (K_U, K_U_v), 
    (M, M_v), (r_ed, r_ed_v), (r_md, r_md_v),
    (V_D, V_D_v), (V_U, V_U_v), (z, z_v)
    ])

eq_C = eq_C_P + eq_C_M + eq_C_D + eq_C_E + eq_C_Em

# Calculate equilibrium value for adsorbed C
e_C_A = sy.Eq(Ka, C_A / (C_D * (MD * z - C_A)))   # Ka = LR / (L * R)
sol_C_A = sy.solve(e_C_A, C_A)[0]

eq_C_A = sol_C_A.subs([(Ka, Ka_v), (C_D, eq_C_D), (MD, MD_v), (z, z_v)])
eq_C2 = eq_C + eq_C_A


#%%
file = open("python_out_Em.txt", "a")

file.write("\n--------------\n" +
           "Time: " + time.strftime('%Y/%m/%d %H:%M:%S') + "\n\n" +
           "Options \n" + 
           "dec_fun: " + str(dec_fun) + " , upt_fun: " + str(upt_fun) +
           "\n\n" + "Solutions" + "\n\n" + 
           "C_P \n" + str(sol_C_P) + "\n\n" + 
           "C_D \n" + str(sol_C_D) + "\n\n" + 
           "C_M \n" + str(sol_C_M) + "\n\n" +
           "C_E \n" + str(sol_C_E) + "\n\n" + 
           "C_Em \n" + str(sol_C_Em) + "\n\n" + "\n\n")
file.close()
