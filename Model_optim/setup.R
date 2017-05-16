setup <- list(
  runinfo = "Description of this run",
  savetxt = "setup", # this apends to output file

  # -------- Model options ----------
  diff_fun  = "hama" ,  # Options: 'hama', 'cubic'
  dec_fun   = "MM" , # One of: 'MM', '2nd', '1st'
  upt_fun   = "1st" , # One of: 'MM', '2nd', '1st'

  # -------- Calibration options ----------
  run.test  = 1 ,  # run model cost once as test?
  run.sens  = 0 ,  # run FME sensitivity analysis?
  run.mfit  = 0 ,  # run modFit for optimization?
  # Observation error: name of column with error values:
  # 'C_R_gm', 'C_R_sdnorm', 'C_R_sd001', 'C_R_sd005', 'C_R_sd01', 'one' or NULL to use weight.
  SRerror  = 'C_R_sd01'  ,
  TRerror  = NULL  ,
  # Weight for cost:  only if error is NULL. One of 'none', 'mean', 'std'.
  SRweight = 'none' ,
  TRweight = 'none' ,
  # Scale variables? TRUE or FALSE
  scalevar = FALSE   ,
  # Choose method for modFit
  mf.method = "Nelder-Mead"     ,
  # Choose cost function
  cost_fun  = "ModCost.R" ,

  # -------- Parameter options ----------
  # csv file with default parameters
  pars.default.file = "parsets/pars_bestfit_man_decMM_upt1st_power.csv" ,
  # csv file with initial valeus for optimized parameters
  pars.optim.file   = "parsets/pars_optim_3.csv"
)
