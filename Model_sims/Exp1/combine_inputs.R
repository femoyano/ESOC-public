# Description:
# A = no leakage, equal EPa and litter
# B = no leakage, dynamic EPa and litter
# C = leakage, dynamic EPa and litter

name <- 'exp1_pars_bestfit_man_decMM_upt1st_power_years50'

# Load first df
all_out$run <- paste0('fixed - ', all_out$ypr, 'mm')
assign(name, all_out)
rm(all_out)

# Load second df
all_out$run <- paste0('dynamic - ', all_out$ypr, 'mm')
assign(name, rbind(get(name), all_out))
rm(all_out)

# Load third df
all_out$run <- paste0('dynamic and leaching - ', all_out$ypr, 'mm')
assign(name, rbind(get(name), all_out))
rm(all_out)

# save
save(list=c(name), file = paste0(name, '.RData'))
