# Description:
# A = no leakage, equal EPa and litter
# B = no leakage, dynamic EPa and litter
# C = leakage, dynamic EPa and litter

name <- 'spin2_pars21_years50'

# Load first df
all_out$run <- paste0('A (', all_out$ypr, ')')
assign(name, all_out)
rm(all_out)

# Load second df
all_out$run <- paste0('B (', all_out$ypr, ')')
assign(name, rbind(get(name), all_out))
rm(all_out)

# Load third df
all_out$run <- paste0('C (', all_out$ypr, ')')
assign(name, rbind(get(name), all_out))
rm(all_out)

# save
save(name, file = paste0(name, '.RData'))
