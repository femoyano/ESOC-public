require("plyr")
require("reshape2")

## clean
rm(list = ls())

##############################################################################################
##### Load and initial process of data
##############################################################################################

mtdata1 <- read.table("data_checked.csv",header=T,sep=",",quote="\"")
# Round moisture to close values
mtdata1$moist_grav <- mtdata1$moist_grav / 100
mtdata1$moist_grav[mtdata1$moist_grav > 0.047 & mtdata1$moist_grav < 0.053] <- 0.05
mtdata1$moist_grav[mtdata1$moist_grav > 0.062 & mtdata1$moist_grav < 0.066] <- 0.065
mtdata1$moist_vol <- mtdata1$moist_grav * 1.8 # 1.8 is the dry soil bulk density
mtdata1$inc_end <- as.POSIXct(strptime(mtdata1$inc_end, format = "%d/%m/%Y %H:%M", tz = ""))
mtdata1$inc_start <- as.POSIXct(strptime(mtdata1$inc_start, format = "%d/%m/%Y %H:%M", tz = ""))
mtdata1$preinc_end <- as.POSIXct(strptime(mtdata1$preinc_end, format = "%d/%m/%Y %H:%M", tz = ""))
mtdata1$preinc_start <- as.POSIXct(strptime(mtdata1$preinc_start, format = "%d/%m/%Y %H:%M", tz = ""))

## ---------------------------------------------------------------------------------------
# End and start times of incubation periods should coincide.
# Differences of more than 1 day are assumed to be errors and are fixed in excel.
# Smaller differences are fixed here.

# First checking for times that differ by more than one day.
if(max(abs(mtdata1$preinc_end - mtdata1$inc_start)) > as.difftime(1, units = "days")) print("Need to correct inc_start and preinc_end dates for mistakes.")
checkdate <- function(df) {
  for (i in 1:(length(df[,1]) - 1)) {
    if (abs(df$inc_end[i] - df$preinc_start[i+1]) > as.difftime(1, units = "days")) print(df$sample[1])
  }
}
print("Samples for which inc_end and preinc_start dates need correction (if any):")
x <- by(mtdata1, mtdata1$sample, checkdate)
rm(x)

# If no more corrections needed, continue by making times equal.
mtdata1$preinc_end <- mtdata1$inc_start

# Calculate times for each incubation step
mtdata1$time_accum <- mtdata1$inc_end - mtdata1$inc_start
mtdata1$time_preinc <- mtdata1$preinc_end - mtdata1$preinc_start

# number the sequence of incubations of each sample
mtdata1 <- mtdata1[order(mtdata1$sample, mtdata1$preinc_start),]
mtdata1 <- ddply(mtdata1, .(sample), function(df) {df$tstage <- seq(1:length(df[,1])); return(df)})

### This entire section should be revised if used. And corrected C respired should come before.
## Calculate respiration rates (accum / time)
units(mtdata1$time_accum) <- "hours"

## ---------------------------------------------------------------------------------------
# Data quality flag
mtdata1$qflag <- logical(length = length(mtdata1[,1]))
# The next was a contaminated flask?
mtdata1$qflag[mtdata1$site == "maize" & mtdata1$temp == 20 & mtdata1$moist_grav == 0.065 & mtdata1$tstage == 4] <- TRUE
mtdata1$CO2[mtdata1$qflag==TRUE] <- NA

## ---------------------------------------------------------------------------------------
# Function to put data in long format
melt_date <- function(df) {
  df <- melt(df, measure.vars = c(5:8), variable.name = "step", na.rm = FALSE,
       value.name = "date", factorsAsStrings = TRUE)
  df$date <- as.POSIXct(df$date, origin = "1970-01-01 00:00.00 CET")
  return(df)
}
mtdata2 <- ddply(mtdata1, .(sample), melt_date)

mtdata2$istep <- rep(0, length(mtdata2[,1]))
mtdata2$istep[mtdata2$step=="preinc_start"] <- 1
mtdata2$istep[mtdata2$step=="preinc_end"] <- 2
mtdata2$istep[mtdata2$step=="inc_start"] <- 3
mtdata2$istep[mtdata2$step=="inc_end"] <- 4

# extract start of first preincubation for each sample
preinc_start <- by(mtdata2, mtdata2$sample, function(df) {min(df$date)}, simplify = TRUE)
preinc_start <- as.POSIXct(as.numeric(preinc_start), origin = "1970-01-01 00:00.00 CET")

# # Remove rows with preincubation times
# mtdata2 <- mtdata2[mtdata2$istep!=1 & mtdata2$istep!=2, ]

# # determine start of first incubation for each sample
# inc_start <- by(mtdata2, mtdata2$sample, function(df) {min(df$date)}, simplify = TRUE)
# inc_start <- as.numeric(inc_start)

# order data
mtdata2 <- mtdata2[order(mtdata2$sample, mtdata2$tstage, mtdata2$istep),]

## ---------------------------------------------------------------------------------------
## simplify, getting rid of unnecessary data ##
fun_normtime <- function(df, preinc_start) {
  s <- unique(df$sample)
  for(i in 1:length(s)) {
    df$sec[df$sample == s[i]] <- df$date[df$sample == s[i]] - preinc_start[i]
  }
  return(df)
}
mtdata2 <- fun_normtime(mtdata2, preinc_start)

mtdata2$hour <- round(mtdata2$sec / 3600)
mtdata2$day <- mtdata2$hour /  24
mtdata2$date <- NULL
mtdata2$sample_old <- NULL
mtdata2$CO2[mtdata2$istep!=4] <- NA
mtdata2$time_accum[mtdata2$istep!=4] <- NA
mtdata2$time_preinc[mtdata2$istep!=2] <- NA
mtdata2$treatment <- as.numeric(factor(interaction(mtdata2$moist_grav, mtdata2$site)))

##############################################################################################
# Select only points where measurements where made (for model data comparison)
##############################################################################################

mtdata3 <- mtdata2[mtdata2$step=="inc_end",]
mtdata3$time_preinc <- NULL
mtdata3$groupid <- NULL
mtdata3$istep <- NULL
mtdata3$step <- NULL
mtdata3$time_accum <- round(as.numeric(mtdata3$time_accum))

# Calculate gC respired per m-3 of soil
bd <- 1.8 # g cm-3
pd <- 2.6 # g cm-3
R <- 8314 # kPa * cm3 / (K * mol)
Cgmol <- 12 # grams of C per mol of CO2
dens_air <- 1.19 / 1000 # g cm-3 at 22.5C
mtdata3$air_vol <- mtdata3$flask_vol - (mtdata3$dry_soil / pd)
mtdata3$air_mol <- mtdata3$air_vol * 100 / (293 * R) # ideal gas law using cm3, kPa and K
mtdata3$C_R <- (mtdata3$CO2 / 1000000) * mtdata3$air_mol * 12 / (mtdata3$dry_soil / 1000) # gC respired per kg soil
mtdata3 <- mtdata3[!is.na(mtdata3$C_R),] # remove missing values

# Calculate means and variance for replicates in each group.
mtdata3$C_R_r <- mtdata3$C_R/mtdata3$time_accum

fun.summary <- function(df) {
  m  <- mean(df$C_R_r)
  n  <- length(df$C_R_r)
  sd <- sd(df$C_R_r)
  cv <- sd/m
  se <- sd/sqrt(n)
  data.frame(hour = round(mean(df$hour)), time_accum = round(mean(df$time_accum)),
             moist_vol = mean(df$moist_vol), temp = mean(df$temp),
             C_R_r = m, C_R_sd = sd, C_R_cv = cv, C_R_se = se, 
             C13 = mean(df$C13), C13sd = sd(df$C13)
             )
}
fun.groupmean <- function(df) {
  df$C_R_gm <- mean(df$C_R_r)
  return(df)
}

mtdata4 <- ddply(mtdata3, .(treatment, site, moist_grav, moist_vol, tstage, temp), fun.summary)
mtdata4 <- ddply(mtdata4, .(treatment), fun.groupmean)

mtdata4 <- mtdata4[order(mtdata4$site, mtdata4$moist_grav),]
rownames(mtdata4) <- NULL
for(s in c("bare_fallow", "maize")) {
  for(t in c(5,20,35)) {
    mtdata4$C_R_cv[mtdata4$site == s & mtdata4$moist_grav == 0.03 & mtdata4$temp == t] <- 
      max(mtdata4$C_R_cv[mtdata4$site == s & mtdata4$moist_grav < 0.07 & mtdata4$temp == t], na.rm = TRUE)
  }
}

mtdata4$C_R_sd[mtdata4$moist_grav == 0.03] <- mtdata4$C_R_r[mtdata4$moist_grav == 0.03] * mtdata4$C_R_cv[mtdata4$moist_grav == 0.03]
mtdata4$C_R_se[mtdata4$moist_grav == 0.03] <- mtdata4$C_R_sd[mtdata4$moist_grav == 0.03]
mtdata4$C_R_sdnorm <- mtdata4$C_R_sd / max(mtdata4$C_R_sd)
mtdata4$C_R_sd001 <- mtdata4$C_R_sdnorm + 0.01
mtdata4$C_R_sd005 <- mtdata4$C_R_sdnorm + 0.05
mtdata4$C_R_sd01 <- mtdata4$C_R_sdnorm + 0.1
mtdata4$one <- 1

##############################################################################################
##### Create dataframe for model input data
##############################################################################################
mtdata5 <- subset(mtdata2, select = c(treatment, site, sample, step, tstage, hour, temp, moist_vol))
mtdata5$moist <- mtdata5$moist_vol
mtdata5$moist_vol <- NULL
rownames(mtdata5) <- NULL

fun2.summary <- function(df) {
  data.frame(hour = round(mean(df$hour, na.rm = TRUE)))
}
mtdata6 <- ddply(mtdata5, .(treatment, site, moist, tstage, step, temp), fun2.summary)

# Define a function that will remove duplicated cases and create an hour of transition
# between temperature changes
cleanInput <- function(df) {
  for (i in 2:(nrow(df)-1)) {
    if (df$hour[i] == df$hour[i-1] & df$temp[i] == df$temp[i-1]) df$remove[i] <- TRUE
    if (df$temp[i] == df$temp[i-1] & df$temp[i] == df$temp[i+1]) df$remove[i] <- TRUE
  }
  df <- df[df$remove == FALSE,]
  for (i in 2:nrow(df)) {
    if (df$hour[i] == df$hour[i-1]) df$hour[i] <- df$hour[i] + 1
  }
  return(df)
}
mtdata6$remove <- FALSE
mtdata7 <- ddply(mtdata6, .(site, moist), cleanInput)
mtdata7$remove <- NULL

mtdata7$litter_met <- 0
mtdata7$litter_str <- 0
mtdata7$temp <- mtdata7$temp + 273  # convert to Kelvin

## Write out data
write.csv(mtdata3, file = "mtdata_co2_all.csv", row.names = FALSE)
write.csv(mtdata4, file = "mtdata_co2.csv", row.names = FALSE)
write.csv(mtdata7, file = "mtdata_model_input.csv", row.names = FALSE)

