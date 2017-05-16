library(plyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

source("AccumCalc.R")

################################################################################
## Data processing
################################################################################

# Get accumulated values to match observations and merge datasets
data.accum <- merge(obs.accum, AccumCalc(mod.out, obs.accum),
                    by.x = c("treatment", "hour"), by.y = c("treatment", "time"))
# get hourly rates [gC kg-1 h-1]
data.accum$C_R_rm <- data.accum$C_R_m / data.accum$time_accum
data.accum$C_Rm_r <- data.accum$C_Rm / data.accum$time_accum
data.accum$C_Rg_r <- data.accum$C_Rg / data.accum$time_accum
data.accum$C_dec_r <- data.accum$C_dec / data.accum$time_accum
data.accum$C_R_ro <- data.accum$C_R_r # Observed data already gC kg-1 h-1
data.accum$C_R <- NULL

dta <- data.accum
# Convert to mg kg-1 h-1
dta$C_R_ro <- data.accum$C_R_ro * 1000
dta$C_R_rm <- data.accum$C_R_rm * 1000
dta$C_Rm_r <- data.accum$C_Rm_r * 1000
dta$C_Rg_r <- data.accum$C_Rg_r * 1000
dta$C_dec_r <- data.accum$C_dec_r * 1000

# Rename and reorder factor levels
dta$soil <- revalue(data.accum$site, c("bare_fallow" = "bare fallow", "maize" = "cropped"))
dta$soil <- factor(dta$soil, levels = c("cropped", "bare fallow"))
dta$temp_group <- data.accum$temp.group

# create a group variable
dta$temp_group <- interaction(dta$soil, dta$temp)

# create a group variable
dta$moist_group <- interaction(dta$soil, dta$moist_vol)


################################################################################
## Temp fits
################################################################################

# Temperature response models -----------
fitExpfun <- function(df, Q10) {
  nls(x ~ a * Q10^(temp/10), start=list(a=0.01, Q10=Q10),
      data = df, na.action = na.exclude)
}

fitArrfun <- function(df, Ea) {
  nls(x ~ a * exp(-Ea/0.008314*(1/(temp+273) - 1/273)),
      start=c(a = 0.1, Ea = Ea), data = df)
}

# # Code for linear fits
# linfitExp <- lm(log(x) ~ temp, data = df, na.action = na.exclude)
# Q10 = exp(10*coef(fitExp)[[2]]); R0 = coef(fitExp)[[1]]
# linfitArr <- lm(log(x) ~ I(1/(temp+273)), data = df, na.action = na.exclude)
# A = coef(fitArr)[[1]]; Ea = coef(fitArr)[[2]] * 0.008314 * (-1)

# Fit a temprature function to each moisture subgroup and plot
FitTemp <- function(df, var, set, trange) {
  texc <- -99; if(trange=='5-20') texc <- 35; if(trange=='20-35') texc <- 5
  df <- df[df$temp!=texc,]
  df$x <- df[[var]]
  fitExp <- try(fitExpfun(df, 5))
  if(class(fitExp)=="try-error") fitExp <- try(fitExpfun(df, 1))
  if(class(fitExp)!="try-error") {
    Q10 = coef(fitExp)[[2]]; R0 = coef(fitExp)[[1]]} else {Q10 = NA; R0 = NA
    }
  fitArr <- try(fitArrfun(df, 120))
  if(class(fitArr)=="try-error") fitArr <- try(fitArrfun(df, 30))
  if(class(fitArr)!="try-error") {

    A = coef(fitArr)[[1]]; Ea = coef(fitArr)[[2]]} else {A = NA; Ea = NA
    }
  out <- data.frame(set = set, trange = trange, soil = df$soil[1], var = var,
                    moist_vol = df$moist_vol[1], Q10 = Q10, R0 = R0, A = A, Ea = Ea)
}

# Fits for respired or decomposed C with temperature for different temperature ranges
fit.temp <- ddply(dta, .(moist_group), .fun = FitTemp, var = "C_R_ro",
                  set = 'obs', trange = '5-35')
fit.temp <- rbind(fit.temp, ddply(dta, .(moist_group), .fun = FitTemp,
                                  var = "C_R_rm", set = 'mod', trange = '5-35'))
fit.temp <- rbind(fit.temp, ddply(dta, .(moist_group), .fun = FitTemp,
                                  var = "C_R_ro", set = 'obs', trange = '5-20'))
fit.temp <- rbind(fit.temp, ddply(dta, .(moist_group), .fun = FitTemp,
                                  var = "C_R_rm", set = 'mod', trange = '5-20'))
fit.temp <- rbind(fit.temp, ddply(dta, .(moist_group), .fun = FitTemp,
                                  var = "C_R_ro", set = 'obs', trange = '20-35'))
fit.temp <- rbind(fit.temp, ddply(dta, .(moist_group), .fun = FitTemp,
                                  var = "C_R_rm", set = 'mod', trange = '20-35'))
fit.temp <- rbind(fit.temp, ddply(dta, .(moist_group), .fun = FitTemp,
                                  var = "C_dec_r", set = 'mod', trange = '5-35'))
fit.temp <- rbind(fit.temp, ddply(dta, .(moist_group), .fun = FitTemp,
                                  var = "C_dec_r", set = 'mod', trange = '5-20'))
fit.temp <- rbind(fit.temp, ddply(dta, .(moist_group), .fun = FitTemp,
                                  var = "C_dec_r", set = 'mod', trange = '20-35'))

################################################################################
## Statistics
################################################################################

# Get RMSE and MAE
res <- dta$C_R_ro - dta$C_R_m
RMSE <- sqrt(mean(res^2))
MAE <- mean(abs(res))
lmmod <- lm(dta$C_R_ro~dta$C_R_rm)
print(RMSE)
print(summary(lmmod))

################################################################################
## Plots
################################################################################

# ------------------------------------------------------------------------------
# Plotting setup  --------------------------------------------------------------
# ------------------------------------------------------------------------------

prefix <- "plot_"
savedir <- file.path("plots")
devname <- "png"
devfun <- png
export <- 1
opar <- par(no.readonly=TRUE)

darkcols <- brewer.pal(8, "Dark2")  # heat.colors(10, alpha = 0.5)
darkcols <- adjustcolor(darkcols, alpha.f = 0.6)
palette(darkcols)
palette(grey.colors(2))

dplot1 <- subset(dta, select = c(moist_vol, C_Rg_r, C_Rm_r, C_dec_r, site, soil,
                                 temp, temp_group))
dplot1 <- melt(data = dplot1, measure.vars = c('C_Rg_r', 'C_Rm_r', 'C_dec_r'))
dplot1$flux <- revalue(dplot1$variable, c('C_Rg_r'='growth',
                                          'C_Rm_r'= 'maintenance',
                                          'C_dec_r'='decomposition'))

dta$temp_group <- factor(dta$temp_group, levels = c("cropped.5", "cropped.20",
                                                    "cropped.35", "bare fallow.5",
                                                    "bare fallow.20", "bare fallow.35"))
dplot2 <- subset(dta, select = c(moist_vol, C_R_ro, C_R_rm, site, soil,
                                 temp, temp_group))
dplot2 <- melt(data = dplot2, measure.vars = c('C_R_ro', 'C_R_rm'))
dplot2$flux <- revalue(dplot2$variable, c('C_R_ro'='measured',
                                          'C_R_rm'= 'modelled'))
col1 <- brewer.pal(3, "Set2")
col2 <- c("#333333B3", "#B3B3B3E6")
col3 <- c("#104E8B77", "#FF8C0077")

# ------------------------------------------------------------------------------
# Plot accumulated model vs data -----------------------------------------------
# ------------------------------------------------------------------------------

plotname <- paste0(prefix, "accum_mod_obs.", devname)
plotfile <- file.path(savedir, plotname)
if(export) devfun(file = plotfile, width = 500, height = 500)
p <- ggplot(data = dta, aes(x = C_R_rm, y = C_R_ro, colour = soil)) +
  xlab(expression(paste("Modeled Accumulated ", CO[2], " (g C)"))) +
  ylab(expression(paste("Measured Accumulated  ", CO[2], " (g C)"))) +
  ylim(c(0,0.6))+
  xlim(c(0,0.6))+
  theme_bw(base_size = 12) +
  geom_point(size = 1) +
  geom_abline(intercept=0, slope=1) +
  scale_color_manual(values = darkcols) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(0.17,0.85))

print(p)

# ------------------------------------------------------------------------------
# Plot mod and obs flux vs moisture values for each temp group -----------------
# ------------------------------------------------------------------------------

plotname <- paste0(prefix, "moist-resp.", devname)
plotfile <- file.path(savedir, plotname)
if(export) devfun(file = plotfile, width = 700, height = 1000) #, width = 5, height = 5)

p <- ggplot(data = dplot2, aes(x = moist_vol, y = value, colour = flux,
                               shape = flux)) +
  geom_smooth(se = FALSE, size = 1.5) +
  geom_point(size = 1.3) +
  scale_color_manual(values = col3) +
  xlab(expression(paste("Soil moisture (", m^3*m^-3,")"))) +
  ylab(expression(paste("Respired ", CO[2], " (",mg~C~kg^-1*soil~h^-1, ")"))) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw(base_size = 12) +
  # scale_color_manual(values = darkcols) +
  theme(axis.line = element_line(colour = "black"),
        legend.position = c(0.68,0.93),
        legend.title = element_blank(),
        legend.key.height=unit(0.5,'cm'), legend.key.width=unit(0.8,'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  ) +
  # facet_wrap("temp_group", scales = "free_y")
  facet_grid(temp~soil, scales = "free")

print(p)

# ------------------------------------------------------------------------------
# Plot modelled fluxes vs moisture values for each temp group ------------------
# ------------------------------------------------------------------------------

plotname <- paste0(prefix, "moist-modflux.", devname)
plotfile <- file.path(savedir, plotname)
if(export) devfun(file = plotfile, width = 1000, height = 700) #, width = 5, height = 5)

p <- ggplot(data = dplot1, aes(x = moist_vol, y = value, group = flux)) +
  geom_smooth(aes(colour = flux))+#, se = FALSE, size = 3) +
  geom_point(aes(colour = flux), size = 1.5) +
  scale_color_manual(values = col1) +
  xlab(expression(paste("Soil moisture (", m^3*m^-3,")"))) +
  ylab(expression(paste("Carbon flux (",mg~C~kg^-1*soil~h^-1, ")"))) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw(base_size = 15) +
  # scale_color_manual(values = darkcols) +
  theme(axis.line = element_line(colour = "black"),
        legend.position = c(0.16,0.9),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        ) +
  # facet_wrap("temp_group", scales = "free")
  facet_grid(soil~temp, scales = "free")

print(p)

# ------------------------------------------------------------------------------
# Plot mod and obs flux vs temperature for each moisture group -----------------
# ------------------------------------------------------------------------------

plotname <- paste0(prefix, "temp-resp.", devname)
plotfile <- file.path(savedir, plotname)
if(export) devfun(file = plotfile, width = 1000, height = 700)
# dta$moist_bin <- dta$moist_bin <- cut(dta$moist_vol,
#   breaks = c(0,0.05,0.1,0.2,0.3,0.4,0.5),
#   labels = c("VWC 0-5 %", "VWC 5-10 %", "VWC 10-20 %",
#              "VWC 20-30 %", "VWC 30-40 %", "VWC 40-50 %"))
dplot2$moist_bin <- cut(dplot2$moist_vol,
                       breaks = c(0,0.1,0.3,0.5),
                       labels = c("VWC 0-10 %", "VWC 10-30 %", "VWC 30-45 %"))
p <- ggplot(data = dplot2, aes(x = temp, y = value, color = flux, shape = flux)) +
  stat_smooth(se = FALSE, size = 2, span = 1.5) +
  geom_point(size = 1.5) +
  xlab(expression(paste("Soil temperature (", degree*C,")"))) +
  ylab(expression(paste("Respired ", CO[2], " (",mg~C~kg^-1*soil~h^-1, ")"))) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = col3) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.position=c(0.12,0.88)) +
  # facet_wrap("moist_bin", scales = "free_y")
facet_grid(soil~moist_bin, scales = "free")

print(p)

# ------------------------------------------------------------------------------
# Plot modelled fluxes vs temperature for each moisture group ------------------
# ------------------------------------------------------------------------------

plotname <- paste0(prefix, "temp-modeflux.", devname)
plotfile <- file.path(savedir, plotname)
if(export) devfun(file = plotfile, width = 1000, height = 700)
dplot1$moist_bin <- cut(dplot1$moist_vol,
                       breaks = c(0,0.1,0.3,0.5),
                       labels = c("VWC 0-10 %", "VWC 10-30 %", "VWC 30-45 %"))
col1 <- brewer.pal(3, "Set2")

p <- ggplot(data = dplot1, aes(x = temp, y = value, group = flux)) +
  geom_smooth(aes(colour = flux))+#, se = FALSE, size = 3) +
  geom_point(aes(colour = flux), size = 1.5) +
  scale_color_manual(values = col1) +
  xlab(expression(paste("Temperature (", degree,"C)"))) +
  ylab(expression(paste("Carbon flux (",mg~C~kg^-1*soil~h^-1, ")"))) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw(base_size = 15) +
  # scale_color_manual(values = darkcols) +
  theme(axis.line = element_line(colour = "black"),
        legend.position = c(0.12,0.9),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  ) +
  # facet_wrap("moist_bin", scales = "free_y")
facet_grid(soil~moist_bin, scales = "free")

print(p)


# ------------------------------------------------------------------------------
## Plot temperature sensitivities ---------------------
# ------------------------------------------------------------------------------

PlotTR <- function(fits, TR, pal) {
  # browser()
  fits <- fits[fits$trange != '5-35',]
  palette(pal)
  fits$group <- paste(fits$soil, fits$trange)
  fits$group <- factor(fits$group, levels = c('cropped 5-20',
                                              'bare fallow 5-20', "cropped 20-35", "bare fallow 20-35"))
  if(TR == 'Ea') ylab <- expression(paste(E[a],
                                          ' (kJ)')) else ylab <- expression(Q[10])
  plotname <- paste0(prefix, TR, '.', devname)
  plotfile <- file.path(savedir, plotname)
  if(export) devfun(file = plotfile, width = 800, height = 800)

  p <- ggplot(data = fits, aes_string(x='moist_vol', y=TR, group='var2',
                                      colour='var2', shape='var2')) +
    scale_linetype_manual(values=c("solid", "longdash")) +
    theme_bw(base_size = 18) +
    scale_color_manual(values=pal) +
    geom_line(size = 1) +
    # geom_smooth(size = 1.5, span=0.3, se=FALSE) +
    geom_point(size = 3) +
    geom_point(size = 1, colour='white') +
    ylab(ylab) +
    xlab(expression(paste("Soil moisture (", m^3*m^-3,")"))) +
    theme(legend.justification=c(1,1), legend.position=c(0.47,0.48),
          legend.box='horizontal',
          # legend.background = element_rect(colour = "grey"),
          legend.title = element_blank(),
          legend.key.height=unit(0.6,'cm'), legend.key.width=unit(1,'cm'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    # facet_wrap(~group, scales = 'free_y')
    facet_grid(trange~soil, scales = 'free_y')#, scales = 'free_y'

  print(p)
}

# Get colors for plots
# Rename values for plotting
fit.temp$var2 <-  revalue(
  fit.temp$var, c(C_R_rm = 'R_mod', C_R_ro = "R_obs", C_dec_r = "D_mod"))
pal <- c('#9ACD32', '#698B22', '#CD8500') # 'olivedrab3', 'olivedrab4', 'orange3'
pal <- c('grey10', 'grey50','grey70')
TR <- "Ea"
PlotTR(fit.temp, TR, pal)

# PlotTR(fit.temp, TR, pal, '5-20')
# PlotTR(fit.temp, TR, pal, '20-35')
# PlotTR(fit.temp, TR, pal, '5-35')

if(export) while(dev.cur() > 1) dev.off()

