# Exp 1 plots

library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(foreach)

dplot <- all_out
dplot <- dplot[dplot$ypr>1000,]
# dplot <- dplot[dplot$ypr == 3000, ]

dplot$month <- ceiling(dplot$time/720)
dplot$year <- dplot$month/12
dplot$C_tot <- dplot$C_P + dplot$C_D + dplot$C_E + dplot$C_M
dplot$group <- interaction(dplot$tempC, dplot$ypr, dplot$mode)


dplot1 <- ddply(.data = dplot, .variables = .(group), function(df) {
  df$C_rel <- df$C_tot/df$C_tot[1]
  df$C_dif <- df$C_tot-df$C_tot[1]
  df$C_P_rel <- df$C_P/df$C_P[1]
  df$C_P_dif <- df$C_P-df$C_P[1]
  return(df)})

ggplot(data = dplot1, aes(group = ypr, colour = ypr)) +
  geom_line(mapping = aes(x = year, y = C_dif)) +
  # geom_smooth(mapping = aes(x = year, y = C_rel)) +
  facet_grid(tempC~mode) + #, scales = "free_y") + #
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.background = element_blank())

dplot2 <- ddply(.data = dplot, .variables = .(group), function(df) {
  df2  <- data.frame(C_rel = tail(df$C_tot, 1)/df$C_tot[1])
  df2$C_dif <- tail(df$C_tot, 1)-df$C_tot[1]
  df2$C_P_rel <- tail(df$C_P, 1)/df$C_P[1]
  df2$C_P_dif <- tail(df$C_P, 1)-df$C_P[1]
  df2$C_P_ini <- df$C_P[1]
  df2$C_tot_ini <- df$C_tot[1]
  df2$Temperature <- as.factor(df$tempC[1])
  df2$VWC <- mean(df$moist)
  df2$mode <- df$mode[1]
  df2$ypr <- df$ypr[1]
  return(df2)})
dplot2$mode <- factor(dplot2$mode, levels = c("fixed", "dynamic"))

ggplot(data = dplot2, aes(x = ypr, y = C_P_rel, colour = Temperature)) +
  scale_color_manual(values = c('grey0', 'grey30', 'grey50', 'grey65', 'grey80', 'grey90'),
                     name = expression(paste(degree,"C"))) +
  geom_line(mapping = aes(), size = 1.2) +
  facet_grid(.~mode, scales = "free_y") + # scales = "free_y"
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # legend.position = c(0.85,0.3),
        axis.title.y = element_text(margin = margin(0,8,0,0)),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="grey")
        ) +
  xlab(expression(paste("Yearly precipitation (mm)"))) +
  ylab("Relative POC change") +
  geom_hline(yintercept = 1, lty = 2)

ggplot(data = dplot2, aes(x = ypr, y = C_P_dif, colour = Temperature)) +
  scale_color_manual(values = c('grey0', 'grey30', 'grey50', 'grey65', 'grey80',
                                'grey90'), name = expression(paste(degree,"C"))) +
  geom_line(size = 1.2) +
  facet_grid(.~mode, scales = "free_y") + # scales = "free_y"
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # legend.position = c(0.85,0.3),
        axis.title.y = element_text(margin = margin(0,8,0,0)),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="grey")
  ) +
  xlab(expression(paste("Yearly precipitation (mm)"))) +
  ylab(expression(paste("POC change (kg  ",m^-2, ")"))) +
  geom_hline(yintercept = 0, lty = 2)
