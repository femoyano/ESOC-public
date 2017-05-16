# Exp 1 plots

library(ggplot2)
library(plyr)
library(RColorBrewer)

dplot <- all_out

dplot$month <- ceiling(dplot$time/720)
dplot$year <- dplot$month/12
dplot$YPE <- as.factor(dplot$nypr)
dplot$tempC <- as.factor(dplot$tempC)
dplot$C_tot <- dplot$C_P + dplot$C_D + dplot$C_E + dplot$C_M
dplot$group <- interaction(dplot$mode, dplot$ypr)

dplot$group <- factor(dplot$group, levels = c('fixed.1000', 'fixed.2000', 'dynamic.1000', 'dynamic.2000'))
# dplot$mode <- factor(dplot$mode, levels = unique(dplot$mode))

dplot <- ddply(.data = dplot, .variables = .(tempC, group),
               function(df) {
                 df$C_rel <- df$C_tot/df$C_tot[1]
                 df$C_dif <- df$C_tot-df$C_tot[1]
                 df$C_P_rel <- df$C_P/df$C_P[1]
                 df$C_P_dif <- df$C_P-df$C_P[1]
                 return(df)})

dplot <- dplot[dplot$tempC != 15 & dplot$year < 50 & dplot$nypr != 3, ]
len <- 8
vars <- data.frame(expand.grid(levels(dplot$tempC)[c(1,3)], levels(dplot$group)[1:4]))
colnames(vars) <- c("tempC", "group")
dat <- data.frame(x = rep(1.5, len), y = rep(1.87, len), vars, labs=c(
  'a','e','b','f','c','g','d','h'))

ggplot(data = dplot, aes(x = year, y = C_P_rel)) +
  geom_line(mapping = aes(colour = YPE)) +
  facet_grid(tempC~group) + #, scales = "free_y"
  scale_y_continuous(limits = c(0.48, 1.92)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # strip.background = element_blank(),
        # strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(0,8,0,0))
        ) +
  ylab("Relative POC change") +
  geom_text(aes(x, y, label=labs, group=NULL), data=dat)
