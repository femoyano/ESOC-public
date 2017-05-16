# Exp 1 plots

library(ggplot2)
library(plyr)
library(RColorBrewer)
library(foreach)

dplot <- all_out

dplot$month <- ceiling(dplot$time/720)
dplot$year <- dplot$month/12
dplot$YPE <- as.factor(dplot$nypr)
dplot$tempC <- as.factor(dplot$tempC)
dplot$C_tot <- dplot$C_P + dplot$C_D + dplot$C_E + dplot$C_M
dplot$group <- interaction(dplot$mode, dplot$ypr)
dplot$group <- factor(dplot$group, levels = c('fixed.1000', 'fixed.2000', 'dynamic.1000', 'dynamic.2000'))
dplot$group2 <- interaction(dplot$ypr, dplot$nypr)

dplot <- foreach(i = unique(dplot$group), .combine = 'rbind') %do% {
  x <- dplot[dplot$group == i,]
  ref <- x[x$tempC == 15 & x$nypr == 90, ]
  ddply(.data = x, .variables = .(tempC, ypr, nypr),
        function(df) {
          df$C_rel <- df$C_tot/df$C_tot[1]
          df$C_dif <- df$C_tot-df$C_tot[1]
          df$C_refrel <- df$C_tot/ref$C_tot
          df$C_refdif <- df$C_tot-ref$C_tot
          df$C_P_rel <- df$C_P/df$C_P[1]
          df$C_P_dif <- df$C_P-df$C_P[1]
          df$C_P_refrel <- df$C_P/ref$C_P
          df$C_P_refdif <- df$C_P-ref$C_P
          return(df)})
}

dplot1 <- dplot[dplot$ypr != "C (1000 mm)" &
                  dplot$ypr != "C (2000 mm)" &
                  dplot$tempC != 15 &
                  dplot$nypr > 3, ]

len <- 8
vars <- data.frame(expand.grid(levels(dplot1$tempC)[c(1,3)], levels(dplot1$group)[1:4]))
colnames(vars) <- c("tempC", "group")
dat <- data.frame(x = c(rep(1.5, 8)),
                  # y = c(0.76, 1.145, 0.76, 1.145, 0.76, 1.145, 0.76, 1.145), # SOC
                  y = c(0.68, 1.17, 0.68, 1.17, 0.68, 1.17, 0.68, 1.17), # POC
                  vars, labs=c('a','e','b','f','c','g','d','h'))

ggplot(data = dplot1, aes(x = year, y = C_P_rel)) +
  geom_line(mapping = aes(colour = YPE)) +
  # geom_smooth() +
  facet_grid(tempC~group, scales = "free_y") + # scales = "free_y"
  ylab("Relative POC change") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  geom_text(aes(x, y, label=labs, group=NULL), data=dat)


dplot2 <- dplot[dplot$ypr != "C (1000 mm)" & dplot$ypr != "C (2000 mm)" &
                  dplot$nypr > 3  &
                  dplot$tempC != 15 &
                  dplot$nypr < 25 &
                  dplot$ypr == 1000, ]

ggplot(data = dplot2, aes(group = group2, fill = group2)) +
  geom_boxplot(mapping = aes(x = YPE, y = C_rel)) +
  # facet_grid(nypr~ypr) + #, scales = "free_y")
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.background = element_blank())
