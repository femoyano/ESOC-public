mtdata3<- mtdata3[order(mtdata3$site, mtdata3$moist_grav, mtdata3$tstage),]
m <- 0.065
View(mtdata3[mtdata3$moist_grav==m,])

tempDat <- mtdata3[mtdata3$moist_grav==m & mtdata3$site=="maize",]
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
#This adds a column of color values
# based on the y values
Col <- rbPal(10)[as.numeric(cut(tempDat$hour, breaks = 10))]
Col <- adjustcolor(Col, alpha.f = 0.7)
with(tempDat, plot(C_R_r ~ temp, pch = 20, cex =2, col = Col))
