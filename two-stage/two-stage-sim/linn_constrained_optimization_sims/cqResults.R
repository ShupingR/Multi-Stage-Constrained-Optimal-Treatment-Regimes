## @file cqResults.R

#library (ggplot2)

setwd ("/Users/kalinnMB/Dropbox/gsl_constrained_code")
#setwd ("/Users/kalinn/Desktop")

sim = read.table ("simResults.txt", header=T)
head (sim)
kappa = seq (12, 20, by=.25)


mxYmnZ = mean (sim$maxYmeanZ)
maxYmn = mean (sim$uncEY)
minZmn = mean (sim$uncEZ[1:26])

pdf ("simResults.pdf", width=7, height=6)
par (mar=c(5,5,1,2))
plot (kappa[4:length (kappa)], sim$EYopt[4:length (kappa)], 
	xlab=expression (paste (kappa)), cex.lab=1.6, cex.axis=1.4, 
	xlim=c(12,20), ylim=c(0, 45), 
	ylab=expression (paste ("Estimated outcomes under ", 
		hat(pi)^opt)), col="slateblue4") 
lines (kappa[4:length (kappa)], sim$EYopt[4:length (kappa)], 
	lwd=2, col="slateblue4")
lines (kappa[4:length (kappa)], sim$EZopt[4:length (kappa)], 
	lwd=2, col="seagreen4")
points (kappa[4:length (kappa)], sim$EZopt[4:length (kappa)], 
	col="seagreen4")
abline (h=maxYmn, col="gray50", lwd=2)
abline (h=minZmn, col="gray50", lwd=2, lty=2)
rect (-10, -5, 12.5, 49, col="darkred")
legend ("topright", c (expression (paste (hat(E)^hat(pi)^opt, "(Y)")), 
	expression (paste (hat(E)^hat (pi)^opt, "(Z)")),  
	expression (paste (hat(E)^hat(pi)^{Y-opt}, "(Y)")),	
	expression (paste (hat(E)^hat(pi)^{Z-opt}, "(Z)"))),	
	col=c("slateblue4", "seagreen4", "gray50", "gray50"), 
	lty=c(1, 1, 1, 2), pch=c(1, 1, -1, -1), seg.len=5, ncol=2,
	lwd=c(2, 2, 2, 2))
dev.off ()

# Grayscale
setEPS ()
postscript ("results250.eps", width=7, height=6)
par (mar=c(5,5,1,2))
plot (kappa[4:length (kappa)], sim$EYopt[4:length (kappa)], 
	xlab=expression (paste (kappa)), cex.lab=1.6, cex.axis=1.4, 
	xlim=c(12,20), ylim=c(10, 45), 
	ylab=expression (paste ("Estimated outcomes under ", 
		hat(pi)^opt))) 
lines (kappa[4:length (kappa)], sim$EYopt[4:length (kappa)], 
	lwd=2, lty=2)
lines (kappa[4:length (kappa)], sim$EZopt[4:length (kappa)], 
	lwd=2)
points (kappa[4:length (kappa)], sim$EZopt[4:length (kappa)], pch=19)
abline (h=maxYmn, col="gray50", lwd=2)
abline (h=minZmn, col="gray50", lwd=2, lty=2)
rect (10, -5, 12.5, 49, col="gray80")
legend ("topright", c (expression (paste (hat(E)^hat(pi)^opt, "(Y)")), 
	expression (paste (hat(E)^hat (pi)^opt, "(Z)")),  
	expression (paste (hat(E)^hat(pi)^{Y-opt}, "(Y)")),	
	expression (paste (hat(E)^hat(pi)^{Z-opt}, "(Z)"))),	
	col=c("black", "black", "gray50", "gray50"), 
	lty=c(2, 1, 1, 2), pch=c(1, 19, -1, -1), seg.len=5, ncol=2,
	lwd=c(2, 2, 2, 2))
dev.off ()















