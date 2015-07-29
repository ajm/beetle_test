args <- commandArgs(trailingOnly = TRUE)
name = args[1]

Omegas <- read.table(file= paste(name, ".results", sep=""))
MeanOmegas <- read.table(file= paste(name, "_means.results", sep=""))
FltrdMeanOmegas <- read.table(file= paste(name, "_fltrd_means.results", sep=""))

Omegas <- Omegas[,1]
FltrdMeanOmegas <- FltrdMeanOmegas[,1]

MeanOmegas <- subset(MeanOmegas, MeanOmegas[,1]<=2)
MeanOmegas <- MeanOmegas[,1]

pdf(paste(name,".pdf", sep=""), width=12, height=5)

par(mfrow=c(1,3))
hist(Omegas, breaks=seq(0, 2, l=101), freq=FALSE, xlab="Omegas", main="Omegas", xlim=c(0, 2), ylim=c(0,10.5))
lines(density(Omegas))
hist(MeanOmegas, breaks=seq(0, 2, l=101), freq=FALSE, col="pink", xlab="Omega means", main="Means of all omegas", xlim=c(0, 2), ylim=c(0, 10.5))
lines(density(MeanOmegas))
hist(FltrdMeanOmegas, breaks=seq(0, 2, l=101), freq=FALSE, col="red", xlab="Filtered omega means", main="Filtered omega means", xlim=c(0,2), ylim=c(0,10.5))
lines(density(FltrdMeanOmegas))

dev.off()


