args <- commandArgs(trailingOnly = TRUE)
name = args[1]

Omegas <- read.table(file= paste(name, ".results", sep=""))
MeanOmegas <- read.table(file= paste(name, "_means.results", sep=""))
FltrdMeanOmegas <- read.table(file= paste(name, "_fltrd_means.results", sep=""))
PamlOmegas <- read.table(file = paste(name, "_paml.results", sep=""))

Omegas <- subset(Omegas, Omegas[,1]<=1)
Omegas <- Omegas[,1]


FltrdMeanOmegas <- subset(FltrdMeanOmegas, FltrdMeanOmegas[,1]<=1)
FltrdMeanOmegas <- FltrdMeanOmegas[,1]

MeanOmegas <- subset(MeanOmegas, MeanOmegas[,1]<=1)
MeanOmegas <- MeanOmegas[,1]

PamlOmegas <- subset(PamlOmegas, PamlOmegas[,1]<=1)
PamlOmegas <- PamlOmegas[,1]

pdf(paste(name,".pdf", sep=""), width=16, height=5)

par(mfrow=c(1,4))
hist(Omegas, breaks=seq(0, 1, l=51), freq=FALSE, xlab="Omegas", main=paste("Omegas, n=", NROW(Omegas), sep=""), xlim=c(0, 1), ylim=c(0,10.5))
lines(density(Omegas))
hist(MeanOmegas, breaks=seq(0, 1, l=51), freq=FALSE, col="pink", xlab="Omega means", main=paste("Means of all omegas, n=", NROW(MeanOmegas), sep=""), xlim=c(0, 1), ylim=c(0, 10.5))
lines(density(MeanOmegas))
hist(FltrdMeanOmegas, breaks=seq(0, 1, l=51), freq=FALSE, col="red", xlab="Filtered omega means", main=paste("Filtered omega means, n=", NROW(FltrdMeanOmegas), sep=""), xlim=c(0,1), ylim=c(0,10.5))
lines(density(FltrdMeanOmegas))
hist(PamlOmegas, breaks=seq(0, 1, l=51), freq=FALSE, col="pink", xlab="Paml", main=paste("Paml omega means, n=", NROW(PamlOmegas), sep=""), xlim=c(0,1), ylim=c(0,10.5))
lines(density(PamlOmegas))


dev.off()


