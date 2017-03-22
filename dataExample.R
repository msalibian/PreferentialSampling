library(RandomFields)
library(geoR)
library(fields)
# sample size
n=100
# choose parameters
# choose preferential parameter (beta=0 is non-preferential)
beta=2
# marginal variance of the field
sigma.sq=1.5
# nugget variance
tau.sq=0.1
# smoothness parameter
kappa=1
# scale
phi=0.15
# mean
mean=4
# set seed
i=1
set.seed(i)
# define grid as 100 by 100 on the unit square
xseq=seq(0,1,length.out=100)
yseq=seq(0,1,length.out=100)
gridFull=expand.grid(xseq,yseq)
# define covariance model for the field S
model <- RMwhittle(nu=kappa, var=sigma.sq, scale=phi)
# generate the raw data for S
rawDat <- RFsimulate(model, x=as.matrix(gridFull))$variable1 + mean
# combine coordinates X with corresponding values for S
obj <- cbind(cbind(gridFull[,1], gridFull[,2]), rawDat)
geodata <- as.geodata(obj, coords.col = 1:2, data.col = 3)
# sample the data according to PP with intensity exp(beta*S(x))
sampData <- sample.geodata(geodata, size = n, prob = exp(beta * geodata$data))
# add nugget variance to Y's
sampData$data <- sampData$data + rnorm(n, mean = 0, sd = sqrt(tau.sq))
# plot the data
image.plot(xseq,yseq,matrix(rawDat, nrow=length(xseq), ncol=length(yseq)),
           xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
points(sampData$coords, pch=19, cex=.5)
