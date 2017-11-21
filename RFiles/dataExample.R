# load required librarys
library(RandomFields)
library(geoR)
library(fields)
# choose sample size
n=100
######################################################################################
# Set field parameters ###############################################################
######################################################################################
# choose preferential parameter (beta=0 is uniform random)
beta=1.5
# marginal variance of the field
sigma.sq=1.5
# nugget variance
tau.sq=0.1
# smoothness parameter
nu=1
# scale (range)
phi=0.15
# mean parameter (constant mean trend)
mean=4
# define grid as 91 by 91 on the unit square
xseq=seq(0,1,length.out=91)
yseq=seq(0,1,length.out=91)
gridFull=expand.grid(xseq,yseq)
# define covariance model for the field S
model <- RMwhittle(nu=nu, var=sigma.sq, scale=phi)
# generate the raw data for S
set.seed(999)
rawDat <- RFsimulate(model, x=as.matrix(gridFull), exactness=TRUE)$variable1 + mean
# combine coordinates X with corresponding values for S
obj <- cbind(cbind(gridFull[,1], gridFull[,2]), rawDat)
geodata <- as.geodata(obj, coords.col = 1:2, data.col = 3)
# sample the data  with intensity exp(beta*S(x))/int(exp(beta*S(u))du)
sampData <- sample.geodata(geodata, size = n, prob = exp(beta * geodata$data))
# add nugget variance to Y's
sampData$data <- sampData$data + rnorm(n, mean = 0, sd = sqrt(tau.sq))
# plot the data
image.plot(xseq,yseq,matrix(rawDat, nrow=length(xseq), ncol=length(yseq)),
           xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
points(sampData$coords, pch=19, cex=.5)
