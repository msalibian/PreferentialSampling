library(RandomFields)
library(geoR)
library(fields)
# load Template Model Builder
library(TMB)
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
# compile and load C++ file
compile("TMBCPP.cpp")
dyn.load(dynlib("TMBCPP"))
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
# estimate parameters ignoring any preferential effects
standardMLE <- likfit(sampData, coords = sampData$coords,
                      data = sampData$data, kappa=kappa, ini=c(0.5, 0.5))
# defined discretisation for TMB
# define grid to simulate Sj's
l=30
Sseq <- seq(0,1,length.out=l)
SGrid <- expand.grid(Sseq,Sseq)
# find closest point in Sj's to data locations
pointer1 <- vector(length=n)
for(i in 1:n){
  nearestPoint <- which.min((SGrid[,1] - sampData$coords[i,1])^2 + (SGrid[,2] - sampData$coords[i,2])^2)
  pointer1[i] <- nearestPoint
}
# calculate distance matrix
distMat = sqrt(outer(SGrid[,1], SGrid[,1],"-")^2 + outer(SGrid[,2], SGrid[,2],"-")^2)
# set up model
data <- list(Y1=sampData$coords[,1], Y2=sampData$coords[,2], Y=sampData$data,
             distMat=distMat, pointer=pointer1)
# parameters for TMB
parameters <- list(
    S=rep(0, nrow(SGrid)),
    mu=standardMLE$beta,
    phi=log(standardMLE$phi),
    sigma=log(sqrt(standardMLE$sigmasq)),
    tau=log(sqrt(standardMLE$tausq)),
    beta=beta
)
# initial paramaters (currently in log scale)
initPar <- c(standardMLE$beta,log(standardMLE$phi),log(sqrt(standardMLE$sigmasq)),
             log(sqrt(standardMLE$tausq+0.0001)),beta)
# construct TMB function and let it integrate out S
obj <- MakeADFun(data,parameters,random=c("S"),DLL="TMBCPP", hessian=FALSE)
# control for optim use
my.control <- list(trace=1,maxit=4000)
# use optim to maximise likelihood (can utilise gradient function depending on
# 'method' argument used)
opt <- optim(initPar,obj$fn,obj$gr)
param <- opt$par

################################################################################
# prediction ###################################################################
# predict at location (0.49, 0.49)
# specificy index of prediction location
locIndex <- 451
# simple kriging for non-preferential parameters
SKDat <- krige.control(obj.model = standardMLE, type.krige = "SK")
predNoPref <- krige.conv(sampData, loc = gridFull, krige = SKDat)
# calculate bias
biasNoPref <-predNoPref$predict[locIndex]-geodata$data[locIndex]
# simple kriging using TMB params
PrefSK <- SKDat
PrefSK$beta=param[1]
PrefSK$cov.pars[1]=exp(param[3])^2
PrefSK$cov.pars[2]=exp(param[2])
PrefSK$nugget=exp(param[4])^2
predPref <- krige.conv(sampData, loc = gridFull, krige = PrefSK)
# calculate TMB bias
biasPref <-predPref$predict[locIndex]-geodata$data[locIndex]

nonPrefParam <- c(standardMLE$beta, standardMLE$phi, standardMLE$sigmasq, standardMLE$tausq)
prefParam <- c(param[1],exp(param[2]),exp(param[3])^2,exp(param[4])^2,param[5])
# predict V2 (mode of laplace)
modePred <- obj$env$last.par.best[1:nrow(SGrid)]

image.plot(Sseq,Sseq,matrix(modePred, nrow=length(Sseq), ncol=length(Sseq)),
           xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
points(sampData$coords, pch=19, cex=.5)
