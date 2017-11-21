library(RandomFields)
library(geoR)
library(fields)
library(prodlim)
# load Template Model Builder
library(TMB)
# load INLA for mesh + sparse matrices
library(INLA)
# sample size
n=100
###########################################
# Set field parameters ####################
###########################################
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
# compile and load C++ file
compile("TMBexample.cpp")
dyn.load(dynlib("TMBexample"))
# define grid as 91 by 91 on the unit square
xseq=seq(0,1,length.out=91)
yseq=seq(0,1,length.out=91)
gridFull=expand.grid(xseq,yseq)
# define covariance model for the field S
model <- RMwhittle(nu=nu, var=sigma.sq, scale=phi)
# generate the raw data for S
set.seed(2)
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
# estimate parameters ignoring any preferential effects
standardMLE <- likfit(sampData, coords = sampData$coords,
                      data = sampData$data, kappa=nu, ini=c(0.5, 0.5))
# define lattice discretisation for TMB
m=31
Sseq <- seq(0,1,length.out=m)
# prediction grid (ie/ lattice)
predGrid <- expand.grid(Sseq,Sseq)
colnames(predGrid) <- c("V1", "V2")
# create larger grid including sampled locations
TMBGrid <- unique(rbind(predGrid, sampData$coords))
# pointer for sampling locations (using C++ indexing)
pointer <- row.match(data.frame(sampData$coords), TMBGrid) -1
## Simple default 15% extension, and refinement based only
## on a minimum angle criterion
mesh <- inla.mesh.create(loc = as.matrix(TMBGrid),
                         extend = T, refine = T)
plot(mesh, asp=1)
points(sampData$coords, col='red')
## Locate the input locations in the output mesh
ii0 <- mesh$idx$loc
# create data frame for TMB -  note indicies using C++ indexing (starts at 0 not 1)
data <- list(Y1=sampData$coords[,1], Y2=sampData$coords[,2], Y=sampData$data,
             pointer=pointer, meshidxloc=mesh$idx$loc-1)
# add elements for sparse precision matrix
data$spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
# indicator (C++ indexing) for lattice points in grid
data$matchedIndic <- row.match(predGrid, TMBGrid) - 1
# Number of points in mesh (including supporting points)
n_s = nrow(data$spde$M0)
########################################################################
# Use TMB to estimate corrected parameters #############################
########################################################################
# parameters for TMB
parameters <- list(
  S=rep(0, n_s),
  mu=standardMLE$beta,
  log_phi=log(standardMLE$phi),
  log_kappa=log(sqrt(1/(4*pi*(standardMLE$phi^-2)*standardMLE$sigmasq))),
  log_tau=log(sqrt(standardMLE$tausq)),
  beta=beta
)
# initial paramaters
initPar <- c(standardMLE$beta, log(standardMLE$phi), log(sqrt(1/(4*pi*(standardMLE$phi^-2)*standardMLE$sigmasq))),
             log(sqrt(standardMLE$tausq+0.0001)), beta)
# construct TMB function and let it integrate out latent field S
obj <- MakeADFun(data,parameters,random=c("S"),DLL="TMBexample", method = "nlminb", hessian=FALSE, silent=FALSE)
# use nlminb to maximise likelihood
opt <- nlminb(initPar,obj$fn,obj$gr)
report_spde <- obj$report()
# obtained preferentially corrected parameters
param <- c(opt$par[1], exp(opt$par[2]), (report_spde$sigma)^2, (exp(opt$par[4])^2), opt$par[5])
########################################################################
# Compare TMB and non-preferential predictions #########################
########################################################################
########################################################################
# Non-preferential predictions through kriging #########################
########################################################################
SKDat <- krige.control(obj.model = standardMLE, type.krige = "SK")
nonPredPref <- krige.conv(sampData, loc = predGrid, krige = SKDat, output=list(signal=T))
########################################################################
# Preferential predictions through TMB #################################
########################################################################
# extract S posterior from TMB
modePredPref <- obj$env$last.par.best[1:nrow(predGrid)]
# match indicies from TMB grid to grid used to generate data
matchedIndic <- row.match(predGrid,gridFull)
# get true field on TMB grid
rawDatSmall <- rawDat[matchedIndic]
# obtain standard errors
sdre <- sdreport(obj)
#
summary(sdre, "fixed")
# prediction variances
predVar <- (summary(sdre, "random")[1:nrow(predGrid),2])^2
# Compare true field with preferential and non-preferential predictions
range1=c(min(c(rawDatSmall,modePredPref,nonPredPref$predict)), max(c(rawDatSmall,modePredPref, nonPredPref$predict)))
# TRUE
image.plot(Sseq,Sseq,matrix(rawDatSmall, nrow=length(Sseq), ncol=length(Sseq)),
           xlab="Longitude", ylab="Latitude", main="True", zlim=range1, col=rev(heat.colors(20)))
points(sampData$coords, pch=19, cex=.5)

range2=c(min(c(predVar, nonPredPref$krige.var)), max(c(predVar, nonPredPref$krige.var)))

par(mfrow=c(2,2))
# Preferential predictions and variances
image.plot(Sseq,Sseq,matrix(modePredPref, nrow=length(Sseq), ncol=length(Sseq)),
           xlab="Longitude", ylab="Latitude", main="Pref Prediction", zlim=range1, col=rev(heat.colors(20)))
points(sampData$coords, pch=19, cex=.5)

image.plot(Sseq,Sseq,matrix(predVar, nrow=length(Sseq), ncol=length(Sseq)),
           xlab="Longitude", ylab="Latitude", main="Pref Variance", zlim=range2, col=rev(heat.colors(20)))
points(sampData$coords, pch=19, cex=.5)
# Non-Preferential prediction and variances
image.plot(Sseq,Sseq,matrix(nonPredPref$predict, nrow=length(Sseq), ncol=length(Sseq)),
           xlab="Longitude", ylab="Latitude", main="NonPref Prediction", zlim=range1, col=rev(heat.colors(20)))
points(sampData$coords, pch=19, cex=.5)

image.plot(Sseq,Sseq,matrix(nonPredPref$krige.var, nrow=length(Sseq), ncol=length(Sseq)),
           xlab="Longitude", ylab="Latitude", main="NonPref Variance", zlim=range2, col=rev(heat.colors(20)))
points(sampData$coords, pch=19, cex=.5)
########################################################################
# Compare Ignorance Scores #############################################
########################################################################
# Ignorance Score Function (see paper)
IGN <- function(pred, act, var) {
  ((pred - act)^2) / var + log(var)
}
IgnScorePref <- IGN(modePredPref, rawDatSmall, predVar)
IgnScoreNonPref <- IGN(nonPredPref$predict, rawDatSmall, nonPredPref$krige.var)
# Compare Mean Ignorance Score (MIGN)
mean(IgnScorePref)
mean(IgnScoreNonPref)
# Plot ignorance scores
range3=c(min(c(IgnScorePref, IgnScoreNonPref)), max(c(IgnScorePref, IgnScoreNonPref)))
par(mfrow=c(1,2))
# Preferential IGN
image.plot(Sseq,Sseq,matrix(IgnScorePref, nrow=length(Sseq), ncol=length(Sseq)),
           xlab="Longitude", ylab="Latitude", main="Pref IGN", zlim=range3, col=rev(heat.colors(20)))
points(sampData$coords, pch=19, cex=.5)
# Non-Preferential IGN
image.plot(Sseq,Sseq,matrix(IgnScoreNonPref, nrow=length(Sseq), ncol=length(Sseq)),
           xlab="Longitude", ylab="Latitude", main="NonPref IGN", zlim=range3, col=rev(heat.colors(20)))
points(sampData$coords, pch=19, cex=.5)
