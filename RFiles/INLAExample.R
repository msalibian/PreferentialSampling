########################################################################
# Thank you to an anonymous reviewer for INLA code  ####################
# adapted in this document                          ####################
########################################################################
library(RandomFields)
library(geoR)
library(fields)
library(prodlim)
# load INLA
library(INLA)
# sample size
n=100
########################################################################
# Set field parameters #################################################
########################################################################
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
# defined discretisation for INLA predictions/[X|S]
m=31
Sseq <- seq(0,1,length.out=m)
# prediction grid (lattice)
predGrid <- expand.grid(Sseq,Sseq)
colnames(predGrid) <- c("V1", "V2")
# full grid (including sampling locations)
INLAgrid <- unique(rbind(predGrid, sampData$coords))
########################################################################
# Use INLA to estimate corrected parameters ############################
########################################################################
y.pref = rep(NA, nrow(INLAgrid))
# find sampling locations in full grid
pointer1 <- row.match(data.frame(sampData$coords), INLAgrid)

pp.pref <- pointer1
y.pref[pp.pref] = sampData$data
# Create INLA mesh
mesh <- inla.mesh.create(loc = as.matrix(INLAgrid),
                         extend = TRUE, refine = TRUE)
plot(mesh, asp=1)
# Locate the input locations in the output mesh
ii0 <- mesh$idx$loc
# Create INLA Matern SPDE model
spde <- inla.spde2.pcmatern(mesh = mesh, constr = FALSE,
                            prior.range = c(.15, 0.01),
                            prior.sigma = c(1.5, 0.01))
# Prepare the data sets for INLA
n2 = nrow(INLAgrid)
ii <- c(ii0, rep(NA, n2))
jj <- c(rep(NA, n2), ii0)
alpha = c(rep(0,n2), rep(1,n2))
mu = c(rep(1,n2), rep(0,n2))
points.pref = rep(0,n2)
points.pref[pp.pref] = 1
# Preferential sampling
yy.pref = matrix(NA,2*n2,2)
yy.pref[1:n2,1] = y.pref
yy.pref[n2+1:n2,2] = points.pref
data.pref.pref = list(yy=yy.pref,mu=mu,ii=ii,jj=jj,alpha=alpha)
# Create data for INLA fitting
data.pref.pref_spde <- list(yy = yy.pref, mu = mu, ii = ii, jj = jj, alpha = alpha)
formula <- yy ~ alpha + mu + f(ii, model = spde) +
  f(jj, copy = "ii", fixed = FALSE, param = c(0, 0.1)) -1
# Fit model
pref.model <-
  inla(formula, family = c("gaussian", "poisson"),
       control.family = list(list(initial = log(1/0.1), fixed = FALSE), list()),
       control.inla = list(strategy = "gaussian", int.strategy = "eb"),
       control.predictor = list(compute = TRUE),
       data = data.pref.pref_spde, verbose = TRUE)
########################################################################
# Compare INLA and non-preferential predictions ########################
########################################################################
########################################################################
# Non-preferential predictions through kriging #########################
########################################################################
SKDat <- krige.control(obj.model = standardMLE, type.krige = "SK")
nonPredPref <- krige.conv(sampData, loc = predGrid, krige = SKDat, output=list(signal=T))
########################################################################
# Preferential predictions through INLA ################################
########################################################################
# Preferential parameters taken from INLA results
prefParam <- c(pref.model$summary.fixed[2,"mean"],
               inla.emarginal(function(x) x / sqrt(8*nu),
                              pref.model$marginals.hyperpar[[2]]),
               inla.emarginal(function(x) x^2, pref.model$marginals.hyperpar[[3]]),
               inla.emarginal(function(x) 1/x, pref.model$marginals.hyperpar[[1]]),
               pref.model$summary.hyperpar[4,"mean"])
# match indicies from INLA grid to grid used to generate data
matchedIndic <- row.match(predGrid,gridFull)
# get true field on INLA grid
rawDatSmall <- rawDat[matchedIndic]
# INLA predictions and prediction variances
predPrefINLA <-  list(predict=pref.model$summary.fitted.values[1:nrow(predGrid), "mean"],
                      variance=pref.model$summary.fitted.values[1:nrow(predGrid), "sd"]^2)
# Compare true field with preferential and non-preferential predictions
range1=c(min(c(rawDatSmall,predPrefINLA$predict,nonPredPref$predict)), max(c(rawDatSmall,predPrefINLA$variance, nonPredPref$predict)))
# TRUE
image.plot(Sseq,Sseq,matrix(rawDatSmall, nrow=length(Sseq), ncol=length(Sseq)),
           xlab="Longitude", ylab="Latitude", main="True", zlim=range1, col=rev(heat.colors(20)))
points(sampData$coords, pch=19, cex=.5)

range2=c(min(c(predPrefINLA$variance, nonPredPref$krige.var)), max(c(predPrefINLA$variance, nonPredPref$krige.var)))
par(mfrow=c(2,2))
# Preferential predictions and variances
image.plot(Sseq,Sseq,matrix(predPrefINLA$predict, nrow=length(Sseq), ncol=length(Sseq)),
           xlab="Longitude", ylab="Latitude", main="Pref Prediction", zlim=range1, col=rev(heat.colors(20)))
points(sampData$coords, pch=19, cex=.5)

image.plot(Sseq,Sseq,matrix(predPrefINLA$variance, nrow=length(Sseq), ncol=length(Sseq)),
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
IgnScorePref <- IGN(predPrefINLA$predict, rawDatSmall, predPrefINLA$variance)
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
