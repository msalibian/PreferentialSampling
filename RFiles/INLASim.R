library(RandomFields)
library(geoR)
library(fields)
# load INLA
library(INLA)
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
y.pref = rep(NA,length(xseq)*length(yseq))
pp.pref = match(sampData$data,rawDat)
# add nugget variance to Y's
sampData$data <- sampData$data + rnorm(n, mean = 0, sd = sqrt(tau.sq))
y.pref[pp.pref] = sampData$data

# plot the data
image.plot(xseq,yseq,matrix(rawDat, nrow=length(xseq), ncol=length(yseq)),
           xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
points(sampData$coords, pch=19, cex=.5)
# estimate parameters ignoring any preferential effects
standardMLE <- likfit(sampData, coords = sampData$coords,
                      data = sampData$data, kappa=kappa, ini=c(0.5, 0.5))
### Prepare the data sets for INLA
n2 = length(xseq)^2
ii=c(1:n2,rep(NA,n2))
jj = c(rep(NA,n2), 1:n2)
alpha = c(rep(0,n2), rep(1,n2))
mu = c(rep(1,n2), rep(0,n2))

points.pref = rep(0,n2)
points.pref[pp.pref] = 1

### preferential sampling
yy.pref = matrix(NA,2*n2,2)
yy.pref[1:n2,1] = y.pref
yy.pref[n2+1:n2,2] = points.pref

data.pref.pref = list(yy=yy.pref,mu=mu,ii=ii,jj=jj,alpha=alpha)
formula.pref1 = yy ~ alpha + mu + f(ii, model = "matern2d", nrow=length(xseq), ncol=length(xseq), nu=1,
                                    initial = c(3, log(10)), fixed=c(FALSE,FALSE),
                                    param = c(1,.1, 1, .1), constr=TRUE) +
  f(jj, copy="ii", fixed=FALSE, param=c(0,0.1)) -1

pref.model.pref1 = inla(formula.pref1, family = c("gaussian", "poisson"),
                        control.family = list(list(initial = log(1/0.1), fixed=FALSE), list()),
                        control.inla= list(strategy = "gaussian"),
                        data = data.pref.pref, verbose = TRUE)
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
# simple kriging using INLA params
PrefSK <- SKDat
PrefSK$beta=pref.model.pref1$summary.fixed[2,1]
PrefSK$cov.pars[1]=1/pref.model.pref1$summary.hyperpar[2,1]
PrefSK$cov.pars[2]=sqrt(8*kappa)/pref.model.pref1$summary.hyperpar[3,1]
PrefSK$nugget=1/pref.model.pref1$summary.hyperpar[1,1]
predPref <- krige.conv(sampData, loc = gridFull, krige = PrefSK)
# calculate INLA bias
biasPref <-predPref$predict[locIndex]-geodata$data[locIndex]
# final parameters
nonPrefParam <- c(standardMLE$beta, standardMLE$phi, standardMLE$sigmasq, standardMLE$tausq)
prefParam <- c(pref.model.pref1$summary.fixed[2,1],sqrt(8*kappa)/(pref.model.pref1$summary.hyperpar[3,1]),
               (1/pref.model.pref1$summary.hyperpar[2,1]),(1/pref.model.pref1$summary.hyperpar[1,1]),
               pref.model.pref1$summary.hyperpar[4,1])

