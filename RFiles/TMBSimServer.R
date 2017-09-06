library(RandomFields)
library(geoR)
library(fields)
library(prodlim)
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

# define grid as 100 by 100 on the unit square
xseq=seq(0,1,length.out=91)
yseq=seq(0,1,length.out=91)
gridFull=expand.grid(xseq,yseq)
# define covariance model for the field S
model <- RMwhittle(nu=kappa, var=sigma.sq, scale=phi)
# generate the raw data for S
prefParams <- NULL
nonPrefParams <- NULL
postBias <- NULL
krigBias <- NULL
nonPrefBias <- NULL
for(k in 1:100){
  set.seed(k)
  rawDat <- RFsimulate(model, x=as.matrix(gridFull))$variable1 + mean
  # combine coordinates X with corresponding values for S
  obj <- cbind(cbind(gridFull[,1], gridFull[,2]), rawDat)
  geodata <- as.geodata(obj, coords.col = 1:2, data.col = 3)
  # sample the data according to PP with intensity exp(beta*S(x))
  sampData <- sample.geodata(geodata, size = n, prob = exp(beta * geodata$data))
  # add nugget variance to Y's
  sampData$data <- sampData$data + rnorm(n, mean = 0, sd = sqrt(tau.sq))
  # plot the data
  # image.plot(xseq,yseq,matrix(rawDat, nrow=length(xseq), ncol=length(yseq)),
  #            xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
  # points(sampData$coords, pch=19, cex=.5)
  # estimate parameters ignoring any preferential effects
  standardMLE <- likfit(sampData, coords = sampData$coords,
                        data = sampData$data, kappa=kappa, ini=c(0.5, 0.5))
  # defined discretisation for TMB
  # define grid to simulate Sj's
  l=31
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
  ########################################################################
  # predict V2 (mode of laplace) #########################################
  ########################################################################
  # extract S posterior from TMB
  modePred <- obj$env$last.par.best[1:nrow(SGrid)]
  # plot prediction
  image.plot(Sseq,Sseq,matrix(modePred, nrow=length(Sseq), ncol=length(Sseq)),
             xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
  points(sampData$coords, pch=19, cex=.5)
  # take best parameters
  param <- obj$env$last.par.best[(nrow(SGrid)+1):(nrow(SGrid)+5)]
  # simple kriging for parameters
  SKDat <- krige.control(obj.model = standardMLE, type.krige = "SK")
  PrefSK <- SKDat
  PrefSK$beta=param[1]
  PrefSK$cov.pars[1]=exp(param[3])^2
  PrefSK$cov.pars[2]=exp(param[2])
  PrefSK$nugget=exp(param[4])^2
  predPref <- krige.conv(sampData, loc = SGrid, krige = PrefSK)
  nonPredPref <- krige.conv(sampData, loc = SGrid, krige = SKDat)
  # plot kriging predictions
  # image.plot(Sseq,Sseq,matrix(predPref$predict, nrow=length(Sseq), ncol=length(Sseq)),
  #            xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
  # points(sampData$coords, pch=19, cex=.5)
  # plot true surface at those points
  #rawDatSmall <- rawDat[seq(1, length(rawDat), 3)]
  # take the subset of the full grid for comparison with predictions from TMB
  matchedIndic <- row.match(SGrid,gridFull)
  rawDatSmall <- rawDat[matchedIndic]
  # image.plot(Sseq,Sseq,matrix(rawDatSmall, nrow=length(Sseq), ncol=length(Sseq)),
  #            xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
  # points(sampData$coords, pch=19, cex=.5)
  # compare real to pred surfaces
  postBias <- rawDatSmall - modePred
  postMSPE <- postBias^2

  krigBias <- rawDatSmall - predPref$predict
  krigMSPE <- krigBias^2

  prefParams <- rbind(prefParams, c(opt$par, opt$convergence))
  nonPrefParams <- rbind(nonPrefParams, c(standardMLE$beta,standardMLE$phi,standardMLE$sigmasq,standardMLE$tausq))
  postBias <- cbind(postBias, rawDatSmall - modePred)
  krigBias <- cbind(krigBias, rawDatSmall - predPref$predict)
  nonPrefBias <- cbind(nonPrefBias, rawDatSmall - nonPredPref$predict)

  write.table(prefParams, "paperPrefParams.txt")
  write.table(nonPrefParams, "paperNonPrefParams.txt")
  write.table(postBias, "paperPostBias.txt")
  write.table(krigBias, "paperKrigBias.txt")
  write.table(nonPrefBias, "paperNonPrefBias.txt")
}
