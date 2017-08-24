library(RandomFields)
library(geoR)
library(fields)
library(mvtnorm)
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

l=30

# compile and load C++ file
i=1
set.seed(i)
# define grid as 100 by 100 on the unit square
xseq=seq(0,1,length.out=l)
yseq=seq(0,1,length.out=l)
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
# defined discretisation for TMB
# define grid to simulate Sj's
CMat <- matrix(0, nrow=n, ncol=l^2)
pointer <- NULL
for(j in 1:n){
  tmp <- which(geodata$data %in% sampData$data[j])
  CMat[j, tmp] <- 1
  pointer <- c(pointer, tmp)
}
# add nugget variance to Y's
sampData$data <- sampData$data + rnorm(n, mean = 0, sd = sqrt(tau.sq))
# plot the data
image.plot(xseq,yseq,matrix(rawDat, nrow=length(xseq), ncol=length(yseq)),
           xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
points(sampData$coords, pch=19, cex=.5)
# estimate parameters ignoring any preferential effects
standardMLE <- likfit(sampData, coords = sampData$coords,
                      data = sampData$data, kappa=kappa, ini=c(0.5, 0.5))
# calculate distance matrix
distMat = sqrt(outer(gridFull[,1], gridFull[,1],"-")^2 + outer(gridFull[,2], gridFull[,2],"-")^2)
covMat <- matern(distMat, phi, kappa)*sigma.sq
#
# set our S*
Sstar <- geodata$data
########################################################################
# [S] ##################################################################
########################################################################
Sden <- dmvnorm(Sstar, mean = rep(mean, l^2), sigma = covMat, log = T)
########################################################################
# [Y|S,X] ##############################################################
########################################################################
meanCond <- Sstar[pointer]
covMatCond <- diag(n)*tau.sq
Ycondden <- dmvnorm(sampData$data, mean = meanCond, sigma = covMatCond, log = T)
########################################################################
# [S|Y,X] ############################################################## THIS IS CURRENTLY S0|Y,X
########################################################################
Sigman <- CMat%*%covMat%*%t(CMat)
SigmanN <- CMat%*%covMat
SigmaNn <- covMat%*%t(CMat)

meanCondS <- SigmaNn%*%solve(Sigman+covMatCond)%*%(sampData$data-mean)
covMatCondS <- covMat - SigmaNn%*%solve(Sigman+covMatCond)%*%SigmanN
Scondden <- dmvnorm(Sstar-mean, mean = meanCondS, sigma = covMatCondS, log = T)
########################################################################
# [X|S] ################################################################
########################################################################
h=(max(xseq)-min(xseq))/l
np<-n
ng<-l^2
area.cell<-h*h
logintegral<-log(area.cell*sum(exp(beta*Sstar)))
Xcondden<-beta*sum(Sstar[pointer]) - np*logintegral
########################################################################
# final ratio ##########################################################
########################################################################
neglogLik <- -((Sden+Xcondden+Ycondden)-Scondden)
neglogLik

diggleRatio <-  -((Sden+Ycondden)-Scondden)
# diggleRatio
# > Sden
# [1] -63.9343
# > Ycondden
# [1] -28.05524
# > Scondden
# [1] 2.245432
# > Xcondden
# [1] 74.45182
#
# > Sden
# [1] -87.37529
# > Ycondden
# [1] -4597.532
# > Scondden
# [1] -4590.672
# > Xcondden
# [1] 74.45182
# > logLik
# [1] -19.78316

Sstar <- geodata$data

Sigman <- CMat%*%covMat%*%t(CMat)
SigmanN <- CMat%*%covMat
SigmaNn <- covMat%*%t(CMat)
########################################################################
# [S0] ##################################################################
########################################################################
Sden <- dmvnorm(Sstar[pointer], mean = rep(mean, n), sigma = Sigman, log = T)
########################################################################
# [Y|S,X] ##############################################################
########################################################################
meanCond <- Sstar[pointer]
covMatCond <- diag(n)*tau.sq
Ycondden <- dmvnorm(sampData$data, mean = meanCond, sigma = covMatCond, log = T)
########################################################################
# [S|Y,X] ##############################################################
########################################################################
meanCondS <- Sigman%*%solve(Sigman+covMatCond)%*%(sampData$data-mean)
covMatCondS <- Sigman - Sigman%*%solve(Sigman+covMatCond)%*%Sigman
Scondden <- dmvnorm(Sstar[pointer]-mean, mean = meanCondS, sigma = covMatCondS, log = T)
########################################################################
# [X|S] ################################################################
########################################################################
h=(max(xseq)-min(xseq))/l
np<-n
ng<-l^2
area.cell<-h*h
logintegral<-log(area.cell*sum(exp(beta*Sstar)))
Xcondden<-beta*sum(Sstar[pointer]) - np*logintegral
########################################################################
# final ratio ##########################################################
########################################################################
neglogLik <- -((Sden+Xcondden+Ycondden)-Scondden)
neglogLik


diggleRatio <-  -((Sden+Ycondden)-Scondden)
# diggleRatio
