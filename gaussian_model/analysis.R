#######################################
###### Load and Explore the Data ######
#######################################

cmaq <- read.csv(file="CMAQ.csv",header=T)
ozone <- read.csv(file="Ozone.csv",header=T)
predlocs <- read.csv(file="PredLocs.csv",header=T)
View(cmaq)
View(ozone)
View(predlocs)

# Merge station measurements of ozone levels with CMAQ Predictions Dataset

# Find the Euclidian Distance between CMAQ predictions and Station Measurements
library(fields)
eucdist <- rdist(cbind(ozone$Latitude,ozone$Longitude),cbind(cmaq$Latitude,cmaq$Longitude))

# Match Station Measurements with closest 1000 CMAQ Predictions
minval <- rep(NA,dim(eucdist)[1])
ozoneexpvar <- matrix(NA,nrow=length(ozone$Latitude),ncol=1000)
longitudeexpvar <- matrix(NA,nrow=length(ozone$Latitude),ncol=1000)
latitudeexpvar <- matrix(NA,nrow=length(ozone$Latitude),ncol=1000)
for(k in 1:dim(ozoneexpvar)[2]){
  for(i in 1:dim(eucdist)[1]){
    minval[i] <- min(eucdist[i,])
    for(j in 1:dim(eucdist)[2]){
      if(minval[i]==eucdist[i,j]) {
        eucdist[i,j] <- 100
        ozoneexpvar[i,k] <- cmaq$CMAQ_O3[j]
        longitudeexpvar[i,k] <- cmaq$Longitude[i]
        latitudeexpvar[i,k] <- cmaq$Latitude[i]
      }
    }
  }
}

# Only keep the first 250 explanatory variables
Ozone <- data.frame(matrix(cbind(ozone$Latitude,ozone$Longitude,ozone$Ozone.8.hr.max.,ozoneexpvar),nrow=800))[,-seq(from=254,to=1003,by=1)]
Ozone.fixed <- Ozone[,-c(1,2)]
View(Ozone)
View(Ozone.fixed)

# Look at the maps of the station observations and the cmaq predictions
library(maps) 
library(LatticeKrig) 

# Plot the measured ozone levels -800 station measurements of max.
# 8 hour O3 on May 22, 2005
quilt.plot(ozone$Longitude,ozone$Latitude,ozone$Ozone.8.hr.max.) 
map('state',add=TRUE)

# Plot the CMAQ Predictions of Ozone levels
# Mathematically simulates (on a fine spatial scale) formulation of ozone based on ground characteristics, temperatures, urban density, etc.
quilt.plot(cmaq$Longitude,cmaq$Latitude,cmaq$CMAQ_O3) 
map('state',add=TRUE)

# Export Maps:
pdf("expmaps.pdf")
par(mfrow=c(1,2))
quilt.plot(ozone$Longitude,ozone$Latitude,ozone$Ozone.8.hr.max.) 
map('state',add=TRUE)
quilt.plot(cmaq$Longitude,cmaq$Latitude,cmaq$CMAQ_O3) 
map('state',add=TRUE)
dev.off()

# Variogram Plot:
library(spatial)
spat.ls <- surf.ls(2,ozone$Longitude,ozone$Latitude,ozone$Ozone.8.hr.max.) # this is the surface
# The variogram should be increasing, then approximately level off at the "sill"
variogram(spat.ls,16)
library(modreg)
test <- variogram(spat.ls,16)

pdf("variogram.pdf")
scatter.smooth(test$x,test$y)
dev.off()

pdf("acf.pdf")
acf(ozone$Ozone.8.hr.max.,main="")
dev.off()

# Research Goals:
#   1. Scientists know that CMAQ is wrong:
#        - Understand the relationship between CMAQ and Station Measurement.
#   2. Want ground-level O3 predictions at lots of locations
#        - Prediction locations provided on website

# Statistical Goals and Issues to consider:
#   1. Estimating relationship between ground level O3 and CMAQ:
#       • Specifying a Model: How to deal with positive support? How are you going to define predictors (what is your X matrix)? How many predictors are you going to use (high dimensional)? Collinearity? Is the relationship linear?
#       • No IID Errors: Station measurements are not independent but are spatially correlated.
#   2. Predicting O3:
#       • Correlation: Ozone is correlated across Lon/Lat.
#       • Prediction Accuracy: How are you going to assess how accurate your predictions are?

# Best Subset Selection
library(MASS)
library(nlme)
library(car)
fit <- gls(log(X3)~X4+X5+X6+X7+X8+X9+X10+X11+X12+X13,correlation=corExp(form=~X2+X1,nugget=TRUE),data=Ozone,method="ML")
vif(fit)
best.subset <- regsubsets(X3~.,Ozone.fixed, nvmax=10)
best.subset.summary <- summary(best.subset)
best.subset.summary$outmat
best.subset.by.adjr2 <- which.max(best.subset.summary$adjr2)
best.subset.by.adjr2
coef(best.subset,6)

####################################################
####### Verify Assumptions of Gaussian Model #######
####################################################

fit.final <- gls(log(X3)~X4+X7+X10+X11+X12+X13,correlation=corExp(form=~X2+X1,nugget=TRUE),data=Ozone,method="ML")
summary(fit.final)

# Compare Fitted values to observed values
pdf("compare.pdf")
par(mfrow=c(1,2))
quilt.plot(ozone$Longitude,ozone$Latitude,ozone$Ozone.8.hr.max.) 
map('state',add=TRUE)
quilt.plot(ozone$Longitude,ozone$Latitude,exp(fitted(fit.final)))
map('state',add=TRUE)
dev.off()

# Decorrelate Regression Model #
nug <- 0.2889362
N <- 800
K <- 2685
R <- diag(N)
R <- (1-nug)*R+nug*diag(N)
L <- t(chol(R))
Linv <- solve(L)
epsilonhat <- Linv%*%(Y-X%*%bhat)
stres <- rep(NA,length(epsilonhat))
for(i in 1:length(epsilonhat)) stres[i] <- (epsilonhat[i]-mean(epsilonhat))/sqrt(var(epsilonhat))

# Check the assumptions of the model with the cholesky residuals:
hist(stres) # normal distribution
yhat <- X%*%bhat
plot(yhat,epsilonhat)
abline(h=0) # constant variance centered around 0

pdf("assump.pdf")
par(mfrow=c(1,2))
hist(stres,freq=FALSE,xlab="Decorrelated Standardized Residuals",main="Histogram of Standardized Residuals")
plot(yhat,epsilonhat,xlab="Predicted Values",ylab="Decorrelated Residuals",main="Fitted vs. Residuals Plot")
abline(h=0)
dev.off()

################################
####### Find Predictions #######
################################

# Find the measurements for the closest to each...
# Find the Euclidian Distance between CMAQ predictions and values to predict
library(fields)
eucdist <- rdist(cbind(predlocs$Latitude,predlocs$Longitude),cbind(cmaq$Latitude,cmaq$Longitude))

# Match Station Measurements with closest 1000 CMAQ Predictions
minval <- rep(NA,dim(eucdist)[1])
ozoneexpvar <- matrix(NA,nrow=length(predlocs$Latitude),ncol=10)
for(k in 1:dim(ozoneexpvar)[2]){
  for(i in 1:dim(eucdist)[1]){
    minval[i] <- min(eucdist[i,])
    for(j in 1:dim(eucdist)[2]){
      if(minval[i]==eucdist[i,j]) {
        eucdist[i,j] <- 100
        ozoneexpvar[i,k] <- cmaq$CMAQ_O3[j]
      }
    }
  }
}

predlocs$cmaq1 <- ozoneexpvar[,1]
predlocs$cmaq4 <- ozoneexpvar[,4]
predlocs$cmaq7 <- ozoneexpvar[,7]
predlocs$cmaq8 <- ozoneexpvar[,8]
predlocs$cmaq9 <- ozoneexpvar[,9]
predlocs$cmaq10 <- ozoneexpvar[,10]

Xstar <- model.matrix(~cmaq1+cmaq4+cmaq7+cmaq8+cmaq9+cmaq10,data=predlocs)
bhat <- cbind(coef(fit.final))
X <- model.matrix(~X4+X7+X10+X11+X12+X13,data=Ozone)
Y <- log(Ozone$X3)

# Calculate the predictions, prediction variance, and the prediction intervals
sig2 <- sigma(fit.final)^2
nug <- 0.2889362
N <- 800
K <- 2685
R <- diag(N+K)
R <- (1-nug)*R+nug*diag(N+K)
pred.mn <- Xstar%*%bhat + R[N+(1:K), (1:N)]%*%solve(R[(1:N),(1:N)])%*%(Y-X%*%bhat) #conditional mean of MVN
pred.var <- sig2*(R[N+(1:K),N+(1:K)]-(R[N+(1:K),(1:N)]%*%solve(R[(1:N),(1:N)])%*%R[(1:N),N+(1:K)])) # conditional variance of MVN
pred.int <- matrix(NA,nrow=K,ncol=2)
for(i in 1:K) pred.int[i,] <- pred.mn[i] + c(-1,1)*qnorm(0.975)*sqrt(diag(pred.var)[i])

# Plot the Predictions
pdf("pred.pdf")
quilt.plot(predlocs$Longitude,predlocs$Latitude,exp(pred.mn)) 
map('state',add=TRUE)
dev.off()

pdf("pred.pdf")
par(mfrow=c(1,3))
quilt.plot(predlocs$Longitude,predlocs$Latitude,exp(pred.int[,1])) 
map('state',add=TRUE)
quilt.plot(predlocs$Longitude,predlocs$Latitude,exp(pred.mn)) 
map('state',add=TRUE)
quilt.plot(predlocs$Longitude,predlocs$Latitude,exp(pred.int[,2])) 
map('state',add=TRUE)
dev.off()

##############################
###### Prediction Power ######
##############################

# R-squared using decorrelated values

rsquared <- 1-(sum(epsilonhat^2)/sum((Y-mean(Y))^2))

# CROSS VALIDATION to find mean bias, RPMSE, coverage, interval width

bias <- vector()
rpmse <- vector()
coverage <- vector()
width <- vector()
for(i in 1:1000){
  x <- sample(dim(Ozone)[1],8)
  test <- Ozone[x,]
  train <- Ozone[-x,]
  testgls  <- gls(log(X3)~X4+X7+X10+X11+X12+X13,correlation=corExp(form=~X2+X1,nugget=TRUE),data=train,method="ML")
  nug <- coef(testgls$modelStruct$corStruct,unconstrained=FALSE)[2]
  sig2 <- testgls$sigma^2
  bhat <- cbind(coef(testgls))
  X <- model.matrix(~X4+X7+X10+X11+X12+X13,data=train)
  Y <- log(train$X3)
  Xstar <- model.matrix(~X4+X7+X10+X11+X12+X13,data=test)
  N <- nrow(Ozone) #Number of observed time periods
  R = diag(N)
  R <- (1-nug)*R+nug*diag(N)
  
  pred.mn <- Xstar%*%bhat + R[x,(1:N)[-x]]%*%solve(R[(1:N)[-x],(1:N)[-x]])%*%(Y-X%*%bhat) #conditional mean of MVN
  pred.var <- sig2*(R[x,x]-(R[x,(1:N)[-x]]%*%solve(R[(1:N)[-x],(1:N)[-x]])%*%R[(1:N)[-x],x])) # conditional variance of MVN
  pred.int <- matrix(NA,nrow=length(x),ncol=2)
  for(j in 1:length(x)) pred.int[j,] <- pred.mn[j] + c(-1,1)*qnorm(0.975)*sqrt(diag(pred.var)[j])
  bias[i] = mean(pred.mn - log(test$X3))
  rpmse[i] = sqrt(mean((pred.mn-log(test$X3))^2))
  coverage[i] = mean(pred.int[,1] < log(test$X3) & 
                    pred.int[,2] > log(test$X3))
  width[i] = mean(pred.int[,2] - pred.int[,1])
}
mean(bias)
mean(rpmse)
mean(coverage)
mean(width)


EPA_Measure <- Ozone$X3
CMAQ1 <- Ozone$X4
CMAQ4 <- Ozone$X7
CMAQ7 <- Ozone$X10
CMAQ8 <- Ozone$X11
CMAQ9 <- Ozone$X12
CMAQ10 <- Ozone$X13
fit.lm <- lm(EPA_Measure~CMAQ1+CMAQ4+CMAQ7+CMAQ8+CMAQ9+CMAQ10)

pdf("avplots.pdf")
avPlots(fit.lm,main="")
dev.off()
plot(Ozone$X4,Ozone$X3)




