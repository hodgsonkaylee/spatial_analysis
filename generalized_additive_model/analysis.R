##############################################
############## STAT 637 Project ##############
##############################################

# Load packages
library(betareg)
library(lme4)
library(nlme)
library(MASS)
library(mgcv)
library(rworldmap)
library(RColorBrewer)
library(car)
library(xtable)

# Load in data
unempdat <- read.csv("unemploydat.csv",header=TRUE)
unempdat$UnemplWB2016 <- unempdat$UnemplWB2016/100

########################
# Exploratory Analysis #
########################

my.colors <- brewer.pal(7, 'YlGnBu')
my.colors <- brewer.pal(9, 'YlGnBu')
cp <- exp(seq(from=log(.001),to=log(.28),length=10))

MyWorldMap <- joinCountryData2Map(unempdat, joinCode = "ISO3", 
                                  nameJoinColumn = "iso3c",
                                  nameCountryColumn = "Country")
pdf("unempplot.pdf")
mapCountryData(MyWorldMap, colourPalette=my.colors,nameColumnToPlot = "UnemplWB2016", catMethod = cp,
               addLegend = FALSE, borderCol="grey", missingCountryCol = "gray",
               mapRegion = "world", mapTitle = "")
dev.off()

pdf("explorplots.pdf",width=15,height=7)
par(mfrow=c(2,6))
plot(unempdat$MonopolyForce,unempdat$UnemplWB2016,
     pch=20,ylab="Unemployment Rates",xlab="State Monopoly on the Use of Force") 
plot(unempdat$PrivPR2016,unempdat$UnemplWB2016,
     pch=20,ylab="Unemployment Rates",xlab="Private Property Rights")
plot(unempdat$SocialSafetyNets,unempdat$UnemplWB2016,
     pch=20,ylab="Unemployment Rates",xlab="Social Safety Nets")
plot(unempdat$YouthRiskFactor2010,unempdat$UnemplWB2016,
     pch=20,ylab="Unemployment Rates",xlab="Youth Risk Factor")
plot(unempdat$MeanYrsSchooling2015,unempdat$UnemplWB2016,
     pch=20,ylab="Unemployment Rates",xlab="Average Years of Schooling")
plot(unempdat$LifeExpect2016,unempdat$UnemplWB2016,
     pch=20,ylab="Unemployment Rates",xlab="Life Expectancy")
plot(unempdat$DscrmViolMin2016,unempdat$UnemplWB2016,
     pch=20,ylab="Unemployment Rates",xlab="Discrimination and Violence against Minorities")
plot(unempdat$GovFrameGE2015,unempdat$UnemplWB2016,
     pch=20,ylab="Unemployment Rates",xlab="Government Framework for Gender Equality")
plot(unempdat$GDP2017,unempdat$UnemplWB2016,
     pch=20,ylab="Unemployment Rates",xlab="GDP per Capita")
plot(as.factor(unempdat$SVS2014),unempdat$UnemplWB2016,
     ylab="Unemployment Rates",xlab="Societal Violence Scale") # factor
plot(as.factor(unempdat$RegimeTypes2017),unempdat$UnemplWB2016,
     ylab="Unemployment Rates",xlab="Regime Types") # factor
dev.off()

##############################################################
# Original Beta and Logistic Models with and without Weights #
##############################################################

binmod <- glm(cbind(round(UnemplWB2016,2)*100,100)~
                 MonopolyForce+as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                 GDP2017+PrivPR2016+SocialSafetyNets+YouthRiskFactor2010+
                 MeanYrsSchooling2015+LifeExpect2016+DscrmViolMin2016+GovFrameGE2015,
               family=binomial,data=unempdat)
binmodw <- glm(cbind(round(UnemplWB2016*Population),Population)~
                 MonopolyForce+as.factor(SVS2014)+as.factor(RegimeTypes2017)+                 
                 GDP2017+PrivPR2016+SocialSafetyNets+YouthRiskFactor2010+
                 MeanYrsSchooling2015+LifeExpect2016+DscrmViolMin2016+GovFrameGE2015,
            family=binomial,data=unempdat)
betmod <- betareg(UnemplWB2016~
                    MonopolyForce+as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                    GDP2017+PrivPR2016+SocialSafetyNets+YouthRiskFactor2010+
                    MeanYrsSchooling2015+LifeExpect2016+DscrmViolMin2016+GovFrameGE2015,
                  data=unempdat)
betmodw <- betareg(UnemplWB2016~
                     MonopolyForce+as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                     GDP2017+PrivPR2016+SocialSafetyNets+YouthRiskFactor2010+
                     MeanYrsSchooling2015+LifeExpect2016+DscrmViolMin2016+GovFrameGE2015,
                  weights=Population,data=unempdat)
# diagnostics
diagnost <- data.frame(rbind(AIC(binmod,binmodw,betmod,betmodw)[,2],
                             BIC(binmod,binmodw,betmod,betmodw)[,2],
                             c(sum(resid(binmod)^2),sum(resid(binmodw)^2),sum(resid(betmod)^2),sum(resid(betmodw)^2)),
                             c(mean(binmod$residuals^2),mean(binmodw$residuals^2),mean(betmod$residuals^2),mean(betmodw$residuals^2))
                             ))
rownames(diagnost) <- c("AIC","BIC","Deviance","MSE")
colnames(diagnost) <- c("Binomial","Binomial (Weighted)","Beta","Beta (Weighted)")
xtable(diagnost)

pdf("residplots.pdf",width=12,height=7)
par(mfrow=c(2,2))
plot(predict(binmod),residuals(binmod),main="Binomial Model (Equal Weights)",
     pch=20,ylab="Standardized Residuals",xlab="Fitted Values")
plot(predict(binmodw),residuals(binmodw,type="response")/sd(residuals(binmodw,type="response")),main="Binomial Model (Weight by Population)",
     pch=20,ylab="Standardized Residuals",xlab="Fitted Values")
plot(predict(betmod),residuals(betmod),main="Beta Model (Equal Weights)",
     pch=20,ylab="Standardized Residuals",xlab="Fitted Values")
plot(predict(betmodw),residuals(betmodw,type="response")/sd(residuals(betmodw,type="response")),main="Beta Model (Weight by Population)",
     pch=20,ylab="Standardized Residuals",xlab="Fitted Values")
dev.off()
# Use the beta regression model without the population weights

######################################
# Beta Model with Spatial Dependence #
######################################

# https://www.casact.org/education/rpm/2010/handouts/PM1-Guszcza.pdf
# https://stats.stackexchange.com/questions/232697/beta-regression-accounting-for-residual-spatial-auto-correlation-in-r
spacmod <- gam(UnemplWB2016~
                 MonopolyForce+as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                 GDP2017+PrivPR2016+SocialSafetyNets+YouthRiskFactor2010+
                 MeanYrsSchooling2015+LifeExpect2016+DscrmViolMin2016+GovFrameGE2015+
                 s(Latitude, Longitude, bs = "gp", m = 2), 
               family = betar(link='logit'), 
               data = unempdat)
# diagnostics
sdiagnost <- data.frame(rbind(AIC(betmod,spacmod)[,2],
                             BIC(betmod,spacmod)[,2],
                             c(sum(resid(betmod)^2),sum(resid(spacmod)^2))
))
rownames(sdiagnost) <- c("AIC","BIC","Deviance")
colnames(sdiagnost) <- c("Beta","Spacial Beta")
xtable(sdiagnost)

pdf("sresidplots.pdf",width=16,height=6)
par(mfrow=c(1,2))
plot(predict(betmod,type="response"),residuals(betmod),main="Beta Model (without Spacial Component)",
     pch=20,ylab="Standardized Residuals",xlab="Fitted Values")
plot(predict(spacmod,type="response"),residuals(spacmod),main="Beta Model (with Spacial Component)",
     pch=20,ylab="Standardized Residuals",xlab="Fitted Values")
dev.off()

#squared correlation between between the linear predictor for the mean and the link-transformed response
plot(density(unempdat$UnemplWB2016,na.rm=T),ylim=c(0,15))
lines(density(predict(betmod,type="response")))
lines(density(predict(spacmod,type="response")))

pdf("modcomp.pdf",width=10,height=7)
hist(unempdat$UnemplWB2016,freq=F,ylim=c(0,14),
     main="Worldwide Unemployment Rates",xlab="Unemployment Rates")
lines(density(unempdat$UnemplWB2016,na.rm=T))
lines(density(predict(spacmod,type="response")),col="blue",lty=2)
lines(density(predict(betmod,type="response")),col="red",lty=2)
lines(density(predict(spacmod.fin,type="response")),col="purple",lwd=2)
legend("topright",legend=c("Observed","Fitted without Spatial Aspect","Fitted with Spatial Aspect","Fitted with Spatial (after VS)"),
       lty=c(1,2,2,1),lwd=c(1,1,1,2),col=c("black","red","blue","purple"))
dev.off()

######################
# Variable Selection #
######################

spacmod <- gam(UnemplWB2016~
                 MonopolyForce+as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                 GDP2017+PrivPR2016+SocialSafetyNets+YouthRiskFactor2010+
                 MeanYrsSchooling2015+LifeExpect2016+DscrmViolMin2016+GovFrameGE2015+
                 s(Latitude, Longitude, bs = "gp", m = 2), 
               family = betar(link='logit'), 
               data = unempdat) # take out MonopolyForce
spacmod1 <- gam(UnemplWB2016~
                 as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                GDP2017+PrivPR2016+SocialSafetyNets+YouthRiskFactor2010+
                 MeanYrsSchooling2015+LifeExpect2016+DscrmViolMin2016+GovFrameGE2015+
                 s(Latitude, Longitude, bs = "gp", m = 2), 
               family = betar(link='logit'), 
               data = unempdat) # take out DscrmViolMin2016
spacmod2 <- gam(UnemplWB2016~
                  as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                  GDP2017+PrivPR2016+SocialSafetyNets+YouthRiskFactor2010+
                  MeanYrsSchooling2015+LifeExpect2016+GovFrameGE2015+
                  s(Latitude, Longitude, bs = "gp", m = 2), 
                family = betar(link='logit'), 
                data = unempdat) # take out PrivPR
spacmod3 <- gam(UnemplWB2016~
                  as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                  GDP2017+SocialSafetyNets+YouthRiskFactor2010+
                  MeanYrsSchooling2015+LifeExpect2016+GovFrameGE2015+
                  s(Latitude, Longitude, bs = "gp", m = 2), 
                family = betar(link='logit'), 
                data = unempdat) # take out LifeExpect2016
spacmod4 <- gam(UnemplWB2016~
                  as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                  GDP2017+SocialSafetyNets+YouthRiskFactor2010+
                  MeanYrsSchooling2015+GovFrameGE2015+
                  s(Latitude, Longitude, bs = "gp", m = 2), 
                family = betar(link='logit'), 
                data = unempdat) # take out GovFrameGE2015
spacmod5 <- gam(UnemplWB2016~
                  as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                  GDP2017+SocialSafetyNets+YouthRiskFactor2010+
                  MeanYrsSchooling2015+
                  s(Latitude, Longitude, bs = "gp", m = 2), 
                family = betar(link='logit'), 
                data = unempdat) # take out MeanYrsSchooling2015
spacmod6 <- gam(UnemplWB2016~
                  as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                  GDP2017+SocialSafetyNets+YouthRiskFactor2010+
                  s(Latitude, Longitude, bs = "gp", m = 2), 
                family = betar(link='logit'), 
                data = unempdat) # take out YouthRiskFactor2010
spacmod7 <- gam(UnemplWB2016~
                  as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                  GDP2017+SocialSafetyNets+
                  s(Latitude, Longitude, bs = "gp", m = 2), 
                family = betar(link='logit'), 
                data = unempdat) # take out SocialSafetyNets
 
aic.vals <- AIC(spacmod,spacmod1,spacmod2,spacmod3,spacmod4,spacmod5,spacmod6,spacmod7)[,2]
bic.vals <- BIC(spacmod,spacmod1,spacmod2,spacmod3,spacmod4,spacmod5,spacmod6,spacmod7)[,2]
dev.vals <- c(sum(resid(spacmod)^2),sum(resid(spacmod1)^2),sum(resid(spacmod2)^2),
  sum(resid(spacmod3)^2),sum(resid(spacmod4)^2),sum(resid(spacmod5)^2),
  sum(resid(spacmod6)^2),sum(resid(spacmod7)^2))
adjr2 <- c(0.661,0.666,0.679,0.684,0.688,0.685,0.695,0.715)

pdf("aicbicplots.pdf")
plot(aic.vals,type="l",lty=1,lwd=2,
     ylab="Diagnostic",xlab="Models from Steps in Backward Selection",
     ylim=c(-610,-300),col="blue",
     main="")
lines(bic.vals,lty=1,lwd=2,col="green")
title(expression(phantom("AIC ") * "and " * phantom("BIC ") * "Diagnostics in Backward Selection"), col.main = "black")
title(expression(phantom("AIC and ") * "BIC " * phantom("Diagnostics in Backward Selection")), col.main = "green")
title(expression("AIC " * phantom("and BIC Diagnostics in Backward Selection")), col.main = "blue")
dev.off()

pdf("devianceplot.pdf",width=7,height=5)
plot(dev.vals,type="l",lty=1,lwd=3,
     ylab="Residual Deviance",xlab="Models from Sequential Steps",
     col="blue",
     main="Deviance in Backward Selection")
dev.off()
pdf("adjr2plot.pdf",width=7,height=5)
plot(adjr2,type="l",lty=1,lwd=3,
     ylab="Adjusted R-Squared",xlab="Models from Sequential Steps",
     col="blue",
     main="Adjusted R-Squared in Backward Selection")
dev.off()

###########################
# Final Model Diagnostics #
###########################

# restrict dataframe to exclude rows with missing values
isna <- apply(unempdat[,c(2,4,5,7,13,14,15)],1,function(x) sum(!is.na(x)))
unempdat.fin <- unempdat[isna==7,]

# final model
spacmod.fin <- gam(UnemplWB2016~
                     as.factor(SVS2014)+as.factor(RegimeTypes2017)+
                     GDP2017+SocialSafetyNets+
                     s(Latitude, Longitude, bs = "gp", m = 2), 
                   family = betar(link='logit'), 
                   data = unempdat.fin)

# add predicted rates to dataframe
unempdat.fin$Predict <- predict(spacmod.fin,type="response")

# plot observed vs. predicted
my.colors <- brewer.pal(9, 'YlGnBu')
MyWorldMap.fin <- joinCountryData2Map(unempdat.fin, joinCode = "ISO3", 
                                  nameJoinColumn = "iso3c",
                                  nameCountryColumn = "Country")
cp <- exp(seq(from=log(.001),to=log(.28),length=10))

pdf("observedplot.pdf")
mapCountryData(MyWorldMap.fin, colourPalette=my.colors,nameColumnToPlot = "UnemplWB2016", catMethod = cp,
               addLegend = FALSE, borderCol="grey", missingCountryCol = "gray",
               mapRegion = "world", mapTitle = "")
dev.off()
pdf("fittedplot.pdf")
mapCountryData(MyWorldMap.fin, colourPalette=my.colors,nameColumnToPlot = "Predict", catMethod = cp,
               addLegend = FALSE, borderCol="grey", missingCountryCol = "gray",
               mapRegion = "world", mapTitle = "")
dev.off()

summary(spacmod.fin)
# 1 ”Pure Autocracy”: both personal autonomy rights and political participation rights be-low the scale midpoint (0.50)
# 2 ”Inclusive Autocracy”: personal autonomy rights below the scalemidpoint, political participation rights above the scale midpoint
# 3 ”Liberal Autocracy”: personal au-tonomy rights above the scale midpoint, political participation rights below
# 4 ”Minimal Democracy”:both personal autonomy rights and political participation rights above the scale midpoint





