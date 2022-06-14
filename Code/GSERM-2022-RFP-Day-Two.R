#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Preliminaries ##################################
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Code for GSERM "Regression for Publishing"
# (June 2022)
#
# Day Two
#
# File created 6/21/2020
#
# File last modified 6/10/2022
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Set working directory first:
#
# setwd("~/Dropbox (Personal)/GSERM/RFP-2022/Slides") # or whatever...

library(readr)    # <-- install these packages 
library(car)      # <-- as necessary...
library(lattice)  # <-- 
library(arm)      # <-- 
library(plotMElm) # <--
library(stargazer)# <-- 

options(scipen = 6) # bias against scientific notation
options(digits = 6) # show fewer decimal places

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Multiplicative interactions... ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Simulations for pictures....
#
# Two dummy predictors:

set.seed(7222009)
N<-400
D1<-rep(c(0,1),times=N/2)
D2<-rep(c(0,0,1,1),times=N/4)
Y <- rnorm(N,(20-10*D2+10*D1+20*D1*D2),2)
df<-data.frame(D1=D1,D2=D2,Y=Y)

pdf("TwoDummyBoxPlotsRD1.pdf",6,6)
par(mar=c(4,4,2,2))
with(df, boxplot(Y~D2+D1,xaxt="n",xlab="Values of D1,D2"))
axis(1,at=c(1,2,3,4),
     labels=c("D1=0, D2=0","D1=0, D2=1",
              "D1=1, D2=0","D1=1, D2=1"))
arrows(1,median(df$Y[which(df$D1==0 & df$D2==0)]),
       3,median(df$Y[which(df$D1==1 & df$D2==0)]),
       lwd=2,length=0.10,col="red")
arrows(2,median(df$Y[which(df$D1==0 & df$D2==1)]),
       4,median(df$Y[which(df$D1==1 & df$D2==1)]),
       lwd=2,length=0.10,col="red")
legend("topleft",bty="n",legend="E(Y) | change in D1",col="red")
dev.off()

pdf("TwoDummyBoxPlotsRD2.pdf",6,6)
par(mar=c(4,4,2,2))
with(df, boxplot(Y~D2+D1,xaxt="n",xlab="Values of D1,D2"))
axis(1,at=c(1,2,3,4),
     labels=c("D1=0, D2=0","D1=0, D2=1",
              "D1=1, D2=0","D1=1, D2=1"))
arrows(1,median(df$Y[which(df$D1==0 & df$D2==0)]),
       2,median(df$Y[which(df$D1==0 & df$D2==1)]),
       lwd=2,length=0.10)
arrows(3,median(df$Y[which(df$D1==1 & df$D2==0)]),
       4,median(df$Y[which(df$D1==1 & df$D2==1)]),
       lwd=2,length=0.10)
legend("topleft",bty="n",legend="E(Y) | change in D2")
dev.off()

# Dummy + continuous:

set.seed(7222009)
N<-200
D<-rep(c(0,1),times=N/2)
X<-rnorm(N,0,5)
color<-ifelse(D==0,"black","red")
df2<-data.frame(D=D,X=X,color=color,
                stringsAsFactors=FALSE)

df2$Y1 <- 50+2*df2$X+3*rnorm(N)
df2$Y2 <- 50+2*df2$X+30*df2$D+3*rnorm(N)
df2$Y3 <- 50+2*df2$X-(4*df2$D*df2$X)+3*rnorm(N)
df2$Y4 <- 50+2*df2$X+30*df2$D-(4*df2$D*df2$X)+3*rnorm(N)

pdf("ScatterInterSameR.pdf",6,6)
par(mar=c(4,4,2,2))
with(df2, plot(X,Y1,pch=D+16,col=color,ylab="Y"))
legend("topleft",bty="n",legend=c("D=0","D=1"),
       pch=c(16,17),col=c("black","red"))
abline(with(df2[df2$D==0,],lm(Y1~X)),col="black",lwd=2)
abline(with(df2[df2$D==1,],lm(Y1~X)),col="red",lwd=2)
abline(v=0,lty=2)
dev.off()

pdf("ScatterInterInterceptR.pdf",6,6)
par(mar=c(4,4,2,2))
with(df2, plot(X,Y2,pch=D+16,col=color,ylab="Y"))
legend("topleft",bty="n",legend=c("D=0","D=1"),
       pch=c(16,17),col=c("black","red"))
abline(with(df2[df2$D==0,],lm(Y2~X)),col="black",lwd=2)
abline(with(df2[df2$D==1,],lm(Y2~X)),col="red",lwd=2)
abline(v=0,lty=2)
dev.off()

pdf("ScatterInterSlopeR.pdf",6,6)
par(mar=c(4,4,2,2))
with(df2, plot(X,Y3,pch=D+16,col=color,ylab="Y"))
legend("topright",bty="n",legend=c("D=0","D=1"),
       pch=c(16,17),col=c("black","red"))
abline(with(df2[df2$D==0,],lm(Y3~X)),col="black",lwd=2)
abline(with(df2[df2$D==1,],lm(Y3~X)),col="red",lwd=2)
abline(v=0,lty=2)
dev.off()

pdf("ScatterInterBothR.pdf",6,6)
par(mar=c(4,4,2,2))
with(df2, plot(X,Y4,pch=D+16,col=color,ylab="Y"))
legend("topright",bty="n",legend=c("D=0","D=1"),
       pch=c(16,17),col=c("black","red"))
abline(with(df2[df2$D==0,],lm(Y4~X)),col="black",lwd=2)
abline(with(df2[df2$D==1,],lm(Y4~X)),col="red",lwd=2)
abline(v=0,lty=2)
dev.off()

# Two continuous: Wireframe plots...

df3<-expand.grid(X1=seq(0,10,1),
                 X2=seq(0,10,1))
df3$YNoInt<-10 + 2*df3$X1 + 2*df3$X2
df3$YInt  <-(10 - 2*df3$X1 - 2*df3$X2 + 4*df3$X1*df3$X2)/5

trellis.par.set("axis.line",list(col="transparent"))

pdf("TwoContinuousNoInteractiveR.pdf",6,6)
par(mar=c(4,4,2,2))
with(df3, wireframe(YNoInt~X1+X2,
                    drape=TRUE,
                    xlab=list("X1",rot=30),
                    ylab=list("X2",rot=-40),
                    zlab=list("Y",rot=90),
                    scales=list(arrows=FALSE,col="black"),
                    zoom=0.85,pretty=TRUE,
                    col.regions=colorRampPalette(c("blue","red"))(100)))
dev.off()

pdf("TwoContinuousInteractiveR.pdf",6,6)
par(mar=c(4,4,2,2))
with(df3, wireframe(YInt~X1+X2,
                    drape=TRUE,
                    xlab=list("X1",rot=30),
                    ylab=list("X2",rot=-40),
                    zlab=list("Y",rot=90),
                    scales=list(arrows=FALSE,col="black"),
                    zoom=0.85,pretty=TRUE,
                    col.regions=colorRampPalette(c("blue","red"))(100)))
dev.off()

# Polynomials...

N<-200
set.seed(7222009)
df4 <- data.frame(X = runif(N,-5,5))
df4$Y2A <- 10 + 1*df4$X - 5*(df4$X^2) + rnorm(N,0,10) # Quad #1
df4$Y2B <- -50 - 1*df4$X + 3*(df4$X^2) + rnorm(N,0,10) # Quad #2
df4$Y3 <- -8 - 6*df4$X + 3*(df4$X^2) + 1*(df4$X^3) + rnorm(N,0,10) # Cubic
df4 <- df4[order(df4$X),]
fitA<-with(df4, lm(Y2A~X+I(X^2)))
fitB<-with(df4, lm(Y2B~X+I(X^2)))
fit3<-with(df4, lm(Y3~X+I(X^2)+I(X^3)))

pdf("TwoQuadraticsR.pdf",7,6)
par(mar=c(4,4,2,2))
with(df4, plot(X,Y2B,pch=16,col="black",ylab="Y",
               ylim=c(min(df4$Y2B),max(df4$Y2A))))
points(df4$X,df4$Y2A,pch=17,col="red")
lines(df4$X,fitted(fitA),lwd=3,col="red")
lines(df4$X,fitted(fitB),lwd=3,col="black")
dev.off()

pdf("CubicR.pdf",7,6)
par(mar=c(4,4,2,2))
with(df4, plot(X,Y3,pch=16,col="black",ylab="Y"))
lines(df4$X,fitted(fit3),lwd=3,col="black")
dev.off()

# Three-way interaction sim:

N <- 100
X <- runif(N,-5,5)
df00<-data.frame(X=X)
df01<-data.frame(X=X)
df10<-data.frame(X=X)
df11<-data.frame(X=X)
set.seed(7222009)
df00$Y<-0.5*df00$X+rnorm(N)
df01$Y<-2*df01$X+rnorm(N)
df10$Y<- -0.5*df10$X+rnorm(N)
df11$Y<- -2*df11$X+rnorm(N)
fit00<-lm(Y~X,data=df00)
fit01<-lm(Y~X,data=df01)
fit10<-lm(Y~X,data=df10)
fit11<-lm(Y~X,data=df11)
hi<-12
lo<- -12

pdf("TwoDummyOneContinuousR.pdf",7,6)
par(mfrow=c(2,2))
par(mar=c(4,4,4,2))
with(df00, plot(X,Y,main="D1=0,D2=0",ylab="Y",
                ylim=c(lo,hi),pch=20))
abline(h=0,lty=2)
abline(reg=fit00,lwd=2)
with(df01, plot(X,Y,main="D1=0,D2=1",ylab="Y",
                ylim=c(lo,hi),pch=20))
abline(h=0,lty=2)
abline(reg=fit01,lwd=2)
with(df10, plot(X,Y,main="D1=1,D2=0",ylab="Y",
                ylim=c(lo,hi),pch=20))
abline(h=0,lty=2)
abline(reg=fit10,lwd=2)
with(df11, plot(X,Y,main="D1=1,D2=1",ylab="Y",
                ylim=c(lo,hi),pch=20))
abline(h=0,lty=2)
abline(reg=fit11,lwd=2)
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Clinton thermometer score example... ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

ClintonTherm<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/ClintonTherm.csv")

summary(ClintonTherm)

fit0<-lm(ClintonTherm~RConserv+GOP,data=ClintonTherm)

# Coefficient / ladder plot:

pdf("ClintonCoefPlot.pdf",7,6)
par(mar=c(2,14,2,2))
coefplot(fit0,main="Estimated Coefficient",cex.var=1.1,
         varnames=c("(Intercept)",
                    "Respondent's Conservatism","GOP"))
dev.off()

# Model w/Interaction:

fit1<-with(ClintonTherm, lm(ClintonTherm~RConserv+GOP+RConserv*GOP))
summary(fit1)

pdf("ClintonInterCoefPlot.pdf",7,6)
par(mar=c(2,16,2,2))
coefplot(fit1,main="Estimated Coefficient",cex.var=1.1,
         varnames=c("(Intercept)",
                    "Respondent's Conservatism","GOP",
                    "Respondent's Conservatism x GOP"))
dev.off()

# Plot of thermometer scores vs. conservatism:

pdf("ClinThermScatterR.pdf",6,6)
scatterplot(ClintonTherm$ClintonTherm~ClintonTherm$RConserv|as.factor(ClintonTherm$GOP),
            legend.plot=FALSE,
            xlab="Respondent Conservatism",
            ylab="Clinton Thermometer Score",
            smooth=FALSE,boxplots=FALSE,
            pch=c(4,16),col=c("red","blue","red","blue"),
            lwd=2,grid=FALSE)
dev.off()

# Separate regressions:

NonReps<-subset(ClintonTherm,GOP==0)
summary(with(NonReps, lm(ClintonTherm~RConserv)))

Reps<-subset(ClintonTherm,GOP==1)
summary(with(Reps, lm(ClintonTherm~RConserv)))


# psi_1:
Psi1<-fit1$coeff[2]+fit1$coeff[4]
Psi1
SPsi1<-sqrt(vcov(fit1)[2,2] + (1)^2*vcov(fit1)[4,4] + 2*1*vcov(fit1)[2,4])
SPsi1
Psi1 / SPsi1 # <-- t-statistic

# psi_2 | RConserv = 1

fit1$coeff[3]+(1 * fit1$coeff[4])

sqrt(vcov(fit1)[3,3] + (1)^2*vcov(fit1)[4,4] + 2*1*vcov(fit1)[3,4])

# Implies t is approximately 2


# psi_2 | RConserv = 7

fit1$coeff[3]+(7 * fit1$coeff[4])

sqrt(vcov(fit1)[3,3] + (7)^2*vcov(fit1)[4,4] + 2*7*vcov(fit1)[3,4])

# t is approximately 11

# Using linearHypothesis:

linearHypothesis(fit1,"RConserv+RConserv:GOP")

# Note: Same as t-test:
sqrt(72.99)

# psi_2 | RConserv = 7:

linearHypothesis(fit1,"GOP+7*RConserv:GOP")

# MFX / psi plots:

ConsSim<-seq(1,7,1)
psis<-fit1$coeff[3]+(ConsSim * fit1$coeff[4])
psis.ses<-sqrt(vcov(fit1)[3,3] + 
                 (ConsSim)^2*vcov(fit1)[4,4] + 2*ConsSim*vcov(fit1)[3,4])

pdf("ClinMFX1.pdf",7,6)
par(mar=c(4,4,2,2))
plot(ConsSim,psis,t="l",lwd=2,xlab="Respondent Conservatism",
     ylab="Estimated Marginal Effect of GOP",ylim=c(-40,0))
lines(ConsSim,psis+(1.96*psis.ses),lty=2,lwd=2)
lines(ConsSim,psis-(1.96*psis.ses),lty=2,lwd=2)
abline(h=0,lwd=1,lty=2)
dev.off()

# Same thing, using plot_me:

pdf("ClinMFX1Alt.pdf",7,6)
plot_me(fit1,"GOP","RConserv",ci_type="fdr")
dev.off()

# Interacting two continuous covariates:

fit2<-with(ClintonTherm,
           lm(ClintonTherm~RConserv+ClintonConserv+RConserv*ClintonConserv))
summary(fit2)

# Hypothesis tests:

fit2$coef[2]+(1*fit2$coef[4])
sqrt(vcov(fit2)[2,2] + (1)^2*vcov(fit2)[4,4] + 2*1*vcov(fit2)[2,4])

linearHypothesis(fit2,"RConserv+1*RConserv:ClintonConserv")

# More hypothesis tests:

# psi_1 | ClintonConserv = mean
fit2$coef[2]+((mean(ClintonTherm$ClintonConserv))*fit2$coef[4])
sqrt(vcov(fit2)[2,2] + (mean(ClintonTherm$ClintonConserv)^2*vcov(fit2)[4,4] +
                          2*(mean(ClintonTherm$ClintonConserv))*vcov(fit2)[2,4]))
pt(((fit2$coef[2]+(2.985*fit2$coef[4])) / sqrt(vcov(fit2)[2,2] + 
                                                 (2.985)^2*vcov(fit2)[4,4] + 2*2.985*vcov(fit2)[2,4])),df=1293)

# psi_2 | RConserv = 1
fit2$coef[3]+(1*fit2$coef[4])

# psi_2 | RConserv = 6
fit2$coef[3]+(6*fit2$coef[4])

# Marginal Effect Plot II:

psis2<-fit2$coef[3]+(ConsSim*fit2$coef[4])
psis2.ses<-sqrt(vcov(fit2)[3,3] + (ConsSim)^2*vcov(fit2)[4,4]
                + 2*ConsSim*vcov(fit2)[3,4])

pdf("ClinMFX2.pdf",6,6)
plot(ConsSim,psis2,t="l",lwd=2,xlab="Respondent's Conservatism",
     ylab="Marginal Effect of Clinton's Conservatism",ylim=c(-10,20))
lines(ConsSim,psis2+(1.96*psis2.ses),lty=2,lwd=2)
lines(ConsSim,psis2-(1.96*psis2.ses),lty=2,lwd=2)
abline(h=0,lty=2,lwd=1,col="red")
dev.off()

# Same thing, using plot_me:

pdf("ClinMFX2Alt.pdf",7,6)
plot_me(fit2,"ClintonConserv","RConserv",ci_type="fdr")
dev.off()

# Contour Plot:

grid<-expand.grid(RConserv=seq(1,7,1),
                  ClintonConserv=seq(1,7,1))
hats<-predict(fit2,newdata=grid)

pdf("ClinContour.pdf",6,6)
levelplot(hats~grid$RConserv*grid$ClintonConserv,
          contour=TRUE,
          cuts=12,pretty=TRUE,xlab="Respondent's Conservatism",
          ylab="Clinton's Conservatism",
          col.regions=heat.colors)
dev.off()

# Wireframe plot:

trellis.par.set("axis.line",list(col="transparent"))

pdf("ClinWireframe.pdf",7,7)
wireframe(hats~grid$RConserv*grid$ClintonConserv,
          drape=TRUE,
          xlab=list("Respondent's Conservatism",rot=30),
          ylab=list("Clinton's Conservatism",
                    rot=-40),zlab=list("Predictions",rot=90),
          scales=list(arrows=FALSE,col="black"),
          zoom=0.85,pretty=TRUE,
          col.regions=colorRampPalette(c("blue","red"))(100))
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Nonlinearity... ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# "Bulging plot":

N<-500
r<-50
set.seed(7222009)
X<-runif(N,-r,r)
TL<-sqrt(r^2-X^2)+rnorm(N)
BL<- -sqrt(r^2-X^2)+rnorm(N)

pdf("BetterBulgingPlot.pdf",10,8)
par(mfrow=c(2,2))
par(mar=c(4,4,2,2))
plot(X[X<0],TL[X<0],pch=20,xlab="X",ylab="Y",
     cex.axis=0.01)
legend("bottomright",bty="n",legend="Log X and/or Square Y")
plot(X[X>0],TL[X>0],pch=20,xlab="X",ylab="Y",
     cex.axis=0.01)
legend("bottomleft",bty="n",legend="Square X and/or Y")
plot(X[X<0],BL[X<0],pch=20,xlab="X",ylab="Y",
     cex.axis=0.01)
legend("topright",bty="n",legend="Log X and/or Y")
plot(X[X>0],BL[X>0],pch=20,xlab="X",ylab="Y",
     cex.axis=0.01)
legend("topleft",bty="n",legend="Square X and/or Log Y")
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Military Spending Example  #   #######
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

Data<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/MilSpend.csv")

Data$MilitarySpending<-Data$MilitarySpending/1000 # Convert
Data$GDP<-Data$GDP/1000      # to thousands of $...

summary(Data)

with(Data, describe(MilitarySpending))
with(Data, describe(GDP))

pdf("NewMilSpendUnivariates.pdf",11,6)
par(mfrow=c(1,2))
par(mar=c(4,4,2,2))
plot(density(Data$MilitarySpending,na.rm=TRUE),main="",
     xlab="Military Spending",lwd=2)
plot(density(Data$GDP,na.rm=TRUE),main="",
     xlab="GDP",lwd=2)
dev.off()

# Ladder Plots in R:

pdf("LadderOfPowersDensitiesGDP-R2.pdf",10,10)
par(mfrow=c(3,3))
with(Data, plot(density(GDP^3,na.rm=TRUE),main="Cubic"))
with(Data, plot(density(GDP^2,na.rm=TRUE),main="Square"))
with(Data, plot(density(GDP,na.rm=TRUE),main="Identity"))
with(Data, plot(density(sqrt(GDP),na.rm=TRUE),main="Square Root"))
with(Data, plot(density(log(GDP),na.rm=TRUE),main="Log"))
with(Data, plot(density(1/sqrt(GDP),na.rm=TRUE),main="1 / Square Root"))
with(Data, plot(density(1/GDP,na.rm=TRUE),main="Inverse"))
with(Data, plot(density(1/GDP^2,na.rm=TRUE),main="1 / Square"))
with(Data, plot(density(1/GDP^3,na.rm=TRUE),main="1 / Cubic"))
dev.off()

pdf("LadderOfPowersDensitiesMilSpend-R2.pdf",10,10)
par(mfrow=c(3,3))
with(Data, plot(density(MilitarySpending^3,na.rm=TRUE),main="Cubic"))
with(Data, plot(density(MilitarySpending^2,na.rm=TRUE),main="Square"))
with(Data, plot(density(MilitarySpending,na.rm=TRUE),main="Identity"))
with(Data, plot(density(sqrt(MilitarySpending),na.rm=TRUE),main="Square Root"))
with(Data, plot(density(log(MilitarySpending+1),na.rm=TRUE),main="Log"))
with(Data, plot(density(1/sqrt(MilitarySpending),na.rm=TRUE),main="1 / Square Root"))
with(Data, plot(density(1/MilitarySpending,na.rm=TRUE),main="Inverse"))
with(Data, plot(density(1/MilitarySpending^2,na.rm=TRUE),main="1 / Square"))
with(Data, plot(density(1/MilitarySpending^3,na.rm=TRUE),main="1 / Cubic"))
dev.off()

# Other plots:

pdf("MilSpendScattersR2.pdf",10,10)
par(mfrow=c(2,2))
par(mar=c(4,4,4,2))
with(Data, plot(GDP,MilitarySpending,pch=20,main="Linear-Linear"))
with(Data, plot(log(GDP),MilitarySpending,pch=20,main="Linear-Log"))
with(Data, plot(GDP,log(MilitarySpending+1),pch=20,main="Log-Linear"))
with(Data, plot(log(GDP),log(MilitarySpending+1),pch=20,main="Log-Log"))
dev.off()

# Regressions:

linlin <- with(Data, lm(MilitarySpending~GDP))
summary(linlin)
linlog <- with(Data, lm(MilitarySpending~log(GDP)))
summary(linlog)
loglin <- with(Data, lm(log(MilitarySpending+1)~GDP))
summary(loglin)
loglog <- with(Data, lm(log(MilitarySpending+1)~log(GDP)))
summary(loglog)

# Nice table:

stargazer(linlin,linlog,loglin,loglog,
          type="latex",
          title="OLS Models of Military Spending",
          dep.var.caption="",
          dep.var.labels=c("Linear Y","Logged Y"),
          column.labels=c("Linear-Linear","Linear-Log",
                           "Log-Linear","Log-Log"),
          digits=2,
          no.space=TRUE,
          intercept.top=TRUE,intercept.bottom=FALSE,
          covariate.labels=c("(Constant)","GDP", "log(GDP)"))

# Residual plot:

pdf("MilSpendGDPResidsDensities-R2.pdf",7,6)
par(mfrow=c(2,2))
par(mar=c(4,4,4,2))
plot(density(linlin$residuals),lwd=2,lty=1,col="black",
     main="Linear-Linear",xlab="Residual Values")
abline(v=0,lty=2,lwd=1)
plot(density(linlog$residuals),lwd=2,lty=2,
     col="blue",main="Linear-Log",xlab="Residual Values")
abline(v=0,lty=2,lwd=1)
plot(density(loglin$residuals),lwd=2,lty=3,col="red",
     main="Log-Linear",xlab="Residual Values")
abline(v=0,lty=2,lwd=1)
plot(density(loglog$residuals),lwd=2,lty=4,col="darkgreen",
     main="Log-Log",xlab="Residual Values")
abline(v=0,lty=2,lwd=1)
dev.off()

# \fin
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=####
