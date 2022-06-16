#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Introduction  #####################
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Code GSERM "Regression for Publishing"
# (June 2022)
#
# Day Four: Event Counts + GLMs
#
# File created 6/20/2021
#
# File last modified 6/11/2022
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# setwd() first:
#
# setwd("~/Dropbox (Personal)/GSERM/Materials 2020") # or whatever...
#
# Packages:

library(readr)    # <-- install all these as necessary
library(stargazer)# <--
library(plyr)     # <-- 
library(car)      # <-- 
library(rms)      # <-- 
library(plm)      # <-- 
library(lmtest)   # <-- 
library(gvlma)    # <-- 
library(dplyr)    # <-- 
library(mfx)      # <--
library(margins)  # <--

# Options:
options(scipen=12)
options(digits=5)

# setwd()...
#
# setwd("~/Dropbox (Personal)/GSERM/RFP-2021")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Event Count Models ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# Various Poisson histograms

set.seed(7222009)
N<-1000
LP05<-rpois(N,0.5)
LP1<-rpois(N,1)
LP5<-rpois(N,5)
LP10<-rpois(N,10)

pdf("PoissonHistogramsR.pdf",7,6)
par(mfrow=c(2,2))
hist(LP05,col="grey",xlim=c(0,25),breaks=seq(0,25,by=1),
     ylim=c(0,1000),xlab="Count",main="Lambda = 0.5")
hist(LP1,col="grey",xlim=c(0,25),breaks=seq(0,25,by=1),
     ylim=c(0,1000),xlab="Count",main="Lambda = 1.0")
hist(LP5,col="grey",xlim=c(0,25),breaks=seq(0,25,by=1),
     ylim=c(0,1000),xlab="Count",main="Lambda = 5")
hist(LP10,col="grey",xlim=c(0,25),breaks=seq(0,25,by=1),
     ylim=c(0,1000),xlab="Count",main="Lambda = 10")
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SCOTUS nullifications data:

Nulls<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/nulls.csv")

# Histogram:

pdf("NullsHist.pdf",6,5)
par(mar=c(4,4,2,2))
with(Nulls, 
     hist(nulls,main="",xlab="Number of Nullifications",
          col="grey"))
dev.off()

# Poisson regression:

nulls.poisson<-glm(nulls~tenure+unified,family="poisson",
                   data=Nulls)
summary(nulls.poisson)

# IRRs:

nulls.poisson.IRR<-poissonirr(nulls~tenure+unified,
                              data=Nulls)
nulls.poisson.IRR

# Predictions:

tenure<-seq(0,20,1)
unified<-1
simdata<-as.data.frame(cbind(tenure,unified))
nullhats<-predict(nulls.poisson,newdata=simdata,se.fit=TRUE)

# NOTE: These are XBs, not predicted counts.
# Transforming:

nullhats$Yhat<-exp(nullhats$fit)
nullhats$UB<-exp(nullhats$fit + 1.96*(nullhats$se.fit))
nullhats$LB<-exp(nullhats$fit - 1.96*(nullhats$se.fit))

# Plot...

pdf("NullsOutOfSampleHatsR.pdf",6,5)
plot(simdata$tenure,nullhats$Yhat,t="l",lwd=3,ylim=c(0,5),ylab=
       "Predicted Count", xlab="Mean Tenure")
lines(simdata$tenure,nullhats$UB,lwd=2,lty=2)
lines(simdata$tenure,nullhats$LB,lwd=2,lty=2)
dev.off()

# Something similar, via -margins-:

pdf("Nulls-Preds-Margin.pdf",6,5)
par(mar=c(4,4,2,2))
cplot(nulls.poisson,"tenure",xlab="Tenure",
      ylab="Predicted Nullifications")
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Offsets with dyadic data...Aggregated counts
# of conflicts between the countries in each
# dyad, 1950-1985...

IR<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/offsetIR.csv")

summary(IR)

cor(IR,use="complete.obs")

IR.fit1<-glm(disputes~allies+openness,data=IR,family="poisson")
summary(IR.fit1)

IR.fit2<-glm(disputes~allies+openness,data=IR,family="poisson",
             offset=log(Ndyads))
summary(IR.fit2)

IR.fit3<-glm(disputes~allies+openness+log(Ndyads),data=IR,
             family="poisson")
summary(IR.fit3)

# z-test:
2*pnorm((0.811-1)/.071)

# Wald test:
wald.test(b=coef(IR.fit3),Sigma=vcov(IR.fit3),Terms=4,H0=1)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Extra-Poisson Variation / Negative Binomial ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# Simulated Poisson / NB example:
#
# Poisson data (N=400):

N<-400
set.seed(7222009)
X <- runif(N,min=0,max=1)
YPois <- rpois(N,exp(0+1*X))           # Poisson
YNB <- rnbinom(N,size=1,mu=exp(0+1*X)) # NB with theta=1.0

describe(cbind(YPois,YNB))

# Density plots:

pdf("PoisNBDensities.pdf",8,6)
par(mar=c(4,4,2,2))
plot(density(YNB),lwd=2,lty=2,col="red",
     ylim=c(0,0.4),main="",xlab="Y")
lines(density(YPois),lwd=2,lty=1,col="black")
legend("topright",bty="n",lty=c(1,2),col=c("black","red"),
       lwd=c(2,2),legend=c("Poisson","Neg. Bin. (theta=1)"))
dev.off()


# Regressions:

summary(glm(YPois~X,family="poisson")) # Poisson
summary(glm.nb(YPois~X))               # NB

# More regressions:

summary(glm(YNB~X,family="poisson"))  # Poisson
summary(glm.nb(YNB~X))                # NB

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Now, a bigger/better sim:

Sims <- 250 # (250 sims each)
theta <- seq(0.1,4,by=0.1) # values of theta
diffs<-matrix(nrow=Sims,ncol=length(theta))

set.seed(7222009)
for(j in 1:length(theta)) {
  for(i in 1:Sims) {
    X<-runif(N,min=0,max=1)
    Y<-rnbinom(N,size=theta[j],mu=exp(0+1*X))
    p<-glm(Y~X,family="poisson")
    nb<-glm.nb(Y~X)
    diffs[i,j]<- ((sqrt(vcov(p))[2,2]) / sqrt(vcov(nb))[2,2])*100
  }
}

Dmeans<-colMeans(diffs)
Dmin<-apply(diffs,2,min)
Dmax<-apply(diffs,2,max)

pdf("PoissonSEsWithNB.pdf",7,6)
par(mar=c(4,4,2,2))
plot(theta,Dmeans,xlim=c(0,4),ylim=c(0,100),
     pch=20,ylab="Percentage",xlab=expression(theta))
segments(theta,Dmin,theta,Dmax)
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # CPB Max(Y) figure (not shown):
# 
# L.CPB <- seq(0.1,10,by=0.1)
# MaxY8 <- (-L.CPB) / (0.8-1)
# MaxY5 <- (-L.CPB) / (0.5-1)
# MaxY2 <- (-L.CPB) / (0.2-1)
# 
# pdf("CPBMaxR.pdf",6,5)
# par(mfrow=c(1,1))
# par(mar=c(4,4,2,2))
# plot(L.CPB,MaxY8,t="l",lwd=3,col="black",xlab="Lambda",
#      ylab="Maximum Value of Y")
# lines(L.CPB,MaxY5,lwd=3,col="red",lty=2)
# lines(L.CPB,MaxY2,lwd=3,col="darkgreen",lty=3)
# legend("topleft",bty="n",lty=c(1,2,3),lwd=3,
#        col=c("black","red","darkgreen"),
#        legend=c("Alpha = 0.8","Alpha = 0.5","Alpha = 0.2"))
# dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SCOTUS amici data:

amici<-read.csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/Amici.csv")

describe(amici)

amici.poisson<-glm(namici~term+civlibs,data=amici,
                   family="poisson")
summary(amici.poisson)

Phats<-fitted.values(amici.poisson)
Uhats<-((amici$namici-Phats)^2 - amici$namici) / (Phats * sqrt(2))
summary(lm(Uhats~Phats))


amici.NB<-glm.nb(namici~term+civlibs,data=amici)
summary(amici.NB)

# alpha:

1 / amici.NB$theta

# Coefficient estimates:

cbind(amici.poisson$coefficients,coef(amici.NB))

# Estimated standard errors:

cbind(diag(sqrt(vcov(amici.poisson))),diag(sqrt(vcov(amici.NB))))

# Plot:

pdf("PoissonNBYHatsR.pdf",6,6)
par(mar=c(4,4,2,2))
plot(amici.poisson$fitted.values,amici.NB$fitted.values,pch=20,
     xlab="Poisson",ylab="Negative Binomial",main="",
     xlim=c(0,3),ylim=c(0,3))
abline(a=0,b=1,lwd=1,lty=2)
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# GLMS!     ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Toy example

X<-c(1,1,2,2,3,3,4,4,5,5)
Y<-c(0,2,1,3,2,4,3,5,4,6)

linmod<-lm(Y~X)
summary(linmod)
linglm<-glm(Y~X,family="gaussian")
summary(linglm)

# 2008 NES Data: 

NES08<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/NES2008.csv")

summary(NES08[,4:16])

pdf("Notes/PolKnowledge.pdf",5,4)
par(mar=c(4,4,2,2))
barplot(xtabs(~NES08$knowledge),pch=19,lcolor="grey",
        ylab="Frequency",ylim=c(0,600),
        xlab="Knowledge Score")
dev.off()

nes08.binom<-glm(cbind(knowledge,4-knowledge)~age+female+white+
                   oftenvote+conservative+prayfreq+heterosexual+married+
                   yrsofschool+income,data=NES08,family=binomial)
summary(nes08.binom)
