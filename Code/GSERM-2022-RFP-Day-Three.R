#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Introduction      ##################
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Code GSERM "Regression for Publishing"
#
# June 2022
#
# Day Three: various things...
#
# File created 6/21/2021
#
# File last modified 6/10/2022
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# setwd() first:
#
# setwd("~/Dropbox (Personal)/GSERM/Materials 2020") # or whatever...

library(readr)    # <-- install all these as necessary
library(stargazer)# <--
library(plyr)     # <-- 
library(car)      # <-- 
library(rms)      # <-- 
library(plm)      # <-- 
library(lmtest)   # <-- 
library(gvlma)    # <-- 
library(dplyr)    # <-- 
library(gmodels)  # <--
library(margins)  # <--
library(ROCR)     # <--
library(pROC)     # <--

# Handy "robust" summary function for lm:

url_robust <- "https://raw.githubusercontent.com/IsidoreBeautrelet/economictheoryblog/master/robust_summary.R"
eval(parse(text = getURL(url_robust, ssl.verifypeer = FALSE)),
     envir=.GlobalEnv)
rm(url_robust)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Outliers and things...   ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

flintstones<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/flintstones.csv")

# No Barney OR Dino:
summary(lm(Y~X,data=subset(flintstones,name!="Dino" & name!="Barney")))

# No Barney (Dino included):
summary(lm(Y~X,data=subset(flintstones,name!="Barney")))

# Variance:

LittleDahl<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/LittleDahl.csv")

library(car)
with(LittleDahl, scatterplotMatrix(~age+tenure+unified+nulls))

Fit<-with(LittleDahl, lm(nulls~age+tenure+unified))
summary(Fit)

FitResid<-with(LittleDahl, (nulls - predict(Fit))) # residuals
FitStandard<-rstandard(Fit) # standardized residuals
FitStudent<-rstudent(Fit) # studentized residuals
FitCooksD<-cooks.distance(Fit) # Cook’s D
FitDFBeta<-dfbeta(Fit) # DFBeta
FitDFBetaS<-dfbetas(Fit) # DFBetaS
FitCOVRATIO<-covratio(Fit) # COVRATIOs

FitStudent[74]
LittleDahl$Congress74<-rep(0,length=104)
LittleDahl$Congress74[74]<-1
summary(with(LittleDahl, lm(nulls~age+tenure+unified+Congress74)))

influencePlot(Fit,id.n=4,labels=LittleDahl$Congress,id.cex=0.8,
              id.col="red",xlab="Leverage")

dfbetasPlots(Fit,id.n=5,id.col="red",main="",pch=19)

plot(FitCOVRATIO~names(FitCOVRATIO),pch=19,xlab="Congress",
     ylab="Value of COVRATIO")
abline(h=1,lty=2)

Outlier<-rep(0,104)
Outlier[74]<-1
Outlier[98]<-1
Outlier[104]<-1
DahlSmall<-LittleDahl[which(Outlier==0),]

summary(lm(nulls~age+tenure+unified,data=DahlSmall))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Africa<-read.csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/africa2001.csv")

summary(Africa)

Fit <- with(Africa, 
            lm(adrate~gdppppd+muslperc+subsaharan+healthexp+
                 literacy+internalwar))
summary(Fit)

# Table:

stargazer(Fit,
          type="latex",
          title="OLS Model of HIV Prevalence Rates, 2001",
          dep.var.caption="",
          dep.var.labels=c("OLS"),
          column.labels=c("Model I"),
          digits=2,
          no.space=TRUE,
          intercept.top=TRUE,intercept.bottom=FALSE,
          covariate.labels=c("(Constant)","GDP Per Capita",
                             "Muslim Percent","Subsaharan Africa",
                             "Health Expenditures","Literacy Rate",
                             "Civil War"))

# What not to do:

Nope <- gvlma(Fit)
display.gvlmatests(Nope)

# Better:

pdf("DefaultLMPlots.pdf",10,8)
par(mfrow=c(2,3))
plot(Fit,which=c(1:6))
dev.off()

# Unpacking that:
#
# #1: Residuals vs. fitted, same as:

pdf("ResidVsFitted.pdf",7,6)
plot(Fit$residuals~Fit$fitted.values,ylab="Residuals",
     xlab="Fitted Values",main="Residuals vs Fitted")
abline(h=0,lty=2)
lines(lowess(Fit$residuals~Fit$fitted.values),lwd=2,
      col="red")
dev.off()

# #2: QQ plot of residuals:

pdf("ResidQQPlot.pdf",6,5)
par(mar=c(4,4,2,2))
plot(Fit,which=2)
dev.off()

# #3: Scale-Location plot:

pdf("ScaleLocationPlot.pdf",6,5)
par(mar=c(4,4,2,2))
plot(Fit,which=3)
dev.off()

# #4: Cook's Distance (D):

pdf("CooksDPlot.pdf",6,5)
par(mar=c(4,4,2,2))
plot(Fit,which=4)
dev.off()


# #5: Residuals vs. Leverage:

pdf("ResidVsLeveragePlot.pdf",6,5)
par(mar=c(4,4,2,2))
plot(Fit,which=5)
dev.off()

# #6: Cook's D vs. Leverage:

pdf("CooksDVsLeveragePlot.pdf",6,5)
par(mar=c(4,4,2,2))
plot(Fit,which=6)
dev.off()

# Which are the outliers?

ASmall<-cbind(Africa[,3],Fit$model)

ASmall[c(36,37,39),]


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Another possibly useful plot (not shown):
# 
# ToPlot<-data.frame(Africa$adrate,Fit$fitted.values,
#                    Fit$residuals,Africa$gdppppd,Africa$muslperc,
#                    Africa$subsaharan,Africa$healthexp,Africa$literacy,
#                    Africa$internalwar)
# 
# pdf("UsefulPlot.pdf",8,7)
# scatterplotMatrix(ToPlot)
# dev.off()


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Variances, "robust" things, etc. ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ANES 2016 pilot study aggregation example...

ANES<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/ANES-pilot-2016.csv")

ANES$ftgay<-ifelse(ANES$ftgay==998,NA,ANES$ftgay)

# Average feeling thermometers about gays and lesbians:

summary(ANES$ftgay)

# States:

ANES$State<-recode(ANES$state,
"1='AL';2='AK';4='AZ';5='AR';6='CA';8='CO';9='CT';
10='DE';11='DC';12='FL';13='GA';15='HI';16='ID';17='IL';
18='IN';19='IA';20='KS';21='KY';22='LA';23='ME';24='MD';
25='MA';26='MI';27='MN';28='MS';29='MO';30='MT';31='NE';
32='NV';33='NH';34='NJ';35='NM';36='NY';37='NC';38='ND';
39='OH';40='OK';41='OR';42='PA';44='RI';45='SC';46='SD';
47='TN';48='TX';49='UT';50='VT';51='VA';53='WA';54='WV';
55='WI';56='WY'")

# Aggregate by state:

ANES$one<-1
StateFT<-ddply(ANES,.(State),summarise,
               Nresp=sum(one),
               meantherm=mean(ftgay,na.rm=TRUE))
summary(StateFT)

respfit<-with(StateFT, lm(meantherm~log(Nresp)))

pdf("StateThermPlot.pdf",6,5)
par(mar=c(4,4,2,2)) 
with(StateFT, plot(Nresp,meantherm,pch=".",col="white",log="x",
                   xlab="ln(N of Respondents)",xlim=c(0.5,200),
                   ylab="Statewide Mean Score"))
with(StateFT, text(Nresp,meantherm,log="x",labels=StateFT$State,
                   cex=0.3*log(StateFT$Nresp+1)))
abline(h=mean(ANES$ftgay,na.rm=TRUE),lwd=2)
abline(h=mean(StateFT$meantherm),lwd=2,lty=2,col="red")
abline(respfit,lwd=3,col="darkgreen")
legend("bottomright",bty="n",legend=c("Disaggregated Mean",
                               "Aggregated Mean",
                               "Aggregated Linear Fit"),
       lwd=c(2,2,3),lty=c(1,2,1),col=c("black","red","darkgreen"))
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# What do "robust" SEs do? A simulation:

set.seed(7222009)
X <- rnorm(10)
Y <- 1 + X + rnorm(10)
df10 <- data.frame(ID=seq(1:10),X=X,Y=Y)

fit10 <- lm(Y~X,data=df10)
summary(fit10)
rob10 <- vcovHC(fit10,type="HC1")
sqrt(diag(rob10))

# "Clone" each observation 100 times

df1K <- df10[rep(seq_len(nrow(df10)), each=100),]
df1K <- pdata.frame(df1K, index="ID")

fit1K <- lm(Y~X,data=df1K)
summary(fit1K)
summary(fit1K, cluster="ID")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Justices data:

Justices<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/Justices.csv")

summary(Justices)

OLSfit<-with(Justices, lm(civrts~score))
summary(OLSfit)

WLSfit<-with(Justices, lm(civrts~score,weights=lnNedit))
summary(WLSfit)

pdf("WLSBubblePlotR.pdf",6,6)
par(mar=c(4,4,2,2))
with(Justices, symbols(score, civrts,circles=Neditorials,
        ylab="Civil Rights Voting",xlab="Segal-Cover Score",
        ylim=c(0,100)))
abline(reg=OLSfit,lwd=2)
abline(reg=WLSfit,lwd=2,lty=2)
with(Justices, points(score,civrts,pch=20))
legend("topleft",bty="n",lty=c(1,2),lwd=2,
       legend=c("OLS","WLS"))
dev.off()

# "Robust"

hccm(OLSfit, type="hc1")

OLSfit2<-ols(civrts~score, x=TRUE, y=TRUE)
RobSEs<-robcov(OLSfit2)
RobSEs

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Binary response models ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# "Toy" example:

set.seed(7222009)
ystar<-rnorm(100)
y<-ifelse(ystar>0,1,0)
x<-ystar+(0.5*rnorm(100))
data<-data.frame(ystar,y,x)
head(data)

pdf("YstarYX-R.pdf",6,5)
par(mar=c(4,4,2,2))
plot(x,ystar,pch=19,ylab="Y* / Y",xlab="X")
points(x,y,pch=4,col="red")
abline(h=0)
legend("topleft",bty="n",pch=c(19,4),col=c("black","red"),
       legend=c("Y*","Y"))
dev.off()

# Probits and logits...

myprobit<-glm(y~x,family=binomial(link="probit"),
              data=data)
summary(myprobit)

mylogit<-glm(y~x,family=binomial(link="logit"),
             data=data)
summary(mylogit)

pdf("LogitProbitHats.pdf",6,5)
plot(mylogit$fitted.values,myprobit$fitted.values,
     pch=20,xlab="Logit Predictions",
     ylab="Probit Predictions")
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# NAFTA example...

NAFTA<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/NAFTA.csv")

summary(NAFTA)

# Probit:

NAFTA.GLM.probit<-glm(vote~democrat+pcthispc+cope93+DemXCOPE,
                   NAFTA,family=binomial(link="probit"))
summary(NAFTA.GLM.probit)

# Logit:

NAFTA.GLM.fit<-glm(vote~democrat+pcthispc+cope93+DemXCOPE,
                   NAFTA,family=binomial)
summary(NAFTA.GLM.fit)

# Add OLS, for fun:

NAFTA.OLS.fit<-lm(vote~democrat+pcthispc+cope93+DemXCOPE,
                  data=NAFTA)

# Table:

stargazer(NAFTA.GLM.probit,NAFTA.GLM.fit,NAFTA.OLS.fit,
          type="latex",
          title="Probit / Logit / OLS Models of the NAFTA Vote",
          dep.var.caption="",
          dep.var.labels=c("NAFTA Vote"),
          column.labels=c("Probit","Logit","OLS"),
          digits=2,
          no.space=TRUE,
          intercept.top=TRUE,intercept.bottom=FALSE,
          covariate.labels=c("(Constant)","Democratic Member",
                             "Hispanic Percent","COPE Score",
                             "Democratic Member x COPE Score"))

# Interactions...

NAFTA.GLM.fit$coeff[4]+NAFTA.GLM.fit$coeff[5]
(NAFTA.GLM.fit$coeff[4]+NAFTA.GLM.fit$coeff[5]) / 
  (sqrt(vcov(NAFTA.GLM.fit)[4,4] + 
          (1)^2*vcov(NAFTA.GLM.fit)[5,5] + 
          2*1*vcov(NAFTA.GLM.fit)[4,5]))

# Same thing, using -car-:

library(car)
linear.hypothesis(NAFTA.GLM.fit,"cope93+DemXCOPE=0")

# Predicted values:

preds<-NAFTA.GLM.fit$fitted.values
hats<-predict(NAFTA.GLM.fit,se.fit=TRUE)

# Plotting in-sample predictions:

XBUB<-hats$fit + (1.96*hats$se.fit) 
XBLB<-hats$fit - (1.96*hats$se.fit)
plotdata<-cbind(as.data.frame(hats),XBUB,XBLB)
plotdata<-data.frame(lapply(plotdata,binomial(link="logit")$linkinv))
par(mfrow=c(1,2))
library(plotrix)
plotCI(cope93[democrat==1],plotdata$fit[democrat==1],ui=plotdata$XBUB[democrat==1],
       li=plotdata$XBLB[democrat==1],pch=20,xlab="COPE Score",ylab="Predicted 
         Pr(Pro-NAFTA Vote)")
plotCI(cope93[democrat==0],plotdata$fit[democrat==0],ui=plotdata$XBUB[democrat==0],
       li=plotdata$XBLB[democrat==0],pch=20,xlab="COPE Score",ylab="Predicted 
         Pr(Pro-NAFTA Vote)")

# Plotting Out-of-sample Predictions:

sim.data<-data.frame(pcthispc=mean(nafta$pcthispc),democrat=rep(0:1,101),
                     cope93=seq(from=0,to=100,length.out=101))
sim.data$DemXCOPE<-sim.data$democrat*sim.data$cope93

OutHats<-predict(NAFTA.GLM.fit,se.fit=TRUE,newdata=sim.data)
OutHatsUB<-OutHats$fit+(1.96*OutHats$se.fit)
OutHatsLB<-OutHats$fit-(1.96*OutHats$se.fit)
OutHats<-cbind(as.data.frame(OutHats),OutHatsUB,OutHatsLB)
OutHats<-data.frame(lapply(OutHats,binomial(link="logit")$linkinv))

par(mfrow=c(1,2))
both<-cbind(sim.data,OutHats)
both<-both[order(both$cope93,both$democrat),]

plot(both$cope93[democrat==1],both$fit[democrat==1],t="l",lwd=2,ylim=c(0,1),
     xlab="COPE Score",ylab="Predicted Pr(Pro-NAFTA Vote)")
lines(both$cope93[democrat==1],both$OutHatsUB[democrat==1],lty=2)
lines(both$cope93[democrat==1],both$OutHatsLB[democrat==1],lty=2)
text(locator(1),label="Democrats")

plot(both$cope93[democrat==0],both$fit[democrat==0],t="l",lwd=2,ylim=c(0,1),
     xlab="COPE Score",ylab="Predicted Pr(Pro-NAFTA Vote)")
lines(both$cope93[democrat==0],both$OutHatsUB[democrat==0],lty=2)
lines(both$cope93[democrat==0],both$OutHatsLB[democrat==0],lty=2)
text(locator(1),label="Republicans")

# Simpler version (for when there aren't interactions,
# and you re interested in a single continuous
# predictor), using -margins-:

pdf("PropHispPredPlot.pdf",7,6)
par(mar=c(4,4,2,2))
cplot(NAFTA.GLM.fit,"pcthispc",xlab="Proportion Hispanic")
dev.off()


# Odds Ratios:

lreg.or <- function(model)
{
  coeffs <- coef(summary(NAFTA.GLM.fit))
  lci <- exp(coeffs[ ,1] - 1.96 * coeffs[ ,2])
  or <- exp(coeffs[ ,1])
  uci <- exp(coeffs[ ,1] + 1.96 * coeffs[ ,2])
  lreg.or <- cbind(lci, or, uci)        
  lreg.or
}

lreg.or(NAFTA.GLM.fit)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Goodness of fit:

table(NAFTA$vote)
table(NAFTA.GLM.fit$fitted.values>0.5,nafta$vote==1)
chisq.test(NAFTA.GLM.fit$fitted.values>0.5,nafta$vote==1)

# PREs, with different cutoffs / taus...
#
# Tau = 0.2:

Hats02<-ifelse(NAFTA.GLM.fit$fitted.values>0.2,1,0)
CrossTable(NAFTA$vote,Hats02,prop.r=FALSE,prop.c=FALSE,
           prop.t=FALSE,prop.chisq=FALSE)

# Tau = 0.8:

Hats08<-ifelse(NAFTA.GLM.fit$fitted.values>0.8,1,0)
CrossTable(NAFTA$vote,Hats08,prop.r=FALSE,prop.c=FALSE,
           prop.t=FALSE,prop.chisq=FALSE)


# ROC curves, plots, etc. (using -ROCR-):

NAFTA.hats<-predict(NAFTA.GLM.fit,type="response")
preds<-ROCR::prediction(NAFTA.hats,NAFTA$vote)

# Plot:

pdf("NAFTA-ROC-plot.pdf",6,5)
plot(performance(preds,"tpr","fpr"),lwd=2,lty=2,
     col="red",xlab="1 - Specificity",ylab="Sensitivity")
abline(a=0,b=1,lwd=3)
dev.off()

# "Bad" model w/only PCTHISPC in X:

NAFTA.bad<-with(NAFTA,
                glm(vote~pcthispc,family=binomial(link="logit")))
NAFTA.bad.hats<-predict(NAFTA.bad,type="response")
bad.preds<-ROCR::prediction(NAFTA.bad.hats,NAFTA$vote)

pdf("BadROCCurve.pdf",6,5)
plot(performance(bad.preds,"tpr","fpr"),lwd=2,lty=2,
     col="red",xlab="1 - Specificity",ylab="Sensitivity")
abline(a=0,b=1,lwd=3)
dev.off()

# Comparing ROCs:

GoodROC<-roc(NAFTA$vote,NAFTA.hats,ci=TRUE)
GoodAUC<-auc(GoodROC)
BadROC<-roc(NAFTA$vote,NAFTA.bad.hats)
BadAUC<-auc(BadROC)

GoodAUC

BadAUC

# Comparison plot:

pdf("TwoROCs.pdf",5,5)
par(mar=c(4,4,2,2))
plot(GoodROC)
lines(BadROC,col="red",lty=2)
dev.off()


# /fin