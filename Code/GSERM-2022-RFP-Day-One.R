# Introduction / Setup    ########################
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Code for GSERM's "Regression for Publishing"
# (June 2021)
#
# Day One: Regression Overview & Review
#
# File created 6/21/2020
#
# File last modified 6/10/2022
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Set working directory first:
#
# setwd("~/Dropbox (Personal)/GSERM/RFP-2022/Slides") # or whatever...
#
# Load relevant R packages:

# install.packages("readr") <- Install readr package, if necessary
library(readr)
# install.packages("arm") <- Install arm package, if necessary
library(arm)
# install.packages("car") <- Install car package, if necessary
library(car)
# install.packages("psych") <- Install psych package, if necessary
library(psych)
# install.packages("lmtest") <- Install lmtest package, if necessary
library(lmtest)
# install.packages("plotrix") <- Install plotrix package, if necessary
library(plotrix)
# install.packages("stargazer") <- Install plotrix package, if necessary
library(stargazer)
# install.packages("dotwhisker") <- Install dotwhisker package
library(dotwhisker)

options(scipen = 6) # bias against scientific notation
options(digits = 4) # show fewer decimal places
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Regression, visually ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# Read data from the github repo:

IR2000<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/CountryData2000.csv")

# I offer no warranty that any of this is the most 
# efficient way to do what I'm doing... all graphics
# are using base R graphics functionality, with
# the exception of the scatterplot matrices, which use
# the -car- package.

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Figure 1

pdf("IMbyLE.pdf",7,7) # make a PDF
par(cex=1.2) # make the things bigger
par(mar=c(4,4,2,2)) # plot margins
with(IR2000, plot(lifeexpectancy,infantmortalityperK, # draw the plot
                  pch="",xlab="Life Expectancy at Birth",
                  ylab="Infant Mortality Per 1000 Live Births"))
par(cex=0.8) # make the things smaller
# Add names:
with(IR2000, text(lifeexpectancy,infantmortalityperK,label=WBcode))
# add a linear-fit line
with(IR2000, abline(lm(infantmortalityperK~lifeexpectancy),lwd=3))
dev.off() # turn off the PDF-maker

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Figure 2 ("residual" plot)

small<-data.frame(infantmortalityperK,lifeexpectancy,WBcode)
small<-small[is.na(small$lifeexpectancy)==FALSE,]
small<-small[is.na(small$infantmortalityperK)==FALSE,]

IMLE.fit<-lm(small$infantmortalityperK~small$lifeexpectancy)
IMLE.res<-resid(IMLE.fit)

pdf("IMbyLEresids.pdf",7,7) # make a PDF
par(cex=1.2)
par(mar=c(4,4,2,2))
with(IR2000, plot(small$lifeexpectancy,IMLE.res,pch="",
                  xlab="Life Expectancy at Birth",
                  ylab="Residuals (Observed minus Expected)"))
par(cex=0.8)
with(IR2000, text(small$lifeexpectancy,IMLE.res,
                  label=small$WBcode))
abline(h=0,lwd=3)
dev.off() # turn off the PDF-maker

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Figure 3

pdf("IMbyFertility.pdf",7,7) # make a PDF
par(mar=c(4,4,2,2)) # margins
par(cex=1.2) #symbol size
with(IR2000, plot(fertility,infantmortalityperK,pch="",xlab="Fertility Rate (Births Per Woman)",
                  ylab="Infant Mortality Per 1000 Live Births"))
par(cex=0.8)
with(IR2000, text(fertility,infantmortalityperK,label=WBcode))
with(IR2000, abline(lm(infantmortalityperK~fertility),lwd=3))
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Figure 4

pdf("IMbyGDP.pdf",7,7) # make a PDF
par(mar=c(4,4,2,2)) # margins
par(cex=1.2)
with(IR2000, plot(NEWrgdpch,infantmortalityperK,pch="",xlab="Real GDP per capita",
                  ylab="Infant Mortality Per 1000 Live Births"))
par(cex=0.8)
with(IR2000, text(NEWrgdpch,infantmortalityperK,label=WBcode))
with(IR2000, abline(lm(infantmortalityperK~NEWrgdpch),lwd=3))
with(IR2000, lines(lowess(infantmortalityperK~NEWrgdpch,span=0.1,iter=1),
                   lwd=3,col="red",lty=2))
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Figure 5

pdf("lnIMbylnGDP.pdf",7,7) # make a PDF
par(mar=c(4,4,2,2)) # margins
par(cex=1.2)
with(IR2000, plot(log(NEWrgdpch),log(infantmortalityperK),pch="",xlab="Logged Real GDP per 
       Capita",ylab="Logged Infant Mortality Per 1000 Live Births"))
par(cex=0.8)
with(IR2000, text(log(NEWrgdpch),log(infantmortalityperK),label=WBcode))
with(IR2000, abline(lm(log(infantmortalityperK)~log(NEWrgdpch)),lwd=3))
with(IR2000, lines(lowess(log(infantmortalityperK)~log(NEWrgdpch)),
                   lwd=3,lty=2,col="red"))
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Figure 6

pdf("IMbyPOLITY.pdf",7,7) # make a PDF
par(mar=c(4,4,2,2)) # margins
par(cex=1.2)
with(IR2000, plot(polity2,infantmortalityperK,pch="",xlab="POLITY IV score",ylab="Infant 
       Mortality Per 1000 Live Births"))
par(cex=0.8)
with(IR2000, text(polity2,infantmortalityperK,label=WBcode))
with(IR2000, abline(lm(infantmortalityperK~polity2),lwd=3))
with(IR2000, lines(lowess(infantmortalityperK~polity2),lwd=3,
                   lty=2,col="red"))
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Figure 7

IR2000$Rich<-as.factor(IR2000$NEWrgdpch>median(IR2000$NEWrgdpch,na.rm=TRUE))

pdf("IMbyPOLITYandGDP.pdf",7,7) # make a PDF
par(mar=c(4,4,2,2)) # margins
par(cex=1.2)
scatterplot(infantmortalityperK~polity2|Rich,IR2000,cex=1.5,
            cex.axis=1.5,cex.lab=1.3,boxplots="",xlab="POLITY IV", 
            ylab="Infant Mortality per 1000 Live Births",
            pch=c(16,17),legend.coords="topleft",grid=FALSE,
            reg.line=FALSE,lwd=c(3,3),col=c(1,2),lty=c(3,4))
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# IM / DPT Regression Example ####
#
# Summary statistics

with(IR2000, describe(infantmortalityperK))
with(IR2000, describe(DPTpct))

# Regression:

IMDPT<-lm(infantmortalityperK~DPTpct,data=IR2000,na.action=na.exclude)
summary.lm(IMDPT)

# ANOVA

anova(IMDPT)

# Scatterplot:

pdf("IMDPT.pdf",6,6) # <- create PDF
with(IR2000, plot(DPTpct,infantmortalityperK,pch=".",
                xlab="DPT Immunization Percentage",
                ylab="Infant Mortality (Deaths per 1000 Births)"))
with(IR2000, text(DPTpct,infantmortalityperK,labels=WBcode))
with(IR2000, abline(v=mean(DPTpct,na.rm=TRUE),lty=2))
with(IR2000, abline(h=mean(infantmortalityperK,na.rm=TRUE),lty=2))
abline(IMDPT,lwd=3)
dev.off()

# Residuals (u):

IR2000$IMDPTres <- with(IR2000, residuals(IMDPT))
describe(IR2000$IMDPTres)

# Residual density plot:

pdf("IMDPTResidualsDensity.pdf",6,6)
with(IR2000, plot(density(IMDPTres,na.rm=TRUE),
                main="Density Plot: Regression Residuals",
                xlab="Residual Value"))
abline(v=0,lty=2,lwd=2)
dev.off()

# Fitted Values:

IR2000$IMDPThat<-fitted.values(IMDPT)
describe(IR2000$IMDPThat)

# Densities plot:

pdf("IMDPTFittedDensity.pdf",6,6)
with(IR2000, plot(density(IMDPThat,na.rm=TRUE),
                main="Density Plot: Actual and Fitted Values",
                xlab="Values of Y"))
with(IR2000, lines(density(infantmortalityperK,na.rm=TRUE),
                 lty=2,col="red"))
with(IR2000, abline(v=mean(infantmortalityperK,na.rm=TRUE),
                  lty=2,lwd=2))
dev.off()

# Correlations:

with(IR2000, cor(infantmortalityperK,DPTpct,use="complete.obs"))
with(IR2000, cor(IMDPTres,infantmortalityperK,use="complete.obs"))
with(IR2000, cor(IMDPTres,DPTpct,use="complete.obs"))
with(IR2000, cor(IMDPThat,infantmortalityperK,use="complete.obs"))
with(IR2000, cor(IMDPThat,DPTpct,use="complete.obs"))
with(IR2000, cor(IMDPTres,IMDPThat,use="complete.obs"))

# Plotting residuals vs. X:

pdf("IMDPTResiduals.pdf",6,6)
with(IR2000, plot(DPTpct,IMDPTres,pch=".",
                xlab="DPT Immunization Percentage",
                ylab="Residuals (Y - Y-hat)"))
with(IR2000, text(DPTpct,IMDPTres,labels=WBcode))
abline(h=0,lwd=2)
dev.off()

# Squared residuals vs. X:

pdf("IMDPTSquaredResiduals.pdf",6,6)
with(IR2000, plot(DPTpct,IMDPTres^2,pch=".",
                xlab="DPT Immunization Percentage",
                ylab="Residuals (Y - Y-hat)"))
with(IR2000, text(DPTpct,IMDPTres^2,labels=WBcode))
dev.off() 


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Multivariate Regression ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# Data on Africa (year = 2001):

Data<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/africa2001.csv")

Data<-with(Data, data.frame(adrate,polity,
                            subsaharan=as.numeric(as.factor(subsaharan))-1,
                            muslperc,literacy))
describe(Data)
cor(Data)

# Scatterplot matrix:

scatterplotMatrix(Data)

# Multivariate regression:

model<-lm(adrate~polity+subsaharan+muslperc+literacy,data=Data)
summary(model)

# Variance-Covariance matrix of betas:

options(digits=4)
vcov(model)

# F-tests...
#
# POLITY = subsaharan = 0:

modelsmall<-lm(adrate~muslperc+literacy,data=Data)
waldtest(model,modelsmall)  # from -lmtest- package

# muslperc = 0.1:

linearHypothesis(model,"muslperc=0.1") # from -car-

# literacy = muslperc:

linearHypothesis(model,"literacy=muslperc")

# Reporting: the output:

summary(model)

# Reporting: The table:

stargazer(model,
          type="latex",
          title="OLS Regression Model of HIV/AIDS Rates in Africa, 2001",
          dep.var.caption="",
          dep.var.labels="Model I",
          digits=2,
          no.space=TRUE,
          intercept.top=TRUE,intercept.bottom=FALSE,
          covariate.labels=c("(Constant)",
                             "POLITY Score",
                             "Subsaharan Africa",
                             "Muslim Percentage of the Population",
                             "Literacy Rate"))

# Multiple models:

model2<-lm(adrate~polity+subsaharan+muslperc,data=Data)

stargazer(model,model2,
          type="latex",
          title="OLS Regression Models of HIV/AIDS Rates in Africa, 2001",
          dep.var.caption="",
          column.labels=c("w/Literacy","w/o Literacy"),
          digits=2,
          no.space=TRUE,
          intercept.top=TRUE,intercept.bottom=FALSE,
          covariate.labels=c("(Constant)",
                             "POLITY Score",
                             "Subsaharan Africa",
                             "Muslim Percentage of the Population",
                             "Literacy Rate"))

# Replacing tables with figures...
#
# Basic model:

pdf("OLS-Ladder-I.pdf",6,4)
dwplot(model,show_intercept=TRUE) + theme_bw()
dev.off()

# Ouch; rescale the predictors...

modelS<-standardize(model)
summary(modelS)

# Redraw the plot:

pdf("OLS-Ladder-II.pdf",6,4)
dwplot(modelS,show_intercept=TRUE) + theme_bw() + 
        geom_vline(xintercept = 0,colour = "black",
                   linetype = 2)
dev.off()

# Add the second model + make pretty:

model2S<-standardize(model2)

pdf("OLS-Ladder-III.pdf",8,4)
m2<-dwplot(list(modelS,model2S),
       show_intercept=TRUE) + 
        geom_vline(xintercept = 0,colour = "black",
                   linetype = 2) 
relabel_predictors(m2,c("(Intercept)"="(Intercept)",
                   z.polity="Standardized POLITY",
                   c.subsaharan="Subsaharan Africa",
                   z.muslperc="Standardized Muslim Percentage",
                   z.literacy="Standardized Literacy Rate")) +
        theme_bw() +
        geom_vline(xintercept = 0,colour = "black",
                   linetype = 2)
dev.off()


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# "Supplemental" things ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# Simple regression (Africa data):

fit<-with(africa, lm(adrate~muslperc))
summary(fit)

# First figure

SEs<-predict(fit,interval="confidence")
Sort<-order(africa$muslperc)

pdf("SRTFig1.pdf",6,5)
plot(africa$muslperc, africa$adrate, 
     xlab="Muslim Percentage of the Population", 
     ylab="HIV Prevalence Rate",pch=16) 
abline(fit,lwd=3)
lines(sort(africa$muslperc),SEs[Sort,2],col="red",lwd=2,lty=2)
lines(sort(africa$muslperc),SEs[Sort,3],col="red",lwd=2,lty=2)
dev.off()

# Add 10:

africa$muslplusten<-africa$muslperc+10
fit2<-with(africa, lm(adrate~muslplusten,data=africa))
summary(fit2)

SEs<-predict(fit2,interval="confidence")

pdf("SRTFig2.pdf",6,5)
plot(africa$muslplusten, africa$adrate, 
     xlab="Muslim Percentage of the Population + 10", 
     ylab="HIV Prevalence Rate",pch=16)
abline(fit2,lwd=3)
lines(sort(africa$muslplusten),SEs[Sort,2],col="red",lwd=2,lty=2)
lines(sort(africa$muslplusten),SEs[Sort,3],col="red",lwd=2,lty=2)
dev.off()

# Multiply Y times -314:

africa$screwyrate<-africa$adrate*(-314)
fit3<-with(africa, lm(screwyrate~muslperc))
summary(fit3)

SEs<-predict(fit3,interval="confidence")

pdf("SRTFig3.pdf",6,5)
plot(africa$muslperc, africa$screwyrate, 
     xlab="Muslim Percentage of the Population", 
     ylab="HIV Prevalence Rate times -314",pch=16,
     ylim=c(-13000,1000))
abline(fit3,lwd=3)
lines(sort(africa$muslperc),SEs[Sort,2],col="red",lwd=2,lty=2)
lines(sort(africa$muslperc),SEs[Sort,3],col="red",lwd=2,lty=2)
dev.off()

# Reversing the scales:

africa$nonmuslimpct <- 100 - africa$muslperc
africa$noninfected <- 100 - africa$adrate
fit4<-lm(noninfected~nonmuslimpct,data=africa)
summary(fit4)

SEs<-predict(fit4,interval="confidence")
Sort2 <- order(africa$nonmuslimpct)

pdf("SRTFig4.pdf",6,5)
plot(africa$nonmuslimpct, africa$noninfected, 
     xlab="Non-Muslim Percentage of the Population", 
     ylab="HIV Non-Prevalence Rate",pch=16,
     ylim=c(60,111))
abline(fit4,lwd=3)
lines(sort(africa$nonmuslimpct),SEs[Sort2,2],col="red",lwd=2,lty=2)
lines(sort(africa$nonmuslimpct),SEs[Sort2,3],col="red",lwd=2,lty=2)
dev.off()

# Centering X:

africa$muslcenter<-africa$muslperc - mean(africa$muslperc, na.rm=TRUE)
fit5<-lm(adrate~muslcenter,data=africa)
summary(fit5)

SEs<-predict(fit5,interval="confidence")

pdf("SRTFig5.pdf",6,5)
plot(africa$muslcenter, africa$adrate, 
     xlab="Centered Muslim Percentage of the Population", 
     ylab="HIV Prevalence Rate",pch=16,
     ylim=c(-10,40))
abline(fit5,lwd=3)
lines(sort(africa$muslcenter),SEs[Sort,2],col="red",lwd=2,lty=2)
lines(sort(africa$muslcenter),SEs[Sort,3],col="red",lwd=2,lty=2)
dev.off()

# Rescaling X for interpretability:

fit6<-lm(adrate~population,data=africa)
summary(fit6)

africa$popmil<-africa$population / 1000000
fit7<-lm(adrate~popmil,data=africa)
summary(fit7)

# Dichotomous X and t-tests:

fit8<-lm(adrate~subsaharan,data=africa)
summary(fit8)

with(africa,
     t.test(adrate~subsaharan, var.equal=TRUE))

################
# Reporting:

fit<-lm(adrate~muslperc, data=africa)
summary.lm(fit)

# Easy LaTeX table using *stargazer*:

stargazer(fit,
          type="latex",
          title="OLS Regression Model of HIV/AIDS Rates in Africa, 2001",
          dep.var.caption="",
          dep.var.labels="Model I",
          digits=2,
          no.space=TRUE,
          intercept.top=TRUE,intercept.bottom=FALSE,
          covariate.labels=c("(Constant)","
                             Muslim Percentage of the Population"))
