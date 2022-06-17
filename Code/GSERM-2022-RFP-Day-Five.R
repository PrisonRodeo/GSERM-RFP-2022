#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Introduction / packages / etc. ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Code GSERM "Regression for Publishing"
# (June 2021)
#
# Day Five: Causal Inference, and
# Panel / TSCS models...
#
# File created 6/16/2022
#
# File last modified 6/16/2022
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Packages:

library(readr)    # <-- install all these as necessary
library(stargazer)# <--
library(plyr)     # <-- 
library(dplyr)    # <-- 
library(car)      # <-- 
library(rms)      # <-- 
library(lmtest)   # <-- 
library(lme4)
library(plm)
library(gtools)
library(texreg)
library(statmod)
library(stargazer)
library(corrplot)
library(rgenoud)
library(MatchIt)
library(Matching)
library(pscl)
library(stargazer)
library(prais)
library(nlme)
library(tseries)
library(pcse)
library(panelView)
library(OrthoPanels)
library(glmmML)
library(bife)
library(censReg)
library(geepack)
library(fixest)
library(pcse)

# Options:
options(scipen=12)
options(digits=4)

# setwd()...

setwd("~/Dropbox (Personal)/GSERM/RFP-2022")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Potential outcomes... ####

Cov0<-data.frame(i=seq(1:6),
                 W=c(rep(0,3),rep(1,3)),
                 Y0=c(8,10,12,8,10,12),
                 Y1=c(10,12,14,10,12,14),
                 Ydiff=rep(2,6),
                 Y=c(8,10,12,10,12,14))

with(Cov0, t.test(Y~W))

PosCov<-data.frame(i=seq(1:6),
                   W=c(rep(0,3),rep(1,3)),
                   Y0=c(8,8,10,10,12,12),
                   Y1=c(10,10,12,12,14,14),
                   Ydiff=rep(2,6),
                   Y=c(8,8,10,12,14,14))

with(PosCov, t.test(Y~W))

NegCov<-data.frame(i=seq(1:6),
                   W=c(rep(0,3),rep(1,3)),
                   Y0=c(12,12,10,10,8,8),
                   Y1=c(14,14,12,12,10,10),
                   Ydiff=rep(2,6),
                   Y=c(12,12,10,12,10,10))

with(NegCov, t.test(Y~W))

# Not run:

Varying<-data.frame(i=seq(1:6),
                    W=c(rep(0,3),rep(1,3)),
                    Y0=c(8,10,12,8,10,12),
                    Y1=c(9,11,14,9,11,14),
                    Ydiff=c(1,2,3,1,2,3),
                    Y=c(8,10,12,9,11,14))

with(Varying, t.test(Y~W))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Bigger simulation (not in slides):

reps <- 1000 # number of sims
N<-1000      # sample size

ATE.U <- numeric(reps) # ATEs w/no confounding
ATE.C <- numeric(reps) # ATEs w/confounding by X
ATE.A <- numeric(reps) # Confounded + adjusted ATEs 

set.seed(7222009)

for(i in 1:reps) {
  W0 <- rep(0,N)
  W1 <- rep(1,N)
  Y0 <- 0 + 1*W0 + rnorm(N) # counterfactual Y for W=0; ATE=1
  Y1 <- 0 + 1*W1 + rnorm(N) # counterfactual Y for W=1; ATE=1
  X <- rnorm(N) # independent thing, for now
  W  <- numeric(N)
  Y  <- numeric(N)
  for(j in 1:N){
    Y[j] <- ifelse(X[j]<0,Y0[j],Y1[j])
    W[j] <- ifelse(X[j]<0,W0[j],W1[j])
  } # == effectively random selection of W, since Cov(X,Y)=0
  
  # boxplot(Y~W,ylim=c(-6,6))
  # t.test(Y~W) # works
  foo<-lm(Y~W) # no problem
  ATE.U[i] <- foo$coefficients[2] # stash the estimate
  # summary(lm(Y~W+X)) # also works; X is unnecessary
  #
  # Now make Y correlated with W **and** X:
  #
  W.cor  <- numeric(N)
  Y.cor  <- numeric(N)
  for(j in 1:N){
    Y.cor[j] <- ifelse(X[j]<0,Y0[j]+X[j],Y1[j]+X[j])
    W.cor[j] <- ifelse(Y.cor[j]==Y0[j]+X[j],W0[j],W1[j])
  } # Selection of W with Cov(X,Y)>0
  #
  # cor(X,W.cor)
  # cor(X,Y.cor)
  # boxplot(Y.cor~W.cor,ylim=c(-6,6))
  # t.test(Y.cor~W.cor) # wrong
  bar<-lm(Y.cor~W.cor) # bad upward bias; should be B~=1
  ATE.C[i] <- bar$coefficients[2] # stash the estimate
  baz<-lm(Y.cor~W.cor+X) # controlling for X fixes things
  ATE.A[i] <- baz$coefficients[2] # stash the estimate
}

# Plot the Betas!

pdf("Notes/ConfoundingSims.pdf",7,6)
par(mar=c(4,4,2,2))
plot(density(ATE.C),t="l",xlim=c(0.5,3),lwd=2,ylim=c(0,7),
     lty=2,col="red",main="",xlab="Estimated ATE")
lines(density(ATE.U),col="black",lwd=2,lty=1)
lines(density(ATE.A),col="darkgreen",lwd=3,lty=5)
abline(v=1,lty=3)
legend("topright",bty="n",lty=c(1,2,5),lwd=c(2,2,3),cex=0.8,
       col=c("black","red","darkgreen"),legend=c("No Confounding",
                                                 "Confounding w/o Adjustment",
                                                 "Confounding with Adjustment"))
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Causal inference... ####
#
# Simulations for balance:

set.seed(7222009)
N <- 10000
W <- c(rep(0,N/2),rep(1,N/2))
Xbal <- rnorm(N)
Xunbal <- rnorm(N,mean=(-1+(2*W)))
dd<-data.frame(W=W,Xbal=Xbal,Xunbal=Xunbal)

pdf("Notes/Balance.pdf",7,5)
par(mar=c(4,4,4,2))
par(mfrow=c(1,2))
# Balanced X:
plot(density(dd[dd$W==0,]$Xbal),lwd=2,main="Balanced X",
     xlab="X")
lines(density(dd[dd$W==1,]$Xbal),lwd=2,lty=2,col="red")
abline(v=mean(dd[dd$W==0,]$Xbal),lty=3)
abline(v=mean(dd[dd$W==1,]$Xbal),lty=3,col="red")
legend("topleft",bty="n",lwd=2,lty=c(1,2),col=c("black","red"),
       legend=c("W=0","W=1"),cex=0.9)
# Unbalanced X:
plot(density(dd[dd$W==0,]$Xunbal),lwd=2,main="Unbalanced X",
     xlab="X",xlim=c(-4,4))
lines(density(dd[dd$W==1,]$Xunbal),lwd=2,lty=2,col="red")
abline(v=mean(dd[dd$W==0,]$Xunbal),lty=3)
abline(v=mean(dd[dd$W==1,]$Xunbal),lty=3,col="red")
dev.off()


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Simulations for overlap:

set.seed(7222009)
N <- 100
X <- rnorm(N,0,1) # Pre-treatment confounder
WO <- rbinom(N,1,pnorm(0.1*X)) # overlapping treatment
WP <- rbinom(N,1,pnorm(1*X)) # partial overlapping treatment
WN <- rbinom(N,1,pnorm(8*X)) # non-overlapping treatment
YO <- 3 + 2*WO + 1*X + rnorm(N)
YP <- 3 + 2*WP + 1*X + rnorm(N)
YN <- 3 + 2*WN + 1*X + rnorm(N)
df<-data.frame(X=X,YO=YO,YP=YP,YN=YN,WO=WO,WP=WP,WN=WN)

pdf("Notes/Overlap.pdf",8,5)
par(mfrow=c(1,3))
par(mar=c(4,4,4,2))
with(df, plot(X,YO,pch=16+WO,col=WO+1,main="Complete Overlap",
              xlab="X",ylab="Y"))
legend("topleft",bty="n",pch=c(16,17),col=c("black","red"),
       legend=c("W=0","W=1"),cex=1.2)
with(df, plot(X,YP,pch=16+WP,col=WP+1,main="Moderate Overlap",
              xlab="X",ylab="Y"))
with(df, plot(X,YN,pch=16+WN,col=WN+1,main="No Overlap",
              xlab="X",ylab="Y"))
dev.off()


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# RDD sim for figure:

N <- 100
set.seed(7222009)
A <- runif(N,0,100)
T <- ifelse(A>50,1,0)
Y <- 8 + 0.1*A + 10*T + rnorm(N,0,4)
drdd<-data.frame(Y=Y,T=T,A=A)
fit0<-lm(Y~A,data=drdd[drdd$T==0,])
fit1<-lm(Y~A,data=drdd[drdd$T==1,])

pdf("Notes/RDD.pdf",5,6)
plot(A,Y,pch=16+T,col=T+1,ylim=c(0,35),
     xlim=c(0,100))
abline(v=50,lty=3)
clip(0,50,0,35)
abline(fit0,lwd=2)
clip(50,100,0,35)
abline(fit1,lwd=2,col="red")
clip(0,100,0,35)
text(20,30,labels=c("T=0"))
text(80,5,labels=c("T=1"),col="red")
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SPORTS data

sports <- read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/HS-Sports.csv")

summary(sports)

# Correlation plot:

pdf("Slides/Sports-CorrPlot.pdf",6,5)
corrplot(cor(sports), method="number",
         number.cex=0.8,tl.col="black")
dev.off()

# Simple OLS-ish models:

Xs <- c("sports","fincome","ses","workage","female",
        "academic","remedial","advanced")
Model <- paste("grades", paste(Xs,collapse=" + "),sep="~")

with(sports, t.test(grades~sports))
summary(lm(Model,data=sports))

# Checking for covariate balance:

set.seed(7222009)
Sports.bal.pre<-MatchBalance(sports~fincome+ses+workage+female+
                               academic+remedial+advanced,
                             data=sports,nboots=1000,digits=3)

# Gather P-values for balance statistics:

PreMatch.Ps<-numeric(7)
for(i in 1:7){
  PreMatch.Ps[i]<-print(Sports.bal.pre$BeforeMatching[[i]]$p.value)
}

# Show them:

pdf("Slides/PreMatch-Balance.pdf",8,5)
par(mar=c(4,7,2,2))
dotchart(PreMatch.Ps,labels=Xs[2:8],xlim=c(0,1),pch=19,
         xlab="Balance: P-Value")
abline(v=0.05,lty=3)
dev.off()

# Now, matching...

# Exact matching:

M.exact <- matchit(sports~fincome+ses+workage+female+academic+
                     remedial+advanced,data=sports,
                   method="exact")
M.exact

# Output matched data:

sports.exact <- match.data(M.exact,group="all")
dim(sports.exact)

# Balance?

set.seed(7222009)
Exact.bal<-MatchBalance(sports~fincome+ses+workage+female+
                          academic+remedial+advanced,
                        data=sports.exact,nboots=1000,digits=3)

# Gather P-values for balance statistics:

Exact.Ps<-numeric(7)
for(i in 1:7){
  Exact.Ps[i]<-print(Exact.bal$BeforeMatching[[i]]$p.value)
}

# Show them:

pdf("Slides/Exact-Balance.pdf",8,5)
par(mar=c(4,7,2,2))
dotchart(PreMatch.Ps,labels=Xs[2:8],xlim=c(0,1),pch=19,
         xlab="Balance: P-Value")
points(Exact.Ps,seq(1:7),pch=17,col="darkgreen")
abline(v=0.05,lty=3)
dev.off()

# Propensity Scores:

PSfit <- glm(sports~fincome+ses+workage+female+academic+remedial+
               advanced,data=sports,family=binomial(link="logit"))
# summary(PSfit)

# Generate scores & check common support:

PS.df <- data.frame(PS = predict(PSfit,type="response"),
                    sports = PSfit$model$sports)

pdf("Slides/PScore-Support.pdf",8,5)
par(mfrow=c(1,2))
with(PS.df[PS.df$sports==0,],
     plot(density(PS),main="Non-Athletes",lwd=2,
          xlab="Propensity Score",xlim=c(0,1)))
with(PS.df[PS.df$sports==1,],
     plot(density(PS),main="Athletes",lwd=2,col="red",
          xlab="Propensity Score",xlim=c(0,1)))
dev.off()

# Matching with propensity scores:

M.prop <- matchit(sports~fincome+ses+workage+female+academic+
                    remedial+advanced,data=sports,
                  method="nearest")
summary(M.prop)

# Get the data:

sports.prop <- match.data(M.prop,group="all")

# Balance:

set.seed(7222009)
Prop.bal<-MatchBalance(sports~fincome+ses+workage+female+
                         academic+remedial+advanced,
                       data=sports.prop,nboots=1000,digits=3)

# Gather P-values for balance statistics:

Prop.Ps<-numeric(7)
for(i in 1:7){
  Prop.Ps[i]<-print(Prop.bal$BeforeMatching[[i]]$p.value)
}

# Show them:

pdf("Notes/Exact-Prop-Balance.pdf",8,5)
par(mar=c(4,7,2,2))
dotchart(PreMatch.Ps,labels=Xs[2:8],xlim=c(0,1),pch=19,
         xlab="Balance: P-Value")
points(Exact.Ps,seq(1:7),pch=17,col="darkgreen")
points(Prop.Ps,seq(1:7),pch=15,col="red")
abline(v=0.05,lty=3)
dev.off()


# Genetic matching (not shown...)

M.genetic <- matchit(sports~fincome+ses+workage+female+academic+
                       remedial+advanced,data=sports,
                     method="genetic")
sports.genetic <- match.data(M.genetic,group="all")

# Re-fit the t-tests:

with(sports, t.test(grades~sports))$statistic # No matching
with(sports.exact, t.test(grades~sports))$statistic # Exact
with(sports.prop, t.test(grades~sports,paired=TRUE))$statistic # PS
with(sports.genetic, t.test(grades~sports))$statistic # Genetic

fit.lm <- lm(Model,data=sports)
fit.exact <- lm(Model,data=sports.exact)
fit.prop <- lm(Model,data=sports.prop)
fit.genetic <- lm(Model,data=sports.genetic)

regs <- texreg(l=list(fit.lm,fit.exact,fit.prop,fit.genetic),
               stars=c(0.05),caption=" ",
               custom.model.names=c("No Matching","Exact",
                                    "Propensity Score","Genetic"),
               custom.coef.names = c("(Intercept)","Plays Sports",
                                     "Family Income",
                                     "Socioeconomic Status",
                                     "Working Age",
                                     "Female","Academic Track",
                                     "Remedial Track",
                                     "Advanced Course"))

regs

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Panel / TSCS Models ####
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# Unit effects, etc.
#
# Cross-country example data, 1945-2014:

Demos<-read_csv("https://github.com/PrisonRodeo/GSERM-RFP-2022/raw/main/Data/DemonstrationsTSCS.csv")

Demos$X<-NULL
Demos$ColdWar <- with(Demos, 
                      ifelse(Year<1990,1,0))

# Create variables; kill duplicates...

Demos$lnGDP <- log(Demos$GDP)
Demos <- Demos[!duplicated(Demos[c("Year", "ccode")]),] # kill duplicates
PDF <- pdata.frame(Demos,index=c("ccode","Year")) # panel data frame

# summary(PDF)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Visualizing panel data...

pdf("PanelDemosViz.pdf",7,5)
panelview(lnDemons~POLITY+Monarch,data=Demos,theme.bw=TRUE,
          outcome.type="continuous",type="outcome",
          by.timing=TRUE,index=c("ccode","Year"),
          main=" ",ylab="log(Demonstrations)",
          legendOff=TRUE)
dev.off()

pdf("PanelMonarchViz.pdf",7,5)
panelview(lnDemons~Monarch,data=Demos,theme.bw=TRUE,
          by.timing=FALSE,index=c("ccode","Year"),
          color=c("white","darkgreen","red"),
          legend.labs=c("Missing","Not Monarchy","Monarchy"),
          main=" ",ylab="Country Code",axis.lab.gap=c(5,5),
          background="white")
dev.off()


# Variation, within and between:

with(Demos, describe(POLITY)) # all variation

pdf("POLITYAll.pdf",6,5)
par(mar=c(4,4,2,2))
plot(density(Demos$POLITY,na.rm=TRUE),
     main="",xlab="POLITY Score",
     lwd=2)
dev.off()

POLITYmeans <- ddply(Demos,.(ccode),summarise,
                     POLITYmean = mean(POLITY))

with(POLITYmeans, describe(POLITYmean)) # "between" variation

pdf("POLITYbetween.pdf",6,5)
par(mar=c(4,4,2,2))
plot(density(POLITYmeans$POLITYmean,na.rm=TRUE),
     main="",xlab="Mean POLITY Scores",
     lwd=2)
dev.off()

Demos <- ddply(Demos, .(ccode), mutate,
               POLITYmean = mean(POLITY))
Demos$POLITYwithin <- with(Demos, POLITY-POLITYmean)

with(Demos, describe(POLITYwithin))

pdf("POLITYwithin.pdf",6,5)
par(mar=c(4,4,2,2))
plot(density(Demos$POLITYwithin,na.rm=TRUE),
     main="",xlab="POLITY: Within-Country Variation",
     lwd=2)
abline(v=0,lty=2)
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Unit Effects models:

# Pooled OLS:

OLS<-plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar, 
         data=PDF,model="pooling")

summary(OLS)

# "Fixed" / within effects:

FE<-plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar, 
        data=PDF,effect="individual",model="within")

summary(FE)

# Make a table:

Table1 <- stargazer(OLS,FE,
                    title="Models of Demonstrations",
                    column.separate=c(1,1),align=TRUE,
                    dep.var.labels.include=FALSE,
                    dep.var.caption="",
                    covariate.labels=c("POLITY","POLITY Squared",
                                       "ln(GDP)","Monarch","Cold War"),
                    header=FALSE,model.names=FALSE,
                    model.numbers=FALSE,multicolumn=FALSE,
                    object.names=TRUE,notes.label="",
                    out="UFX1.tex")

# Between effects:

BE<-plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar, 
        data=PDF, effect="individual",model="between")

summary(BE)

Table2 <- stargazer(OLS,FE,BE,
                    title="Models of Demonstrations",
                    column.separate=c(1,1,1),align=TRUE,
                    dep.var.labels.include=FALSE,
                    dep.var.caption="",
                    covariate.labels=c("POLITY","POLITY Squared",
                                       "ln(GDP)","Monarch","Cold War"),
                    header=FALSE,model.names=FALSE,
                    model.numbers=FALSE,multicolumn=FALSE,
                    object.names=TRUE,notes.label="",
                    out="UFX2.tex")

# REs:

RE<-plm(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar,
        data=PDF, effect="individual", model="random")

summary(RE)

AltRE<-lmer(lnDemons~POLITY+I(POLITY^2)+lnGDP+Monarch+ColdWar+
              (1|ccode), data=Demos)

summary(AltRE)

Table3 <- stargazer(OLS,FE,BE,RE,
                    title="Models of Demonstrations",
                    column.separate=c(1,1,1,1),align=TRUE,
                    dep.var.labels.include=FALSE,
                    dep.var.caption="",
                    covariate.labels=c("POLITY","POLITY Squared",
                                       "ln(GDP)","Monarch","Cold War"),
                    header=FALSE,model.names=FALSE,
                    model.numbers=FALSE,multicolumn=FALSE,
                    object.names=TRUE,notes.label="",
                    out="UFX3.tex")


