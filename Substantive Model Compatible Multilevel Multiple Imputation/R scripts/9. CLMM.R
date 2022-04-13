##############################################
####  CLMM Scenario. SMC-JOMO paper.      ####
##############################################

# Load libraries

library(ordinal)
library(jomo)
library(mitml)
library(lme4)
library(lmerTest)
library(readstata13)
library(MASS)
library(mice)
library(boot)
library(ordinalimputation)
library(miceadds)

# Set working directory

# setwd("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Missing data/SMC-jomo/Tutorial+Simulations paper/R script + Data")

# Set seed

set.seed(1)

# Load data

kemri<-read.dta13("orig_data.dta")

# Estimate parameters for data generatin mechanism using complete records:

kemri$completion<-kemri$completion/100
kemri$completord<-ifelse(kemri$completion<0.5,1,ifelse(kemri$completion<0.8,2,
                                                       ifelse(kemri$completion<0.9,3,4)))
kemri$scaled.age<-scale(kemri$age)
kemri$completord<-factor(kemri$completord)
kemri$child_gender<-factor(kemri$child_gender)
analysis.fml<-as.formula(completord~child_gender+scaled.age+(1|clinician))
fit.CR<-clmm(analysis.fml, data=kemri, family = binomial )
summary(fit.CR)

Y<-kemri[,c("child_gender","scaled.age")]
clus<-kemri[,c("clinician")]
jointmodel<-jomo.MCMCchain(Y=Y, clus=clus, nburn=1000)

# Set simulation parameters

n.sim<-1000 # Number of simulations
n.obs<-1500 # Number of observations
n.clus<-20 # Number of clusters
p.miss<-0.4
n.imp<-10
n.burn<-500
n.between<-500
mu<-apply(jointmodel$collectbeta, c(1,2),mean)
omega<-apply(jointmodel$collectomega, c(1,2),mean)
psi<-apply(jointmodel$collectcovu,c(1,2),mean)
beta.smc<-coef(summary(fit.CR))[,1]
reff.smc<-as.numeric(VarCorr(fit.CR))/5

#Initialise matrices of results:

est.FD<-matrix(NA,n.sim,length(beta.smc))
est.CR<-matrix(NA,n.sim,length(beta.smc))
est.MI.hom<-matrix(NA,n.sim,length(beta.smc))
est.MI.het<-matrix(NA,n.sim,length(beta.smc))
est.MI.SMC<-matrix(NA,n.sim,length(beta.smc))
se.FD<-matrix(NA,n.sim,length(beta.smc))
se.CR<-matrix(NA,n.sim,length(beta.smc))
se.MI.hom<-matrix(NA,n.sim,length(beta.smc))
se.MI.het<-matrix(NA,n.sim,length(beta.smc))
se.MI.SMC<-matrix(NA,n.sim,length(beta.smc))
cov.FD<-matrix(NA,n.sim,length(beta.smc))
cov.CR<-matrix(NA,n.sim,length(beta.smc))
cov.MI.hom<-matrix(NA,n.sim,length(beta.smc))
cov.MI.het<-matrix(NA,n.sim,length(beta.smc))
cov.MI.SMC<-matrix(NA,n.sim,length(beta.smc))
rev.FD<-rep(NA,n.sim)
rev.CR<-rep(NA,n.sim)
rev.MI.hom<-rep(NA,n.sim)
rev.MI.het<-rep(NA,n.sim)
rev.MI.SMC<-rep(NA,n.sim)
analysis.model<-as.formula(Y~child_gender+scaled.age+(1|clus))


for (i in 1:n.sim) {
  
  clus<-sample(1:n.clus,n.obs,replace = TRUE)
  v<-mvrnorm(n.clus,c(0,0),psi)
  Xvar<-mvrnorm(n.obs,mu,omega)+v[clus,]
  X<-data.frame(as.numeric(Xvar[,2]<0), Xvar[,1])
  
  colnames(X)<-c("child_gender","scaled.age")
  
  u<-rnorm(n.clus,0,sqrt(reff.smc))
  pY1<-inv.logit(beta.smc[1]-as.matrix(X)%*%beta.smc[4:5]-u[clus])
  pY2<-inv.logit(beta.smc[2]-as.matrix(X)%*%beta.smc[4:5]-u[clus])-pY1
  pY3<-inv.logit(beta.smc[3]-as.matrix(X)%*%beta.smc[4:5]-u[clus])-pY2-pY1
  pY4<-1-pY1-pY2-pY3
  Y<-rep(NA,n.obs)
  for ( jjj in 1:n.obs) {
    Y[jjj]<-which(rmultinom(1,1,prob=c(pY1[jjj],pY2[jjj],pY3[jjj],pY4[jjj]))==1)
  }
  data.sim<-data.frame(Y,X,clus)
  data.sim$Y<-factor(data.sim$Y)
  data.sim$child_gender<-factor(data.sim$child_gender)
  
  # Fit model on full data
  
  fit.FD<-clmm(analysis.model, data = data.sim)
  est.FD[i,]<-coef(summary(fit.FD))[,1]
  se.FD[i,]<-coef(summary(fit.FD))[,2]
  cov.FD[i,]<-((est.FD[i,]-se.FD[i,]*qnorm(0.975)<beta.smc)&(est.FD[i,]+se.FD[i,]*qnorm(0.975)>beta.smc))
  rev.FD[i]<-as.numeric(VarCorr(fit.FD))
  
  # Introduce MAR in all variables using ampute function in mice:
  
  data.miss<-data.sim[,c("Y","child_gender","scaled.age")]
  p.miss.cg.MAR<-ifelse(data.miss$Y==1,0.1,ifelse(data.miss$Y==2,0.15,
                                                  ifelse(data.miss$Y==3,0.2,0.25)))
  p.miss.age.MAR<-ifelse(data.miss$Y==1,0.1,ifelse(data.miss$Y==2,0.15,
                                                   ifelse(data.miss$Y==3,0.2,0.25)))
  
  for (j in 1:n.obs) {
    if (runif(1,0,1)<p.miss.cg.MAR[j]) data.miss$child_gender[j]<-NA
    if (runif(1,0,1)<p.miss.age.MAR[j]) data.miss$scaled.age[j]<-NA
  }
  data.miss$clus<-data.sim$clus
  data.miss$child_gender<-factor(data.miss$child_gender)
  data.miss$Y<-factor(data.miss$Y)
  
  fit.CR<-clmm(analysis.model, data = data.miss)
  est.CR[i,]<-coef(summary(fit.CR))[,1]
  se.CR[i,]<-coef(summary(fit.CR))[,2]
  cov.CR[i,]<-((est.CR[i,]-se.CR[i,]*qnorm(0.975)<beta.smc)&(est.CR[i,]+se.CR[i,]*qnorm(0.975)>beta.smc))
  rev.CR[i]<-as.numeric(VarCorr(fit.CR))
  
  # Homoscedastic JM MI
  
  imp.model<-as.formula(Y+scaled.age+child_gender~(1|clus))
  imps<-jomoImpute(data.miss, formula = imp.model, n.burn=n.burn,
                   n.iter=n.between, m=n.imp)
  implist <- mitmlComplete(imps, print=1:n.imp)
  
  fits<-matrix(NA,n.imp,length(beta.smc))
  ses<-matrix(NA, n.imp, length(beta.smc))
  ran.var<-rep(NA,n.sim)
  for (k in 1:n.imp) {
    fit.m<-clmm(Y~child_gender+scaled.age+(1|clus), data = implist[[k]])
    fits[k,]<-coef(summary(fit.m))[,1]
    ses[k,]<-coef(summary(fit.m))[,2]
    ran.var[k]<-as.numeric(VarCorr(fit.m))
  }
  
  est.MI.hom[i,]<-apply(fits,2,mean)
  se.MI.hom[i,]<-sqrt(apply(ses,2,mean)^2+(1+1/n.imp)*apply(fits,2,var))
  cov.MI.hom[i,]<-((est.MI.hom[i,]-se.MI.hom[i,]*qnorm(0.975)<beta.smc)&(est.MI.hom[i,]+se.MI.hom[i,]*qnorm(0.975)>beta.smc))
  rev.MI.hom[i]<- mean(ran.var)
  
  # Heteroscedastic JM MI
  
  imps2<-jomoImpute(data.miss, formula = imp.model, n.burn=n.burn,
                    n.iter=n.between, m=n.imp, random.L1 = "mean")
  implist2 <- mitmlComplete(imps2, print=1:n.imp)
  
  fits2<-matrix(NA,n.imp,length(beta.smc))
  ses2<-matrix(NA, n.imp, length(beta.smc))
  ran.var2<-rep(NA,n.sim)
  for (k in 1:n.imp) {
    fit.m2<-clmm(Y~child_gender+scaled.age+(1|clus), data = implist2[[k]])
    fits2[k,]<-coef(summary(fit.m2))[,1]
    ses2[k,]<-coef(summary(fit.m2))[,2]
    ran.var2[k]<-as.numeric(VarCorr(fit.m2))
  }
  
  est.MI.het[i,]<-apply(fits2,2,mean)
  se.MI.het[i,]<-sqrt(apply(ses2,2,mean)^2+(1+1/n.imp)*apply(fits2,2,var))
  cov.MI.het[i,]<-((est.MI.het[i,]-se.MI.het[i,]*qnorm(0.975)<beta.smc)&(est.MI.het[i,]+se.MI.het[i,]*qnorm(0.975)>beta.smc))
  rev.MI.het[i]<- mean(ran.var2)
  
  # SMC JM MI
  
  imps3<-jomo.clmm(data.miss, formula = analysis.model, nburn=n.burn,
                   nbetween=n.between, nimp=n.imp)
  implist3 <- jomo2mitml.list(imps3)
  
  fits3<-matrix(NA,n.imp,length(beta.smc))
  ses3<-matrix(NA, n.imp, length(beta.smc))
  ran.var3<-rep(NA,n.sim)
  for (k in 1:n.imp) {
    fit.m3<-clmm(Y~child_gender+scaled.age+(1|clus), data = implist3[[k]])
    fits3[k,]<-coef(summary(fit.m3))[,1]
    ses3[k,]<-coef(summary(fit.m3))[,2]
    ran.var3[k]<-as.numeric(VarCorr(fit.m3))
  }
  
  est.MI.SMC[i,]<-apply(fits3,2,mean)
  se.MI.SMC[i,]<-sqrt(apply(ses3,2,mean)^2+(1+1/n.imp)*apply(fits3,2,var))
  cov.MI.SMC[i,]<-((est.MI.SMC[i,]-se.MI.SMC[i,]*qnorm(0.975)<beta.smc)&(est.MI.SMC[i,]+se.MI.SMC[i,]*qnorm(0.975)>beta.smc))
  rev.MI.SMC[i]<- mean(ran.var3)
  
  cat("Simulation ", i," completed\n")
}

apply(est.FD,2,mean)
apply(est.CR,2,mean)
apply(est.MI.hom,2,mean)
apply(est.MI.het,2,mean)
apply(est.MI.SMC,2,mean)
apply(se.FD,2,mean)
apply(se.CR,2,mean)
apply(se.MI.hom,2,mean)
apply(se.MI.het,2,mean)
apply(se.MI.SMC,2,mean)
apply(se.FD,2,sd)
apply(se.CR,2,sd)
apply(se.MI.hom,2,sd)
apply(se.MI.het,2,sd)
apply(se.MI.SMC,2,sd)
apply(cov.FD,2,mean)
apply(cov.CR,2,mean)
apply(cov.MI.hom,2,mean)
apply(cov.MI.het,2,mean)
apply(cov.MI.SMC,2,mean)

# save.image( file="Results_clmm.RData")
