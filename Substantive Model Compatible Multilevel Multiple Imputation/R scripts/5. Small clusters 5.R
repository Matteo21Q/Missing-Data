##############################################
####  Sensitivity to small (~5) clusters. ####
##############################################

# Load libraries

library(jomo)
library(mitml)
library(lme4)
library(lmerTest)
library(readstata13)
library(MASS)
library(mice)

# Set working directory

# setwd("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Missing data/SMC-jomo/Tutorial+Simulations paper/R script + Data")

# Set seed

set.seed(1)

# Load data

kemri<-read.dta13("orig_data.dta")

# Estimate parameters for data generatin mechanism using complete records:

kemri$completion<-kemri$completion/100
kemri$scaled.age<-scale(kemri$age)
kemri$par<-factor(kemri$par)
analysis.fml<-as.formula(completion~par+scaled.age+(1+par|clinician))
fit.CR<-lmer(analysis.fml, data=kemri )
summary(fit.CR)

Y<-kemri[,c("par","scaled.age")]
clus<-kemri[,c("clinician")]
jointmodel<-jomo.MCMCchain(Y=Y, clus=clus, nburn=1000)

# Set simulation parameters

n.sim<-1000 # Number of simulations
n.obs<-1500 # Number of observations
n.clus<-300 # Number of clusters
p.miss<-0.4
n.imp<-10
n.burn<-500
n.between<-500
mu<-apply(jointmodel$collectbeta, c(1,2),mean)
omega<-apply(jointmodel$collectomega, c(1,2),mean)
psi<-apply(jointmodel$collectcovu,c(1,2),mean)
beta.smc<-coef(summary(fit.CR))[,1]
res.smc<-attr(VarCorr(fit.CR),"sc")
reff.smc<-matrix(as.numeric(VarCorr(fit.CR)$clinician), 2,2)

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
rev.FD<-array(NA,c(n.sim,nrow(reff.smc),ncol(reff.smc)))
rev.CR<-array(NA,c(n.sim,nrow(reff.smc),ncol(reff.smc)))
rev.MI.hom<-array(NA,c(n.sim,nrow(reff.smc),ncol(reff.smc)))
rev.MI.het<-array(NA,c(n.sim,nrow(reff.smc),ncol(reff.smc)))
rev.MI.SMC<-array(NA,c(n.sim,nrow(reff.smc),ncol(reff.smc)))
analysis.model<-as.formula(Y~par+scaled.age+(1+par|clus))


for (i in 1:n.sim) {
  
  clus<-sample(1:n.clus,n.obs,replace = TRUE)
  v<-mvrnorm(n.clus,c(0,0),psi)
  Xvar<-mvrnorm(n.obs,mu,omega)+v[clus,]
  X<-data.frame(1, as.numeric(Xvar[,2]<0), Xvar[,1])
  Z<-data.frame(1, as.numeric(Xvar[,2]<0))
  
  colnames(X)<-c("cons","par","scaled.age")
  colnames(Z)<-c("cons","par")
  
  u<-mvrnorm(n.clus,c(0,0),reff.smc)
  Y<-rnorm(n.obs,as.matrix(X)%*%beta.smc,res.smc)+apply(as.matrix(Z)*u[clus,],1,sum)
  data.sim<-data.frame(Y,X,clus)
  data.sim$par<-factor(data.sim$par)
  
  # Fit model on full data
  
  fit.FD<-lmer(analysis.model, data = data.sim)
  est.FD[i,]<-coef(summary(fit.FD))[,1]
  se.FD[i,]<-coef(summary(fit.FD))[,2]
  cov.FD[i,]<-((est.FD[i,]-se.FD[i,]*qnorm(0.975)<beta.smc)&(est.FD[i,]+se.FD[i,]*qnorm(0.975)>beta.smc))
  rev.FD[i,,]<-matrix(as.numeric(VarCorr(fit.FD)$clus), 2,2)
  
  # Introduce MAR in all variables using ampute function in mice:

  data.miss<-data.sim[,c("Y","par","scaled.age")]
  p.miss.par.MAR<-(1+exp(2-1.5*data.miss$Y))^(-1)
  p.miss.age.MAR<-(1+exp(2-1.5*data.miss$Y))^(-1)
  
  for (j in 1:n.obs) {
    if (runif(1,0,1)<p.miss.par.MAR[j]) data.miss$par[j]<-NA
    if (runif(1,0,1)<p.miss.age.MAR[j]) data.miss$scaled.age[j]<-NA
  }
  data.miss$clus<-data.sim$clus
  data.miss$par<-factor(data.miss$par)
  
  fit.CR<-lmer(analysis.model, data = data.miss)
  est.CR[i,]<-coef(summary(fit.CR))[,1]
  se.CR[i,]<-coef(summary(fit.CR))[,2]
  cov.CR[i,]<-((est.CR[i,]-se.CR[i,]*qnorm(0.975)<beta.smc)&(est.CR[i,]+se.CR[i,]*qnorm(0.975)>beta.smc))
  rev.CR[i,,]<-matrix(as.numeric(VarCorr(fit.CR)$clus), 2,2)

  # Homoscedastic JM MI
  
  imp.model<-as.formula(Y+scaled.age+par~(1|clus))
  imps<-jomoImpute(data.miss, formula = imp.model, n.burn=n.burn, 
                   n.iter=n.between, m=n.imp)
  implist <- mitmlComplete(imps, print=1:n.imp)
  
  fit.m<-with(implist, lmer(Y~par+scaled.age+(1+par|clus)))
  fit.JM.hom<-testEstimates(fit.m, var.comp=T)
  est.MI.hom[i,]<-fit.JM.hom$estimates[,1]
  se.MI.hom[i,]<-fit.JM.hom$estimates[,2]
  cov.MI.hom[i,]<-((est.MI.hom[i,]-se.MI.hom[i,]*qnorm(0.975)<beta.smc)&(est.MI.hom[i,]+se.MI.hom[i,]*qnorm(0.975)>beta.smc))
  rev.MI.hom[i,,]<- matrix(c(fit.JM.hom$var.comp[1], fit.JM.hom$var.comp[2], fit.JM.hom$var.comp[2],fit.JM.hom$var.comp[3]),2,2)
  
  # Heteroscedastic JM MI
  
  imps2<-jomoImpute(data.miss, formula = imp.model, n.burn=n.burn, 
                   n.iter=n.between, m=n.imp, random.L1 = "mean")
  implist2 <- mitmlComplete(imps2, print=1:n.imp)
  
  fit.m2<-with(implist2, lmer(Y~par+scaled.age+(1+par|clus)))
  fit.JM.het<-testEstimates(fit.m2, var.comp=T)
  est.MI.het[i,]<-fit.JM.het$estimates[,1]
  se.MI.het[i,]<-fit.JM.het$estimates[,2]
  cov.MI.het[i,]<-((est.MI.het[i,]-se.MI.het[i,]*qnorm(0.975)<beta.smc)&(est.MI.het[i,]+se.MI.het[i,]*qnorm(0.975)>beta.smc))
  rev.MI.het[i,,]<- matrix(c(fit.JM.het$var.comp[1], fit.JM.het$var.comp[2], fit.JM.het$var.comp[2],fit.JM.het$var.comp[3]),2,2)
  
  # SMC JM MI
  
  imps3<-jomo.lmer(data.miss, formula = analysis.model, nburn=n.burn, 
                    nbetween=n.between, nimp=n.imp)
  implist3 <- jomo2mitml.list(imps3)
  
  fit.m3<-with(implist3, lmer(Y~par+scaled.age+(1+par|clus)))
  fit.JM.SMC<-testEstimates(fit.m3, var.comp=T)
  est.MI.SMC[i,]<-fit.JM.SMC$estimates[,1]
  se.MI.SMC[i,]<-fit.JM.SMC$estimates[,2]
  cov.MI.SMC[i,]<-((est.MI.SMC[i,]-se.MI.SMC[i,]*qnorm(0.975)<beta.smc)&(est.MI.SMC[i,]+se.MI.SMC[i,]*qnorm(0.975)>beta.smc))
  rev.MI.SMC[i,,]<- matrix(c(fit.JM.SMC$var.comp[1], fit.JM.SMC$var.comp[2], fit.JM.SMC$var.comp[2],fit.JM.SMC$var.comp[3]),2,2)
  
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

# save.image( file="Results_smallclusters 5.RData")
