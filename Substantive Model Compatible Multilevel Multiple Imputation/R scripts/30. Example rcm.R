##############################################
####  CRT Data analysis.  SMC-JOMO paper. ####
##############################################

# Load libraries

library(jomo)
library(mitml)
library(lme4)
library(lmerTest)
library(readstata13)

# Set working directory

setwd("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Missing data/SMC-jomo/Tutorial+Simulations paper/R script + Data")

# Set seed

set.seed(1)

# Load data

kemri<-read.dta13("orig_data.dta")

##########################################
#### Fit complete records:          ######
##########################################

analysis.fml<-as.formula(completion~par+child_gender+gender+yr_exp+I(yr_exp^2)+
                           hospstatus+kemri_cme+parXkemri_cme+(1+par|hosp:clinician)+
                           (1|hosp) )
fit.CR<-lmer(analysis.fml, data=kemri )
summary(fit.CR)

###########################################
#### Check convergence of MCMC       ######
###########################################

kemri$by<-kemri$hospstatus
kemri[(kemri$hospstatus==0) & (kemri$par==1),"by"]<-2
kemri[(kemri$hospstatus==1) & (kemri$par==1),"by"]<-3
kemri$hosp<-factor(kemri$hosp)
kemri$gender<-factor(kemri$gender)
kemri$child_gender<-factor(kemri$child_gender)
kemri$kemri_cme<-factor(kemri$kemri_cme)


imp.fml14<-list(as.formula(completion+child_gender~hosp+(1|clinician)),
                as.formula(gender+yr_exp+kemri_cme~hosp))
imp.fml58<-list(as.formula(completion+child_gender~hosp+(1|clinician)),
                as.formula(gender+yr_exp~hosp))

kemri14<-kemri[kemri$hosp%in%1:4,]
kemri58<-kemri[kemri$hosp%in%5:8,]
kemri14$hosp<-factor(kemri14$hosp)
kemri58$hosp<-factor(kemri58$hosp)

# imp.list14<-jomoImpute(data=kemri14, formula=imp.fml14, n.burn=10000, n.iter=2000, group="par")
# imp.list14b<-jomoImpute(data=kemri14, formula=imp.fml14, n.burn=10000, n.iter=2000, group="par",random.L1 = "mean")
# 
# imp.list58<-jomoImpute(data=kemri58, formula=imp.fml58, n.burn = 1000, n.iter = 1000, group="par")
# imp.list58b<-jomoImpute(data=kemri58, formula=imp.fml58, n.burn= 12000, n.iter= 2000, group="par",random.L1 = "mean")


#################################################
####  Impute with homoscedastic model        ####
#################################################

kemri58$const<-1
khosps<-data.frame(model.matrix(~kemri58$hosp))
colnames(khosps)<-c("hosp58","hosp6","hosp7","hosp8")
kemri58<-data.frame(kemri58,khosps)
kemri1<-kemri58[kemri58$par==1,]
kemri0<-kemri58[kemri58$par==0,]

Y<-kemri1[,c("completion", "child_gender")]
X<-kemri1[,c("hosp58","hosp6","hosp7","hosp8")]
Y2<-kemri1[,c("yr_exp", "gender"), drop=F]
X2<-X
clus<-kemri1$clinician
setwd("C:/Users/rmjlmqu/Documents")
imps<-jomo(Y=Y, Y2=Y2, X=X, X2=X2, clus=clus, nburn = 12000, nbetween=2000, nimp=10)

Y<-kemri0[,c("completion", "child_gender")]
X<-kemri0[,c("hosp58","hosp6","hosp7","hosp8")]
Y2<-kemri0[,c("yr_exp", "gender"), drop=F]
X2<-X
clus<-kemri0$clinician
imps2<-jomo(Y=Y, Y2=Y2, X=X, X2=X2, clus=clus, nburn = 12000, nbetween=2000, nimp=10)

kemri14$const<-1
khosps14<-data.frame(model.matrix(~kemri14$hosp))
colnames(khosps14)<-c("hosp14","hosp2","hosp3","hosp4")
kemri14<-data.frame(kemri14,khosps14)
kemri1<-kemri14[kemri14$par==1,]
kemri0<-kemri14[kemri14$par==0,]

Y<-kemri1[,c("completion", "child_gender")]
X<-kemri1[,c("hosp14","hosp2","hosp3","hosp4")]
Y2<-kemri1[,c("yr_exp", "gender", "kemri_cme"), drop=F]
X2<-X
clus<-kemri1$clinician
imps3<-jomo(Y=Y, Y2=Y2, X=X, X2=X2, clus=clus, nburn = 2000, nbetween=1000, nimp=10)

Y<-kemri0[,c("completion", "child_gender")]
X<-kemri0[,c("hosp14","hosp2","hosp3","hosp4")]
Y2<-kemri0[,c("yr_exp", "gender", "kemri_cme"), drop=F]
X2<-X
clus<-kemri0$clinician
imps4<-jomo(Y=Y, Y2=Y2, X=X, X2=X2, clus=clus, nburn = 2000, nbetween=1000, nimp=10)

imps$hospstatus<-0
imps2$hospstatus<-0
imps3$hospstatus<-1
imps4$hospstatus<-1
imps$par<-1
imps2$par<-0
imps3$par<-1
imps4$par<-0
imps$hosp<-ifelse(imps$hosp6==1,6,ifelse(imps$hosp7==1,7,
                                         ifelse(imps$hosp8==1,8,5)))
imps2$hosp<-ifelse(imps2$hosp6==1,6,ifelse(imps2$hosp7==1,7,
                                           ifelse(imps2$hosp8==1,8,5)))
imps3$hosp<-ifelse(imps3$hosp2==1,2,ifelse(imps3$hosp3==1,3,
                                           ifelse(imps3$hosp4==1,4,1)))
imps4$hosp<-ifelse(imps4$hosp2==1,2,ifelse(imps4$hosp3==1,3,
                                           ifelse(imps4$hosp4==1,4,1)))
imps<-imps[,-which(colnames(imps)%in%c("hosp58","hosp6","hosp7","hosp8"))]
imps2<-imps2[,-which(colnames(imps2)%in%c("hosp58","hosp6","hosp7","hosp8"))]
imps3<-imps3[,-which(colnames(imps3)%in%c("hosp14","hosp2","hosp3","hosp4"))]
imps4<-imps4[,-which(colnames(imps4)%in%c("hosp14","hosp2","hosp3","hosp4"))]

imps$kemri_cme<-0
imps2$kemri_cme<-0

imputs<-data.frame(rbind(imps,imps2,imps3,imps4))
imputs<-imputs[order(imputs$Imputation),]
imputs$id<-1:nrow(kemri)

impList <- jomo2mitml.list(imputs)

fit.m<-with.mitml.list(impList, lmer(as.formula(completion~par+child_gender+gender+yr_exp+I(yr_exp^2)+
                                                  hospstatus+kemri_cme+par:kemri_cme+(1+par|hosp:clus)+
                                                  (1|hosp) )))
fit.MI<-testEstimates(fit.m, var.comp = T)
summary(fit.MI)

fit.MI$estimates[,c(1,2,5)]
coef(summary(fit.CR))[,c(1,2,5)]

#################################################
####  Impute with heteroscedastic model      ####
#################################################

kemri58$const<-1
khosps<-data.frame(model.matrix(~kemri58$hosp))
colnames(khosps)<-c("hosp58","hosp6","hosp7","hosp8")
kemri58<-data.frame(kemri58,khosps)
kemri1<-kemri58[kemri58$par==1,]
kemri0<-kemri58[kemri58$par==0,]

Y<-kemri1[,c("completion", "child_gender")]
X<-kemri1[,c("hosp58","hosp6","hosp7","hosp8")]
Y2<-kemri1[,c("yr_exp", "gender"), drop=F]
X2<-X
clus<-kemri1$clinician
setwd("C:/Users/rmjlmqu/Documents")
imps<-jomo(Y=Y, Y2=Y2, X=X, X2=X2, clus=clus, nburn = 12000, nbetween=2000, nimp=10, meth="random") 

Y<-kemri0[,c("completion", "child_gender")]
X<-kemri0[,c("hosp58","hosp6","hosp7","hosp8")]
Y2<-kemri0[,c("yr_exp", "gender"), drop=F]
X2<-X
clus<-kemri0$clinician
imps2<-jomo(Y=Y, Y2=Y2, X=X, X2=X2, clus=clus, nburn = 12000, nbetween=2000, nimp=10, meth="random")

kemri14$const<-1
khosps14<-data.frame(model.matrix(~kemri14$hosp))
colnames(khosps14)<-c("hosp14","hosp2","hosp3","hosp4")
kemri14<-data.frame(kemri14,khosps14)
kemri1<-kemri14[kemri14$par==1,]
kemri0<-kemri14[kemri14$par==0,]

Y<-kemri1[,c("completion", "child_gender")]
X<-kemri1[,c("hosp14","hosp2","hosp3","hosp4")]
Y2<-kemri1[,c("yr_exp", "gender", "kemri_cme"), drop=F]
X2<-X
clus<-kemri1$clinician
imps3<-jomo(Y=Y, Y2=Y2, X=X, X2=X2, clus=clus, nburn = 12000, nbetween=2000, nimp=10, meth = "random")

Y<-kemri0[,c("completion", "child_gender")]
X<-kemri0[,c("hosp14","hosp2","hosp3","hosp4")]
Y2<-kemri0[,c("yr_exp", "gender", "kemri_cme"), drop=F]
X2<-X
clus<-kemri0$clinician
imps4<-jomo(Y=Y, Y2=Y2, X=X, X2=X2, clus=clus, nburn = 12000, nbetween=2000, nimp=10, meth = "random")

imps$hospstatus<-0
imps2$hospstatus<-0
imps3$hospstatus<-1
imps4$hospstatus<-1
imps$par<-1
imps2$par<-0
imps3$par<-1
imps4$par<-0
imps$hosp<-ifelse(imps$hosp6==1,6,ifelse(imps$hosp7==1,7,
                                         ifelse(imps$hosp8==1,8,5)))
imps2$hosp<-ifelse(imps2$hosp6==1,6,ifelse(imps2$hosp7==1,7,
                                           ifelse(imps2$hosp8==1,8,5)))
imps3$hosp<-ifelse(imps3$hosp2==1,2,ifelse(imps3$hosp3==1,3,
                                           ifelse(imps3$hosp4==1,4,1)))
imps4$hosp<-ifelse(imps4$hosp2==1,2,ifelse(imps4$hosp3==1,3,
                                           ifelse(imps4$hosp4==1,4,1)))
imps<-imps[,-which(colnames(imps)%in%c("hosp58","hosp6","hosp7","hosp8"))]
imps2<-imps2[,-which(colnames(imps2)%in%c("hosp58","hosp6","hosp7","hosp8"))]
imps3<-imps3[,-which(colnames(imps3)%in%c("hosp14","hosp2","hosp3","hosp4"))]
imps4<-imps4[,-which(colnames(imps4)%in%c("hosp14","hosp2","hosp3","hosp4"))]

imps$kemri_cme<-0
imps2$kemri_cme<-0

imputs2<-data.frame(rbind(imps,imps2,imps3,imps4))
imputs2<-imputs2[order(imputs2$Imputation),]
imputs2$id<-1:nrow(kemri)

impList2 <- jomo2mitml.list(imputs2)

fit.m2<-with.mitml.list(impList2, lmer(as.formula(completion~par+child_gender+gender+yr_exp+I(yr_exp^2)+
                                                    hospstatus+kemri_cme+par:kemri_cme+(1+par|hosp:clus)+
                                                    (1|hosp) )))
fit.MI.het<-testEstimates(fit.m2, var.comp = T)
summary(fit.MI.het)

fit.MI.het$estimates[,c(1,2,5)]
coef(summary(fit.CR))[,c(1,2,5)]

##############################################
####  SMC-JOMO                            ####
##############################################


smcimp.fml<-as.formula(completion~par+child_gender+gender+yr_exp+I(yr_exp^2)+
                         hospstatus+kemri_cme+par:kemri_cme
                       +(1+par|clinician))

kemri.sub<-kemri[,c("completion","par","child_gender","gender","yr_exp",
                    "hospstatus","kemri_cme","hosp",
                    "clinician")]
levels<-c(1,1,1,2,2,2,2,2,1)
imp.SMC<-jomo.lmer(smcimp.fml,kemri.sub,levels, nburn = 2000, nbetween=1000, nimp=10)

imp.SMC$hosp<-kemri$hosp


implist3<-jomo2mitml.list(imp.SMC)
fit.m3<-with.mitml.list(implist3, lmer(as.formula(completion~par+child_gender+gender+yr_exp+I(yr_exp^2)+
                                                    hospstatus+kemri_cme+par:kemri_cme+(1+par|hosp:clus)+
                                                    (1|hosp) )))
fit.MI.SMC<-testEstimates(fit.m3, var.comp = T)
summary(fit.MI.SMC)

coef(summary(fit.CR))[,c(1,2,5)]
fit.MI$estimates[,c(1,2,5)]
fit.MI.het$estimates[,c(1,2,5)]
fit.MI.SMC$estimates[,c(1,2,5)]

# save(fit.CR, fit.MI, fit.MI.het, fit.MI.SMC, file="Results_Example_rcm.RData")
