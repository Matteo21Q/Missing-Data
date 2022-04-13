##############################################################
####  Figure 3. Zip plot                                  ####
##############################################################

# Set working directory and load rsimsum:
setwd("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Missing data/SMC-jomo/Tutorial+Simulations paper/R script + Data")
library(rsimsum)
library(ggplot2)

meths<-c("MI.hom", "MI.het", "MI.SMC")

#####################
####Load data:  ####
#####################

load("Results_Level2.RData")
res<-data.frame(matrix(NA,3000,4))
True<-beta.smc[3]

for ( i in 1:3) {
  res[((i-1)*1000+1):(i*1000),1]<-1:1000
  res[((i-1)*1000+1):(i*1000),2]<-meths[i]
  meanz<-get(paste("est",meths[i],sep="."))[,3]
  res[((i-1)*1000+1):(i*1000),3]<-meanz
  sez<-get(paste("se",meths[i],sep="."))[,3]
  res[((i-1)*1000+1):(i*1000),4]<-sez
  
} 
colnames(res)<-c("dataset", "method", "estimate","se")
View(res)

# COnvert to rsimsum format:
s <- simsum(data = res, estvarname = "estimate", true = True, 
            se = "se", methodvar = "method", ref = "MI.hom", x=T)
autoplot(s, type = "zip", zoom = 0.5)
