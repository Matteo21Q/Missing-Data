##############################################################
####  Generate table of results. GLMER and CLMM           ####
##############################################################

# Set working directory:
setwd("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Missing data/SMC-jomo/Tutorial+Simulations paper/R script + Data")

# Initialise table of results
Table.of.results<-data.frame(matrix(NA,12,17))
Table.of.results[c(1,7),1]<-c("True value")
Table.of.results[c(2,8),1]<-c("Full Data")
Table.of.results[c(3,9),1]<-c("Complete Records")
Table.of.results[c(4,10),1]<-c("JM-Hom (5)")
Table.of.results[c(5,11),1]<-c("JM-Het (6)")
Table.of.results[c(6,12),1]<-c("SMC-JM (7)")

# set methods acronyms:
meths<-c("FD","CR","MI.hom","MI.het","MI.SMC")

##############################################
#### Load results GLMER                   ####

load("Results_GLMER.RData")

# Set true values:
Table.of.results[1,c(2,7,12)]<-beta.smc
Table.of.results[1,c(6,11,6)]<-95.0
Table.of.results[1,17]<-reff.smc


# Add simulation results:
for ( i in 1: length(meths)) {
  
  Table.of.results[1+i,c(2,7,12)]<-apply(get(paste("est",meths[i],sep=".")),2,mean)
  Table.of.results[1+i,c(3,8,13)]<-abs((Table.of.results[1+i,c(2,7,12)]-beta.smc)/beta.smc)*100
  Table.of.results[1+i,c(4,9,14)]<-apply(get(paste("se",meths[i],sep=".")),2,mean)
  Table.of.results[1+i,c(5,10,15)]<-apply(get(paste("est",meths[i],sep=".")),2,sd)
  Table.of.results[1+i,c(6,11,16)]<-apply(get(paste("cov",meths[i],sep=".")),2,mean)*100
  Table.of.results[1+i,c(17)]<-mean(get(paste("rev",meths[i],sep=".")))
  
} 

##############################################
#### Load results CLMM                    ####

load("Results_clmm.RData")

# Set true values:
Table.of.results[7,c(2,3,4,7,12)]<-beta.smc

# Add simulation results:
for ( i in 1: length(meths)) {
  
  Table.of.results[7+i,c(7,12)]<-apply(get(paste("est",meths[i],sep=".")),2,mean, na.rm=T)[4:5]
  Table.of.results[7+i,c(8,13)]<-abs((Table.of.results[7+i,c(7,12)]-beta.smc[4:5])/beta.smc[4:5])*100
  Table.of.results[7+i,c(9,14)]<-apply(get(paste("se",meths[i],sep=".")),2,mean, na.rm=T)[4:5]
  Table.of.results[7+i,c(10,15)]<-apply(get(paste("est",meths[i],sep=".")),2,sd)[4:5]
  Table.of.results[7+i,c(11,16)]<-apply(get(paste("cov",meths[i],sep=".")),2,mean, na.rm=T)[4:5]*100
  Table.of.results[7+i,c(2:4)]<-apply(get(paste("est",meths[i],sep=".")),2,mean)[1:3]
  Table.of.results[1+i,c(17)]<-mean(get(paste("rev",meths[i],sep=".")))
  
} 

Table.of.results[,2:17]<-round(Table.of.results[,2:17],3)
View(Table.of.results)
colnames(Table.of.results)<-c("Method","Mean", "% rel. bias", "Model SE", "Emp SE","% Cov","Mean ", "% rel. bias ", "Model SE ", "Emp SE ","% Cov "," Mean", " % rel. bias", " Model SE", " Emp SE"," % Cov", " Mean ")
row.names(Table.of.results)<-NULL
write.table(Table.of.results, "Table_i.txt", sep=",")