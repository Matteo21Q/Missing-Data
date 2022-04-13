##############################################################
####  Generate table of results. Auxiliary                ####
##############################################################

# Set working directory:
setwd("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Missing data/SMC-jomo/Tutorial+Simulations paper/R script + Data")

# Initialise table of results
Table.of.results<-data.frame(matrix(NA,12,19))
Table.of.results[c(1,7),1]<-c("True value")
Table.of.results[c(2,8),1]<-c("Full Data")
Table.of.results[c(3,9),1]<-c("Complete Records")
Table.of.results[c(4,10),1]<-c("JM-Hom (5)")
Table.of.results[c(5,11),1]<-c("JM-Het (6)")
Table.of.results[c(6,12),1]<-c("SMC-JM (7)")

# set methods acronyms:
meths<-c("FD","CR","MI.hom","MI.het","MI.SMC")

##############################################
#### Load results auxiliary 1             ####

load("Results_Auxiliary1.RData")

# Set true values:
Table.of.results[1,c(2,7,12)]<-beta.smc
Table.of.results[1,c(6,11,6)]<-95.0
Table.of.results[1,17]<-reff.smc[1,1]
Table.of.results[1,18]<-reff.smc[1,2]
Table.of.results[1,19]<-reff.smc[2,2]


# Add simulation results:
for ( i in 1: length(meths)) {
  
  Table.of.results[1+i,c(2,7,12)]<-apply(get(paste("est",meths[i],sep=".")),2,mean)
  Table.of.results[1+i,c(3,8,13)]<-abs((Table.of.results[1+i,c(2,7,12)]-beta.smc)/beta.smc)*100
  Table.of.results[1+i,c(4,9,14)]<-apply(get(paste("se",meths[i],sep=".")),2,mean)
  Table.of.results[1+i,c(5,10,15)]<-apply(get(paste("est",meths[i],sep=".")),2,sd)
  Table.of.results[1+i,c(6,11,16)]<-apply(get(paste("cov",meths[i],sep=".")),2,mean)*100
  Table.of.results[1+i,c(17)]<-mean(get(paste("rev",meths[i],sep="."))[,1,1])
  Table.of.results[1+i,c(18)]<-mean(get(paste("rev",meths[i],sep="."))[,1,2])
  Table.of.results[1+i,c(19)]<-mean(get(paste("rev",meths[i],sep="."))[,2,2])
  
} 

##############################################
#### Load results auxiliary 2           ####

load("Results_Auxiliary2.RData")

# Set true values:
Table.of.results[7,c(2,7,12)]<-beta.smc
Table.of.results[7,c(6,11,6)]<-95.0
Table.of.results[7,17]<-reff.smc[1,1]
Table.of.results[7,18]<-reff.smc[1,2]
Table.of.results[7,19]<-reff.smc[2,2]


# Add simulation results:
for ( i in 1: length(meths)) {
  
  Table.of.results[7+i,c(2,7,12)]<-apply(get(paste("est",meths[i],sep=".")),2,mean)
  Table.of.results[7+i,c(3,8,13)]<-abs((Table.of.results[7+i,c(2,7,12)]-beta.smc)/beta.smc)*100
  Table.of.results[7+i,c(4,9,14)]<-apply(get(paste("se",meths[i],sep=".")),2,mean)
  Table.of.results[7+i,c(5,10,15)]<-apply(get(paste("est",meths[i],sep=".")),2,sd)
  Table.of.results[7+i,c(6,11,16)]<-apply(get(paste("cov",meths[i],sep=".")),2,mean)*100
  Table.of.results[7+i,c(17)]<-mean(get(paste("rev",meths[i],sep="."))[,1,1])
  Table.of.results[7+i,c(18)]<-mean(get(paste("rev",meths[i],sep="."))[,1,2])
  Table.of.results[7+i,c(19)]<-mean(get(paste("rev",meths[i],sep="."))[,2,2])
  
} 

Table.of.results[,2:19]<-round(Table.of.results[,2:19],3)
View(Table.of.results)
colnames(Table.of.results)<-c("Method","Mean", "% rel. bias", "Model SE", "Emp SE","% Cov","Mean ", "% rel. bias ", "Model SE ", "Emp SE ","% Cov "," Mean", " % rel. bias", " Model SE", " Emp SE"," % Cov", " Mean ", "  Mean", "Mean  ")
row.names(Table.of.results)<-NULL
write.table(Table.of.results, "Table_l.txt", sep=",")