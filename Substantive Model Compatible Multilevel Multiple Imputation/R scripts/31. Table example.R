##############################################################
####  Generate table of results. Random int and base-case ####
##############################################################

# Set working directory:
setwd("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Missing data/SMC-jomo/Tutorial+Simulations paper/R script + Data")

# Initialise table of results
Table.of.results<-data.frame(matrix(NA,31,13))
Table.of.results[c(1,15),2]<-c("Complete Records")
Table.of.results[c(1,15),5]<-c("JM-Hom")
Table.of.results[c(1,15),8]<-c("JM-Het")
Table.of.results[c(1,15),11]<-c("SMC-JM")

Table.of.results[c(2,16),1]<-"Variable"
Table.of.results[c(2,16),c(2,5,8,11)]<-c("Est.")
Table.of.results[c(2,16),c(3,6,9,12)]<-c("SE")
Table.of.results[c(2,16),c(4,7,10,13)]<-c("p-value")

Table.of.results[c(3:10),1]<-c("Intercept", "PAR use", "female child", "female clinician",
                            "years of experience", "interventional hospital", "CME attendance",
                            "PAR x CME interaction")
Table.of.results[c(17:25),1]<-c("Intercept", "PAR use", "female child", "female clinician",
                                     "years of experience", "years of exp squared", "interventional hospital", "CME attendance",
                                     "PAR x CME interaction")

Table.of.results[c(11,26),1]<-"Variance COmponents"
Table.of.results[c(11,26),c(2,5,8,11)]<-"Est."

Table.of.results[12:14,1]<-c("Hospital", "Clinican", "Child")
Table.of.results[27:31,1]<-c("Hospital", "Clinican Int.", "Clinician slope", "Clinician corr", "Child")

##############################################
#### Load results random intercept model  ####

load("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Missing data/SMC-jomo/Tutorial+Simulations paper/R script + Data/Results_Example_rim.RData")

# Complete records column:
Table.of.results[3:10,2:4]<-round(coef(summary(fit.CR))[,c(1,2,5)],3)
Table.of.results[12:14,2]<-round(c(as.numeric(VarCorr(fit.CR))[c(2,1)], attr(VarCorr(fit.CR), "sc")^2),3)

# JM-Hom column:
Table.of.results[3:10,5:7]<-round(fit.MI$estimates[,c(1,2,5)],3)
Table.of.results[12:14,5]<-round(fit.MI$var.comp[c(2,1,3)],3)

# JM-Het column:
Table.of.results[3:10,8:10]<-round(fit.MI.het$estimates[,c(1,2,5)],3)
Table.of.results[12:14,8]<-round(fit.MI.het$var.comp[c(2,1,3)],3)

# SMC-JM column:
Table.of.results[3:10,11:13]<-round(fit.MI.SMC$estimates[,c(1,2,5)],3)
Table.of.results[12:14,11]<-round(fit.MI.SMC$var.comp[c(2,1,3)],3)

################################################
#### Load results random coefficient model  ####

load("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Missing data/SMC-jomo/Tutorial+Simulations paper/R script + Data/Results_Example_rcm.RData")

# Complete records column:
Table.of.results[17:25,2:4]<-round(coef(summary(fit.CR))[,c(1,2,5)],3)
Table.of.results[27:31,2]<-round(c( attr(VarCorr(fit.CR)$hosp, "stddev")^2, 
                                    attr(VarCorr(fit.CR)$"hosp:clinician", "stddev")^2,
                                    attr(VarCorr(fit.CR)$"hosp:clinician", "correlation")[2,1]*
                                      attr(VarCorr(fit.CR)$"hosp:clinician", "stddev")[1]*
                                      attr(VarCorr(fit.CR)$"hosp:clinician", "stddev")[2],
                                    attr(VarCorr(fit.CR), "sc")^2),3)

# JM-Hom column:
Table.of.results[17:25,5:7]<-round(fit.MI$estimates[,c(1,2,5)],3)
Table.of.results[27:31,5]<-round(fit.MI$var.comp[c(4,1,3,2,5)],3)

# JM-Het column:
Table.of.results[17:25,8:10]<-round(fit.MI.het$estimates[,c(1,2,5)],3)
Table.of.results[27:31,8]<-round(fit.MI.het$var.comp[c(4,1,3,2,5)],3)

# SMC-JM column:
Table.of.results[17:25,11:13]<-round(fit.MI.SMC$estimates[,c(1,2,5)],3)
Table.of.results[27:31,11]<-round(fit.MI.SMC$var.comp[c(4,1,3,2,5)],3)

View(Table.of.results)
row.names(Table.of.results)<-NULL

write.table(Table.of.results, "Table3.txt", sep=",")