##############################################################
####  Figure 1. Random slopes sensitivity                 ####
##############################################################

# Set working directory:
setwd("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Missing data/SMC-jomo/Tutorial+Simulations paper/R script + Data")

#Initialise data.frame of results:
res<-data.frame(matrix(NA,27,5))
meths<-c("MI.hom","MI.het","MI.SMC")
res[,1]<-c(1:3)
res[,2]<-rep(c(1,2,5),each=3)

#####################
####Load data1:  ####
#####################

load("Results_Quadratic.RData")

#True values
True1<-beta.smc[2:4]

# Values estimated for various methods:

for ( i in 1:3) {
  meanz<-apply(get(paste("est",meths[i],sep=".")),2,mean)[2:4]
  res[c(1,2,3),i+2]<-((meanz-True1)/True1)*100
} 


#####################
####Load data2:  ####
#####################

load("Results_Quadratic_sens2.RData")

#True values
True1<-beta.smc[2:4]

# Values estimated for various methods:

for ( i in 1:3) {
  meanz<-apply(get(paste("est",meths[i],sep=".")),2,mean)[2:4]
  res[c(4,5,6),i+2]<-((meanz-True1)/True1)*100
} 

#####################
####Load data3:  ####
#####################

load("Results_Quadratic_sens5.RData")

#True values
True1<-beta.smc[2:4]

# Values estimated for various methods:

for ( i in 1:3) {
  meanz<-apply(get(paste("est",meths[i],sep=".")),2,mean)[2:4]
  res[c(7,8,9),i+2]<-((meanz-True1)/True1)*100
} 

#####################
####Load data4:  ####
#####################

load("Results_Cubic.RData")

#True values
True1<-beta.smc[2:4]

# Values estimated for various methods:

for ( i in 1:3) {
  meanz<-apply(get(paste("est",meths[i],sep=".")),2,mean)[2:4]
  res[c(10:12),i+2]<-((meanz-True1)/True1)*100
} 

#####################
####Load data5:  ####
#####################

load("Results_Cubic_sens2.RData")

#True values
True1<-beta.smc[2:4]

# Values estimated for various methods:

for ( i in 1:3) {
  meanz<-apply(get(paste("est",meths[i],sep=".")),2,mean)[2:4]
  res[c(13:15),i+2]<-((meanz-True1)/True1)*100
} 

#####################
####Load data6:  ####
#####################

load("Results_Cubic_sens5.RData")

#True values
True1<-beta.smc[2:4]

# Values estimated for various methods:

for ( i in 1:3) {
  meanz<-apply(get(paste("est",meths[i],sep=".")),2,mean)[2:4]
  res[c(16:18),i+2]<-((meanz-True1)/True1)*100
} 

#####################
####Load data7:  ####
#####################

load("Results_Interaction.RData")

#True values
True1<-beta.smc[2:4]

# Values estimated for various methods:

for ( i in 1:3) {
  meanz<-apply(get(paste("est",meths[i],sep=".")),2,mean)[2:4]
  res[c(19:21),i+2]<-((meanz-True1)/True1)*100
} 

#####################
####Load dat8:  ####
#####################

load("Results_Interaction_sens2.RData")

#True values
True1<-beta.smc[2:4]

# Values estimated for various methods:

for ( i in 1:3) {
  meanz<-apply(get(paste("est",meths[i],sep=".")),2,mean)[2:4]
  res[c(22:24),i+2]<-((meanz-True1)/True1)*100
} 

#####################
####Load data9:  ####
#####################

load("Results_Interaction_sens5.RData")

#True values
True1<-beta.smc[2:4]

# Values estimated for various methods:

for ( i in 1:3) {
  meanz<-apply(get(paste("est",meths[i],sep=".")),2,mean)[2:4]
  res[c(25:27),i+2]<-((meanz-True1)/True1)*100
} 



pdf("PolyInter.pdf", ,width=7,height=7)

plot(rep(1,nrow(res)),res[,3], xaxt="n", xlim=c(0,4), ylim=c(-100,100), pch=19,
     xlab="Imputation Method", ylab="Relative bias (%)", 
     main="Relative bias for fixed-effect parameter estimates")
points(rep(2, nrow(res)),res[,4], pch=17, col="blue")
points(rep(3, nrow(res)),res[,5], pch=15, col="darkgreen")
abline(h=0, col="red")
axis(1,at=c(1,2,3), labels = c("MI-Hom", "MI-Het", "MI-SMC"))
dev.off()
