##############################################################
####  Figure 1. Random slopes sensitivity                 ####
##############################################################

# Set working directory:
setwd("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Missing data/SMC-jomo/Tutorial+Simulations paper/R script + Data")

#Initialise data.frame of results:
res<-data.frame(matrix(NA,15,10))
meths<-c("FD","CR","MI.hom","MI.het","MI.SMC")
res[,1]<-meths
res[,2]<-rep(c(1,2,5),each=5)

#####################
####Load data1:  ####
#####################

load("Results_Base_case.RData")

#True values
True1<-beta.smc
True1v<-reff.smc[2,2]

# Values estimated for various methods:

for ( i in 1:5) {
  
  meanz<-apply(get(paste("est",meths[i],sep=".")),2,mean)
  res[i,3]<-((meanz-True1)/True1)[3]*100
  res[i,5]<-sqrt(sum((get(paste("est",meths[i],sep="."))[3]-meanz[3])^2)/(1000*999))
  meanu<-mean(get(paste("rev",meths[i],sep="."))[,2,2])
  res[i,4]<-((meanu-True1v)/True1v)*100
  res[i,6]<-sqrt(sum((get(paste("rev",meths[i],sep="."))[,2,2]-meanu)^2)/(1000*999))
  res[i,7]<-((meanz[3]-res[i,5]*qnorm(0.975)-True1[3])/True1[3])*100
  res[i,8]<-((meanz[3]+res[i,5]*qnorm(0.975)-True1[3])/True1[3])*100
  res[i,9]<-((meanu-res[i,6]*qnorm(0.975)-True1v)/True1v)*100
  res[i,10]<-((meanu+res[i,6]*qnorm(0.975)-True1v)/True1v)*100
  
} 


#####################
####Load data2:  ####
#####################

load("Results_Base_case_sensitivity_2.RData")

#True values
True1<-beta.smc
True1v<-reff.smc[2,2]

# Values estimated for various methods:

for ( i in 1:5) {
  
  meanz<-apply(get(paste("est",meths[i],sep=".")),2,mean)
  res[i+5,3]<-((meanz-True1)/True1)[3]*100
  res[i+5,5]<-sqrt(sum((get(paste("est",meths[i],sep="."))[3]-meanz[3])^2)/(1000*999))
  meanu<-mean(get(paste("rev",meths[i],sep="."))[,2,2])
  res[i+5,4]<-((meanu-True1v)/True1v)*100
  res[i+5,6]<-sqrt(sum((get(paste("rev",meths[i],sep="."))[,2,2]-meanu)^2)/(1000*999))
  res[i+5,7]<-((meanz[3]-res[i+5,5]*qnorm(0.975)-True1[3])/True1[3])*100
  res[i+5,8]<-((meanz[3]+res[i+5,5]*qnorm(0.975)-True1[3])/True1[3])*100
  res[i+5,9]<-((meanu-res[i+5,6]*qnorm(0.975)-True1v)/True1v)*100
  res[i+5,10]<-((meanu+res[i+5,6]*qnorm(0.975)-True1v)/True1v)*100
  
} 

#####################
####Load data2:  ####
#####################

load("Results_Base_case_sensitivity_5.RData")

#True values
True1<-beta.smc
True1v<-reff.smc[2,2]

# Values estimated for various methods:

for ( i in 1:5) {
  
  meanz<-apply(get(paste("est",meths[i],sep=".")),2,mean)
  res[i+10,3]<-((meanz-True1)/True1)[3]*100
  res[i+10,5]<-sqrt(sum((get(paste("est",meths[i],sep="."))[3]-meanz[3])^2)/(1000*999))
  meanu<-mean(get(paste("rev",meths[i],sep="."))[,2,2])
  res[i+10,4]<-((meanu-True1v)/True1v)*100
  res[i+10,6]<-sqrt(sum((get(paste("rev",meths[i],sep="."))[,2,2]-meanu)^2)/(1000*999))
  res[i+10,7]<-((meanz[3]-res[i+10,5]*qnorm(0.975)-True1[3])/True1[3])*100
  res[i+10,8]<-((meanz[3]+res[i+10,5]*qnorm(0.975)-True1[3])/True1[3])*100
  res[i+10,9]<-((meanu-res[i+10,6]*qnorm(0.975)-True1v)/True1v)*100
  res[i+10,10]<-((meanu+res[i+10,6]*qnorm(0.975)-True1v)/True1v)*100
  
} 
View(res)

colnames(res)<-c("method","magnitude","estimate","var","se_est","se_var","low_est","up_est", "low_var", "up_var")


pdf("Random slopes plot.pdf", ,width=14,height=7)
m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.85,0.15))

plot(res[res$method=="MI.hom","magnitude"], res[res$method=="MI.hom","estimate"],
     main="PAR fixed-effect estimate", xlab="Magnitude of effect", ylab="Relative bias (%)",
     xaxt="n", yaxt="n", ylim=c(-30,30), type="b")
polygon(x=c(1,2,5,5,2,1), y= c(res[res$method=="MI.hom","low_est"],
        rev(res[res$method=="MI.hom","up_est"])), col="grey79", border="NA")
mycol <- rgb(0, 0, 255, max = 255, alpha = 25, names = "blue50")
polygon(x=c(1,2,5,5,2,1), y= c(res[res$method=="MI.het","low_est"],
                               rev(res[res$method=="MI.het","up_est"])), col=mycol, border="NA")
mycol2 <- rgb(0, 255, 0, max = 255, alpha = 25, names = "green50")
polygon(x=c(1,2,5,5,2,1), y= c(res[res$method=="MI.SMC","low_est"],
                               rev(res[res$method=="MI.SMC","up_est"])), col=mycol2, border="NA")
abline(h=0, col="red")
lines(res[res$method=="MI.hom","magnitude"], res[res$method=="MI.hom","estimate"],
      type="b", pch=19)
lines(res[res$method=="MI.het","magnitude"], res[res$method=="MI.het","estimate"],
      type="b", pch=17, lty=2, col="blue")
lines(res[res$method=="MI.SMC","magnitude"], res[res$method=="MI.SMC","estimate"],
      type="b", pch=15, lty=3, col="darkgreen")
axis(1, at=c(1,2,5), labels=c("1X","2X","5X"))
axis(2, at=seq(-30,30,10), labels = paste(seq(-30,30,10),"%", sep=""), las=2)

plot(res[res$method=="MI.hom","magnitude"], res[res$method=="MI.hom","var"],
     main="Random slope variance", xlab="Magnitude of variance", ylab="Relative bias (%)",
     xaxt="n", yaxt="n", ylim=c(-10,10), type="b")
polygon(x=c(1,2,5,5,2,1), y= c(res[res$method=="MI.hom","low_var"],
                               rev(res[res$method=="MI.hom","up_var"])), col="grey79", border="NA")
mycol <- rgb(0, 0, 255, max = 255, alpha = 25, names = "blue50")
polygon(x=c(1,2,5,5,2,1), y= c(res[res$method=="MI.het","low_var"],
                               rev(res[res$method=="MI.het","up_var"])), col=mycol, border="NA")
mycol2 <- rgb(0, 255, 0, max = 255, alpha = 25, names = "green50")
polygon(x=c(1,2,5,5,2,1), y= c(res[res$method=="MI.SMC","low_var"],
                               rev(res[res$method=="MI.SMC","up_var"])), col=mycol2, border="NA")
abline(h=0, col="red")
lines(res[res$method=="MI.hom","magnitude"], res[res$method=="MI.hom","var"],
      type="b", pch=19)
lines(res[res$method=="MI.het","magnitude"], res[res$method=="MI.het","var"],
      type="b", pch=17, lty=2, col="blue")
lines(res[res$method=="MI.SMC","magnitude"], res[res$method=="MI.SMC","var"],
      type="b", pch=15, lty=3, col="darkgreen")
axis(1, at=c(1,2,5), labels=c("1X","2X","5X"))
axis(2, at=seq(-10,10,5), labels = paste(seq(-10,10,5),"%", sep=""), las=2)

# Now draw the legend:
par(mar = c(0.4,0.4,0.4,0.4))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("MI-Hom", "MI-Het", "MI-SMC"), 
       col=c("black","blue", "darkgreen"),  pch=c(19,17,15),lty=c(1,2,3), cex=1, ncol=3)
par(mar = c(5.1, 4.1, 4.1, 2.1))



dev.off()
