#-------------------------------------------------------------------------
# Plot for object mcmc (trace, density, acf)
#-------------------------------------------------------------------------
library(TSA) #acf

plot.mcmc <- function(A, ylab0=colnames(A), cex.lab0=1.2, mar0=c(2,4,1.9,1.9))
{
par(mfrow=c(ncol(A),3), mar=mar0,cex.lab=cex.lab0)
 for(i in 1:ncol(A)){
 ts.plot(A[,i], ylab=ylab0[i],main="")
 plot(density(A[,i]), ylab=ylab0[i], main="")
 plot(acf(A[,i], lag=100,plot=F),main='',xlab='Defasagem')#autocorr.plot(out[,i],auto.layout=F) 
 }
}


#-------------------------------------------------------------------------
plot.mcmc2 <- function(A, ylab0=colnames(A),breaks0=20, cex.lab0=1.2, mar0=c(2,4,1.9,1.9))
{
par(mfrow=c(ncol(A),3), mar=mar0, cex.lab=cex.lab0)
 for(i in 1:ncol(A)){
 ts.plot(A[,i], ylab=ylab0[i],main="")
 hist(A[,i],ylab=ylab0[i],main='',breaks=breaks0); abline(v=mean(A[,i]),lwd=1.5,lty=2,col=4)
 plot(acf(A[,i], lag=100,plot=F),main='',xlab='Defasagem') #acf(A[,i], lag=100, ylab=ylab0[i],main="")  
 }
}


