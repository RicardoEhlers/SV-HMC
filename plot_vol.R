
plot.vol <- function(y,x.est,col0=c(1,4,4),lwd0=1.5) {
#  Plot absolute returns y and estimated volatilities exp(x.est/2).
   x.est =as.matrix(x.est)
   ts.plot(abs(y),type='h',ylab="|retornos|",col='grey',xlab='Tempo')
   #lines(exp(x.est[,1]/2),col=col0[1])
   for(i in 1:ncol(x.est))lines(exp(x.est[,i]/2),col=col0[i],lwd=lwd0)
}


plot.vol2 <- function(y,x.est,rgb0 = c(30,144,255),alpha0=85,ylab0="|Returns|",xlab0='Time') {
#  Plot absolute returns y and estimated volatilities exp(x.est/2).
   Teste = data.frame(exp(x.est/2))
   means = Teste[,1]; mins = Teste[,2]; maxs = Teste[,3] 
   xcoords = 1:nrow(Teste)
   par(mfrow=c(1, 1), mar = c(3.2, 3.2, 1., .5), mgp = c(2, .6, 0))
   ts.plot(abs(y),type='h',ylab=ylab0,col='grey',xlab=xlab0)
   rangecolor <- rgb(rgb0[1],rgb0[2],rgb0[3],alpha=alpha0,maxColorValue=255)
   polygon(x=c(xcoords,rev(xcoords)),y=c(maxs,rev(means)),col=rangecolor,border=NA)
   polygon(x=c(xcoords,rev(xcoords)),y=c(mins,rev(means)),col=rangecolor,border=NA)
   lines(x=xcoords,y=means,col=1)
}

# for rgb values of colors see http://en.wikipedia.org/wiki/Web_colors)


plot.vol3 <- function(y,x.est,ylab0="|Returns|",xlab0='Time') {
  #  Plot absolute returns y and estimated volatilities exp(x.est/2).
  Teste=data.frame(100*exp(x.est/2))
  means = Teste[,1]; mins <- Teste[,2]; maxs <- Teste[,3] 
  xcoords = 1:nrow(Teste)
  par(mfrow=c(1, 1), mar = c(3.2, 3.2, 1., .5), mgp = c(2, .6, 0))
  ts.plot(abs(100*y),type='h',ylab=ylab0,col='grey',xlab=xlab0)
  rangecolor <- rgb(30,144,255,alpha=85,maxColorValue=255)
  polygon(x=c(xcoords,rev(xcoords)),y=c(maxs,rev(means)),col=rangecolor,border=NA)
  polygon(x=c(xcoords,rev(xcoords)),y=c(mins,rev(means)),col=rangecolor,border=NA)
  lines(x=xcoords,y=means,col=1)
}

