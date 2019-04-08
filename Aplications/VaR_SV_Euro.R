#======================================================================#
#             Códigos dissertação Valor em Risco EURO/USD              #
#======================================================================#
#----- Bibliotecas:
library(stochvol)
library(rstan)
library(coda)
library(fBasics)
library(forecast)
library(mcmcplots)
library(tsbugs) # data pound
library(loo)
library(TSA) # for acf
library(sn)
library(fGarch)

#---- Carregando funções externas:
source('/home/ddias/Dropbox/r_code/SV_model_rstan.R') # códigos dos modelos SV
source('/home/ddias/Dropbox/r_code/plot_mcmc.R')
source('/home/ddias/Dropbox/r_code/plot_vol.R')

#-----------------------------------------------------------------------
# Cálculo do VaR para a distribuição Normal
#-----------------------------------------------------------------------
#----- Dados:
data(exrates)
#euro = exrates[-c(1:509,2816:nrow(exrates)),] # tirando os anos de (2000,2001(509),2011,2012) 2002(765)
euro = exrates[-c(1:1019),]
y = euro$USD 
r  <- diff(log(euro$USD))

hh = 2076:2120
hh = paste('h[',hh,']',sep='')

VaRt = matrix(NA,ncol=3,nrow=length(hh))
for(i in 1:length(hh)){ #nrow(dados)
  
  dados = euro[1:((nrow(euro)-length(hh))+i), ]
  Pr <- dados[,23]
  
  #------------------------
  # Log - Retornos
  #------------------------
  T <- length(Pr)
  Pt <- Pr[2:T]
  Pt_1 <- Pr[1:(T - 1)]
  rt <- log(Pt/Pt_1)
  
  #--- Série de log-retornos depreciada:
  ret = rt - mean(rt)
  
  #-----------------------------------------------------------------------
  # Apliação ao conjunto de dados - rstan modelo t-Student
  #-----------------------------------------------------------------------
  nsamples = 10000
  
  #--- Valores iniciais:
  phi.0 =  0.99               
  sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
  nu.0 = 6
  
  data = list(T =length(ret), y= ret)
  
  inits <-list(list(mu = -10, phiT=(phi.0 + 1)/2, s2=sigma.0^2, nu=nu.0))
  fit.t <- stan(model_code=stan_code_t, data=data, init=inits, chains=1, iter=nsamples)
  
  #----------- Extração de valores necessários:
  v.t = extract(fit.t, pars = c('mu','phi','sigma','nu',hh[i]))
  theta.t = cbind(v.t$mu,v.t$phi,v.t$sigma,v.t$nu)
  colnames(theta.t) = c('mu','phi','sigma','nu')
  
  hT = v.t[[5]]
  
  #-----------------------------------------------------------------------
  # Valor em Risco um passo a frente
  #-----------------------------------------------------------------------
  
  M=1000
  h_T = NULL; r_T = matrix(NA, ncol=M, nrow=nrow(theta.t))
  for(l in 1:M){
    for(s in 1:nrow(theta.t))
    {
      h_T[s] = rnorm(1, theta.t[s,1] + theta.t[s,2]*(hT[s] - theta.t[s,1]), theta.t[s,3]) 
      r_T[s,l] = rt(1, theta.t[s,4])*exp(h_T[s]/2) 
    }
  }
  VaRt[i,] = quantile(r_T,c(0.01, 0.05, 0.1))
  cat("Iteration = ", i, "\n") 
}


VaRt

par(mar = c(4, 4, 1, 1))
plot(100*ret[2076:2120],ylim=c(-5,5),xlab='Tempo',ylab='Retornos %',pch=19)
lines(VaRt[,1],col=2)
lines(VaRt[,2],col=3)
lines(VaRt[,3],col=4)
legend('topleft',col=2:4,title=expression(VaR[t-Student]),legend=c('(1%)','(5%)','(10%)'),pch=15,bty='n')


#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo GED
#-----------------------------------------------------------------------

VaRged = matrix(NA,ncol=3,nrow=length(hh))
for(i in 1:length(hh)){ #nrow(dados)
  
  dados = euro[1:((nrow(euro)-length(hh))+i), ]
  Pr <- dados[,23]
  
  #------------------------
  # Log - Retornos
  #------------------------
  T <- length(Pr)
  Pt <- Pr[2:T]
  Pt_1 <- Pr[1:(T - 1)]
  rt <- log(Pt/Pt_1)
  
  #--- Série de log-retornos depreciada:
  ret = rt - mean(rt)
  
  #-----------------------------------------------------------------------
  # Apliação ao conjunto de dados - rstan modelo normal
  #-----------------------------------------------------------------------
  nsamples = 10000
  
  #--- Valores iniciais:
  phi.0 =  0.985               
  sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
  nu.0 = 1.5
  
  data = list(T =length(ret), y= ret)
  
  inits <-list(list(mu = -10, phiT=(phi.0 + 1)/2, s2=sigma.0^2, nu=nu.0))
  fit.ged <- stan(model_code=stan_code_ged, data=data, init=inits, chains=1, iter=nsamples)
  
  #----------- Extração de valores necessários:
  v.ged = extract(fit.ged, pars = c('mu','phi','sigma','nu',hh[i]))
  theta.ged = cbind(v.ged$mu,v.ged$phi,v.ged$sigma,v.ged$nu)
  colnames(theta.ged) = c('mu','phi','sigma','nu')
  
  hT = v.ged[[5]]
  
  #-----------------------------------------------------------------------
  # Valor em Risco um passo a frente
  #-----------------------------------------------------------------------
  
  M=1000
  h_T = NULL; r_T = matrix(NA, ncol=M, nrow=nrow(theta.ged))
  for(l in 1:M){
    for(s in 1:nrow(theta.ged))
    {
      h_T[s] = rnorm(1, theta.ged[s,1] + theta.ged[s,2]*(hT[s] - theta.ged[s,1]), theta.ged[s,3]) 
      r_T[s,l] = rged(1, nu=theta.ged[s,4])*exp(h_T[s]/2) 
    }
  }
  VaRged[i,] = quantile(r_T,c(0.01, 0.05, 0.1))
  cat("Iteration = ", i, "\n") 
}

VaRged

par(mar = c(4, 4, 1, 1))
plot(100*ret[2076:2120],ylim=c(-5,5),xlab='Tempo',ylab='Retornos %',pch=19)
lines(VaRged[,1],col=2)
lines(VaRged[,2],col=3)
lines(VaRged[,3],col=4)
legend('topleft',col=2:4,title=expression(VaR[GED]),legend=c('(1%)','(5%)','(10%)'),pch=15,bty='n')


#save.image("/home/ddias/Documentos/VaRIBOVESPA_N.RData")

