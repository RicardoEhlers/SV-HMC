#======================================================================#
#               Códigos dissertação Valor em Risco                     #
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

#---- Carregando funções externas:
source('/home/ddias/Dropbox/r_code/SV_model_rstan.R') # códigos dos modelos SV
source('/home/ddias/Dropbox/r_code/plot_mcmc.R')
source('/home/ddias/Dropbox/r_code/plot_vol.R')

#-----------------------------------------------------------------------
# Cálculo do VaR para a distribuição Normal
#-----------------------------------------------------------------------
#----- Dados:
#http://www.professores.uff.br/joel/lib/exe/fetch.php?media=disciplinas:retornos.pdf
dataxy <- read.csv("/home/ddias/Dropbox/Documentos/Mestrado ICMC-USP/Séries Temporais/Series_financeiras/conjunto de dados/ibovespa.csv", header = TRUE, sep = ",")
dados = dataxy[c(3254:3294),]

dataxy = dataxy[-c(1:1236),]

hh = 2016:2057
hh = paste('h[',hh,']',sep='')

VaR = matrix(NA,ncol=3,nrow=length(hh))
for(i in 1:length(hh)){ #nrow(dados)
  
  dados = dataxy[1:((nrow(dataxy)-length(hh))+i), ]
  Pr <- dados[,5]
  
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
  phi.0 =  0.98               
  sigma.0 = 0.10 
  data = list(T = length(ret), y = ret)
  
  inits <-list(list(mu = -9, phiT=(phi.0 + 1)/2, s2=sigma.0^2))
  fit.n <- stan(model_code=stan_code_n, data=data, init=inits, chains=1, iter=nsamples)
  
  #----------- Extração de valores necessários:
  v.n = extract(fit.n, pars = c('mu','phi','sigma',hh[i]))   
  theta.n = cbind(v.n$mu,v.n$phi,v.n$sigma)
  colnames(theta.n) = c('mu','phi','sigma')
  
  hT = v.n[[4]]
  
  #-----------------------------------------------------------------------
  # Valor em Risco um passo a frente
  #-----------------------------------------------------------------------
  
  M=1000
  h_T = NULL; r_T = matrix(NA, ncol=M, nrow=nrow(theta.n))
  for(l in 1:M){
    for(s in 1:nrow(theta.n))
    {
      h_T[s] = rnorm(1, theta.n[s,1] + theta.n[s,2]*(hT[s] - theta.n[s,1]), theta.n[2,3]) 
      r_T[s,l] = rnorm(1,0,exp(h_T[s]/2))
    }
  }
  VaR[i,] = quantile(r_T,c(0.01, 0.05, 0.1))
  cat("Iteration = ", i, "\n") 
}


VaRN = 100*VaR

par(mar = c(4, 4, 1, 1))
plot(100*ret[2017:2057],ylim=c(-5,5),xlab='Tempo',ylab='Retornos %',pch=19)
lines(VaRN[,1],col=2)
lines(VaRN[,2],col=3)
lines(VaRN[,3],col=4)
legend('topleft',col=2:4,title=expression(VaR[Gaussiano]),legend=c('(1%)','(5%)','(10%)'),pch=15,bty='n')

#save.image("/home/ddias/Documentos/VaRIBOVESPA_N.RData")


#=================================================================================================================================
#-----------------------------------------------------------------------
# Cálculo do VaR para a distribuição Skew-Normal
#-----------------------------------------------------------------------
#----- Dados:
#http://www.professores.uff.br/joel/lib/exe/fetch.php?media=disciplinas:retornos.pdf
dataxy <- read.csv("/home/ddias/Dropbox/Documentos/Mestrado ICMC-USP/Séries Temporais/Series_financeiras/conjunto de dados/ibovespa.csv", header = TRUE, sep = ",")
dados = dataxy[c(3254:3294),]

dataxy = dataxy[-c(1:1236),]

hh = 2016:2057
hh = paste('h[',hh,']',sep='')


VaR = matrix(NA,ncol=3,nrow=length(hh))
for(i in 1:length(hh)){ #nrow(dados)
  
  dados = dataxy[1:((nrow(dataxy)-length(hh))+i), ]
  Pr <- dados[,5]
  
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
  phi.0 =  0.99               
  sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
  alpha.0 = 0.01

  data = list(T =length(ret), y= ret)

  inits <-list(list(mu = -10, phiT=(phi.0 + 1)/2, s2=sigma.0^2, alpha=alpha.0))
  fit.sn <- stan(model_code=stan_code_sn, data=data, init=inits, chains=1, iter=nsamples)
  
  #----------- Extração de valores necessários:
  v.sn = extract(fit.sn, pars = c('mu','phi','sigma','alpha',hh[i]))
  theta.sn = cbind(v.sn$mu,v.sn$phi,v.sn$sigma,v.sn$alpha)
  colnames(theta.sn) = c('mu','phi','sigma','alpha')

  hT = v.sn[[5]]
  
  #-----------------------------------------------------------------------
  # Valor em Risco um passo a frente
  #-----------------------------------------------------------------------
  
  M=1000
  h_T = NULL; r_T = matrix(NA, ncol=M, nrow=nrow(theta.sn))
  for(l in 1:M){
    for(s in 1:nrow(theta.sn))
    {
      h_T[s] = rnorm(1, theta.sn[s,1] + theta.sn[s,2]*(hT[s] - theta.sn[s,1]), theta.sn[s,3]) 
      r_T[s,l] = rsn(1,xi=0, omega=exp(h_T[s]/2), alpha = theta.sn[s,4]) 
    }
  }
  VaR[i,] = quantile(r_T,c(0.01, 0.05, 0.1))
  cat("Iteration = ", i, "\n") 
}

VaRSN = 100*VaR

par(mar = c(4, 4, 1, 1))
plot(100*ret[2017:2057],ylim=c(-5,5),xlab='Tempo',ylab='Retornos %',pch=19)
lines(VaRSN[,1],col=2)
lines(VaRSN[,2],col=3)
lines(VaRSN[,3],col=4)
legend('topleft',col=2:4,title=expression(VaR[Skew-Normal]),legend=c('(1%)','(5%)','(10%)'),pch=15,bty='n')




par(mar = c(4, 4, 1, 1))
plot(100*ret[2017:2057],ylim=c(-5,5),xlab='Tempo',ylab='Retornos %',pch=19)
lines(VaRN[,2],col=2)
lines(VaRSN[,2],col=4)
legend('topleft',col=c(2,4),title='VaR(5%)',legend=c('Gaussiano','Skew-Normal'),pch=15,bty='n')

