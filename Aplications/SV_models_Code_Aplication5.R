#======================================================================#
#                Códigos Dissertação - Modelos SV                      #
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

#---- Carregando funções externas:
source('/home/ddias/Dropbox/r_code/SV_model_rstan.R') # códigos dos modelos SV
source('/home/ddias/Dropbox/r_code/plot_mcmc.R')
source('/home/ddias/Dropbox/r_code/plot_vol.R')

#-----------------------------------------------------------------------
# Análise descritiva do conjunto de dados Bitcoin/USD
#-----------------------------------------------------------------------
#----- Dados:
#https://br.investing.com/currencies/btc-usd-chart
btc <- read.csv("/home/ddias/Dropbox/Documentos/Mestrado ICMC-USP/Séries Temporais/Series_financeiras/conjunto de dados/BTC_USD.csv", header = TRUE, sep = ",")

btc$date = as.Date(btc$Data, "%d.%m.%Y")
btc$rdate<-format(btc$date, format = "%Y")

btc= subset(btc, rdate > 2016)

Pr <- btc[,2]

titulo = "Série de Preços 1 BTC em USD"
par(mar = c(4, 4, 1, 1))
plot(Pr, type = "l", xaxt = "n", xlab = "Dias", ylab = "Preços 1 BTC em USD")
mth <- unique(btc$rdate)
axis(1, at = match(mth, btc$rdate), labels = mth, cex.axis = 1)

r  <- diff(log(Pr))

#--- Série Depreciada:
ret <- logret(Pr, demean=TRUE)

basicStats(ret)

par(mar = c(4, 4, 1, 1))
plot(100*ret, type = "l", xaxt = "n", xlab = "Dias", ylab = "Log-retornos da série Euro/USD")
# add x-axis
qtr <- mth[seq(1,length(mth),1)] #mth[seq(1,length(mth),2)]
axis(1, at = match(mth, btc$rdate), labels = mth, cex.axis = 1)

hist(100*ret,breaks=50,main='',ylab='Frequência',xlab='Log-retornos da série Euro/Dólar')

plot(acf(ret,lag=100),main='',ylim=c(-0.2,0.2),xlab='Defasagem',ylab='ACF Log-retornos')

plot(acf(ret^2,lag=100),main='',ylim=c(-0.20,0.30),xlab='Defasagem',ylab='ACF quadrado Log-retornos')

#=========================================== Modelo 1 =============================================
nsamples = 20000
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo normal
#-----------------------------------------------------------------------
#--- Initial Values:
phi.0 =  0.99               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
data = list(T = length(ret), y = ret)

inits <-list(list(mu = -10, phiT=(phi.0 + 1)/2, s2=sigma.0^2))
fit.n <- stan(model_code=stan_code_n, data=data, init=inits, chains=1, iter=nsamples)

summary(fit.n)$summary[c(1,(length(ret)+4),(length(ret)+5)), ] 

#----------- Uso do Coda:
v.n = extract(fit.n, pars = c('mu','phi','sigma'))
theta.n = cbind(v.n$mu,v.n$phi,v.n$sigma)
colnames(theta.n) = c('mu','phi','sigma')

m.n = window(as.mcmc(theta.n),start=1,thin=1)
summary(m.n)
effectiveSize(m.n)

x.est.n = summary(fit.n)$summary[c(4:(length(ret)+3)),c(1,4,8)]
ht.n = as.vector(x.est.n[,1])

#--- Cálculo do DIC, LOO e WAIC:
logl.n = extract_log_lik(fit.n, parameter_name = "loglik1")   #logp(yi|tt)
lik.n = extract_log_lik(fit.n, parameter_name = "loglik")     #logp(y|tt(S))

Loo.n = loo(logl.n)  
Loo.n
Waic.n = waic(logl.n) 
Waic.n

lik_hat.n =  sum(dnorm(ret,0,exp(ht.n/2),log=T)) # log(p(y|theta_hat))
#lik_hat.n = sum(log(dn(ret, exp(ht.n/2))))
DIC.n = dic(lik_hat.n,lik.n)
DIC.n

#---- Gráficos de Convergência: 
#mcmcplot(m)
pars_names.n = c(expression(mu), expression(phi), expression(sigma))
#plot.mcmc(theta.n,pars_names, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))
plot.mcmc2(theta.n,pars_names.n, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))







