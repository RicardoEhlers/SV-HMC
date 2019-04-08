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
# Análise descritiva do conjunto de dados Ibovespa
#-----------------------------------------------------------------------
#----- Dados:
#http://www.professores.uff.br/joel/lib/exe/fetch.php?media=disciplinas:retornos.pdf
dataxy <- read.csv("/home/ddias/Dropbox/Documentos/Mestrado ICMC-USP/Séries Temporais/Series_financeiras/conjunto de dados/ibovespa.csv", header = TRUE, sep = ",")
dataxy = dataxy[-c(1:1236,3254:3294),]

dataxy$date = as.Date(dataxy$Date, "%m/%d/%Y")

Pr <- dataxy[,5]

titulo = "Série de Preços de Fechamento da IBOVESPA"
par(mar = c(4, 4, 1, 1))
plot(Pr, type = "l", xaxt = "n", xlab = "Dias", ylab = "Preços de Fechamento do IBOVESPA")
dataxy$rdate<-format(dataxy$date, format = "%Y")
mth <- unique(dataxy$rdate)
axis(1, at = match(mth, dataxy$rdate), labels = mth, cex.axis = 1)

ts.plot(Pr, main='',xlab='Dias',ylab=expression(P[t]))
acf(Pr,lag=100,main='',xlab='Defasagem')

nf <- layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,3))
par(mar=c(5.1, 4.1, 1.1, 2.1))
boxplot(Pr, horizontal=TRUE,  outline=FALSE)
hist(Pr)

#------------------------
# Retornos aritméticos:
#------------------------
T <- length(Pr)
Pt <- Pr[2:T]
Pt_1 <- Pr[1:(T - 1)]
Rt <- (Pt - Pt_1)/Pt_1

# identificando e tratando o outliers
out <- which(Rt > mean(Rt) + 6*sd(Rt) | Rt < mean(Rt) - 6*sd(Rt))
#Rt[out] <- mean(Rt)

layout(1:2)
plot(ts(Rt), main = "Serie de Retornos")
acf(ts(Rt), main = "Autocorrelacao dos Retornos")

#------------------------
# Log - Retornos
#------------------------
rt <- log(Pt/Pt_1)

basicStats(rt)

par(mar = c(4, 4, 1, 1))
plot(rt, type = "l", xaxt = "n", xlab = "Dias", ylab = "Log-retornos da série IBOVESPA")
axis(1, at = match(mth, dataxy$rdate), labels = mth, cex.axis = 1)

hist(100*rt,freq=F,breaks=50,ylab='Densidade',xlab=expression(r[t]),main='',xlim=c(-15,15))

#Autocorrelacao dos Log-Retornos
plot(acf(ts(rt),lag=100,plot=F),main='',ylim=c(-0.07,0.07),xlab='Defasagem',ylab='ACF Log-retornos')

plot(acf(rt^2,lag=100,plot=F),main='',ylim=c(-0.05,0.45),xlab='Defasagem',ylab='ACF quadrado Log-retornos')

#--- Série de log-retornos depreciada:
ret = rt - mean(rt)
basicStats(ret)

#=========================================== Modelo 1 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo normal
#-----------------------------------------------------------------------
nsamples = 10000

#--- Initial Values:
x.0 = est.x(ret)/2
phi.0 =  0.98               
sigma.0 = 0.10 #var(x.0)*(1-phi.0^2)
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


#--- Gráfico de ajuste das volatilidades
#plot.vol(ret,x.est.n,col0=c(1,4,4))
plot.vol3(ret,x.est.n)

write.table(x.est.n,'/home/ddias/Dropbox/Images_Paper/x.est_ibovespa_n.txt')

#plot.ts(cbind(abs(ret),vol[,1],x.est.n[,1])) # comparar com o garch

#=========================================== Modelo 2 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo t-student
#-----------------------------------------------------------------------
#--- Initial Values
x.0 = est.x(ret)/2
phi.0 =  0.98               
sigma.0 = 0.10 #var(x.0)*(1-phi.0^2)
nu.0 = 6

data = list(T =length(ret), y= ret)

#h=x.0
inits <-list(list(mu = -10, phiT=(phi.0 + 1)/2, s2=sigma.0^2, nu=nu.0))
fit.t <- stan(model_code=stan_code_t, data=data, init=inits, chains=1, iter=nsamples)

summary(fit.t)$summary[c(1,4,(length(ret)+5),(length(ret)+6)), ] 

#----------- Uso do Coda:
v.t = extract(fit.t, pars = c('mu','phi','sigma','nu'))
theta.t = cbind(v.t$mu,v.t$phi,v.t$sigma,v.t$nu)
colnames(theta.t) = c('mu','phi','sigma','nu')

m.t = window(as.mcmc(theta.t),start=1,thin=1)
summary(m.t)
effectiveSize(m.t)

x.est.t = summary(fit.t)$summary[c(5:(length(ret)+4)),c(1,4,8)]
ht.t = as.vector(x.est.t[,1])
nu_hat.t = summary(m.t)$statistics[4]

#--- Cálculo do DIC, LOO e WAIC:
logl.t = extract_log_lik(fit.t, parameter_name = "loglik1")   #logp(yi|tt)
lik.t = extract_log_lik(fit.t, parameter_name = "loglik")     #logp(y|tt(S))

Loo.t = loo(logl.t)
Loo.t
Waic.t = waic(logl.t) 
Waic.t

lik_hat.t =  sum(log(dt.st(ret, nu_hat.t, exp(ht.t/2))))
DIC.t = dic(lik_hat.t,lik.t)
DIC.t

#---- Gráficos de Convergência: 
pars_names= c(expression(mu), expression(phi), expression(sigma), expression(nu))
plot.mcmc2(theta.t,pars_names, cex.lab0=1.4, mar0 = c(1.8,4.5,1.9,1.9))

x = extract(fit.t, pars=c('h[1]','h[15]','h[50]','h[100]','h[245]','h[375]','h[450]','h[900]','h[1030]','h[1570]'))
h.t = cbind(x$'h[1]',x$'h[15]',x$'h[50]',x$'h[100]',x$'h[245]',x$'h[375]',x$'h[450]',x$'h[900]',x$'h[1030]',x$'h[1570]')
pars_ht = c('h[1]','h[15]','h[50]','h[100]','h[245]','h[375]','h[450]','h[900]','h[1030]','h[1570]')
plot.mcmc2(h.t, pars_ht, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))

#--- Gráfico
plot.vol2(ret,x.est.t)
plot.vol3(ret,x.est.t)

#=========================================== Modelo 3 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo skew-normal
#-----------------------------------------------------------------------
#--- Initial Values
x.0 = est.x(ret)/2
phi.0 =  0.99               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
alpha.0 = 0.01

data = list(T =length(ret), y= ret)

#h=x.0
inits <-list(list(mu = -10, phiT=(phi.0 + 1)/2, s2=sigma.0^2, alpha=alpha.0))
fit.sn <- stan(model_code=stan_code_sn, data=data, init=inits, chains=1, iter=nsamples)

summary(fit.sn)$summary[c(1,4,(length(ret)+5),(length(ret)+6)), ]

#----------- Uso do Coda:
v.sn = extract(fit.sn, pars = c('mu','phi','sigma','alpha'))
theta.sn = cbind(v.sn$mu,v.sn$phi,v.sn$sigma,v.sn$alpha)
colnames(theta.sn) = c('mu','phi','sigma','alpha')

m.sn = window(as.mcmc(theta.sn),start=1,thin=1)
summary(m.sn)
effectiveSize(m.sn)

x.est.sn = summary(fit.sn)$summary[c(5:(length(ret)+4)),c(1,4,8)]
ht.sn = as.vector(x.est.sn[,1])
alpha_hat.sn = summary(m.sn)$statistics[4]

#--- Cálculo do DIC, LOO e WAIC:
logl.sn = extract_log_lik(fit.sn, parameter_name = "loglik1")   #logp(yi|tt)
lik.sn = extract_log_lik(fit.sn, parameter_name = "loglik")     #logp(y|tt(S))

Loo.sn = loo(logl.sn)  
Loo.sn
Waic.sn = waic(logl.sn) 
Waic.sn

lik_hat.sn = sum(log(dsn(ret, alpha_hat.sn, exp(ht.sn/2))))
DIC.sn = dic(lik_hat.sn,lik.sn)
DIC.sn

#---- Gráficos de Convergência 
pars_names= c(expression(mu), expression(phi), expression(sigma), expression(nu))
plot.mcmc2(theta.sn,pars_names, cex.lab0=1.4, mar0 = c(1.8,4.5,1.9,1.9))

#--- Gráfico
plot.vol2(ret,x.est.sn)
plot.vol3(ret,x.est.sn)

write.table(x.est.sn,'/home/ddias/Dropbox/Images_Paper/x.est_ibovespa_sn.txt')


#=========================================== Modelo 4 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo GED
#-----------------------------------------------------------------------
#--- Initial Values
x.0 = est.x(ret)/2
phi.0 =  0.985               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
nu.0 = 1.5

data = list(T =length(ret), y= ret)

#h=x.0
inits <-list(list(mu = -10, phiT=(phi.0 + 1)/2, s2=sigma.0^2, nu=nu.0))
fit.ged <- stan(model_code=stan_code_ged, data=data, init=inits, chains=1, iter=nsamples)

summary(fit.ged)$summary[c(1,4,(length(ret)+5),(length(ret)+6)), ]

#----------- Uso do Coda:
v.ged = extract(fit.ged, pars = c('mu','phi','sigma','nu'))
theta.ged = cbind(v.ged$mu,v.ged$phi,v.ged$sigma,v.ged$nu)
colnames(theta.ged) = c('mu','phi','sigma','nu')

m.ged = window(as.mcmc(theta.ged),start=1,thin=1)
summary(m.ged)
effectiveSize(m.ged)

x.est.ged = summary(fit.ged)$summary[c(5:(length(ret)+4)),c(1,4,8)]
ht.ged = as.vector(x.est.ged[,1])
nu_hat.ged = summary(m.ged)$statistics[4]

#--- Cálculo do DIC, LOO e WAIC:
logl.ged = extract_log_lik(fit.ged, parameter_name = "loglik1")   #logp(yi|tt)
lik.ged = extract_log_lik(fit.ged, parameter_name = "loglik")     #logp(y|tt(S))

Loo.ged = loo(logl.ged)   
Loo.ged
Waic.ged = waic(logl.ged)  
Waic.ged

lik_hat.ged = sum(log(dged(ret, nu_hat.ged, exp(ht.ged/2))))
DIC.ged = dic(lik_hat.ged,lik.ged)
DIC.ged

#---- Gráficos de Convergência 
pars_names= c(expression(mu), expression(phi), expression(sigma), expression(nu))
plot.mcmc2(theta.ged,pars_names, cex.lab0=1.4, mar0 = c(1.8,4.5,1.9,1.9))


Rhat = summary(fit.ged)$summary[,10]
plot(Rhat,ylim=c(1,1.5))

#--- Gráfico
plot.vol2(ret,x.est.ged)
plot.vol3(ret,x.est.ged)

#save.image("/home/ddias/Documentos/SVIbovespa.RData")
