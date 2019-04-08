#======================================================================#
#                Códigos da Aplicação 1 Qualificação                   #
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
#library(sn)  # for skew-normal

#install.packages("packageName")
#---- Carregando funções externas:
source('/home/ddias/Dropbox/r_code/SV_model_rstan.R') # códigos dos modelos SV
source('/home/ddias/Dropbox/r_code/plot_mcmc.R')
source('/home/ddias/Dropbox/r_code/plot_vol.R')


#-----------------------------------------------------------------------
# Análise descritiva do conjunto de dados Pound/USD
#-----------------------------------------------------------------------
# Pound/Dollar data
dd = scan("/home/ddias/Dropbox/Projeto_David/MALA/sv.txt") # data(svpdx)
pound = dd - mean(dd)                                      # série depreciada

data(svpdx)
pound = svpdx$pdx - mean(svpdx$pdx)
svpdx$pound = pound
n = length(pound)

par(mar = c(4, 4, 1, 1))
plot(svpdx$pound, type = "l", xaxt = "n", xlab = "Dias", ylab = "Log-retornos da série Libras/Dólar")

# add x-axis
svpdx$rdate<-format(svpdx$date, format = "%Y")
mth <- unique(svpdx$rdate)
axis(1, at = match(mth[-1], svpdx$rdate), labels = mth[-1], cex.axis = 1)

# usando o formato com mês
#svpdx$rdate<-format(svpdx$date, format = "%b %Y")
#mth <- unique(svpdx$rdate)
#qtr <- mth[seq(1,length(mth),2)]
#axis(1, at = match(qtr, svpdx$rdate), labels = qtr, cex.axis = 1)
#axis(1, at = match(mth, svpdx$rdate), labels = FALSE, tcl = -0.2)

basicStats(pound)
hist(pound,breaks=50,xlim=c(-4,4),main='',ylab='Frequência',xlab='Log-retornos da série Libras/Dólar')

plot(acf(pound,lag=100),main='',ylim=c(-0.1,0.1),xlab='Defasagem',ylab='ACF Log-retornos')

plot(acf(pound^2,lag=100),main='',ylim=c(-0.06,0.30),xlab='Defasagem',ylab='ACF quadrado Log-retornos')

#=========================================== Modelo 1 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo normal
#-----------------------------------------------------------------------
nsamples = 10000
#--- Initial Values
x.0 = est.x(pound)/2
phi.0 =  0.985               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
data = list(T =length(pound), y= pound)

inits <-list(list(mu = -10, phiT=(phi.0 + 1)/2, s2=sigma.0^2))
fit.n <- stan(model_code=stan_code_n, data=data, init=inits, chains=1, iter=nsamples)

summary(fit.n)$summary[c(1,949,950), ] 

#----------- Uso do Coda:
v.n = extract(fit.n, pars = c('mu','phi','sigma'))
theta.n = cbind(v.n$mu,v.n$phi,v.n$sigma)
colnames(theta.n) = c('mu','phi','sigma')

m.n = window(as.mcmc(theta.n),start=1,thin=1)
summary(m.n)
effectiveSize(m.n)

x.est.n = summary(fit.n)$summary[c(4:948),c(1,4,8)]
ht.n = as.vector(x.est.n[,1])

#--- Cálculo do DIC, LOO e WAIC:
logl.n = extract_log_lik(fit.n, parameter_name = "loglik1")   #logp(yi|tt)
lik.n = extract_log_lik(fit.n, parameter_name = "loglik")     #logp(y|tt(S))

                   # Prior Gamma                #Prior IG                      #Prior Inv-X²
Loo.n = loo(logl.n)   #1809.3 (53.8), 1811.1 (54.1) # 1814.0 (54.5), 1809.6 (54.0) #  1811.7 (54.0) 1813.0 (54.3)
Loo.n
Waic.n = waic(logl.n) #1805.0 (53.1), 1807.1 (53.4) # 1810.3 (53.8), 1805.8 (53.4) #  1808.0 (53.4) 1809.0 (53.5)
Waic.n

lik_hat.n =  sum(dnorm(pound,0,exp(ht.n/2),log=T)) # log(p(y|theta_hat))
#lik_hat = sum(log(dn(pound,exp(ht/2))))
DIC.n = dic(lik_hat.n,lik.n)
DIC.n

#---- Gráficos de Convergência 
#mcmcplot(m)
pars_names.n = c(expression(mu), expression(phi), expression(sigma))
#plot.mcmc(theta.n,pars_names, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))
plot.mcmc2(theta.n,pars_names.n, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))


#--- Gráfico de ajuste das volatilidades
#plot.vol(pound,x.est.n,col0=c(1,4,4))
plot.vol2(pound,x.est.n)

#write.table(x.est.sn,'/home/ddias/Dropbox/Images_Paper/x.est.n.txt')

#save.image("/home/david/Documentos/SVPound_n1")

#=========================================== Modelo 2 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo t-student
#-----------------------------------------------------------------------
#--- Initial Values
x.0 = est.x(pound)/2
phi.0 =  0.985               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
nu.0 = 6

data = list(T =length(pound), y= pound)

#h=x.0
inits <-list(list(mu = -10, phiT=(phi.0 + 1)/2, s2=sigma.0^2, nu=nu.0))
fit.t <- stan(model_code=stan_code_t, data=data, init=inits, chains=1, iter=nsamples)

summary(fit.t)$summary[c(1,4,950,951), ] 

#----------- Uso do Coda:
v.t = extract(fit.t, pars = c('mu','phi','sigma','nu'))
theta.t = cbind(v.t$mu,v.t$phi,v.t$sigma,v.t$nu)
colnames(theta.t) = c('mu','phi','sigma','nu')

m.t = window(as.mcmc(theta.t),start=1,thin=1)
summary(m.t)
effectiveSize(m.t)

x.est.t = summary(fit.t)$summary[c(5:949),c(1,4,8)]
ht.t = as.vector(x.est.t[,1])
nu_hat.t = summary(m.t)$statistics[4]

#--- Cálculo do DIC, LOO e WAIC:
logl.t = extract_log_lik(fit.t, parameter_name = "loglik1")   #logp(yi|tt)
lik.t = extract_log_lik(fit.t, parameter_name = "loglik")     #logp(y|tt(S))
                   # prior Gamma                 #Prior IG                     #Prior Inv-X²
Loo.t = loo(logl.t)   #1813.4 (52.9), 1811.5 (52.9) # 1812.3 (52.9), 1814.2 (52.9) # 1814.1 (52.9), 1815.7 (53.0)
Loo.t
Waic.t = waic(logl.t) #1812.8 (52.8), 1810.7 (52.9) # 1811.5 (52.8), 1813.4 (52.8) # 1813.5 (52.8), 1815.1 (52.9)
Waic.t

lik_hat.t =  sum(log(dt.st(pound, nu_hat.t, exp(ht.t/2))))
DIC.t = dic(lik_hat.t,lik.t)
DIC.t

#---- Gráficos de Convergência 
pars_names= c(expression(mu), expression(phi), expression(sigma), expression(nu))
plot.mcmc2(theta.t,pars_names, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))

#--- Gráfico
plot.vol2(pound,x.est.t)

#save.image("/home/david/Documentos/Quali_Aplication/SVmodels/Prior_InvX/SVPound_5mil_t2")


#=========================================== Modelo 3 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo skew-normal
#-----------------------------------------------------------------------
#--- Initial Values
x.0 = est.x(pound)/2
phi.0 =  0.985               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
alpha.0 = 0.01

data = list(T =length(pound), y= pound)

#h=x.0
inits <-list(list(mu = -10, phiT=(phi.0 + 1)/2, s2=sigma.0^2, alpha=alpha.0))
fit.sn <- stan(model_code=stan_code_sn, data=data, init=inits, chains=1, iter=nsamples)

summary(fit.sn)$summary[c(1,4,950,951), ] 

#----------- Uso do Coda:
v.sn = extract(fit.sn, pars = c('mu','phi','sigma','alpha'))
theta.sn = cbind(v.sn$mu,v.sn$phi,v.sn$sigma,v.sn$alpha)
colnames(theta.sn) = c('mu','phi','sigma','alpha')

m.sn = window(as.mcmc(theta.sn),start=1,thin=1)
summary(m.sn)
effectiveSize(m.sn)

x.est.sn = summary(fit.sn)$summary[c(5:949),c(1,4,8)]
ht.sn = as.vector(x.est.sn[,1])
alpha_hat.sn = summary(m.sn)$statistics[4]

#--- Cálculo do DIC, LOO e WAIC:
logl.sn = extract_log_lik(fit.sn, parameter_name = "loglik1")   #logp(yi|tt)
lik.sn = extract_log_lik(fit.sn, parameter_name = "loglik")     #logp(y|tt(S))
                   # prior Gamma                #Prior IG                      #Prior Inv-X²
Loo.sn = loo(logl.sn)   #1812.7 (54.2), 1812.1 (54.0) # 1812.0 (54.0), 1811.4 (53.9) # 1811.6 (53.9), 1811.2 (54.0)
Loo.sn
Waic.sn = waic(logl.sn) #1808.8 (53.5), 1808.2 (53.3) # 1808.3 (53.4), 1807.7 (53.3) # 1808.5 (53.4), 1807.6 (53.3)
Waic.sn

#lik_hat = sum(dsn(pound, xi=0, omega=exp(ht/2), alpha=alpha_hat, log=T) )   
lik_hat.sn = sum(log(dsn(pound, alpha_hat.sn, exp(ht.sn/2))))
DIC.sn = dic(lik_hat.sn,lik.sn)
DIC.sn

#---- Gráficos de Convergência 
pars_names= c(expression(mu), expression(phi), expression(sigma), expression(nu))
plot.mcmc2(theta.sn,pars_names, cex.lab0=1.4, mar0 = c(1.8,4.5,1.9,1.9))

#--- Gráfico
plot.vol2(pound,x.est.sn)

#write.table(x.est.sn,'/home/ddias/Dropbox/Images_Paper/x.est.sn.txt')

#save.image("/home/david/Documentos/Quali_Aplication/SVmodels/Prior_InvX/SVPound_5mil_sn2")


#=========================================== Modelo 4 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo ExpModNorm
#-----------------------------------------------------------------------
#--- Initial Values
x.0 = est.x(pound)/2
phi.0 =  0.985               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
lambda.0 = 5

data = list(T =length(pound), y= pound)

#h=x.0
inits <-list(list(mu = -10, phiT=(phi.0 + 1)/2, s2=sigma.0^2, lambda=lambda.0))
fit.emg <- stan(model_code=stan_code_emn, data=data, init=inits, chains=1, iter=nsamples)

summary(fit.emg)$summary[c(1,4,950,951), ] 

#----------- Uso do Coda:
v.emg = extract(fit.emg, pars = c('mu','phi','sigma','lambda'))
theta.emg = cbind(v.emg$mu,v.emg$phi,v.emg$sigma,v.emg$lambda)
colnames(theta.emg) = c('mu','phi','sigma','lambda')

m.emg = window(as.mcmc(theta.emg),start=1,thin=1)
summary(m.emg)
effectiveSize(m.emg)

x.est.emg = summary(fit.emg)$summary[c(5:949),c(1,4,8)]
ht.emg = as.vector(x.est.emg[,1])
lambda_hat.emg = summary(m.emg)$statistics[4]

#--- Cálculo do DIC, LOO e WAIC:
logl.emg = extract_log_lik(fit.emg, parameter_name = "loglik1")   #logp(yi|tt)
lik.emg = extract_log_lik(fit.emg, parameter_name = "loglik")     #logp(y|tt(S))
					       #prior Gamma                  #prior IG                     #Prior Inv-X²
Loo.emg = loo(logl.emg)   # 1849.1 (54.3), 1850.4 (54.3) | 1849.4 (54.1), 1848.9 (54.3) #1850.8 (54.3) 1851.0 (54.4)
Loo.emg
Waic.emg = waic(logl.emg) # 1845.5 (53.6), 1847.0 (53.7) | 1846.7 (53.6), 1846.0 (53.8) #1847.8 (53.7) 1848.3 (53.9)
Waic.emg

lik_hat.emg = sum(log(demn(pound, lambda_hat.emg, exp(ht.emg/2))))
DIC.emg = dic(lik_hat.emg,lik.emg)
DIC.emg

#---- Gráficos de Convergência 
pars_names= c(expression(mu), expression(phi), expression(sigma), expression(nu))
plot.mcmc2(theta.emg,pars_names, cex.lab0=1.4, mar0 = c(1.8,4.5,1.9,1.9))

#--- Gráfico
plot.vol2(pound,x.est.emg)

#save.image("/home/david/Documentos/Quali_Aplication/SVmodels/Prior_InvX/SVPound_5mil_emg2")

#=========================================== Modelo 5 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo GED
#-----------------------------------------------------------------------
#--- Initial Values
x.0 = est.x(pound)/2
phi.0 =  0.985               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
nu.0 = 1.5

data = list(T =length(pound), y= pound)

#h=x.0
inits <-list(list(mu = -10, phiT=(phi.0 + 1)/2, s2=sigma.0^2, nu=nu.0))
fit.ged <- stan(model_code=stan_code_ged, data=data, init=inits, chains=1, iter=nsamples)

summary(fit.ged)$summary[c(1,4,950,951), ] 

#----------- Uso do Coda:
v.ged = extract(fit.ged, pars = c('mu','phi','sigma','nu'))
theta.ged = cbind(v.ged$mu,v.ged$phi,v.ged$sigma,v.ged$nu)
colnames(theta.ged) = c('mu','phi','sigma','nu')

m.ged = window(as.mcmc(theta.ged),start=1,thin=1)
summary(m.ged)
effectiveSize(m.ged)

x.est.ged = summary(fit.ged)$summary[c(5:949),c(1,4,8)]
ht.ged = as.vector(x.est.ged[,1])
nu_hat.ged = summary(m.ged)$statistics[4]

#--- Cálculo do DIC, LOO e WAIC:
logl.ged = extract_log_lik(fit.ged, parameter_name = "loglik1")   #logp(yi|tt)
lik.ged = extract_log_lik(fit.ged, parameter_name = "loglik")     #logp(y|tt(S))
                  # prior Inv-X²                # prior Gamma                  # Prior IG
Loo.ged = loo(logl.ged)   #1813.9 (53.5), 1817.0 (53.6) # 1812.0 (53.4), 1815.0 (53.7) # 1813.2 (53.0), 1815.5 (53.6)
Loo.ged
Waic.ged = waic(logl.ged) #1812.0 (53.2), 1815.5 (53.3) # 1809.5 (53.0), 1812.2 (53.2) # 1810.4 (53.2), 1813.3 (53.2) 
Waic.ged

lik_hat.ged = sum(log(dged(pound, nu_hat.ged, exp(ht.ged/2))))
DIC.ged = dic(lik_hat.ged,lik.ged)
DIC.ged

#---- Gráficos de Convergência 
pars_names= c(expression(mu), expression(phi), expression(sigma), expression(nu))
plot.mcmc2(theta.ged,pars_names, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))

#--- Gráfico
plot.vol2(pound,x.est.ged)

#save.image("/home/david/Documentos/SVPound")




#as.Date(datas, "%d/%m/%Y")
#results[[1]]$summary
#plot(stoc$para)

#save.image("/home/david/Dropbox/SVPound")

#load("/home/david/Dropbox/SVPound")
#App = load('/home/david/Dropbox/Documentos/Mestrado ICMC-USP/Tópicos de Pesquisa I/Aplication1')




