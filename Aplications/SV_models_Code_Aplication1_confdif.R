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

#---- Carregando funções externas:
source('/home/david/Dropbox/r_code/SV_model_rstan_confdif.R') # códigos dos modelos SV
source('/home/david/Dropbox/r_code/plot_mcmc.R')
source('/home/david/Dropbox/r_code/plot_vol.R')


#-----------------------------------------------------------------------
# Análise descritiva do conjunto de dados Pound/USD
#-----------------------------------------------------------------------
# Pound/Dollar data
dd = scan("/home/david/Dropbox/Projeto_David/MALA/sv.txt") # data(svpdx)
pound = dd - mean(dd)                                      # série depreciada
n = length(dd)


data(svpdx)
pound = svpdx$pdx - mean(svpdx$pdx)
svpdx$pound = pound
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
#--- Initial Values
x.0 = est.x(pound)/2
phi.0 =  0.985               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
data = list(T =length(pound), y= pound)
beta.0 = sqrt(var(pound)/exp(0.5*sigma.0^2/(1-phi.0^2)))

inits <-list(list(beta = beta.0, phiT=(phi.0 + 1)/2, s2=sigma.0^2))
fit_n <- stan(model_code=stan_code_n, data=data, init=inits, chains=1, iter=5000)

summary(fit_n)$summary[c(1,949,950), ] 

#----------- Uso do Coda:
v = extract(fit_n, pars = c('beta','phi','sigma'))
theta_n = cbind(v$beta,v$phi,v$sigma)
colnames(theta_n) = c('beta','phi','sigma')

m = window(as.mcmc(theta_n),start=1,thin=1)
summary(m)
effectiveSize(m)

x.est.hmc_n = summary(fit_n)$summary[c(4:948),c(1,4,8)]
ht = as.vector(x.est.hmc_n[,1])
beta_hat = summary(m)$statistics[1]

#--- Cálculo do DIC, LOO e WAIC:
logl = extract_log_lik(fit_n, parameter_name = "loglik1")   #logp(yi|tt)
lik = extract_log_lik(fit_n, parameter_name = "loglik")     #logp(y|tt(S))
                    # prior IG                  #Prior Gamma                   # Prior Inv-X²
Loo = loo(logl)   #1811.3 (54.0), 1809.4 (53.8) # 1810.9 (54.3), 1810.3 (54.0) # 1812.8 (54.1), 1811.4 (54.1)
Loo

Waic = waic(logl) #1806.3 (53.2), 1803.9 (52.9) # 1803.0 (52.9), 1805.1 (53.1) # 1809.4 (53.5), 1805.5 (53.1)
Waic

lik_hat =  sum(dnorm(pound,0,beta_hat*exp(ht/2),log=T)) # log(p(y|theta_hat))
DIC = dic(lik_hat,lik)
DIC

#---- Gráficos de Convergência 
#mcmcplot(m)
pars_names = c(expression(beta), expression(phi), expression(sigma))
#plot.mcmc(theta_n,pars_names, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))
plot.mcmc2(theta_n,pars_names, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))


#--- Gráfico de ajuste das volatilidades
#plot.vol(pound,x.est.hmc_n,col0=c(1,4,4))
plot.vol2(pound,x.est.hmc_n)

#save.image("/home/david/Documentos/Quali_Aplication/Modified_SVmodels/Prior_InvX/SVPound_5mil_n2")


#=========================================== Modelo 2 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo t-student
#-----------------------------------------------------------------------
#--- Initial Values
x.0 = est.x(pound)/2
phi.0 =  0.985               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
nu.0 = 6
beta.0 = sqrt(var(pound)/exp(0.5*sigma.0^2/(1-phi.0^2)))

data = list(T =length(pound), y= pound)

#h=x.0
inits <-list(list(beta = beta.0, phiT=(phi.0 + 1)/2, s2=sigma.0^2, nu=nu.0))
fit_t <- stan(model_code=stan_code_t, data=data, init=inits, chains=1, iter=5000)

summary(fit_t)$summary[c(1,4,950,951), ] 

#----------- Uso do Coda:
v = extract(fit_t, pars = c('beta','phi','sigma','nu'))
theta_t = cbind(v$beta,v$phi,v$sigma,v$nu)
colnames(theta_t) = c('beta','phi','sigma','nu')

m= window(as.mcmc(theta_t),start=1,thin=1)
summary(m)
effectiveSize(m)

x.est.hmc_t = summary(fit_t)$summary[c(5:949),c(1,4,8)]
beta_hat = summary(m)$statistics[1]
nu_hat = summary(m)$statistics[4]
ht = as.vector(x.est.hmc_t[,1])

#--- Cálculo do DIC, LOO e WAIC:
logl = extract_log_lik(fit_t, parameter_name = "loglik1")   #logp(yi|tt)
lik = extract_log_lik(fit_t, parameter_name = "loglik")     #logp(y|tt(S))
                   #prior IG                     # Prior Gamma                 # Prior Inv-X²
Loo = loo(logl)   #1812.0 (52.9), 1813.5 (53.0) # 1812.2 (52.9), 1811.8 (53.0) # 1813.2 (53.0), 1812.1 (52.9)
Loo

Waic = waic(logl) #1810.9 (52.8), 1812.6 (52.9) # 1811.2 (52.7), 1810.7 (52.8) # 1812.5 (52.9), 1811.2 (52.8)
Waic

lik_hat = sum(log(dt.st(pound,nu_hat,beta_hat*exp(ht/2)))) # log(p(y|theta_hat))
DIC = dic(lik_hat,lik)
DIC

#---- Gráficos de Convergência 
pars_names= c(expression(beta), expression(phi), expression(sigma), expression(nu))
plot.mcmc2(theta_t,pars_names, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))

#--- Gráfico
plot.vol2(pound,x.est.hmc_t)

#save.image("/home/david/Documentos/Quali_Aplication/Modified_SVmodels/Prior_InvX/SVPound_5mil_t2")


#=========================================== Modelo 3 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo skew-normal
#-----------------------------------------------------------------------
#--- Initial Values
x.0 = est.x(pound)/2
phi.0 =  0.985               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
alpha.0 = 0.01
beta.0 = sqrt(var(pound)/exp(0.5*sigma.0^2/(1-phi.0^2)))

data = list(T =length(pound), y= pound)

#h=x.0
inits <-list(list(beta = beta.0, phiT=(phi.0 + 1)/2, s2=sigma.0^2, alpha=alpha.0))
fit_sn <- stan(model_code=stan_code_sn, data=data, init=inits, chains=1, iter=5000)

summary(fit_sn)$summary[c(1,4,950,951), ] 

#----------- Uso do Coda:
v = extract(fit_sn, pars = c('beta','phi','sigma','alpha'))
theta_sn = cbind(v$beta,v$phi,v$sigma,v$alpha)
colnames(theta_sn) = c('beta','phi','sigma','alpha')

m= window(as.mcmc(theta_sn),start=1,thin=1)
summary(m)
effectiveSize(m)

x.est.hmc_sn = summary(fit_sn)$summary[c(5:949),c(1,4,8)]
ht = as.vector(x.est.hmc_sn[,1])
beta_hat = summary(m)$statistics[1]
alpha_hat = summary(m)$statistics[4]

#--- Cálculo do DIC, LOO e WAIC:
logl = extract_log_lik(fit_sn, parameter_name = "loglik1")   #logp(yi|tt)
lik = extract_log_lik(fit_sn, parameter_name = "loglik")     #logp(y|tt(S))
                   #prior iG                    #Prior Gamma                   # Prior Inv-X²
Loo = loo(logl)   #1809.2 (53.7), 1809.2 (53.7) # 1812.3 (54.0), 1811.7 (54.0) # 1810.6 (53.7), 1811.5 (54.1)
Loo
Waic = waic(logl) #1804.3 (52.9), 1803.0 (52.7) # 1808.1 (53.3), 1805.8 (53.0) # 1806.3 (53.0), 1806.3 (53.2)
Waic

lik_hat = sum(log(dsn(pound, alpha_hat, beta_hat*exp(ht/2))))
DIC = dic(lik_hat,lik) # 1801.438
DIC

#---- Gráficos de Convergência 
pars_names= c(expression(beta), expression(phi), expression(sigma), expression(alpha))
plot.mcmc2(theta_sn,pars_names, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))

#--- Gráfico
plot.vol2(pound,x.est.hmc_sn)

#save.image("/home/david/Documentos/Quali_Aplication/Modified_SVmodels/Prior_InvX/SVPound_5mil_sn2")

#=========================================== Modelo 4 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo ExpModNorm
#-----------------------------------------------------------------------
#--- Initial Values
x.0 = est.x(pound)/2
phi.0 =  0.985               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
lambda.0 = 5
beta.0 = sqrt(var(pound)/exp(0.5*sigma.0^2/(1-phi.0^2)))

data = list(T =length(pound), y= pound)

#h=x.0
inits <-list(list(beta=beta.0, phiT=(phi.0 + 1)/2, s2=sigma.0^2, lambda=lambda.0))
fit_emn <- stan(model_code=stan_code_emn, data=data, init=inits, chains=1, iter=5000)

summary(fit_emn)$summary[c(1,4,950,951), ] 

#----------- Uso do Coda:
v = extract(fit_emn, pars = c('beta','phi','sigma','lambda'))
theta_emn = cbind(v$beta,v$phi,v$sigma,v$lambda)
colnames(theta_emn) = c('beta','phi','sigma','lambda')

m= window(as.mcmc(theta_emn),start=1,thin=1)
summary(m)
effectiveSize(m)

x.est.hmc_emn = summary(fit_emn)$summary[c(5:949),c(1,4,8)]
ht = as.vector(x.est.hmc_emn[,1])
beta_hat = summary(m)$statistics[1]
lambda_hat = summary(m)$statistics[4]

#--- Cálculo do DIC, LOO e WAIC:
logl = extract_log_lik(fit_emn, parameter_name = "loglik1")   #logp(yi|tt)
lik = extract_log_lik(fit_emn, parameter_name = "loglik")     #logp(y|tt(S))
                    #prior IG                   # Prior Gamma                  # Prior Inv-X²
Loo = loo(logl)   #1849.5 (54.3), 1847.6 (53.7) # 1846.2 (53.5), 1848.8 (54.0) # 1847.6 (53.8), 1848.1 (54.0)
Loo
Waic = waic(logl) #1844.5 (53.5), 1844.1 (53.2) # 1841.4 (52.9), 1843.2 (53.1) # 1842.9 (53.1), 1842.9 (53.1)
Waic

lik_hat = sum(log(demn(pound, lambda_hat, beta_hat*exp(ht/2))))  # E(logp(y|theta)) D.bar  
DIC = dic(lik_hat,lik)
DIC

#---- Gráficos de Convergência 
pars_names= c(expression(beta), expression(phi), expression(sigma), expression(lambda))
plot.mcmc2(theta_emn,pars_names, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))

#--- Gráfico
plot.vol2(pound,x.est.hmc_emn)

#save.image("/home/david/Documentos/Quali_Aplication/Modified_SVmodels/Prior_InvX/SVPound_5mil_emg2")


#=========================================== Modelo 5 =============================================
#-----------------------------------------------------------------------
# Apliação ao conjunto de dados - rstan modelo GED
#-----------------------------------------------------------------------
#--- Initial Values
x.0 = est.x(pound)/2
phi.0 =  0.985               
sigma.0 = 0.075 #var(x.0)*(1-phi.0^2)
nu.0 = 1.5
beta.0 = sqrt(var(pound)/exp(0.5*sigma.0^2/(1-phi.0^2)))

data = list(T =length(pound), y= pound)

#h=x.0
inits <-list(list(beta = beta.0, phiT=(phi.0 + 1)/2, s2=sigma.0^2, nu=nu.0))
fit_ged <- stan(model_code=stan_code_ged, data=data, init=inits, chains=1, iter=5000)

summary(fit_ged)$summary[c(1,4,950,951), ] 

#----------- Uso do Coda:
v = extract(fit_ged, pars = c('beta','phi','sigma','nu'))
theta_ged = cbind(v$beta,v$phi,v$sigma,v$nu)
colnames(theta_ged) = c('beta','phi','sigma','nu')

m= window(as.mcmc(theta_ged),start=1,thin=1)
summary(m)
effectiveSize(m)

x.est.hmc_ged = summary(fit_ged)$summary[c(5:949),c(1,4,8)]
ht = as.vector(x.est.hmc_ged[,1])
beta_hat = summary(m)$statistics[1]
nu_hat = summary(m)$statistics[4]

#--- Cálculo do DIC, LOO e WAIC:
logl = extract_log_lik(fit_ged, parameter_name = "loglik1")   #logp(yi|tt)
lik = extract_log_lik(fit_ged, parameter_name = "loglik")     #logp(y|tt(S))
                  #prior IG                     #Prior Gamma                   # Prior Inv-X²
Loo = loo(logl)   #1812.5 (53.5), 1813.5 (53.5) # 1810.8 (53.3), 1814.3 (53.5) # 1811.1 (53.3), 1812.0 (53.4)
Loo
Waic = waic(logl) #1808.5 (52.9), 1810.3 (53.0) # 1807.2 (52.7), 1811.4 (53.1) # 1808.3 (52.9), 1808.1 (52.8)
Waic

lik_hat = sum(log( dged(pound, nu_hat, beta_hat*exp(ht/2))))
DIC = dic(lik_hat,lik) #Prior Inv-X² 1806.948, 1807.03
DIC

#---- Gráficos de Convergência 
pars_names= c(expression(beta), expression(phi), expression(sigma), expression(nu))
plot.mcmc2(theta_ged,pars_names, cex.lab0=1.4, mar0 = c(2,4.5,1.9,1.9))

#--- Gráfico
plot.vol2(pound,x.est.hmc_ged)

#save.image("/home/david/Documentos/Quali_Aplication/Modified_SVmodels/Prior_InvX/SVPound_5mil_ged2")




