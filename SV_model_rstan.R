#--------------------------------------------------------------------------
# SV model:
# yt = exp(ht/2)
# ht = mu + phi*(h[t-1] - mu) + etat,     etat~N(0,sigeta)
# ht ~ N(mu, sigeta/(1-phi²))
#--------------------------------------------------------------------------
#-------------------------------------------------------------------------
# Stochastic Volatility models normal with rstan 
#-------------------------------------------------------------------------
stan_code_n = "
data{
  int<lower=0> T;   // time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters{
  real mu;                      // mean log volatility
  real<lower=-1,upper=1> phiT;  // persistence of volatility
  real<lower=0> s2;             // white noise shock scale
  vector[T] h;                  // log volatility at time t
}
transformed parameters{
 real<lower=-1,upper=1> phi;
 real<lower=0> sigma;
 phi = (2*phiT - 1);
 sigma = pow(s2,2);
}
model{
  mu ~ normal(-10, 1);
  phiT ~ beta(20, 1.5);
  s2 ~ gamma(0.5, 1/(0.1*2)); //s2 ~ inv_gamma(2.5,0.025); // s2 ~ scaled_inv_chi_square(10, 0.05);            // perguntar s2 ~ cauchy(0,5);

  //--- Sampling volatilitys:
  h[1] ~ normal(mu, sigma/sqrt(1 - phi*phi));
  for(t in 2:T)
     {
      h[t] ~ normal(mu + phi*(h[t - 1] - mu), sigma);
     }
  
  //--- Sampling observations:  
  for(t in 1:T)
     {
      y[t] ~ normal(0, exp(h[t]/2));
     }
}
//--- To calculate the DIC, WAIC and LOO:
generated quantities{
 real loglik;
 vector[T] loglik1;
 loglik = 0;
 for(t in 1:T)
    {
     loglik = loglik + normal_lpdf(y[t]| 0, exp(h[t]/2));      //dic
     loglik1[t] = normal_lpdf(y[t]| 0, exp(h[t]/2));           //waic
    }
}
"

#-------------------------------------------------------------------------
# Stochastic Volatility models t-student with rstan 
#-------------------------------------------------------------------------
stan_code_t = '
data{
  int<lower=0> T;   // time points (equally spaced)
  vector[T] y;        // mean corrected return at time t
}
parameters{
 real mu;                       // mean log volatility
 real<lower=-1,upper=1> phiT;   // persistence of volatility
 real<lower=0> s2;              // white noise shock scale
 real<lower=4> nu;              // parameter of the t-student distribution
 vector[T] h;                   // log volatility at time t
}
transformed parameters{
 real<lower=-1,upper=1> phi;
 real<lower=0> sigma;
 phi = (2*phiT - 1);
 sigma = pow(s2,2);
}
model{
 mu ~ normal(-10, 1);
 phiT ~ beta(20, 1.5);
 s2 ~ gamma(0.5, 1/(0.1*2)); //s2 ~ inv_gamma(2.5,0.025); //s2 ~ scaled_inv_chi_square(10, 0.05);         
 nu ~ exponential(0.333333);

  //--- Sampling volatilitys:
  h[1] ~ normal(mu, sigma/sqrt(1 - phi*phi));
  for(t in 2:T)
     {
      h[t] ~ normal(mu + phi*(h[t - 1] - mu), sigma);
     }
 
  //--- Sampling observations:  
  for(t in 1:T)
     {
      y[t] ~ student_t(nu,0, exp(h[t]/2));
     }
}
//--- To calculate the DIC and WAIC:
generated quantities{
 real loglik;
 vector[T] loglik1;
 loglik = 0;
 for(t in 1:T)
    {
     loglik = loglik + student_t_lpdf(y[t]| nu, 0, exp(h[t]/2));
     loglik1[t] = student_t_lpdf(y[t]| nu, 0, exp(h[t]/2));           
    }
}
'

#-------------------------------------------------------------------------
# Stochastic Volatility models skew-normal with rstan 
#-------------------------------------------------------------------------
stan_code_sn = "
data{
  int<lower=0> T;   // time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters{
  real mu;                      // mean log volatility
  real<lower=-1,upper=1> phiT;  // persistence of volatility
  real<lower=0> s2;             // white noise shock scale
  real alpha;                    // shape skew-normal distribution  
  vector[T] h;                  // log volatility at time t
}
transformed parameters{
 real<lower=-1,upper=1> phi;
 real<lower=0> sigma;
 phi = (2*phiT - 1);
 sigma = pow(s2,2);
}
model{
  mu ~ normal(-10, 1);
  phiT ~ beta(20, 1.5);
  s2 ~ gamma(0.5, 1/(0.1*2)); // s2 ~ inv_gamma(2.5,0.025); //s2 ~ scaled_inv_chi_square(10, 0.05);                
  alpha ~ normal(0,1);

  //--- Sampling volatilitys:
  h[1] ~ normal(mu, sigma/sqrt(1 - phi*phi));
  for(t in 2:T)
     {
      h[t] ~ normal(mu + phi*(h[t - 1] - mu), sigma);
     }
  
  //--- Sampling observations:  
  for(t in 1:T)
     {
      y[t] ~ skew_normal(0, exp(h[t]/2), alpha);
     }
}
//--- To calculate the DIC, WAIC and LOO:
generated quantities{
 real loglik;
 vector[T] loglik1;
 loglik = 0;
 for(t in 1:T)
    {
     loglik = loglik + skew_normal_lpdf(y[t]| 0, exp(h[t]/2), alpha);
     loglik1[t] = skew_normal_lpdf(y[t]| 0, exp(h[t]/2), alpha);           
    }
}
"

#-------------------------------------------------------------------------
# Stochastic Volatility models GED with rstan 
#-------------------------------------------------------------------------
stan_code_ged = "
functions{
  real ged_lpdf(real y, real h, real n)
  {
   real lb2;
   real lpdf;
   lb2 = (2^(-2/n)) * tgamma(1/n) / tgamma(3/n); 
   lpdf = - (h/2) - (1+(1/n))*log(2) + log(n) - log(sqrt(lb2)) - lgamma(1/n) - (0.5)*(fabs(y/(exp(h/2)*(sqrt(lb2))))^n);                              
   return lpdf;
  }
}
data{
 int<lower=0> T;                // time points (equally spaced)
 vector[T] y;                   // mean corrected return at time t
}
parameters{
 real mu;                       // mean log volatility
 real<lower=-1,upper=1> phiT;   // persistence of volatility
 real<lower=0> s2;              // white noise shock scale
 real<lower=0> nu;              // parameter of the GED distribution
 vector[T] h;                   // log volatility at time t
}
transformed parameters{
 real<lower=-1,upper=1> phi;
 real<lower=0> sigma;
 phi = (2*phiT - 1);
 sigma = pow(s2,2);
}
model{
 mu ~ normal(-10, 1);
 phiT ~ beta(20, 1.5);
 s2 ~ gamma(0.5, 1/(0.1*2)); // s2 ~ inv_gamma(2.5,0.025); // s2 ~ scaled_inv_chi_square(10, 0.05);                    
 nu ~ scaled_inv_chi_square(10, 0.05);            // nu ~ inv_gamma(2.5,0.025);           
 
 //--- Sampling volatilitys:
  h[1] ~ normal(mu, sigma/sqrt(1 - phi*phi));
  for(t in 2:T)
     {
      h[t] ~ normal(mu + phi*(h[t - 1] - mu), sigma);
     }
  
 //--- Sampling observations:  
 for(t in 1:T)
    {
     target += ged_lpdf(y[t]|h[t], nu);
    }
}
//--- To calculate the DIC and WAIC:
generated quantities{
 real loglik;
 vector[T] loglik1;
 loglik = 0;
 for(t in 1:T)
   {
    loglik = loglik + ged_lpdf(y[t]| h[t], nu);
    loglik1[t] = ged_lpdf(y[t]| h[t], nu);           
   }
}
"

#======================================= End Models ===========================================
#-------------------------------------------------------------------------
#                           Aux. Functions
#-------------------------------------------------------------------------

# error function
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

# complementary error function
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)

#Initial values
est.x <- function(y) {
  # Initial estimate of the latent variable.
  log(y^2)-mean(log(y^2))
}

dn = function(y, sigma) (1/(sqrt(2*pi)*sigma))*exp(-((y^2)/(2*(sigma^2))))

dt.st = function(y, nu, sigma) (1/(sqrt(pi*(nu-2))*sigma)) * (gamma((nu+1)/2)/gamma(nu/2)) * ((1+((y^2)/((nu-2)*(sigma^2))))^(-(nu+1)/2)) 

dsn = function(y, nu, sigma) (1/(sqrt(2*pi)*sigma)) * exp(-((y^2)/(2*(sigma^2)))) * (1+erf(nu*(y/(sqrt(2)*sigma)))) 

demn = function(y, nu, sigma) (nu/2) * exp((nu/2) * (nu*(sigma^2) - 2*y)) * erfc((nu*(sigma^2) - y)/(sqrt(2)*sigma))

dged = function(y, nu, sigma){
       lambda = sqrt((2^(-2/nu)) * gamma(1/nu) / gamma(3/nu)) 
       p = (nu/( sigma*lambda*(2^(1+(1/nu)))*gamma(1/nu))) * exp(-(1/(2*((lambda)^nu)))*(abs(y/sigma)^nu)) 
       return(p)
}


# deviance information criterion:
dic = function(lik_hat, lik)
      {
        pD = 2*(lik_hat - mean(lik))      # mean(lik) = E(logp(y|theta)) D.bar
        DIC = -2*lik_hat + 2*pD  
        return(DIC)
      }                










