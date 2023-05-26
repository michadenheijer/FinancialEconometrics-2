rm(list=ls())    
# install.packages("here")
setwd(here::here())

######## 1
# Read data files
msft = read.table("/Users/micha/Documents/GitHub/FinancialEconometrics-2/MSFT.txt",header=FALSE)
xom = read.table("/Users/micha/Documents/GitHub/FinancialEconometrics-2/XOM.txt",header=FALSE)
bac = read.table("/Users/micha/Documents/GitHub/FinancialEconometrics-2/BAC.txt",header=FALSE)
market = read.table("/Users/micha/Documents/GitHub/FinancialEconometrics-2/market.txt",header=FALSE)

plot(market[,1], type="l")
plot(xom[,1], type="l")



# Get log returns
msft_log_returns = diff(log(msft[,1]) * 100)
xom_log_returns = diff(log(xom[,1]) * 100)
bac_log_returns = diff(log(bac[,1]) * 100)
market_log_returns = diff(log(market[,1]) * 100)

# Plot log returns
plot(msft_log_returns, type="l", xlab="t", ylab="MSFT Single Day Return (%)")
plot(xom_log_returns, type="l", xlab="t", ylab="XOM Single Day Return (%)")
plot(bac_log_returns, type="l", xlab="t", ylab="BAC Single Day Return (%)")
plot(market_log_returns, type="l", xlab="t", ylab="Market Single Day Return (%)")

# Estimate beta, using OLS (without intercept)
est_beta <- function(market_log_returns, stock_log_returns) {
  return (solve(t(market_log_returns) %*% market_log_returns) %*% t(market_log_returns) %*% stock_log_returns)
}

# MSFT beta
est_beta(market_log_returns, msft_log_returns)

# XOM beta
est_beta(market_log_returns, xom_log_returns)

# BAC beta
est_beta(market_log_returns, bac_log_returns)

####### 3

# Import observation driven model
source("/Users/micha/Documents/GitHub/FinancialEconometrics-2/llik_OD_regression.R") 

## This code is for finding the parameters and plotting beta_t
estObservationDrivenModelAndPlot <- function(xt, yt, title) {
  a <- 0.2/sd(xt*yt) # initial value for alpha
  phi <- 0.9  # initial value for beta
  omega <- (cov(xt,yt)/var(xt))*(1-phi) # initial value for omega
  sig2 <- var(yt)
  
  # Add all init values to a vector
  par_ini <- c(omega,log(phi/(1-phi)),log(a),log(sig2))
  
  # Optimize the log likelyhood function
  est <- optim(par=par_ini,fn=function(par)-llik_OD_regression(yt,xt,par), method = "BFGS")
  
  # Get the estimated values
  omega_hat <- est$par[1]
  phi_hat <- exp(est$par[2])/(1+exp(est$par[2]))
  alpha_hat <- exp(est$par[3])
  sigma2_hat <- exp(est$par[4])
  
  # Add all to single vector
  theta_hat <- c(omega_hat,phi_hat,alpha_hat,sigma2_hat)
  
  ##### PLOT
  n <- length(xt)
  beta <- rep(0,n)
  beta[1] <- omega_hat/(1-phi_hat) # initialize betat at unconditional expectation
  
  # Create the plot
  for(t in 2:n){
    beta[t] <- omega_hat+phi_hat*beta[t-1]+alpha_hat*(yt[t-1]-beta[t-1]*xt[t-1])*xt[t-1];
  }
  
  # Show the plot
  plot(beta,type="l", col=2, main= paste(title, " beta_t"), ylab="",xlab="")
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted") # Grid
  
  return(list(theta_hat = theta_hat, beta_T = beta[length(beta)], x_T = xt[length(xt)], y_T = yt[length(yt)]))
}

## MSFT
msft_par <- estObservationDrivenModelAndPlot(market_log_returns, msft_log_returns, "MSFT")
msft_par

## XOM
xom_par <- estObservationDrivenModelAndPlot(market_log_returns, xom_log_returns, "XOM")
xom_par

## BAC
bac_par <- estObservationDrivenModelAndPlot(market_log_returns, bac_log_returns, "BAC")
bac_par

###### 4
get_beta_t_plus_1 <- function(theta_hat, x_T, y_T, beta_T){
  omega_hat <- theta_hat[1]
  phi_hat <- theta_hat[2]
  alpha_hat <- theta_hat[3]
  sigma2_hat <- theta_hat[4]
  
  beta_T_plus_1 <- omega_hat + phi_hat * beta_T + alpha_hat * (y_T - beta_T * x_T)*x_T
  
  return(beta_T_plus_1)
}

# MSFT
get_beta_t_plus_1(msft_par$theta_hat, msft_par$x_T, msft_par$y_T, msft_par$beta_T)

# XOM
get_beta_t_plus_1(xom_par$theta_hat, xom_par$x_T, xom_par$y_T, xom_par$beta_T)

# BAC
get_beta_t_plus_1(bac_par$theta_hat, bac_par$x_T, bac_par$y_T, bac_par$beta_T)



####### 5

# Import GARCH likelyhood func
source("/Users/micha/Documents/GitHub/FinancialEconometrics-2/llik_fun_GARCH.R")

estBeta_t_CCC <- function(xt, yt){
  # estimating GARCH of xt
  alpha_ini <- 0.2 
  beta_ini <- 0.6  
  omega_ini <- var(xt)*(1-alpha_ini-beta_ini) 
  par_ini <- c(log(omega_ini),log(alpha_ini/(1-alpha_ini)),log(beta_ini/(1-beta_ini)))
  
  est1 <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par,xt), method = "BFGS")
  
  
  (omega_hat1 <- exp(est1$par[1]))
  (alpha_hat1 <- exp(est1$par[2])/(1+exp(est1$par[2])))
  (beta_hat1 <- exp(est1$par[3])/(1+exp(est1$par[3])))
  
  # estimating GARCH of yt
  alpha_ini <- 0.2 
  beta_ini <- 0.6  
  omega_ini <- var(yt)*(1-alpha_ini-beta_ini) 
  par_ini <- c(log(omega_ini),log(alpha_ini/(1-alpha_ini)),log(beta_ini/(1-beta_ini)))
  
  est2 <- optim(par=par_ini,fn=function(par)-llik_fun_GARCH(par, yt), method = "BFGS")
  
  (omega_hat2 <- exp(est2$par[1]))
  (alpha_hat2 <- exp(est2$par[2])/(1+exp(est2$par[2])))
  (beta_hat2 <- exp(est2$par[3])/(1+exp(est2$par[3])))
  
  # get corr
  n <- length(xt)
  
  s1 <- rep(0,n)
  s2 <- rep(0,n)
  
  s1[1] <- var(xt)
  s2[1] <- var(yt)
  
  for(t in 2:n){
    s1[t] <- omega_hat1 + alpha_hat1*xt[t-1]^2 + beta_hat1*s1[t-1]
    s2[t] <- omega_hat2 + alpha_hat2*yt[t-1]^2 + beta_hat2*s2[t-1]
  }
  
  e1 <- xt/sqrt(s1)
  e2 <- yt/sqrt(s2)
  
  r <- cor(e1,e2)
  
  s12=r*sqrt(s1)*sqrt(s2)
  
  ## s1 is market variance, s2 is stock variance
  beta <- s12/s1
  plot(beta, type="l")
  
  # w1, a1, b1 from variance equation 1st series:
  theta_hat1 <- c(omega_hat1,alpha_hat1,beta_hat1)
  
  # w2, a2, b2 from variance equation 2nd series:
  theta_hat2 <- c(omega_hat2,alpha_hat2,beta_hat2)
}

# MSFT
estBeta_t_CCC(market_log_returns, msft_log_returns)

# XOM
estBeta_t_CCC(market_log_returns, xom_log_returns)

# BAC
estBeta_t_CCC(market_log_returns, bac_log_returns)


####### 6
source("/Users/micha/Documents/GitHub/FinancialEconometrics-2/sim_m_REG.R") 

estimateParModel <- function(xt, yt) {
  # set initial variables
  hb <- cov(yt,xt)/var(xt)
  yr <- yt-hb*xt
  
  xy <- yr*xt
  acvfxy <- acf(xy, lag.max=15, type ="covariance", plot=F)$acf[-1]
  sample_m <- c(var(yr),hb,acvfxy)
  
  n <- length(xt)
  M <- 20
  H <- M*n
  
  # start simulation
  set.seed(234)
  eta <- rnorm(H)
  eps <- rnorm(H)
  e <- cbind(eta,eps)
  
  a1 <- 0.95 
  a0 <- cov(xt,yt)/var(xt)*(1-a1)
  s_eta <- 0.2
  s_eps <- var(yr)
  
  par_ini <- c(a0, log(a1/(1-a1)), log(s_eta), log(s_eps))
  
  # Optimize from simulation for estimation
  est <- optim(par=par_ini,fn=function(par) mean((sim_m_REG(e,xt,par)-sample_m)^2), method = "Nelder-Mead", control = list(maxit = 100, reltol=0)) # Use different optimization method, because other fails
  
  # Obtain parameter extimate using the link functions
  a0_hat <- est$par[1]
  a1_hat <- exp(est$par[2])/(1+exp(est$par[2]))
  s_eta <- exp(est$par[3])
  s_eps <- exp(est$par[4])
  
  theta_hat <- c(a0_hat,a1_hat,s_eta,s_eps)
  return(theta_hat)
}

# MSFT
estimateParModel(market_log_returns, msft_log_returns)

# XOM
estimateParModel(market_log_returns, xom_log_returns)

# BAC
estimateParModel(market_log_returns, bac_log_returns)
