rm(list=ls())    
# install.packages("here")
setwd(here::here())

######## 1
# Read data files
msft = read.table("MSFT.txt",header=FALSE)
xom = read.table("XOM.txt",header=FALSE)
bac = read.table("BAC.txt",header=FALSE)
market = read.table("market.txt",header=FALSE)


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
source("llik_OD_regression.R") 

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
  
  return(list(theta_hat = theta_hat, beta_T = beta[length(beta)]))
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
get_beta_t_plus_1 <- function(theta_hat, x_t, beta_t){
  omega_hat <- theta_hat[1]
  phi_hat <- theta_hat[2]
  alpha_hat <- theta_hat[3]
  sigma2_hat <- theta_hat[4]
  
  
}
