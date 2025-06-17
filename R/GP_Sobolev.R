# Implements the proposed plug-in Gaussian Process method with a Sobolev kernel
# under the assumption of homoscedastic (homogeneous) regression error.


library(Rsolnp)
library(MASS)
library(RandomFieldsUtils)


# Sobolev kernel function.
# Computes the Sobolev kernel matrix between x1 and x2.
# Inputs:
# - x1, x2: Numeric vectors of input locations.
# Output:
# - A matrix of kernel evaluations.
sob <- function(x1, x2){
  n1 <- length(x1)
  n2 <- length(x2)
  XX <- t(matrix(rep(x1, n2), ncol = n2))
  xx <- matrix(rep(x2, n1), ncol = n1)
  return(1 + XX * xx + pmin(XX, xx)^2 * (3 * pmax(XX, xx) - pmin(XX, xx)) / 6)
}


# First derivative of the Sobolev kernel with respect to the second argument.
# Inputs:
# - x1, x2: Numeric vectors of input locations.
# Output:
# - Matrix of first derivative kernel evaluations.
sob_prime <- function(x1, x2){
  n1 <- length(x1)
  n2 <- length(x2)
  XX <- t(matrix(rep(x1, n2), ncol = n2))
  xx <- matrix(rep(x2, n1), ncol = n1)
  return(XX + pmin(XX, xx) * XX - pmin(XX, xx)^2 / 2)
}


# Second derivative of the Sobolev kernel with respect to the second argument.
# Inputs:
# - x1, x2: Numeric vectors of input locations.
# Output:
# - Matrix of second derivative kernel evaluations.
sob_2prime <- function(x1, x2){
  n1 <- length(x1)
  n2 <- length(x2)
  XX <- t(matrix(rep(x1, n2), ncol = n2))
  xx <- matrix(rep(x2, n1), ncol = n1)
  return(1 + pmin(XX, xx))
}


# Estimates the regularization parameter (lambda) using empirical Bayes
# by maximizing the marginal likelihood for the Sobolev kernel.
# Inputs:
# - x: Design points.
# - y: Observed responses.
# Output:
# - Estimated lambda (scalar).
get_lambda_sobolev <- function(x, y){
  n <- length(x)
  K_XX <- sob(x, x)
  eigendecomposition <- eigen(K_XX)
  eigenvector <- eigendecomposition$vector
  eigenvalue <- eigendecomposition$values
  
  log_post_sobolev <- function(theta){
    diag.new <- eigenvalue + n * theta
    pred.comp <- c(eigenvector %*% (c(crossprod(eigenvector, y)) / diag.new))
    
    log_p <- n * log(mean(y * pred.comp)) + sum(log(diag.new))
    log_p}
  
  powell <- solnp(pars = c(0.01), fun = log_post_sobolev, LB = 2e-3, UB = 1)
  lambda <- powell$pars
  lambda}


# Computes the leave-one-out cross-validation error for the GP with Sobolev kernel.
# Inputs:
# - x, y: Design points and observed responses.
# - lambda: Regularization parameter.
# Output:
# - Scalar value of LOO mean squared error.
get_LOO_sobolev <- function(x, y, lambda){
  mse <- rep(0, n)
  for(i in 1:n){
    testy <- y[i]
    testx <- x[i]
    trainy <- y[-i]
    trainx <- x[-i]
    
    K_XX <- sob(trainx, trainx)
    K_Xx <- sob(trainx, testx)
    eigendecomposition <- eigen(K_XX)
    eigenvector <- eigendecomposition$vector
    eigenvalue <- eigendecomposition$values
    
    diag.new <- diag(1/(eigenvalue + (n-1) * lambda))
    pred.comp <- c(eigenvector %*% diag.new %*% c(crossprod(eigenvector, trainy)))
    
    f_hat <- c(K_Xx %*% pred.comp)
    
    mse[i] <- (f_hat - testy)^2
  }
  return(mean(mse))
}


# Main function: Fits a plug-in Gaussian Process model with Sobolev kernel under
# the assumption of homoscedastic (constant) noise and returns posterior estimates.
# 
# Inputs:
# - x, y: Design points and observed responses.
# - x_new: New input locations where predictions are desired.
# - alpha: Significance level for credible bands (default = 0.05).
# - LOO: Whether to return leave-one-out error (default = FALSE).
# - CI: Whether to return pointwise and simultaneous credible bands (default = FALSE).
# - RMSE: Whether to return RMSE compared to known truth f0_new and f0_new_prime (default = FALSE).
#
# Output (a list containing the following elements):
# - loo_error: Leave-one-out cross-validation error.
# - rmse: A numeric vector containing RMSEs for the estimated regression function and its first derivative.
# - f_hat: Estimated values of the regression function at x_new.
# - f_prime_hat: Estimated values of the first derivative at x_new.
# - f_hat_lb_pw: Lower bounds of the pointwise credible bands for the regression function.
# - f_hat_ub_pw: Upper bounds of the pointwise credible bands for the regression function.
# - f_hat_lb_st: Lower bounds of the simultaneous L-infinity credible bands for the regression function.
# - f_hat_ub_st: Upper bounds of the simultaneous L-infinity credible bands for the regression function.
# - f_prime_lb_pw: Lower bounds of the pointwise credible bands for the first derivative.
# - f_prime_ub_pw: Upper bounds of the pointwise credible bands for the first derivative.
# - f_prime_hat_lb_st: Lower bounds of the simultaneous L-infinity credible bands for the first derivative.
# - f_prime_hat_ub_st: Upper bounds of the simultaneous L-infinity credible bands for the first derivative.
get_GPR_sobolev <- function(x, y, x_new, alpha = 0.05, LOO = FALSE, CI = FALSE, RMSE = FALSE){
  n <- length(x)
  n_new <- length(x_new)
  
  lambda <- get_lambda_sobolev(x, y)
  loo_error <- NULL
  
  if(LOO == TRUE){
    loo_error <- get_LOO_sobolev(x, y, lambda)
  }

  K_XX <- sob(x, x)
  K_Xx <- sob(x, x_new)
  K_xx <- sob(x_new, x_new)
  K10_Xx <- sob_prime(x, x_new)
  K11_xx <- sob_2prime(x_new, x_new)
  
  eigendecomposition <- eigen(K_XX)
  eigenvector <- eigendecomposition$vector
  eigenvalue <- eigendecomposition$values
  
  diag.new <- diag(1/(eigenvalue + n * lambda))
  pred.comp <- c(eigenvector %*% diag.new %*% c(crossprod(eigenvector, y)))
  K_Xx.eigenvector <- K_Xx %*% eigenvector
  K10_Xx.eigenvector <- K10_Xx %*% eigenvector
  
  f_hat <- c(K_Xx %*% pred.comp)
  f_prime_hat <- c(K10_Xx %*% pred.comp)
  
  rmse <- NULL
  if(RMSE == TRUE){
    rmse = c(sqrt(mean((f_hat - f0_new)^2)), sqrt(mean((f_prime_hat - f0_new_prime)^2)))
  }
  
  f_hat_lb_pw = f_hat_ub_pw = f_hat_lb_st = f_hat_ub_st <- NULL
  f_prime_lb_pw = f_prime_ub_pw = f_prime_lb_st = f_prime_ub_st <- NULL
  
  if(CI == TRUE){
    sigma_hat <- sqrt(lambda * sum(y * pred.comp))
    
    cov <- sigma_hat^2 / (n * lambda) * 
      (K_xx - K_Xx.eigenvector %*% diag.new %*% (t(K_Xx.eigenvector)))
    
    sds <- sqrt(diag(cov))
    f_hat_lb_pw <- f_hat - qnorm(1-alpha/2) * sds
    f_hat_ub_pw <- f_hat + qnorm(1-alpha/2) * sds
    
    # Compute the L-infinity simultaneous credible band for mean
    pf_hat_sample <- mvrnorm(n = 1000, mu = rep(0, n_new), Sigma = cov)
    pr_hat_max <- apply(abs(pf_hat_sample), 1, max)
    radius <- quantile(pr_hat_max, 1 - alpha)
    f_hat_lb_st <- f_hat - radius
    f_hat_ub_st <- f_hat + radius
    
    cov_11 <- sigma_hat^2 / (n * lambda) * 
      (K11_xx - K10_Xx.eigenvector %*% diag.new %*% (t(K10_Xx.eigenvector)))
    
    sds_11 <- sqrt(diag(cov_11))
    f_prime_lb_pw <- f_prime_hat - qnorm(1-alpha/2) * sds_11
    f_prime_ub_pw <- f_prime_hat + qnorm(1-alpha/2) * sds_11
    
    # Compute the L-infinity simultaneous credible band for derivative of mean
    pf_prime_sample <- mvrnorm(n = 1000, mu = rep(0, n_new), Sigma = cov_11)
    pr_prime_max <- apply(abs(pf_prime_sample), 1, max)
    radius_prime <- quantile(pr_prime_max, 1 - alpha)
    f_prime_lb_st <- f_prime_hat - radius_prime
    f_prime_ub_st <- f_prime_hat + radius_prime
  }
  
  list(loo_error = loo_error,
       rmse = rmse,
       f_hat = f_hat,
       f_prime_hat = f_prime_hat,
       f_hat_lb_pw = f_hat_lb_pw,
       f_hat_ub_pw = f_hat_ub_pw,
       f_hat_lb_st = f_hat_lb_st,
       f_hat_ub_st = f_hat_ub_st,
       f_prime_lb_pw = f_prime_lb_pw, 
       f_prime_ub_pw = f_prime_ub_pw,
       f_prime_lb_st = f_prime_lb_st,
       f_prime_ub_st = f_prime_ub_st)}
