library(Rsolnp)
library(MASS)
library(RandomFieldsUtils)

# Matern kernel
mat <- function(x1, x2, nu){
  matrix(matern(abs(outer(x2, x1, "-")), nu, derivative = 0), ncol = length(x1))}

mat_prime <- function(x1, x2, nu){
  XX <- outer(x2, x1, "-")
  out <- matrix(matern(abs(XX), nu, derivative = 1), ncol = length(x1)) * sign(XX)
  out}

mat_2prime <- function(x1, x2, nu){
  XX <- outer(x2, x1, "-")
  out <- -matrix(matern(abs(XX), nu, derivative = 2), ncol = length(x1))
  out}

get_lambda_matern <- function(x, y, nu){
  n <- length(x)
  K_XX <- mat(x, x, nu)
  eigendecomposition <- eigen(K_XX)
  eigenvector <- eigendecomposition$vector
  eigenvalue <- eigendecomposition$values
  
  log_post_matern <- function(theta){
    
    diag.new <- eigenvalue + n * theta
    pred.comp <- c(eigenvector %*% (c(crossprod(eigenvector, y)) / diag.new))
    
    log_p <- n * log(mean(y * pred.comp)) + sum(log(diag.new))
    log_p}
  
  powell <- solnp(pars = c(0.01), fun = log_post_matern, LB = 2e-3, UB = 1)

  lambda <- powell$pars
  lambda}

get_LOO_matern <- function(x, y, nu, lambda){
  mse <- rep(0, n)
  for(i in 1:n){
    testy <- y[i]
    testx <- x[i]
    trainy <- y[-i]
    trainx <- x[-i]
    
    K_XX <- mat(trainx, trainx, nu)
    K_Xx <- mat(trainx, testx, nu)
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


get_GPR_matern <- function(x, y, x_new, nu, alpha = 0.05, LOO = FALSE, CI = FALSE, RMSE = FALSE){
  n <- length(x)
  n_new <- length(x_new)
  
  lambda <- get_lambda_matern(x, y, nu)
  
  loo_error <- NULL
  if(LOO == TRUE){
    loo_error <- get_LOO_matern(x, y, nu, lambda)
  }

  K_XX <- mat(x, x, nu)
  K_Xx <- mat(x, x_new, nu)
  K_xx <- mat(x_new, x_new, nu)
  K10_Xx <- mat_prime(x, x_new, nu)
  K11_xx <- mat_2prime(x_new, x_new, nu)

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
    
    # Compute the L-infinity simultaneous credible band for derivative of posterior mean
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
    
    # Compute the L-infinity simultaneous credible band for derivative of posterior mean
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

