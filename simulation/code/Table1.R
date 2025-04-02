# remotes::install_github("cran/RandomFieldsUtils")
library(foreach)
library(doSNOW)
library(doParallel)

source("GP_Matern.R")
source("GP_SE.R")
source("GP_Sobolev.R")
source("GP_Bspline.R")

# True regression function
f0 <- function(x){sum(sqrt(2) * I^(-4) * sin(I) * cos((I - 0.5) * pi * x))}

# Derivative of true regression function
f0_prime <- function(x){sum(sqrt(2) * I^(-4) * sin(I) * (0.5 - I) * pi * sin((I - 0.5) * pi * x))}

num_sim = 100
cl = makeCluster(num_sim)
registerDoSNOW(cl)
pb <- txtProgressBar(max = num_sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


########################################
#
# n = 100
#
########################################

res_100 = foreach(i = 1:num_sim, .options.snow = opts, .combine = rbind) %dopar% {
  
  library(RandomFieldsUtils)
  library(Rsolnp)
  library(splines)
  
  I <- 1:9999
  set.seed(i)
  n <- 100
  x <- sort(runif(n))
  y0 <- sapply(x, f0) # true value at samples
  sigma0 <- sqrt(0.1) # true regression standard deviation
  y <- y0 + sigma0 * rnorm(n)
  
  n_new <- 100 # number of grid points
  x_new <- seq(0, 1, length = n_new) # grid points
  f0_new <- sapply(x_new, f0) # true value at grid points
  f0_new_prime <- sapply(x_new, f0_prime) # true derivative at grid points
  
  # Matern kernel
  nu_seq <- seq(3, 10, 0.5)
  loo_seq <- c()
  
  for (i in 1:length(nu_seq)) {
    res_Matern <- get_GPR_matern(x, y, x_new, nu_seq[i], alpha = 0.05, LOO = TRUE)
    loo_seq <- c(loo_seq, res_Matern$loo_error)
  }
  
  best_nu <- nu_seq[which.min(loo_seq)]
  res_Matern <- get_GPR_matern(x, y, x_new, best_nu, alpha = 0.05, LOO = TRUE, RMSE = TRUE)
  
  # Squared exponential kernel
  res_SE <- get_GPR_SE(x, y, x_new, alpha = 0.05, LOO = TRUE, RMSE = TRUE)
  
  # Sobolev kernel
  res_Sobolev <- get_GPR_sobolev(x, y, x_new, alpha = 0.05, LOO = TRUE, RMSE = TRUE)
  
  # B-splines
  res_Bspline <- get_Bspline(x, y, x_new, alpha = 0.05)
  
  if(res_Matern$loo_error == min(res_Matern$loo_error, res_SE$loo_error, res_Sobolev$loo_error)){
    cv_mse <- res_Matern$rmse
  }else if(res_SE$loo_error == min(res_Matern$loo_error, res_SE$loo_error, res_Sobolev$loo_error)){
    cv_mse <- res_SE$rmse
  }else{
    cv_mse <- res_Sobolev$rmse
  }
  
  return(c(res_Matern$rmse, res_SE$rmse, res_Sobolev$rmse, cv_mse, res_Bspline$rmse, res_Bspline$N))
}

save(res_100, file = "table1_data_n100.RData")

# Column 1 - Mean RMSE of estimating f0
round(matrix(colMeans(res_100[,-11]), ncol = 2, byrow = TRUE), 3)[,1]

# Column 1 - Mean RMSE of estimating f0
round(matrix(apply(res_100[,-11], 2, sd), ncol = 2, byrow = TRUE), 3)[,1]

# Column 4 - Mean RMSE of estimating derivative of f0
round(matrix(colMeans(res_100[,-11]), ncol = 2, byrow = TRUE), 2)[,2]

# Column 4 - Mean RMSE of estimating derivative of f0
round(matrix(apply(res_100[,-11], 2, sd), ncol = 2, byrow = TRUE), 3)[,2]

# Median of Bspline RMSE
round(matrix(apply(res_100[,9:10], 2, median), ncol = 2, byrow = TRUE), 3)


########################################
#
# n = 500
#
########################################

res_500 = foreach(i = 1:num_sim, .options.snow = opts, .combine = rbind) %dopar% {
  
  library(RandomFieldsUtils)
  library(Rsolnp)
  library(splines)
  
  I <- 1:9999
  set.seed(i)
  n <- 500
  x <- sort(runif(n))
  y0 <- sapply(x, f0) # true value at samples
  sigma0 <- sqrt(0.1) # true regression standard deviation
  y <- y0 + sigma0 * rnorm(n)
  
  n_new <- 100 # number of grid points
  x_new <- seq(0, 1, length = n_new) # grid points
  f0_new <- sapply(x_new, f0) # true value at grid points
  f0_new_prime <- sapply(x_new, f0_prime) # true derivative at grid points
  
  # Matern kernel
  nu_seq <- seq(3, 10, 0.5)
  loo_seq <- c()
  
  for (i in 1:length(nu_seq)) {
    res_Matern <- get_GPR_matern(x, y, x_new, nu_seq[i], alpha = 0.05, LOO = TRUE)
    loo_seq <- c(loo_seq, res_Matern$loo_error)
  }
  
  best_nu <- nu_seq[which.min(loo_seq)]
  res_Matern <- get_GPR_matern(x, y, x_new, best_nu, alpha = 0.05, LOO = TRUE, RMSE = TRUE)
  
  # Squared exponential kernel
  res_SE <- get_GPR_SE(x, y, x_new, alpha = 0.05, LOO = TRUE, RMSE = TRUE)
  
  # Sobolev kernel
  res_Sobolev <- get_GPR_sobolev(x, y, x_new, alpha = 0.05, LOO = TRUE, RMSE = TRUE)
  
  # B-splines
  res_Bspline <- get_Bspline(x, y, x_new, alpha = 0.05)
  
  if(res_Matern$loo_error == min(res_Matern$loo_error, res_SE$loo_error, res_Sobolev$loo_error)){
    cv_mse <- res_Matern$rmse
  }else if(res_SE$loo_error == min(res_Matern$loo_error, res_SE$loo_error, res_Sobolev$loo_error)){
    cv_mse <- res_SE$rmse
  }else{
    cv_mse <- res_Sobolev$rmse
  }
  
  return(c(res_Matern$rmse, res_SE$rmse, res_Sobolev$rmse, cv_mse, res_Bspline$rmse, res_Bspline$N))
}

save(res_500, file = "table1_data_n500.RData")

# Column 2 - Mean RMSE of estimating f0
round(matrix(colMeans(res_500[,-11]), ncol = 2, byrow = TRUE), 3)[,1]

# Column 2 - Mean RMSE of estimating f0
round(matrix(apply(res_500[,-11], 2, sd), ncol = 2, byrow = TRUE), 4)[,1]

# Column 5 - Mean RMSE of estimating derivative of f0
round(matrix(colMeans(res_500[,-11]), ncol = 2, byrow = TRUE), 2)[,2]

# Column 5 - Mean RMSE of estimating derivative of f0
round(matrix(apply(res_500[,-11], 2, sd), ncol = 2, byrow = TRUE), 3)[,2]

# Median of Bspline RMSE
round(matrix(apply(res_500[,9:10], 2, median), ncol = 2, byrow = TRUE), 3)


########################################
#
# n = 1000
#
########################################

res_1000 = foreach(i = 1:num_sim, .options.snow = opts, .combine = rbind) %dopar% {
  
  library(RandomFieldsUtils)
  library(Rsolnp)
  library(splines)
  
  I <- 1:9999
  set.seed(i)
  n <- 1000
  x <- sort(runif(n))
  y0 <- sapply(x, f0) # true value at samples
  sigma0 <- sqrt(0.1) # true regression standard deviation
  y <- y0 + sigma0 * rnorm(n)
  
  n_new <- 100 # number of grid points
  x_new <- seq(0, 1, length = n_new) # grid points
  f0_new <- sapply(x_new, f0) # true value at grid points
  f0_new_prime <- sapply(x_new, f0_prime) # true derivative at grid points
  
  # Matern kernel
  nu_seq <- seq(3, 10, 0.5)
  loo_seq <- c()
  
  for (i in 1:length(nu_seq)) {
    res_Matern <- get_GPR_matern(x, y, x_new, nu_seq[i], alpha = 0.05, LOO = TRUE)
    loo_seq <- c(loo_seq, res_Matern$loo_error)
  }

  best_nu <- nu_seq[which.min(loo_seq)]
  res_Matern <- get_GPR_matern(x, y, x_new, best_nu, alpha = 0.05, LOO = TRUE, RMSE = TRUE)
  
  # Squared exponential kernel
  res_SE <- get_GPR_SE(x, y, x_new, alpha = 0.05, LOO = TRUE, RMSE = TRUE)
  
  # Sobolev kernel
  res_Sobolev <- get_GPR_sobolev(x, y, x_new, alpha = 0.05, LOO = TRUE, RMSE = TRUE)
  
  # B-splines
  res_Bspline <- get_Bspline(x, y, x_new, alpha = 0.05)
  
  if(res_Matern$loo_error == min(res_Matern$loo_error, res_SE$loo_error, res_Sobolev$loo_error)){
    cv_mse <- res_Matern$rmse
  }else if(res_SE$loo_error == min(res_Matern$loo_error, res_SE$loo_error, res_Sobolev$loo_error)){
    cv_mse <- res_SE$rmse
  }else{
    cv_mse <- res_Sobolev$rmse
  }
  
  return(c(res_Matern$rmse, res_SE$rmse, res_Sobolev$rmse, cv_mse, res_Bspline$rmse, res_Bspline$N))
}

save(res_1000, file = "table1_data_n1000.RData")

# Column 3 - Mean RMSE of estimating f0
round(matrix(colMeans(res_1000[,-11]), ncol = 2, byrow = TRUE), 3)[,1]

# Column 3 - Mean RMSE of estimating f0
round(matrix(apply(res_1000[,-11], 2, sd), ncol = 2, byrow = TRUE), 4)[,1]

# Column 6 - Mean RMSE of estimating derivative of f0
round(matrix(colMeans(res_1000[,-11]), ncol = 2, byrow = TRUE), 2)[,2]

# Column 6 - Mean RMSE of estimating derivative of f0
round(matrix(apply(res_1000[,-11], 2, sd), ncol = 2, byrow = TRUE), 3)[,2]

# Median of Bspline RMSE
round(matrix(apply(res_1000[,9:10], 2, median), ncol = 2, byrow = TRUE), 3)
