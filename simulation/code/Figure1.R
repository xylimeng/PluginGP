# remotes::install_github("cran/RandomFieldsUtils")
library(foreach)
library(doSNOW)
library(doParallel)

par(mai=c(0.4, 0.4, 0.1, 0.1))

source("GP_Matern.R")
source("GP_SE.R")
source("GP_Sobolev.R")
source("GP_Bspline.R")

# True regression function
f0 <- function(x){sum(sqrt(2) * I^(-4) * sin(I) * cos((I - 0.5) * pi * x))}

# Derivative of true regression function
f0_prime <- function(x){sum(sqrt(2) * I^(-4) * sin(I) * (0.5 - I) * pi * sin((I - 0.5) * pi * x))}

plot_function <- function(res) {
  plot(x_new, res$f_hat, type = "l", lwd = 3, xlab = "", ylab = "", ylim = c(-0.1, 1.5))
  lines(x_new, f0_new, lty = 3, lwd = 3)
  lines(x_new, res$f_hat_lb_st, lty = 5, lwd = 3)
  lines(x_new, res$f_hat_ub_st, lty = 5, lwd = 3)
}

plot_derivative <- function(res) {
  plot(x_new, res$f_prime_hat, type = "l", lwd = 3, xlab = "", ylab = "", ylim = c(-2.5, 1.2))
  lines(x_new, f0_new_prime, lty = 3, lwd = 3)
  lines(x_new, res$f_prime_lb_st, lty = 5, lwd = 3)
  lines(x_new, res$f_prime_ub_st, lty = 5, lwd = 3)
}

########################################
#
# Example 1
#
########################################

I <- 1:9999
set.seed(1000)
n <- 1000
x <- sort(runif(n))
y0 <- sapply(x, f0) # true value at samples
sigma0 <- sqrt(0.1) # true regression standard deviation
y <- y0 + sigma0 * rnorm(n)  #generate some data

n_new <- 100 # number of grid points
x_new <- seq(0, 1, length = n_new) # grid points
f0_new <- sapply(x_new, f0) # true value at grid points
f0_new_prime <- sapply(x_new, f0_prime) # true derivative at grid points

# Matern kernel
nu_seq <- seq(3, 10, 0.5)
loo_seq <- c()

num_sim <- length(nu_seq)
cl <- makeCluster(num_sim)
registerDoSNOW(cl)
pb <- txtProgressBar(max = num_sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

loo_seq <- foreach(i = 1:num_sim, .options.snow = opts, .combine = rbind) %dopar% {

  library(RandomFieldsUtils)
  library(Rsolnp)
  library(splines)
  
  res_Matern <- get_GPR_matern(x, y, x_new, nu_seq[i], alpha = 0.05, LOO = TRUE)
  
  return(res_Matern$loo_error)
}

# The best nu selected is 3
best_nu <- nu_seq[which.min(loo_seq)]
res_Matern <- get_GPR_matern(x, y, x_new, best_nu, alpha = 0.05, CI = TRUE)
plot_function(res_Matern)
plot_derivative(res_Matern)

# SE kernel
res_SE <- get_GPR_SE(x, y, x_new, alpha = 0.05, CI = TRUE)
plot_function(res_SE)
plot_derivative(res_SE)

# Sobolev kernel
res_Sobolev <- get_GPR_sobolev(x, y, x_new, alpha = 0.05, CI = TRUE)
plot_function(res_Sobolev)
plot_derivative(res_Sobolev)

# B-splines
res_Bspline <- get_Bspline(x, y, x_new, alpha = 0.05, CI = TRUE)
plot_function(res_Bspline)
plot_derivative(res_Bspline)


########################################
#
# Example 2
#
########################################

I <- 1:9999
set.seed(20)
n <- 1000
x <- sort(runif(n))
y0 <- sapply(x, f0) # true value at samples
sigma0 <- sqrt(0.1) # true regression standard deviation
y <- y0 + sigma0 * rnorm(n)  #generate some data

n_new <- 100 # number of grid points
x_new <- seq(0, 1, length = n_new) # grid points
f0_new <- sapply(x_new, f0) # true value at grid points
f0_new_prime <- sapply(x_new, f0_prime) # true derivative at grid points

# Matern kernel
nu_seq <- seq(3, 10, 0.5)
loo_seq <- c()

num_sim <- length(nu_seq)
cl <- makeCluster(num_sim)
registerDoSNOW(cl)
pb <- txtProgressBar(max = num_sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

loo_seq <- foreach(i = 1:num_sim, .options.snow = opts, .combine = rbind) %dopar% {
  .libPaths("~/R/x86_64-redhat-linux-gnu-library/4.0")
  library(RandomFieldsUtils)
  library(Rsolnp)
  library(splines)
  
  res_Matern <- get_GPR_matern(x, y, x_new, nu = nu_seq[i], alpha = 0.05, LOO = TRUE)
  
  return(res_Matern$loo_error)
}

# The best nu selected is 3
best_nu <- nu_seq[which.min(loo_seq)]
res_Matern <- get_GPR_matern(x, y, x_new, nu = best_nu, alpha = 0.05, CI = TRUE)
plot_function(res_Matern)
plot_derivative(res_Matern)

# SE kernel
res_SE <- get_GPR_SE(x, y, x_new, alpha = 0.05, CI = TRUE)
plot_function(res_SE)
plot_derivative(res_SE)

# Sobolev kernel
res_Sobolev <- get_GPR_sobolev(x, y, x_new, alpha = 0.05, CI = TRUE)
plot_function(res_Sobolev)
plot_derivative(res_Sobolev)

# B-splines
res_Bspline <- get_Bspline(x, y, x_new, alpha = 0.05, CI = TRUE)
plot_function(res_Bspline)
plot_derivative(res_Bspline)
