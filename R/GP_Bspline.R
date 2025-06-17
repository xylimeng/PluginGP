# Implements the B-spline method proposed by Yoo and Ghosal (2016), 
# which serves as a benchmark method in the associated paper.
# The implementation is based on code provided by the original authors.
#
# Main function: get_Bspline()
#
# Inputs:
# - x: A numeric vector of design points.
# - Y: A numeric vector of observed responses corresponding to x.
# - x_new: A numeric vector of evaluation points for prediction.
# - Alpha: Significance level for credible bands (default is 0.05).
# - CI: Logical flag indicating whether to include credible bands in the output (default is FALSE).
#
# Output (a list containing the following elements):
# - N: Number of interior knots selected via cross-validation.
# - rmse: A numeric vector with RMSEs for the estimated regression function and its first derivative.
# - f_hat: Estimated regression function values at x_new.
# - f_prime_hat: Estimated first derivative values at x_new.
# - f_hat_lb_st: Lower bound of the simultaneous L-infinity credible band for the regression function.
# - f_hat_ub_st: Upper bound of the simultaneous L-infinity credible band for the regression function.
# - f_prime_hat_lb_st: Lower bound of the simultaneous L-infinity credible band for the first derivative.
# - f_prime_hat_ub_st: Upper bound of the simultaneous L-infinity credible band for the first derivative.

library(splines)

#calculate sigma2hat for empirical Bayes
sigma2emp <- function(y, cmiddle, B, eta){
  ydiff <- y - B %*% eta
  sigma2 = (crossprod(ydiff) - crossprod(forwardsolve(cmiddle, crossprod(B, ydiff)))) / n
  return(sigma2)
}

#returns posterior mean and variance
pfmeanf <- function(y, cmiddle, B, b, Omegainv, eta){
  ans = forwardsolve(cmiddle, crossprod(B, y) + Omegainv %*% eta)
  pmean = crossprod(t(b), backsolve(t(cmiddle), ans))  #(3.4), (3.5), r = 0
  return(pmean)
}

pfvarf <- function(cmiddle, b){ 
  pSigma = crossprod(forwardsolve(cmiddle, t(b)))  #(3.6), r = 0
  return(pSigma)
}

#determine the best number of splines
optimal <- function(y, x, a, Omegainv, eta){
  mse <- rep(0, n)
  for(k in 1:n){
    testy <- y[k]
    testx <- x[k]
    trainy <- y[-k]
    trainx <- x[-k]
    
    q <- 4
    B <- bs(trainx, knots = a, degree = q - 1, intercept = TRUE, Boundary.knots = c(0, 1))
    BB = crossprod(B)
    V = BB + Omegainv
    cM <- t(chol(V))
    bnewx <- predict(B, testx)
    
    pfmean = pfmeanf(trainy, cmiddle = cM, B = B, Omegainv = Omegainv, eta = eta, b = bnewx)
    mse[k] <- (pfmean - testy) ^ 2  #1st col Bayes
  }
  
  return(mean(mse))
}


bsprime <- function(x, derivs, df = NULL, knots = NULL, degree, intercept = TRUE, Boundary.knots = range(x))
{
  nx <- names(x)
  x <- as.vector(x)
  Boundary.knots <- sort(Boundary.knots)
  ord <- 1L + (degree <- as.integer(degree))
  if(!is.null(df) && is.null(knots)){
    nIknots <- df - ord + (1L - intercept)
    knots <- if(nIknots > 0L) {
      knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
                         2L)[-c(1L, nIknots + 2L)]
      stats::quantile(x, knots)
    }
  }
  Aknots <- sort(c(rep(Boundary.knots, ord), knots))
  
  basis <- splineDesign(Aknots, x, ord, derivs)
  n.col <- ncol(basis)
  dimnames(basis) <- list(nx, 1L:n.col)
  basis
}

get_Bspline <- function(x, y, x_new, alpha = 0.05, CI = FALSE){
  n <- length(x)
  n_new <- length(x_new)
  
  q = 4  #cubic splines
  degree = q - 1
  Nmax = 10
  mseBF = rep(0, Nmax)
  tau2 <- 1
  tau2inv <- 1 / tau2
  
  for(j in 1:Nmax){
    N <- j  #optimal in terms of mse
    J <- N + q
    
    eta <- rep(0, J)
    Omegainv <- tau2inv * diag(rep(1, J))
    
    a = seq(from = 0, to = 1, length.out = N + 2)[-c(1, N + 2)]
    mseBF[j] <- optimal(y = y, x = x, a = a, Omegainv = Omegainv, eta = eta)
  }
  
  #-----------------------
  
  N = which(mseBF == min(mseBF))
  J <- N + q

  eta <- rep(0, J)
  Omegainv <- (1 / tau2) * diag(rep(1, J))  #sigma^-2 * Omega^-1
  
  a = seq(from = 0, to = 1, length.out = N + 2)[-c(1, N + 2)]  #intex of knots from 0 to N+1
  
  B <- bs(x, knots = a, degree = degree, intercept = TRUE, Boundary.knots = c(0,1))
  BB <- crossprod(B)
  V = BB + Omegainv
  cmiddle <- t(chol(V))
  
  bnew <- predict(B, x_new)
  bnew_prime = bsprime(x_new, rep(1, length(x_new)), knots = a, degree = degree, intercept = TRUE, Boundary.knots = c(0,1))
  
  f_hat = pfmeanf(y = y, cmiddle = cmiddle, B = B, Omegainv = Omegainv, eta = eta, b = bnew)
  f_prime_hat = pfmeanf(y = y, cmiddle = cmiddle, B = B, Omegainv = Omegainv, eta = eta, b = bnew_prime)
  
  f_hat_lb_st = f_hat_ub_st = f_prime_lb_st = f_prime_ub_st <- NULL
  
  if(CI == TRUE){
    pfvar <- bnew %*% solve(V, t(bnew))
    pfhat <- mvrnorm(n = 1000, mu = rep(0, n_new), Sigma = pfvar)
    prmax = apply(abs(pfhat), 1, max)
    prmaxquan = quantile(prmax, 1 - alpha) 
    psigma2 <- as.numeric(sigma2emp(y = y, cmiddle = cmiddle, B = B, eta = eta))
    radiusmax <- sqrt(psigma2) * prmaxquan / 2
    f_hat_lb_st <- f_hat - radiusmax
    f_hat_ub_st <- f_hat + radiusmax
    
    pfvar_prime <- bnew_prime %*% solve(V, t(bnew_prime))
    pfhat_prime <- mvrnorm(n = 1000, mu = rep(0, n_new), Sigma = pfvar_prime)
    prmax_prime = apply(abs(pfhat_prime), 1, max)
    prmaxquan_prime = quantile(prmax_prime, 1 - alpha) 
    radiusmax_prime <- sqrt(psigma2) * prmaxquan_prime / 2
    f_prime_lb_st <- f_prime_hat - radiusmax_prime
    f_prime_ub_st <- f_prime_hat + radiusmax_prime

  }
  

  list(N = N,
       rmse = c(sqrt(mean((f_hat - f0_new)^2)), sqrt(mean((f_prime_hat - f0_new_prime)^2))),
       f_hat = f_hat,
       f_prime_hat = f_prime_hat,
       f_hat_lb_st = f_hat_lb_st,
       f_hat_ub_st = f_hat_ub_st,
       f_prime_lb_st = f_prime_lb_st,
       f_prime_ub_st = f_prime_ub_st)}

