# This R script generates Figure 2(a) in the paper.


# remotes::install_github("cran/RandomFieldsUtils")
source("GP_Matern.R")
source("GP_SE.R")
source("GP_Sobolev.R")

par(mai=c(0.4, 0.4, 0.1, 0.1))


# Read data
dat = read.csv("CSIRO_Recons_gmsl_yr_2015.txt", header = F, sep = "")
range = max(dat$V1) - min(dat$V1)
x = dat$V1[1:130]/range
y = dat$V2[1:130]/1000
n = length(x)

# Number of grid points
n_new = 100
# Grid points
x_new = seq(1880/range, 2010/range,, n_new)

# Matern kernel
nu_seq <- seq(2, 10, 0.5)
loo_seq <- c()

for (i in 1:length(nu_seq)) {
  res_Matern <- get_GPR_matern(x, y, x_new, nu_seq[i], alpha = 0.1, LOO = TRUE)
  loo_seq <- c(loo_seq, res_Matern$loo_error)
}

# The best nu selected is 3
best_nu <- nu_seq[which.min(loo_seq)]
res_Matern <- get_GPR_matern(x, y, x_new, 2, alpha = 0.1, LOO = TRUE, CI = TRUE)

# SE kernel
res_SE <- get_GPR_SE(x, y, x_new, alpha = 0.1, LOO = TRUE, CI = TRUE)

# Sobolev kernel
res_Sobolev <- get_GPR_sobolev(x, y, x_new, alpha = 0.1, LOO = TRUE, CI = TRUE)

# Sobolev kernel has the least LOO error
res_Matern$loo_error
res_SE$loo_error
res_Sobolev$loo_error

# Plot
plot(x_new*range, res_Sobolev$f_prime_hat/range*1000,
     lwd = 3, col = "red", type = "l", ylim = c(1, 2),
     ylab = "Global Mean Sea-Level Rise (mm/yr)", xlab = "Year AD",
     xaxp = c(1880, 2010, 13))

polygon(c(rev(x_new*range), x_new*range),
        c(rev(res_Sobolev$f_prime_lb_pw/range*1000), res_Sobolev$f_prime_ub_pw/range*1000),
        col = rgb(0, 0, 0, 0.5), border = NA)

polygon(c(rev(x_new*range), x_new*range),
        c(rev(res_Sobolev$f_prime_lb_st/range*1000), res_Sobolev$f_prime_ub_st/range*1000),
        col = rgb(0, 0, 0, 0.1), border = NA)

lines(x_new*range, res_Sobolev$f_prime_hat/range*1000,
      lwd = 3, col = "black", type = "l", ylim = c(0, 3))
abline(h = 1.7, lty = 2)
text(1890, 1.75, "1.7 mm/yr")

