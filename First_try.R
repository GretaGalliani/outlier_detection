library(MASS)
source("main.R")

mu_1 <- c(0,0)
Sigma_1 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)

data_1 <- mvrnorm(n = 5, mu_1, Sigma_1)

mu_2 <- c(2,2)
Sigma_2 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)

data_2 <- mvrnorm(n = 5, mu_1, Sigma_1)

data = rbind(data_1, data_2)



algorithm <- function(Y, S_init, sigma_init, theta_init, beta_init, xi_mu, xi_cov, Q_param, P_param, n_iter )