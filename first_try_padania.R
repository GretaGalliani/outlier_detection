library(MASS)

mu_1 <- c(0,0)
Sigma_1 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)

data_1 <- mvrnorm(n = 5, mu_1, Sigma_1)

mu_2 <- c(2,2)
Sigma_2 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)

data_2 <- mvrnorm(n = 5, mu_1, Sigma_1)

data = rbind(data_1, data_2)

S_init = c(1,2,3,4,5,6,7,8,9,10)

# NB NON INIZIALIZZARE BETA A 1 ASSOLUTAMENTE
beta_init <- 0.9999999
sigma_init <- 0.5
theta_init <- 0.5

xi_mu <- list()
xi_cov <- list()
mu = c(0,0)
cov = matrix(c(1,0,0,1), nrow = 2, ncol = 2)
for (i in 1:10){
  xi_mu <- append(xi_mu, list(mu))
  xi_cov <- append(xi_cov, list(cov))
}

Q_param <- list("nu_0"= 3, "mu_0" = c(0,0), "lambda_0" = matrix(c(1,0,0,1), nrow = 2, ncol = 2), "k_0" = 1)
P_param <- list("nu_0"= 3, "mu_0" = c(0,0), "lambda_0" = matrix(c(1,0,0,1), nrow = 2, ncol = 2), "k_0" = 1)


source("main_letresicilie.R")
result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, xi_mu, xi_cov, Q_param, P_param, 1)
warnings()