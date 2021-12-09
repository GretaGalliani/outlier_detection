library(MASS)

mu_1 <- c(0,0)
Sigma_1 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)

data_1 <- mvrnorm(n = 5, mu_1, Sigma_1)

mu_2 <- c(2,2)
Sigma_2 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)

data_2 <- mvrnorm(n = 5, mu_1, Sigma_1)

data = rbind(data_1, data_2)

