library(MASS)
source("main.R")

set.seed(21081997)

Q_param = list()
P_param = list()
d=2


Q_param$k_0 = 2
Q_param$mu_0 = c(0,0)
Q_param$nu_0 = 3 # it must be > (p-1)
Q_param$lambda_0 = matrix(c(1,0,0,1), nrow = 2, ncol = 2)

# computation of degrees of freedom
df = Q_param$nu_0-d+1


### GROUP 1

# sigma of group 1
xi_sigma_1 = LaplacesDemon::rinvwishart(Q_param$nu_0, as.inverse(Q_param$lambda_0))

# Sampling of mu 
xi_mu_1 = MASS::mvrnorm(1, Q_param$mu_0, xi_sigma_1/Q_param$k_0)


# Sampling 5 data from group 1
data_1 <- mvrnorm(n = 5, xi_mu_1, xi_sigma_1)


### GROUP 2
# sigma of group 2
xi_sigma_2 = LaplacesDemon::rinvwishart(Q_param$nu_0, as.inverse(Q_param$lambda_0))

# Sampling of mu 2
xi_mu_2 = MASS::mvrnorm(1, Q_param$mu_0, xi_sigma_2/Q_param$k_0)


# Sampling 5 data from group 2
data_2 <- mvrnorm(n = 5, xi_mu_2, xi_sigma_2)




# SAMPLING FROM THE CONTAMINATED PART
P_param$k_0 = 200
P_param$mu_0 = c(100,100)
P_param$nu_0 = 300 # it must be > (p-1)
P_param$lambda_0 = matrix(c(100,10,10,500), nrow = 2, ncol = 2)

# computation of degrees of freedom
df = P_param$nu_0-d+1

# sigma of group 0
xi_sigma_0 = LaplacesDemon::rinvwishart(P_param$nu_0, as.inverse(P_param$lambda_0))

# Sampling of mu 0
xi_mu_0 = MASS::mvrnorm(1, P_param$mu_0, xi_sigma_0/P_param$k_0)

# Sampling a contaminated data

data_3 <- mvrnorm(n = 1, xi_mu_0, xi_sigma_0)

data = rbind(data_1, data_2, data_3)


plot(data[-11,])

S_init = c(1,2,3,4,5,6,7,8,9,10,11)
beta_init <- 0.5
sigma_init <- 0.5
theta_init <- 0.5

a_beta = 1
b_beta = 1
a_sigma = 1
b_sigma = 1

xi_mu <- list()
xi_cov <- list()


for (i in 1:5){
  xi_mu <- append(xi_mu, list(xi_mu_1))
  xi_cov <- append(xi_cov, list(xi_sigma_1))
}

for (i in 6:10){
  xi_mu <- append(xi_mu, list(xi_mu_2))
  xi_cov <- append(xi_cov, list(xi_sigma_2))
}


xi_mu <- append(xi_mu, list(xi_mu_0))
xi_cov <- append(xi_cov, list(xi_sigma_0))


source("main.R")
result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, a_beta, b_beta, a_sigma, b_sigma, xi_mu, xi_cov, Q_param, P_param, 3)
