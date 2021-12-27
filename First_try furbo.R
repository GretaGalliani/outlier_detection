library(MASS)
source("main.R")

set.seed(24091998)

Q_param = list()
P_param = list()
d=2


Q_param$k_0 = 0.01
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
outlier <- NULL
outlier <- cbind(data_3[1],data_3[2])
data = rbind(data_1, data_2, outlier)


plot(data[-11,])

S_init = c(1,2,3,4,5,6,7,8,9,10,3)
beta_init <- 0.5
sigma_init <- 0.2
theta_init <- 0.2

beta_param = list()
sigma_param = list()
theta_param = list()
beta_param$a = 1
beta_param$b = 1
sigma_param$a = 1
sigma_param$b = 1
theta_param$a = 1
theta_param$b = 1

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
result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, beta_param, sigma_param, theta_param, xi_mu, xi_cov, Q_param, P_param, 15000, 5000, 10)

result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, beta_param, sigma_param, theta_param, xi_mu, xi_cov, Q_param, P_param, 10, 0, 1)



tail(result$sigma, 20)
plot(result$sigma, type='l')

tail(result$theta, 20)
plot(result$theta, type='l')

result$S[1:30,]
plot(data, col=result$S[150,])

max <- c()
for (i in 1:1000){
  max <- c(max, max(result$S[i,]))
}

mean(max) # mean number of clusters by the algorithm 


# IMPLEMENTING MIN BINDER LOSS
library(mcclust)

# These functions needs to have indexes of the groups >=1
aux = result$S + 1

# For a sample of clusterings of the same objects the proportion of clusterings 
# in which observation $i$ and $j$ are together in a cluster is computed and 
# a matrix containing all proportions is given out. 
psm <- comp.psm(aux)

# finds the clustering that minimizes the posterior expectation of Binders loss function
min_bind <-  minbinder(psm, cls.draw = NULL, method = c("avg", "comp", "draws", 
                                           "laugreen","all"), max.k = NULL, include.lg = FALSE, 
          start.cl = NULL, tol = 0.001)

par(mfrow=c(1,2)) 

# best cluster according to binder loss (without outlier)
plot(data[-11,], col=min_bind$cl, pch = 19)

# real cluster
real <- c(1,1,1,1,1,2,2,2,2,2)
plot(data[-11,], col=real, pch = 19)

# IMPLEMENTING MIN VARIATION OF INFORMATION
# devtools::install_github("sarawade/mcclust.ext")
library(mcclust.ext)

aux <- result$S[,-11]
psm2 <- comp.psm(aux)

# finds the clustering that minimizes  the lower bound to the posterior expected Variation of Information from Jensen's Inequality
min_vi <- minVI(psm2, cls.draw=NULL, method=c("avg","comp","draws","greedy","all"), 
      max.k=NULL, include.greedy=FALSE, start.cl=NULL, maxiter=NULL,
      l=NULL, suppress.comment=TRUE)

par(mfrow=c(1,2))

# best cluster according to iv loss (without outlier)
plot(data[-11,], col=min_vi$cl, pch = 19)

# real cluster
plot(data[-11,], col=real, pch = 19)
