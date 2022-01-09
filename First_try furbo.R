library(MASS)
library(RColorBrewer)
source("main.R")

set.seed(357846)

Q_param = list()
P_param = list()
d=2


Q_param$k_0 = 0.1
Q_param$mu_0 = c(0,0)
Q_param$nu_0 = 3 # it must be > (p-1)
Q_param$lambda_0 = matrix(c(3,0,0,3), nrow = 2, ncol = 2)

# computation of degrees of freedom
df = Q_param$nu_0-d+1


### GROUP 1

# sigma of group 1
# xi_sigma_1 = LaplacesDemon::rinvwishart(Q_param$nu_0, as.inverse(Q_param$lambda_0))
xi_sigma_1 = LaplacesDemon::rinvwishart(Q_param$nu_0, Q_param$lambda_0)

# Sampling of mu 
xi_mu_1 = MASS::mvrnorm(1, Q_param$mu_0, xi_sigma_1/Q_param$k_0)


# Sampling 10 data from group 1
data_1 <- mvrnorm(n = 10, xi_mu_1, xi_sigma_1)


### GROUP 2
# sigma of group 2
# xi_sigma_2 = LaplacesDemon::rinvwishart(Q_param$nu_0, as.inverse(Q_param$lambda_0))
xi_sigma_2 = LaplacesDemon::rinvwishart(Q_param$nu_0, Q_param$lambda_0)

# Sampling of mu 2
xi_mu_2 = MASS::mvrnorm(1, Q_param$mu_0, xi_sigma_2/Q_param$k_0)


# Sampling 5 data from group 2
data_2 <- mvrnorm(n = 10, xi_mu_2, xi_sigma_2)




# SAMPLING FROM THE CONTAMINATED PART
P_param$k_0 = 0.01
P_param$mu_0 = c(0,0)
P_param$nu_0 = 3 # it must be > (p-1)
P_param$lambda_0 = matrix(c(3,0,0,3), nrow = 2, ncol = 2)

# computation of degrees of freedom
df = P_param$nu_0-d+1

# sigma of group 0
# xi_sigma_0 = LaplacesDemon::rinvwishart(P_param$nu_0, as.inverse(P_param$lambda_0))
xi_sigma_0 = LaplacesDemon::rinvwishart(P_param$nu_0, P_param$lambda_0)

# Sampling of mu 0
xi_mu_0 = MASS::mvrnorm(1, P_param$mu_0, xi_sigma_0/P_param$k_0)

# Sampling a contaminated data

data_3 <- mvrnorm(n = 1, xi_mu_0, xi_sigma_0)
outlier1 <- NULL
outlier1 <- cbind(data_3[1],data_3[2])

# Repeating other two times
# sigma of group 0
# xi_sigma_0 = LaplacesDemon::rinvwishart(P_param$nu_0, as.inverse(P_param$lambda_0))
xi_sigma_0 = LaplacesDemon::rinvwishart(P_param$nu_0, P_param$lambda_0)

# Sampling of mu 0
xi_mu_0 = MASS::mvrnorm(1, P_param$mu_0, xi_sigma_0/P_param$k_0)

# Sampling a contaminated data

data_3 <- mvrnorm(n = 1, xi_mu_0, xi_sigma_0)
outlier2 <- NULL
outlier2 <- cbind(data_3[1],data_3[2])

# sigma of group 0
# xi_sigma_0 = LaplacesDemon::rinvwishart(P_param$nu_0, as.inverse(P_param$lambda_0))
xi_sigma_0 = LaplacesDemon::rinvwishart(P_param$nu_0, P_param$lambda_0)

# Sampling of mu 0
xi_mu_0 = MASS::mvrnorm(1, P_param$mu_0, xi_sigma_0/P_param$k_0)

# Sampling a contaminated data

data_3 <- mvrnorm(n = 1, xi_mu_0, xi_sigma_0)
outlier3 <- NULL
outlier3 <- cbind(data_3[1],data_3[2])

data = rbind(data_1, data_2, outlier1, outlier2, outlier3)


# Plotting the final data
pal = brewer.pal(n = 5, name = "Set2")
col_real = c(rep(pal[1],10), rep(pal[2],10), pal[3], pal[4], pal[5])
plot(data, col = col_real, pch = 19, main = "Real data")

S_init = rep(1,23)
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

# for (i in 1:10)
# {
#   xi_mu <- append(xi_mu, list(xi_mu_0))
#   xi_cov <- append(xi_cov, list(xi_sigma_0))
# }

init_mu <- colMeans(data)
init_var <- cov(data)

for (i in 1:dim(data)[1]){
  xi_mu <- append(xi_mu, list(init_mu))
  xi_cov <- append(xi_cov, list(init_var))
}

# S_init = c(rep(1,10), rep(2,10), rep(0,3))

source("main.R")
result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, beta_param, sigma_param, theta_param, xi_mu, xi_cov, Q_param, P_param, 15000, 1000, 1)

# result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, beta_param, sigma_param, theta_param, xi_mu, xi_cov, Q_param, P_param, 10, 0, 1)



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
aux = result$S

# These functions needs to have indexes of the groups >=1
for (i in 1:dim(aux)[1]){
  for (j in 1:dim(aux)[2]){
    if (aux[i,j]==0){
      aux[i,j] = max(aux[i,]) + 1
    }
  }
}

# For a sample of clusterings of the same objects the proportion of clusterings 
# in which observation $i$ and $j$ are together in a cluster is computed and 
# a matrix containing all proportions is given out. 
psm <- comp.psm(aux)

# finds the clustering that minimizes the posterior expectation of Binders loss function
min_bind <-  minbinder(psm, cls.draw = NULL, method = c("avg", "comp", "draws", 
                                           "laugreen","all"), max.k = NULL, include.lg = FALSE, 
          start.cl = NULL, tol = 0.001)

par(mfrow=c(1,3)) 
# real cluster
plot(data, col = col_real, pch = 19, main = "Real data")


# best cluster according to binder loss 
pal = brewer.pal(n = max(min_bind$cl), name = "Set2")
col_min_bind = rep(0,dim(data)[1])

for (i in 1:length(col_min_bind))
{
  col_min_bind[i] = pal[min_bind$cl[i]]
}


plot(data, col=col_min_bind, pch = 19, main = "Partition minimizing Binder Loss")



# IMPLEMENTING MIN VARIATION OF INFORMATION
# devtools::install_github("sarawade/mcclust.ext")
library(mcclust.ext)


# finds the clustering that minimizes  the lower bound to the posterior expected Variation of Information from Jensen's Inequality
min_vi <- minVI(psm, cls.draw=NULL, method=c("avg","comp","draws","greedy","all"), 
      max.k=NULL, include.greedy=FALSE, start.cl=NULL, maxiter=NULL,
      l=NULL, suppress.comment=TRUE)

pal = brewer.pal(n = max(min_bind$cl), name = "Set2")
col_min_vi = rep(0,dim(data)[1])

for (i in 1:length(col_min_bind))
{
  col_min_vi[i] = pal[min_vi$cl[i]]
}


plot(data, col=col_min_vi, pch = 19, main = "Partition minimizing VI Loss")
