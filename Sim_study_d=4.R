set.seed(04021997)
library(MASS)
library(RColorBrewer)

#### SAMPLING FROM THE DENSITY ####
# BASE MEASURE
d = 4 #dimension
mean_a = rep(-3,d) #vector of means
sigma_b = diag(d) #sigma
m = 90 #m number of samples (outliers excluded)
m1 = rbinom(1, size=m, prob = 0.5) #number of samples coming from the first gaussian 
m2 = m-m1 # from the second one

val1 = mvrnorm (m1,mu = mean_a, Sigma = sigma_b) #samples from first multivariate function
val2 = mvrnorm (m2,mu = -mean_a, Sigma = sigma_b) #samples from second multivariate function
allval = rbind(val1,val2) #combine

# CONTAMINATED MEASURE
s=10 #number of outliers
# sampling
i=0 
while(i<s){ # cycle to find s outliers
  value=mvrnorm(1,mu= rep(0,d),Sigma= 3^2*diag(d)) # sampling from a multivariate normal distribution
  module = norm(as.matrix(value), type="2")
  chi=qchisq(0.9, df = d)
  # If we are sampling from the over-disperse truncated Gaussian distribution,
  # We keep only points verifying this condition
  if(module^2>9*chi) 
  {
    i=i+1
    allval = rbind(allval,value)
  }
}

# Constructing the dataset
rownames(allval)=NULL
allval
data <- allval

# Plotting the data with the original groups
pal = brewer.pal(n = 9, name = "Set1")
col_real = c(rep(pal[1], m1), rep(pal[2], m-m1))
for (i in 1:s){
  col_real = c(col_real, pal[(i+2)%%9])
}

pairs(data, col = col_real, pch = 19)


#### INITIALIZATION ####
Q_param = list()
P_param = list()

d = dim(data)[2]

Q_param$k_0 = 1
Q_param$mu_0 = rep(0,d)
Q_param$nu_0 = d + 3 # it must be > (p-1)
Q_param$lambda_0 = diag(diag(cov(data)))

P_param$k_0 = 0.25
P_param$mu_0 = rep(0,d)
P_param$nu_0 = d + 3 # it must be > (p-1)
P_param$lambda_0 = diag(diag(cov(data)))

n = dim(data)[1]
S_init = rep(1,n)

beta_init <- 0.5
sigma_init <- 0.5
theta_init <- 1

beta_param = list()
sigma_param = list()
theta_param = list()
beta_param$a = 1
beta_param$b = 1
sigma_param$a = 1
sigma_param$b = 1
theta_param$a = 2
theta_param$b = 1

xi_mu <- list()
xi_cov <- list()

init_mu <- colMeans(data)
init_var <- cov(data)

for (i in 1:n){
  xi_mu <- append(xi_mu, list(init_mu))
  xi_cov <- append(xi_cov, list(init_var))
}

#### RUNNING THE ALGORITHM ####
source("main.R")
result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, beta_param, sigma_param, theta_param, xi_mu, xi_cov, Q_param, P_param, 15000, 10000, 1)

#### PARAMETER ANALYSIS ####
x11()
par(mfrow=c(1,3))
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(result$sigma,xlab="iteration",ylab=expression(sigma), type = 'l')  
plot(result$theta,xlab="iteration",ylab=expression(theta), type = 'l')  
plot(result$beta,xlab="iteration",ylab=expression(beta), type = 'l')  

dev.off()

## marginal traceplots
x11()
par(mfrow=c(1,3))
plot(ts(result$sigma),xlab="iteration",ylab=expression(sigma))
plot(ts(result$theta),xlab="iteration",ylab=expression(theta))
plot(ts(result$beta),xlab="iteration",ylab=expression(beta))


library(coda)

# Plot of AUTOCORRELATION
par(mfrow=c(1,3))
tmp1 <- acf(result$sigma, main='Autocorrelation of sigma')
tmp2 <- acf(result$theta, main='Autocorrelation of theta')
tmp3 <- acf(result$beta, main='Autocorrelation of beta')

# ESS
effectiveSize(result$sigma)
effectiveSize(result$theta)
effectiveSize(result$beta)

#### CLUSTER ANALYSIS ####
max <- c()
for (i in 1:1000){
  max <- c(max, max(result$S[i,]))
}

mean(max) # mean number of clusters by the algorithm 


# WE COUNT THE NUMBER OF SINGLETONS
source("auxiliary_functions/auxiliary_functions.R")

n_singletons <- c()
for (i in 1:dim(result$S)[1]){
  n_singletons <- c(n_singletons, m1(result$S[i,]))
}

n_singletons

# IMPLEMENTING MIN BINDER LOSS
library(mcclust)

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


par(mfrow=c(1,2)) 

# best cluster according to binder loss (without outlier)
pairs(data, col=min_bind$cl, pch = 19)

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
pairs(data, col=min_vi$cl, pch = 19, main = "Our algorithm")

v = as.vector(3:(s+3))
col = c(rep(1,m1), rep(2,m-m1), v)

# real cluster
plot(data, col=col, pch = 19, main = "Real data")
