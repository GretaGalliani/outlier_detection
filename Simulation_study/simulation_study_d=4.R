#This function implement the algorithm for the simulation study:
#Generate a set of n data from a mixture of two Gaussian distribution 
#Generate a set of S outliers from an over-disperse truncated gaussian
#dimension = 4

#We used this script to developed different models tuning some hyperparameters:
# c in {1,1.25,1.5}
#K0_P in {0.5,0.25}
#K0_Q in {0.5,1}
#Number of samples without outliers {90,240}

set.seed(04021997)
library(MASS)
library(RColorBrewer)

# SAMPLING FROM THE DENSITY ####
# BASE MEASURE
d = 4 #dimension
mean_a = rep(-3,d) #vector of means
sigma_b = diag(d) #sigma
m = 240 #m number of samples (outliers excluded) {90, 240}
m1 = rbinom(1, size=m, prob = 0.5) #number of samples coming from the first gaussian 
m2 = m-m1 # from the second one

val1 = mvrnorm (m1,mu = mean_a, Sigma = sigma_b) #samples from first multivariate function
val2 = mvrnorm (m2,mu = -mean_a, Sigma = sigma_b) #samples from second multivariate function
allval = rbind(val1,val2) #combine

# CONTAMINATED MEASURE
s=10 #number of outliers
c = 1.5 # augmenting factor
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
    allval = rbind(allval,c*value) 
    #C has the role to shrink or expand the nuisance observations towards the origin
  }
}

# Constructing the dataset
rownames(allval)=NULL
allval
data <- allval

# Plotting the data with the original groups
pal = brewer.pal(n = 9, name = "Set1")
col_real = c(rep(pal[1], m1), rep(pal[3], m-m1), rep(pal[2],s))
pc = c(rep(16,m), rep(17,s))

pairs(data, col = col_real, pch = pc, main = "Real data")


# INITIALIZATION for P0 and Q0 ####
Q_param = list()
P_param = list()

d = dim(data)[2]

Q_param$k_0 = 0.5 #{1, 0.5}
Q_param$mu_0 = rep(0,d)
Q_param$nu_0 = d + 3 # it must be > (p-1)
Q_param$lambda_0 = diag(diag(cov(data)))

P_param$k_0 = 0.5 #{0.25,0.5}
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

# RUNNING THE ALGORITHM ####
source("algorithm_v1/main.R")
result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, beta_param, sigma_param, theta_param, xi_mu, xi_cov, Q_param, P_param, 15000, 1000, 10)

# PARAMETER ANALYSIS ####
x11()
par(mfrow=c(1,3))
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(result$sigma,xlab="iteration",ylab=expression(sigma), type = 'l')  
plot(result$theta,xlab="iteration",ylab=expression(theta), type = 'l')  
plot(result$beta,xlab="iteration",ylab=expression(beta), type = 'l')  

dev.off()

# Marginal traceplots
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


# NUMBER OF SINGLETONS
source("auxiliary_functions/auxiliary_functions.R")

n_singletons <- c()
for (i in 1:dim(result$S)[1]){
  n_singletons <- c(n_singletons, m1(result$S[i,]))
}

n_singletons

## LOSS FUNCTION ####
#Binder loss function

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



# best cluster according to binder loss 
# NOTE: We report only the graph for CPY1, n=100, c=1
tab_bind <- table(min_bind$cl)

pal = brewer.pal(n = 9, name = "Set1")
col_bind = c(rep(pal[1], tab[[1]]), rep(pal[3], tab[[2]]), rep(pal[4], tab[[3]]), rep(pal[5], tab[[4]]),
             rep(pal[6], tab[[5]]), rep(pal[7], tab[[6]]), rep(pal[8], tab[[7]]), rep(pal[2],12))


pc = c(rep(16,90), rep(17,10))

pairs(data, col = col_bind, pch = pc, main = "Clustering minimizing Binder loss")

pairs(data, col = min_bind$cl, pch = pc, main = "Clustering minimizing Binder loss")

# IMPLEMENTING MIN VARIATION OF INFORMATION
# devtools::install_github("sarawade/mcclust.ext")
library(mcclust.ext)

psm2 <- comp.psm(aux)

# finds the clustering that minimizes  the lower bound to the posterior expected Variation of Information from Jensen's Inequality
min_vi <- minVI(psm2, cls.draw=NULL, method=c("avg","comp","draws","greedy","all"), 
                max.k=NULL, include.greedy=FALSE, start.cl=NULL, maxiter=NULL,
                l=NULL, suppress.comment=TRUE)


# best cluster according to iv loss 
# NOTE: We report only the graph for CPY1, n=100, c=1
tab_vi <- table(min_vi$cl)

pal = brewer.pal(n = 9, name = "Set1")
col_vi = c(rep(pal[1], tab[[1]]), rep(pal[3], tab[[2]]), rep(pal[2],10))


pc = c(rep(16,90), rep(17,10))

pairs(data, col = col_vi, pch = pc, main = "Clustering minimizing VI loss")

pairs(data, col = min_vi$cl, pch = pc, main = "Clustering minimizing VI loss")


# RELEVANT PARAMETERS ####
#Beta (mean and sd)
mean(result$beta)
sd(result$beta)

#Theta (mean and sd)
mean(result$theta)
sd(result$theta)

#Sigma (mean and sd)
mean(result$sigma)
sd(result$sigma)

#Number of clusters (mean and sd)
mean(max)
sd(max)

#Numbero fo singletons (mean and sd)
mean(n_singletons)
sd(n_singletons)

#Number of singletons estimated by Binder loss function 
length(which(tab_bind==1))
#Number of clusters estimated by Binder loss function 
length(which(tab_bind>1))

#Number of singletons estimated by Variation of information
length(which(tab_vi==1))
#Number of clusters estimated by Variation of information 
length(which(tab_vi>1))

#Acceptance rate for sigma e theta 
result$acc_theta
result$acc_sigma