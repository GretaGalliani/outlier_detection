library(MASS)
library(RColorBrewer)
library(robustbase)
library(TSstudio)

# Dataset location:
# https://finance.yahoo.com/quote/DOGE-USD/history/

# Import dataset
doge = read.csv('DOGE-USD.csv')

# Calculating the log return variable from Adj.Close
doge$LogReturn = rep(0, dim(doge)[1])
for (i in 2:dim(doge)[1]){
  doge$LogReturn[i] = log(doge$Adj.Close[i]/doge$Adj.Close[i-1])*100
}

# Generate data structure to be given as input to the algorithm
data = as.matrix(doge$LogReturn)
data=data[1:365,]
data= scale(data)
n = dim(data)[1]
d = dim(data)[2]

# Initialization of the parameters for the priors
Q_param = list()
P_param = list()

# Contaminated component
Q_param$k_0 = 1
Q_param$mu_0 = mean(data)
Q_param$nu_0 = d+5
Q_param$lambda_0 = var(data)/6

# Contaminant diffuse component
P_param$k_0 = 0.01
P_param$mu_0 = mean(data)
P_param$nu_0 = d+5
P_param$lambda_0 = var(data)/6

# Initialization of the parameters for the Pitman-Yor and initial partition
S_init = rep(1, n)
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
theta_param$a = 2
theta_param$b = 1

xi_mu <- list()
xi_cov <- list()

init_mu <- mean(data)
init_var <- var(data)


for (i in 1:dim(data)[1]){
  xi_mu <- append(xi_mu, list(init_mu))
  xi_cov <- append(xi_cov, list(init_var))
}

# Import algorithm framework
source("main.R")
result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, beta_param, sigma_param, theta_param, xi_mu, xi_cov, Q_param, P_param, 12000, 2000, 10)

aux = result$S

# MCMC Diagnostics

# Sigma

#Traceplot
plot(result$sigma, type='l', xlim=c(0,1000), ylim=c(0,1), xlab='iter', ylab=expression(sigma))
abline(h = seq(0, 1, 0.1), lty = 1, col = "gray", lwd=0.5)
abline(v = seq(0, 1000, 100), lty = 1, col = "gray", lwd=0.5)
lines(result$sigma, type='l')
lines(rep(mean(result$sigma),length(result$sigma)), type='l', col='red' )
dev.off()

# Posterior summary
result$acc_sigma
mean(result$sigma)
sd(result$sigma)

# Posterior density
plot(density(result$sigma), ylim = c(0,4), type='l', 
     xlab = expression(sigma), col = 'blue', 
     main = expression(paste('Posterior density of ', sigma)))
abline(v = seq(0, 1.1, 0.1), lty = 1, col = "gray", lwd=0.5)
abline(h = seq(0, 6, 0.5), lty = 1, col = "gray", lwd=0.5)
abline(v = mean(result$sigma), lty = 2, col = "blue", lwd=1)
polygon(density(result$sigma), col = rgb(0, 0, 1, alpha = 0.5))
lines(density(result$sigma), type='l')
dev.off()

# Theta

# Traceplot
plot(result$theta, type='l', xlab='iter', ylab=expression(theta))
abline(h = seq(0, 13, 1), lty = 1, col = "gray", lwd=0.5)
abline(v = seq(0, 1000, 100), lty = 1, col = "gray", lwd=0.5) 
lines(result$theta, type='l')
lines(rep(mean(result$theta),length(result$theta)), type='l', col='red' )
dev.off()

# Posterior summary
result$acc_theta
mean(result$theta)
sd(result$theta)

# Posterior density
plot(density(result$theta), xlim = c(0,13), ylim = c(0,0.35), type='l', 
     xlab = expression(theta), col = 'blue', 
     main = expression(paste('Posterior density of ', theta)))
abline(h = seq(0, 0.4, 0.05), lty = 1, col = "gray", lwd=0.5)
abline(v = seq(0, 12, 1), lty = 1, col = "gray", lwd=0.5)
abline(v = mean(result$theta), lty = 2, col = "blue", lwd=1)
polygon(density(result$theta), col = rgb(0, 0, 1, alpha = 0.5))
lines(density(result$theta), type='l')
dev.off()

# Beta

# Traceplot
plot(result$beta, type='l', ylim=c(0,1), xlab='iter', ylab=expression(beta))
abline(h = seq(0, 1, 0.1), lty = 1, col = "gray", lwd=0.5)
abline(v = seq(0, 1000, 100), lty = 1, col = "gray", lwd=0.5)
lines(result$beta, type='l')
lines(rep(mean(result$beta),length(result$beta)), type='l', col='red' )
dev.off()

# Posterior summary
mean(result$beta)
sd(result$beta)

# Posterior density
plot(density(result$beta), type='l', xlab = expression(beta), col = 'blue', 
     main = expression(paste('Posterior density of ', beta)))
abline(h = seq(0, 11, 1), lty = 1, col = "gray", lwd=0.5)
abline(v = seq(0.6, 1, 0.025), lty = 1, col = "gray", lwd=0.5)
abline(v = mean(result$beta), lty = 2, col = "blue", lwd=1)
polygon(density(result$beta), col = rgb(0, 0, 1, alpha = 0.5))
lines(density(result$beta), type='l')
dev.off()

# M1 bar

# M1 bar computation for all the iterations
result_m1_bar = {}
for(i in 1:dim(result$S)[1])
  result_m1_bar = append(result_m1_bar, m1_bar(result$S[i,]))

# Traceplot
plot(result_m1_bar, type='l', xlab='iter', ylab=expression(bar(m[1])))
abline(h = seq(50, 150, 10), lty = 1, col = "gray", lwd=0.5)
abline(v = seq(0, 1000, 100), lty = 1, col = "gray", lwd=0.5)
lines(result_m1_bar, type='l')
lines(rep(mean(result_m1_bar),length(result_m1_bar)), type='l', col='red' )
dev.off()

# Posterior summary
mean(result_m1_bar)
sd(result_m1_bar)

# Posterior density
plot(NULL, xlim = c(55,140), ylim = c(0,0.05), main = expression(paste('Posterior density of ', bar(m[1]))),
     xlab = expression(bar(m[1])), ylab = 'Density')
abline(h = seq(0, 0.08, 0.005), lty = 1, col = "gray", lwd=0.5)
abline(v = seq(55, 140, 6), lty = 1, col = "gray", lwd=0.5)
abline(v = mean(result_m1_bar), lty = 2, col = "blue", lwd=1)
hist(result_m1_bar, col = rgb(0, 0, 1, alpha = 0.5), freq=FALSE, add=TRUE, 
     xlim = c(55,140),ylim = c(0,0.06), breaks = length(unique(result_m1_bar)))
dev.off()

# Plots of autocorrelation
library(coda)
par(mfrow=c(1,2))
tmp1 <- acf(result$sigma, main='Autocorrelation of sigma')
tmp2 <- acf(result$theta, main='Autocorrelation of theta')
graphics.off()

# Effective Sample Size
effectiveSize(result$sigma)
effectiveSize(result$theta)
effectiveSize(result$beta)
effectiveSize(result_m1_bar)

# IMPLEMENTING MIN BINDER LOSS
library(mcclust)

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
image(psm)

# Finds the clustering that minimizes the posterior expectation of Binders loss function
min_bind <-  minbinder(psm, cls.draw = NULL, method = c("avg", "comp", "draws", 
                                                        "laugreen","all"), max.k = NULL, include.lg = FALSE, 
                       start.cl = NULL, tol = 0.001)

# best cluster according to binder loss 
pal = rainbow(max(min_bind$cl))
col_min_bind = rep(0,dim(data)[1])

bind_tab = table(min_bind$cl)
bind_pch = rep(17,n)
bind_pch[min_bind$cl %in% which(bind_tab>1)]=16

cl = which(bind_tab>1)
for (i in 1:length(col_min_bind))
{
  if(min_bind$cl[i] %in% which(bind_tab==1))
    col_min_bind[i] = '#000000'
  else{
    col_min_bind[i] = pal[which(cl==min_bind$cl[i])]
  }
}

plot(doge$LogReturn[1:365], col=col_min_bind, pch = bind_pch, main = "Partition minimizing Binder Loss")
dev.off()


# IMPLEMENTING MIN VARIATION OF INFORMATION
# devtools::install_github("sarawade/mcclust.ext")
library(mcclust.ext)

# finds the clustering that minimizes  the lower bound to the posterior expected Variation of Information from Jensen's Inequality
min_vi <- minVI(psm, cls.draw=NULL, method=c("avg","comp","draws","greedy","all"), 
                max.k=NULL, include.greedy=FALSE, start.cl=NULL, maxiter=NULL,
                l=NULL, suppress.comment=TRUE)

pal = rainbow(max(min_vi$cl))
col_min_vi = rep(0,dim(data)[1])

vi_tab = table(min_vi$cl)
vi_pch = rep(17,n)
vi_pch[min_vi$cl %in% which(vi_tab>1)]=16

cl = which(vi_tab>1)
for (i in 1:length(col_min_vi))
{
  if(min_vi$cl[i] %in% which(vi_tab==1))
    col_min_vi[i] = '#000000'
  else{
    col_min_vi[i] = pal[which(cl==min_vi$cl[i])]
  }
}

plot(doge$LogReturn[1:365], col=col_min_vi, pch = vi_pch, main = "Partition minimizing VI Loss", ylab="Log return", xlab="Day")
dev.off()

# Print outliers and clusters
out_i = as.numeric(which(vi_tab==1))
cl_i = as.numeric(which(vi_tab>1))

print(paste0('Number of outliers: ', length(out_i)))
print(paste0('Number of clusters: ', length(cl_i)))
