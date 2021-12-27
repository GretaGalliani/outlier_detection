### LORO COSTRUISCONO I DATI 

set.seed(24091998)

# INIZIALIZZAZIONE - P0 DIVERSO DA Q0
Q_param = list()
P_param = list()

d = 2

Q_param$k_0 = 1
Q_param$mu_0 = c(0,0)
Q_param$nu_0 = d + 3 # it must be > (p-1)
Q_param$lambda_0 = diag(diag(cov(data)))

P_param$k_0 = 0.25
P_param$mu_0 = c(0,0)
P_param$nu_0 = d + 3 # it must be > (p-1)
P_param$lambda_0 = diag(diag(cov(data)))

n = dim(data)[1]
S_init = rep(1,n)

beta_param = list()
sigma_param = list()
theta_param = list()
beta_param$a = 1
beta_param$b = 1
sigma_param$a = 1
sigma_param$b = 1
theta_param$a = 2
theta_param$b = 0.02

xi_mu <- list()
xi_cov <- list()

init_mu <- colMeans(data)
init_var <- cov(data)

for (i in 1:n){
  xi_mu <- append(xi_mu, list(init_mu))
  xi_cov <- append(xi_cov, list(init_var))
}

# RUNNING THE ALGORITHM 
source("main.R")
result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, beta_param, sigma_param, theta_param, xi_mu, xi_cov, Q_param, P_param, 15000, 5000, 10)

# PARAMETER ANALYSIS
x11()
par(mfrow=c(1,3))
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(result$sigma,xlab="iteration",ylab=expression(sigma))  
plot(result$theta,xlab="iteration",ylab=expression(theta))  
plot(result$beta,xlab="iteration",ylab=expression(beta))  

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

# CLUSTER ANALYSIS