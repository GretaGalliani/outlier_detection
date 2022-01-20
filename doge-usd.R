library(MASS)
library(RColorBrewer)
library(robustbase)
library(TSstudio)

# Import dataset DOGE-USD
doge = read.csv('DOGE-USD.csv')
date = as.Date.character(doge$Date)
doge$Date = NULL
doge$LogReturn = rep(0, dim(doge)[1])

for (i in 2:dim(doge)[1]){
  doge$LogReturn[i] = log(doge$Adj.Close[i]/doge$Adj.Close[i-1])*100
}

doge$Adj.Close = NULL


ts_plot(data.frame(date = date, LogReturn=doge$LogReturn))

data = as.matrix(doge$LogReturn)
n = dim(data)[1]
d = dim(data)[2]

# Initialization of the parameters for the priors
Q_param = list()
P_param = list()

# Contaminated component
Q_param$k_0 = 1.25
Q_param$mu_0 = mean(data)
Q_param$nu_0 = d+3
Q_param$lambda_0 = var(data)

# Contaminant diffuse component
P_param$k_0 = 1
P_param$mu_0 = mean(data)
P_param$nu_0 = d+3
P_param$lambda_0 = var(data)

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
theta_param$a = 1
theta_param$b = 1

xi_mu <- list()
xi_cov <- list()

init_mu <- colMeans(data)
init_var <- cov(data)


for (i in 1:dim(data)[1]){
  xi_mu <- append(xi_mu, list(init_mu))
  xi_cov <- append(xi_cov, list(init_var))
}

source("main.R")
result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, beta_param, sigma_param, theta_param, xi_mu, xi_cov, Q_param, P_param, 2500, 100, 1)



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

par(mfrow=c(1,2)) 


# best cluster according to binder loss 
pal = brewer.pal(n = max(min_bind$cl), name = "Set1")
col_min_bind = rep(0,dim(data)[1])

bind_tab = table(min_bind$cl)
bind_pch = rep(17,n)
bind_pch[min_bind$cl %in% which(bind_tab>1)]=16

for (i in 1:length(col_min_bind))
{
  col_min_bind[i] = pal[min_bind$cl[i]]
}


plot(data, col=col_min_bind, pch = bind_pch, main = "Partition minimizing Binder Loss")



# IMPLEMENTING MIN VARIATION OF INFORMATION
# devtools::install_github("sarawade/mcclust.ext")
library(mcclust.ext)


# finds the clustering that minimizes  the lower bound to the posterior expected Variation of Information from Jensen's Inequality
min_vi <- minVI(psm, cls.draw=NULL, method=c("avg","comp","draws","greedy","all"), 
                max.k=NULL, include.greedy=FALSE, start.cl=NULL, maxiter=NULL,
                l=NULL, suppress.comment=TRUE)

pal = brewer.pal(n = max(min_vi$cl), name = "Set1")
col_min_vi = rep(0,dim(data)[1])

for (i in 1:length(col_min_bind))
{
  col_min_vi[i] = pal[min_vi$cl[i]]
}

vi_tab = table(min_vi$cl)
vi_pch = rep(17,n)
vi_pch[min_vi$cl %in% which(vi_tab>1)]=16


plot(data, col=col_min_vi, pch = vi_pch, main = "Partition minimizing VI Loss")

