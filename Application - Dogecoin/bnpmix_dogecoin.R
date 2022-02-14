library(BNPmix)

library(MASS)
library(RColorBrewer)
library(robustbase)
library(TSstudio)

# Dataset location:
# https://finance.yahoo.com/quote/DOGE-USD/history/

# Import dataset
doge = read.csv('Application - Dogecoin/Dogecoin_dataset.csv')

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

result <- PYdensity(data, mcmc = list(niter = 6000, nburn = 1000, hyper = FALSE)
                    , prior = list(m0 = colMeans(data), k0 = 0.9, n0 = d+5, 
                                   Sigma0 = diag(d)/6))

library(mcclust.ext)
library(mcclust)

aux = result$clust

# For a sample of clusterings of the same objects the proportion of clusterings 
# in which observation $i$ and $j$ are together in a cluster is computed and 
# a matrix containing all proportions is given out. 
psm <- comp.psm(aux)


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
