library(MASS)
library(RColorBrewer)
library(robustbase)
library(TSstudio)
library(rworldmap)
library(countrycode)
library(BNPmix)

# Dataset location:
# https://www.kaggle.com/amritachatterjee09/clustering-help-international-ngo-aid/data

# Covariates:
# country:	Name of the country
# child_mort:	Death of children under 5 years of age per 1000 live births
# exports:	Exports of goods and services per capita. Given as %age of the GDP 
#           per capita
# health:	Total health spending per capita. Given as %age of GDP per capita
# imports:	Imports of goods and services per capita. Given as %age of the GDP 
#           per capita
# Income:	Net income per person
# Inflation:	The measurement of the annual growth rate of the Total GDP
# life_expec:	The average number of years a new born child would live if the 
#             current mortality patterns are to remain the same
# total_fer:	The number of children that would be born to each woman if the 
#             current age-fertility rates remain the same.
# gdpp:	The GDP per capita. Calculated as the Total GDP divided by the total 
#       population.

# Import dataset
df = read.csv('countries_dataset.csv')

# Change the namecode for Micronesia, needed when generating maps
df[102,]$country='Federated States of Micronesia'

# Shuffle rows
set.seed(1997)
rows = sample(nrow(df))
df = df[rows,]

# Data visualization
plot(df)

# Apply logtransformations where needed
df$life_expec=log(df$life_expec)
df$gdpp=log(df$gdpp)
df$child_mort=log(df$child_mort)
df$income=log(df$income)
df$total_fer=log(df$total_fer)

# Hold out country names
country = df$country
data = df
data$country = NULL

# Generate data structure to be given as input to the algorithm
data = as.matrix(scale(data))
n = dim(data)[1]
d = dim(data)[2]

# Initialization of the parameters for the priors
Q_param = list()
P_param = list()

# Contaminated component
Q_param$k_0 = 0.7
Q_param$mu_0 = colMeans(data)
Q_param$nu_0 = d+4
Q_param$lambda_0 = diag(d)*0.6

# Contaminant diffuse component
P_param$k_0 = 0.05
P_param$mu_0 = colMeans(data)
P_param$nu_0 = d+4
P_param$lambda_0 = diag(d)*0.6

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

init_mu <- colMeans(data)
init_var <- cov(data)

for (i in 1:dim(data)[1]){
  xi_mu <- append(xi_mu, list(init_mu))
  xi_cov <- append(xi_cov, list(init_var))
}

# Import algorithm framework
source("algorithm_v1/main.R")

# Run the MCMC
result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, 
                    beta_param, sigma_param, theta_param, xi_mu, xi_cov, 
                    Q_param, P_param, 12000, 2000, 10)

# Import the results for the best model
load('countries_result.RData')

aux = result$S

# MCMC Diagnostics

# Sigma

# Traceplot
plot(result$sigma, type='l', ylim=c(0,1), xlab='iter', ylab=expression(sigma))
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
plot(density(result$sigma), xlim = c(0.3,1), ylim = c(0,6), type='l', 
     xlab = expression(sigma), col = 'blue', 
     main = expression(paste('Posterior density of ', sigma)))
abline(v = seq(0.3, 1, 0.05), lty = 1, col = "gray", lwd=0.5)
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
plot(density(result$theta), xlim = c(0,15), ylim = c(0,0.3), type='l', 
     xlab = expression(theta), col = 'blue', 
     main = expression(paste('Posterior density of ', theta)))
abline(h = seq(0, 0.3, 0.025), lty = 1, col = "gray", lwd=0.5)
abline(v = seq(0, 15, 1), lty = 1, col = "gray", lwd=0.5)
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
abline(h = seq(5, 35, 2.5), lty = 1, col = "gray", lwd=0.5)
abline(v = seq(0, 1000, 100), lty = 1, col = "gray", lwd=0.5)
lines(result_m1_bar, type='l')
lines(rep(mean(result_m1_bar),length(result_m1_bar)), type='l', col='red' )
dev.off()

# Posterior summary
mean(result_m1_bar)
sd(result_m1_bar)

# Posterior density
plot(NULL, xlim = c(0,40), ylim = c(0,0.12), main = expression(paste('Posterior density of ', bar(m[1]))),
     xlab = expression(bar(m[1])), ylab = 'Density')
abline(h = seq(0, 0.12, 0.02), lty = 1, col = "gray", lwd=0.5)
abline(v = seq(0, 40, 2), lty = 1, col = "gray", lwd=0.5)
abline(v = mean(result_m1_bar), lty = 2, col = "blue", lwd=1)
hist(result_m1_bar, col = rgb(0, 0, 1, alpha = 0.5), freq=FALSE, add=TRUE, 
     xlim = c(0,40),ylim = c(0,0.12), breaks = length(unique(result_m1_bar)))
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

# Compute the probability of a single country to be an outlier 
# based on relative frequency

freq = rep(0,n)
for (i in 1:n){
  freq[i] = length(which(aux[,i]==0))
}
rel_freq = freq/dim(aux)[1]

# Implementing Binder loss
library(mcclust)
# devtools::install_github("sarawade/mcclust.ext")
library(mcclust.ext)

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
plotpsm(psm)

# Find the clustering that minimizes the posterior expectation of Binders loss function
min_bind <-  minbinder(psm, cls.draw = NULL, method = c("avg", "comp", "draws", 
                                                        "laugreen","all"), max.k = NULL, include.lg = FALSE, 
                       start.cl = NULL, tol = 0.001)

# Best clustering according to binder loss 
pal = rainbow(max(min_bind$cl))
col_min_bind = rep(0,dim(data)[1])

bind_tab = table(min_bind$cl)
bind_pch = rep(17,n)
bind_pch[min_bind$cl %in% which(bind_tab>1)]=16

for (i in 1:length(col_min_bind))
{
  col_min_bind[i] = pal[min_bind$cl[i]]
}

plot(df$gdpp, df$life_expec, col=col_min_bind, pch = bind_pch, main = "Partition minimizing Binder Loss")


# Implementing Variation of Information loss

# Find the clustering that minimizes  the lower bound to the posterior expected Variation of Information from Jensen's Inequality
min_vi <- minVI(psm, cls.draw=NULL, method=c("avg","comp","draws","greedy","all"), 
                max.k=NULL, include.greedy=FALSE, start.cl=NULL, maxiter=NULL,
                l=NULL, suppress.comment=TRUE)

pal = rainbow(max(min_vi$cl))
col_min_vi = rep(0,dim(data)[1])

for (i in 1:length(col_min_bind))
{
  col_min_vi[i] = pal[min_vi$cl[i]]
}

vi_tab = table(min_vi$cl)
vi_pch = rep(17,n)
vi_pch[min_vi$cl %in% which(vi_tab>1)]=16

plot(df$gdpp, df$life_expec, col=col_min_vi, pch = vi_pch, main = "Partition minimizing VI Loss")

# Print outliers and clusters
out_i = as.numeric(which(vi_tab==1))
cl_i = as.numeric(which(vi_tab>1))
out_names = country[which(min_vi$cl %in% out_i)]

print(paste0('Number of outliers: ', length(out_i)))
print(paste0('Number of clusters: ', length(cl_i)))

print("Outliers: ")
print(out_names)

count = 1
for (i in 1:length(cl_i)){
  print(paste0("Cluster number ", count))
  print(country[which(min_vi$cl == cl_i[i])])
  count = count+1
}

# Scale outliers in order to be labeled with 0
clust_map = min_vi$cl 
clust_map[which(clust_map %in% out_i)] = 0

uniq = 1:length(unique(clust_map))
levels = as.numeric(levels(factor(clust_map)))
for (i in uniq){
  if(levels[i]>uniq[i]-1){
    clust_map[which(clust_map==levels[i])]=i-1
  }
}

# Produce the map highlighting the outliers

out_map = data.frame(country = out_names, cluster = rep(0, length(out_names)))
o_map = joinCountryData2Map(out_map, joinCode = "NAME", nameJoinColumn = "country", 
                          nameCountryColumn = "country", suggestForFailedCodes = TRUE, 
                          mapResolution = "medium", projection = NA, verbose = FALSE)

out_mapParams = mapCountryData(o_map, nameColumnToPlot = 'cluster', catMethod = 'categorical', 
                           oceanCol = "azure2", missingCountryCol = 'white', 
                           colourPalette = 'red', mapTitle = '', 
                           addLegend = FALSE)
dev.off()

# Focus on United States and Haiti

# Function marking with a red cross the analyzed outlier among the data cloud
red_mark = function(name, countries){
  col = rep('#000000', length(countries))
  pch = rep(20, length(countries))
  col[which(country==name)] = '#FF0000'
  pch[which(country==name)] = 8
  return(list(col=col, pch=pch))
}

US = red_mark('United States', country)
Ha = red_mark('Haiti', country)

# Select relevant covariates to be plotted
p = df[,c(2,4,7,8,10)]

plot(p, col=US$col, pch=US$pch, main = 'Focus on United States')
dev.off()

plot(p, col=Ha$col, pch=Ha$pch, main = 'Focus on Haiti')
dev.off()

# Produce the map of the clustering structure

map_data = data.frame(country = country, cluster = clust_map)

map = joinCountryData2Map(map_data, joinCode = "NAME", nameJoinColumn = "country", 
                          nameCountryColumn = "country", suggestForFailedCodes = TRUE, 
                          mapResolution = "medium", projection = NA, verbose = FALSE)

col = c('black', 'blue', 'red', 'wheat', 'yellow', 'chocolate', 'orange', 
        'lightblue', 'goldenrod1', 'darkgreen', 'green2', 'magenta')

mapParams = mapCountryData(map, nameColumnToPlot = 'cluster', catMethod = 'categorical', 
                           oceanCol = "azure2", missingCountryCol = 'white', 
                           colourPalette = col, mapTitle = '', 
                           addLegend = FALSE)

do.call(addMapLegendBoxes, c(mapParams,
                             x = 'bottom',
                             title = "Clusters (outliers as 0)",
                             horiz = TRUE,
                             bg = "transparent",
                             bty = "n"))
dev.off()


