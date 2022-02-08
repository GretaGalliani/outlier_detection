library(MASS)
library(RColorBrewer)
library(robustbase)
library(TSstudio)
library(rworldmap)
library(countrycode)

# Dataset location:
# https://www.kaggle.com/amritachatterjee09/clustering-help-international-ngo-aid/data
# Import dataset about countries development
df = read.csv('Country-data.csv')
df[102,]$country='Federated States of Micronesia'

set.seed(1997)
rows = sample(nrow(df))
df = df[rows,]

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

#plot(df)
#plot(df$gdpp, df$life_expec)

df$life_expec=log(df$life_expec)
df$gdpp=log(df$gdpp)
df$child_mort=log(df$child_mort)
df$income=log(df$income)
df$total_fer=log(df$total_fer)

country = df$country

data = df
data$country = NULL

data = as.matrix(scale(data))
n = dim(data)[1]
d = dim(data)[2]

# Initialization of the parameters for the priors
Q_param = list()
P_param = list()

# Contaminated component
Q_param$k_0 = 0.8
Q_param$mu_0 = colMeans(data)
Q_param$nu_0 = d+4
Q_param$lambda_0 = diag(d)*0.5

# Contaminant diffuse component

P_param$k_0 = 0.05
P_param$mu_0 = colMeans(data)
P_param$nu_0 = d+4
P_param$lambda_0 = diag(d)*0.5

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

source("main.R")
result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, 
                    beta_param, s2igma_param, theta_param, xi_mu, xi_cov, 
                    Q_param, P_param, 6000, 1500, 5)
save(result, file='country_Q08_P005_CovDiv_05_.RData')

# Diagnostics

# Sigma
jpeg("sigma.jpg", width = 500, height = 500)
plot(result$sigma, type='l', ylim=c(0,1), xlab='iter', ylab=expression(sigma))
lines(rep(mean(result$sigma),length(result$sigma)), type='l', col='red' )
grid(nx=NULL, ny=NULL, lty = 1, col = "gray", lwd = 1)
dev.off()

plot(density(result$sigma))
result$acc_sigma
mean(result$sigma)
sd(result$sigma)

# Theta
jpeg("theta.jpg", width = 500, height = 500)
plot(result$theta, type='l', xlab='iter', ylab=expression(theta))
lines(rep(mean(result$theta),length(result$theta)), type='l', col='red' )
grid(nx=NULL, ny=NULL, lty = 1, col = "gray", lwd = 1) 
dev.off()

plot(density(result$theta))
result$acc_theta
mean(result$theta)
sd(result$theta)

# Beta
jpeg("beta.jpg", width = 500, height = 500)
plot(result$beta, type='l', ylim=c(0,1), xlab='iter', ylab=expression(beta))
lines(rep(mean(result$beta),length(result$beta)), type='l', col='red' )
grid(nx=NULL, ny=NULL, lty = 1, col = "gray", lwd = 1)
dev.off()

plot(density(result$beta))
mean(result$beta)
sd(result$beta)

# M1 bar
result_m1_bar = {}
for(i in 1:dim(result$S)[1])
  result_m1_bar = append(result_m1_bar, m1_bar(result$S[i,]))

jpeg("m1_bar.jpg", width = 500, height = 500)
plot(result_m1_bar, type='l', xlab='iter', ylab=expression(bar(m[1])))
lines(rep(mean(result_m1_bar),length(result_m1_bar)), type='l', col='red' )
grid(nx=NULL, ny=NULL, lty = 1, col = "gray", lwd = 1)
dev.off()

plot(density(result$beta))
mean(result$beta)
sd(result$beta)





# Plot of AUTOCORRELATION
library(coda)
par(mfrow=c(1,2))
tmp1 <- acf(result$sigma, main='Autocorrelation of sigma')
tmp2 <- acf(result$theta, main='Autocorrelation of theta')

# ESS
effectiveSize(result$sigma)
effectiveSize(result$theta)

max <- c()
for (i in 1:1000){
  max <- c(max, max(result$S[i,]))
}

mean(max) # mean number of clusters by the algorithm 


# IMPLEMENTING MIN BINDER LOSS
library(mcclust)

aux = result$S

# Compute the probability of a single country to be an outlier 
# based on relative frequency

freq = rep(0,n)
for (i in 1:n){
  freq[i] = length(which(aux[,i]==0))
}
rel_freq = freq/dim(aux)[1]

# Take the 5 most probable outliers by relative frequencies
country[order(rel_freq, decreasing = TRUE)[1:10]]



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

# finds the clustering that minimizes the posterior expectation of Binders loss function
min_bind <-  minbinder(psm, cls.draw = NULL, method = c("avg", "comp", "draws", 
                                                        "laugreen","all"), max.k = NULL, include.lg = FALSE, 
                       start.cl = NULL, tol = 0.001)

par(mfrow=c(1,2)) 


# best cluster according to binder loss 
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


# IMPLEMENTING MIN VARIATION OF INFORMATION
# devtools::install_github("sarawade/mcclust.ext")
library(mcclust.ext)


# finds the clustering that minimizes  the lower bound to the posterior expected Variation of Information from Jensen's Inequality
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

plot(df, col=col_min_vi, pch=vi_pch)

# Try a logarithmic transformation
plot(log(df[,c(2,10)]), col=col_min_vi, pch=vi_pch)

# Clusters and outliers
out_i = as.numeric(which(vi_tab==1))
cl_i = as.numeric(which(vi_tab>1))
out_names = country[which(min_vi$cl %in% out_i)]

# NUMBER OF OUTLIERS AND CLUSTERS
print(paste0('Number of outliers: ', length(out_i)))
print(paste0('Number of clusters: ', length(cl_i)))

factor(min_vi$cl) 
table(min_vi$cl)

# Scale outliers
# clust_map = min_bind$cl
# out_i = as.numeric(which(bind_tab==1))
# cl_i = as.numeric(which(bind_tab>1))
# out_names = country[which(min_bind$cl %in% out_i)]
clust_map = min_vi$cl 
clust_map[which(clust_map %in% out_i)] = 0

uniq = 1:length(unique(clust_map))
levels = as.numeric(levels(factor(clust_map)))
for (i in uniq){
  if(levels[i]>uniq[i]-1){
    clust_map[which(clust_map==levels[i])]=i-1
  }
}


print("Outliers: ")
print(out_names)

count = 1
for (i in 1:length(cl_i)){
  print(paste0("Cluster number ", count))
  print(country[which(min_vi$cl == cl_i[i])])
  count = count+1
}

# Map plot
map_data = data.frame(country = country, cluster = clust_map)

map = joinCountryData2Map(map_data, joinCode = "NAME", nameJoinColumn = "country", 
                    nameCountryColumn = "country", suggestForFailedCodes = TRUE, 
                    mapResolution = "coarse", projection = NA, verbose = FALSE)

#jpeg("map_Q07_P005_CovDiv_08.jpg", width = 700, height = 700)
mapCountryData(map, nameColumnToPlot = 'cluster', catMethod = 'categorical', 
               colourPalette = 'rainbow', mapTitle = 'Country clusters')
#dev.off()


# Focus on United States, Luxembourg, Haiti and Burundi

red_mark = function(name, countries){
  col = rep('#000000', length(countries))
  pch = rep(20, length(countries))
  col[which(country==name)] = '#FF0000'
  pch[which(country==name)] = 8
  return(list(col=col, pch=pch))
}


p = df[,c(3,4,7,8,10)]


for (name in out_names){
  mark = red_mark(name, country)
  x11()
  plot(p, col=mark$col, pch=mark$pch, main = paste0('Focus on ',name))
}

US = red_mark('United States', country)
Lux = red_mark('Luxembourg', country)
Ha = red_mark('Haiti', country)
Bur = red_mark('Burundi', country)


plot(p, col=US$col, pch=US$pch, main = 'Focus on United States')
plot(p, col=Lux$col, pch=Lux$pch, main = 'Focus on Luxembourg')
plot(p, col=Ha$col, pch=Ha$pch, main = 'Focus on Haiti')
plot(p, col=Bur$col, pch=Bur$pch, main = 'Focus on Burundi')

graphics.off()

# DA QUI IN POI INUTILE
# Compare outliers value with clusters means
# Consider the covariate gdpp
data_ns = df
data_ns$country = NULL
data_ns = as.matrix(data_ns)

cl_means <- list()
for (i in 1:length(cl_i)-1){
  cl_data = data_ns[which(clust_map==i),]
  m = colMeans(cl_data)
  cl_means = append(cl_means, list(m))
}

gdpp_cl =rep(0,n)
for (i in 1:n){
  if (clust_map[i] == 0)
    gdpp_cl[i]=data_ns[i,9]
  else
    gdpp_cl[i]=as.numeric(cl_means[[clust_map[i]]][9])
}

map_data_gdpp = data.frame(country = country, gdpp = data_ns[,9])

map_gdpp = joinCountryData2Map(map_data_gdpp, joinCode = "NAME", nameJoinColumn = "country", 
                          nameCountryColumn = "country", suggestForFailedCodes = TRUE, 
                          mapResolution = "coarse", projection = NA, verbose = FALSE)

mapCountryData(map_gdpp, nameColumnToPlot = 'gdpp', catMethod = 'fixedWidth', 
               colourPalette = 'heat', mapTitle = 'Country clusters')

length(unique(result$S[500,]))

