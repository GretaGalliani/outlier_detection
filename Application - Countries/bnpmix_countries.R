library(BNPmix)

library(MASS)
library(RColorBrewer)
library(robustbase)
library(TSstudio)
library(rworldmap)
library(countrycode)

# Import dataset about countries development
df = read.csv('Application - Countries/countries_dataset.csv')
df[102,]$country='Federated States of Micronesia'

set.seed(1997)
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

## PYdensity
result <- PYdensity(data, mcmc = list(niter = 6000, nburn = 1000, hyper = FALSE)
                    , prior = list(m0 = colMeans(data), k0 = 0.7, n0 = d+3, 
                                   Sigma0 = diag(d)*0.6))

library(mcclust.ext)
library(mcclust)

aux = result$clust

# For a sample of clusterings of the same objects the proportion of clusterings 
# in which observation $i$ and $j$ are together in a cluster is computed and 
# a matrix containing all proportions is given out. 
psm <- comp.psm(aux)


# Find the clustering that minimizes  the lower bound to the posterior expected Variation of Information from Jensen's Inequality
min_vi <- minVI(psm, cls.draw=NULL, method=c("avg","comp","draws","greedy","all"), 
                max.k=NULL, include.greedy=FALSE, start.cl=NULL, maxiter=NULL,
                l=NULL, suppress.comment=TRUE)

vi_tab = table(min_vi$cl)
vi_pch = rep(17,n)
vi_pch[min_vi$cl %in% which(vi_tab>1)]=16

# Clusters and outliers
out_i = as.numeric(which(vi_tab==1))
cl_i = as.numeric(which(vi_tab>1))
out_names = country[which(min_vi$cl %in% out_i)]

# Map plot
map_data = data.frame(country = country, cluster = min_vi$cl)

map = joinCountryData2Map(map_data, joinCode = "NAME", nameJoinColumn = "country", 
                          nameCountryColumn = "country", suggestForFailedCodes = TRUE, 
                          mapResolution = "medium", projection = NA, verbose = FALSE)

col = c('red', 'blue', 'green2', 'yellow', 'magenta')

mapParams = mapCountryData(map, nameColumnToPlot = 'cluster', catMethod = 'categorical', 
                           oceanCol = "azure2", missingCountryCol = 'white', 
                           colourPalette = col, mapTitle = '', 
                           addLegend = FALSE)

do.call(addMapLegendBoxes, c(mapParams,
                             x = 'bottom',
                             title = "Non-contaminated model clustering structure",
                             horiz = TRUE,
                             bg = "transparent",
                             bty = "n"))
dev.off()
