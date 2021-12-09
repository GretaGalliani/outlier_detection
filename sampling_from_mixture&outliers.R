library(MASS)

#y_test = seq (-1,1,length.out=10) 
#f_data = 0.5 * dmvt(y_test,delta = mean_a, sigma =sigma_b) + 0.5 * dmvt(y_test,delta = -mean_a, sigma =sigma_b) 
#n = 10

d = 2 #dimension
mean_a = rep(-3,d)#vector of means
sigma_b = diag(d)#sigma
m = 25 #m+s total number of samples (outliers included)
m1 = rbinom(1, size=m, prob = 0.5) #number of samples coming from the first gaussian 
m2 = m-m1 # from the second one

val1 = mvrnorm (m1,mu = mean_a, Sigma = sigma_b) #samples from first multivariate function
val2 = mvrnorm (m2,mu = -mean_a, Sigma = sigma_b)#samples from second multivariate function
allval = rbind(val1,val2)#combine
allval = allval[sample(m,m),] #randomizing rows

s=5#number of outliers
vector_sample = c(sample(1:m,s,replace = FALSE)) #take s=5 samples from m to make them outliers

allval[vector_sample,] = allval[vector_sample,]*100 #make them outliers

allval #matrix containing all the samples (outliers included)
#qchisq(0.9, df = 2)#the value of the quantile 0.9 of a chisquared = 4.60517
#s = 5 
#rtmvnorm(s,mean= c(0,0),sigma= 3^2*diag(2),lower= )