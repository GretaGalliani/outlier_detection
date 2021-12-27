library(MASS)

#y_test = seq (-1,1,length.out=10) 
#f_data = 0.5 * dmvt(y_test,delta = mean_a, sigma =sigma_b) + 0.5 * dmvt(y_test,delta = -mean_a, sigma =sigma_b) 
#n = 10

d = 2 #dimension
mean_a = rep(-3,d)#vector of means
sigma_b = diag(d)#sigma
m = 20 #m number of samples (outliers excluded)
m1 = rbinom(1, size=m, prob = 0.5) #number of samples coming from the first gaussian 
m2 = m-m1 # from the second one

val1 = mvrnorm (m1,mu = mean_a, Sigma = sigma_b) #samples from first multivariate function
val2 = mvrnorm (m2,mu = -mean_a, Sigma = sigma_b) #samples from second multivariate function
allval = rbind(val1,val2) #combine


s=5 #number of outliers
i=0 
while(i<s){ #cycle to find s outliers
  value=mvrnorm(1,mu= rep(0,d),Sigma= 3^2*diag(d)) #sampling from a multivariate normal distribution
  module = norm(as.matrix(value), type="2")
  chi=qchisq(0.9, df = d)
  if(module^2>3*sqrt(chi)) #If we are sampling from the over-disperse truncated Gaussian distribution
    {
    i=i+1
    allval = rbind(allval,val)
    }
}
allval = allval[sample(m+s,m+s),] #randomizing rows
allval