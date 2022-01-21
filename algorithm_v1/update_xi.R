##### UPDATE XIs
# Library to sample from an inverse wishart
library(LaplacesDemon)

# Function to update the groups' parameters (mu and cov)
# INPUT: Y -> data 
#        S -> vector of labels Si i=1,...,n
#        k -> number of distinct groups
#        Q_param -> parameters of the prior for Q

# OUTPUT: xi_mu_star -> list xi_mu_starj j=1,...,k which contains for each group the  value of mu 
#         xi_cov_star -> list xi_cov_starj j=1,...,k which contains for each group the value of cov
update_xi <- function(Y, S, k, Q_param){
  
  # Creation of the lists for mu and cov
  xi_mu = list()
  xi_sigma = list()
  
  # For each group
  for (j in 1:k){
    # Select data belonging to cluster j
    Yj <- as.matrix(Y[which(S==j),])
    
    # If there is only one data in the group, sample_cov is set to 0
    if (length(which(S==j))==1){
      n = 1
      d = length(Yj)
      sample_cov = 0
      sample_mean = rowMeans(as.matrix(Yj))
    }
    # else, compute sample_cov
    else{
      n = dim(Yj)[1]
      d = dim(Yj)[2]
      if (d==1) 
        sample_cov = var(Yj)
      else
        sample_cov = cov(Yj)
      sample_mean = colMeans(Yj)
      
    }
    
    # Compute the parameters for the marginals
    k_n = Q_param$k_0 + n
    mu_n = Q_param$k_0/(k_n)*Q_param$mu_0 + n/k_n*sample_mean
    nu_n = Q_param$nu_0 + n
    lambda_n = Q_param$lambda_0 + sample_cov + Q_param$k_0*n/k_n*(sample_mean-Q_param$mu_0)%*%t((sample_mean-Q_param$mu_0))
    
    # computation of degrees of freedom
    df = nu_n-d+1
    
    # Sampling of mu (from a multivariate t-student)
    xi_mu_j = LaplacesDemon::rmvt(1, mu_n, lambda_n/k_n/df, df)
    
    
    # Sampling of cov (from an inverse-wishart), based on cholesky decomposition
    # xi_sigma_j = LaplacesDemon::rinvwishart(nu_n, as.inverse(lambda_n))
    xi_sigma_j = LaplacesDemon::rinvwishart(nu_n, lambda_n)
    
    # Appending of the sampled values
    xi_mu <- append(xi_mu, list(xi_mu_j))
    xi_sigma <- append(xi_sigma, list(xi_sigma_j))
  }
  
  return(list("xi_mu_star" = xi_mu, "xi_cov_star" = xi_sigma))
}






