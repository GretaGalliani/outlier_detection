##### UPDATE XIs
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
    Yj <- Y[which(S==j),]
    print(which(S==j))
    print(Yj)
    
    Yj <- as.matrix(Y[which(S==j),])
    print(Yj)
    
    # Extraction of the number of data points
    n = dim(Yj)[1]
    print(n)
    
    # Extraction of the dimension of data
    d = dim(Yj)[2]
    
    # Sample covariance
    if (n==1)
      sample_cov = 0
    else
      sample_cov = cov(Yj)
    
    # Compute the parameters for the marginals
    k_n = Q_param$k_0 + n
    mu_n = (Q_param$k_0*Q_param$mu_0+n*mean(Yj))/k_n
    nu_n = Q_param$nu_0 + n
    lambda_n = Q_param$lambda_0 + sample_cov + Q_param$k_0*n/k_n*(mean(Yj)-Q_param$mu_0)%*%t((mean(Yj)-Q_param$mu_0))
    
    # Application of inversion based on cholesky decomposition
    lambda_n_chol <- chol(lambda_n)
    inv_lambda_n_chol <- chol2inv(lambda_n_chol)
    
    # computation of degrees of freedom
    df = nu_n-d+1
    
    # Sampling of mu (from a multivariate t-student)
    xi_mu_j = rmvt(1, mu_n, lambda_n/k_n/df, df)
    
    # Sampling of cov (from an inverse-wishart), based on cholesky decomposition
    xi_sigma_j = rinvwishartc(nu_n, chol(inv_lambda_n_chol))
    
    # Appending of the sampled values
    xi_mu <- append(xi_mu, list(xi_mu_j))
    xi_sigma <- append(xi_sigma, list(xi_sigma_j))
  }
  
  return(list("xi_mu_star" = xi_mu, "xi_cov_star" = xi_sigma))
}






