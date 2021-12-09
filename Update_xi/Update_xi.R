##### UPDATE XIs


update_xi <- function(Y, S, k, Q_param){
  
  mu_0 = Q_param$mu_0
  k_0 = Q_param$k_0
  lambda_0 = Q_param$lambda_0
  nu_0 = Q_param$nu_0
  
  xi_mu = list()
  xi_sigma = list()
  
  for (j in 1:k){
    Yj = Y[S==j]    # select data belonging to cluster j
    sample_cov = cov(Yj)
    n = dim(Yj)[1]
    d = dim(Yj)[2]
    # Compute the parameters for the marginals
    k_n = k0 + n
    mu_n = (k0*mu0+n*mean(Yj))/k_n
    nu_n = nu0 + n
    lambda_n = lambda_0 + sample_cov + k0*n/k_n*(mean(Yj)-mu0)%*%t((mean(Yj)-mu0))
    
    df = nu_n-d+1
    xi_mu_j = rmvt(1, mu_n, lambda_n/k_n/df, df)
    xi_sigma_j = rinvwishart(nu_n, inv(lambda_n))
    
    xi_mu <- append(xi_mu, list(xi_mu_j))
    xi_sigma <- append(xi_sigma, list(xi_sigma_j))
  }
  
  return(list("xi_mean" = xi_mu, "xi_sigma" = xi_sigma))
}

# How the created data structure works:
# list of lists where each inner list contains a vector and a matrix
# sigma = diag(2)
# l = list()
# l = append(l, list(sigma))
# l


rwish<-function(n,mu0,S)
{
  sS <- chol(S)
  S<-array( dim=c( dim(S),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(mu0 * dim(S)[1]), mu0, dim(S)[1]) %*% sS
    S[,,i]<- t(Z)%*%Z
  }
  return (S[,,1:n])
}
