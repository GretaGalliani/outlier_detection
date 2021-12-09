##### UPDATE XIs
library(LaplacesDemon)

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
    
    lambda_n_chol <- chol(lambda_n)
    inv_lambda_n_chol <- chol2inv(lambda_n_chol)
    
    df = nu_n-d+1
    xi_mu_j = rmvt(1, mu_n, lambda_n/k_n/df, df)
    
    xi_sigma_j = rinvwishartc(nu_n, chol(inv_lambda_n_chol))
    
    xi_mu <- append(xi_mu, list(xi_mu_j))
    xi_sigma <- append(xi_sigma, list(xi_sigma_j))
  }
  
  return(list("xi_mean" = xi_mu, "xi_sigma" = xi_sigma))
}






