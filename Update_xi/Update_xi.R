##### UPDATE XIs


update_xi <- function(Y, S, k, xi_old_mean, xi_old_sd, mu0, k0, lambda0, nu0){
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
    
    xi_sd = rinvwishart(nu_n, inv(lambda_n))
    
    df = nu_n-d+1
    xi_mu = rmvt(1, mu_n, lambda_n/k_n/df, df)
    
  }
  
  
  
}