##### UPDATE XIs


update_xi <- function(Y, S, k, xi_old_mean, xi_old_sd, mu0, k0, nu0 ){
  for (j in [1:k]){
    Yj = Y[S==j]    # select data belonging to cluster j
    S0 = diag(cov(Yj))
    n = len(Yj)
    mu_n = (k0*mu0+n*mean(Yj))/(k0+n)
    k_n = k0 + n
    nu_n = nu0 + n
    Sn = S0
    
    
  }
  
  
  
}