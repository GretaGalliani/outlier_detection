

algorithm <- function(Y, S_init, sigma_init, theta_init, beta_init, xi_mu, xi_cov, Q_param, P_param, n_iter ){
  
  source("update_clusters/update_clusters.R")
  source("update_xi/update_xi.R")
  #PASSO A
  #Input:
  # Y <- data.frame
  # xi_mu_star <- lista di n vettori --> struttura che contiene medi gruppi
  # xi_cov_star <- lista di n matrici --> struttura che contiene matrici cov dei gruppi
  # beta_old <- numero
  # theta_old <- numero
  # sigma_old <- numero
  # S_old <- vettore cluster passo precedente
  # k_old <- numero
  #nu_0_P & nu_0_Q
  #mu_0_P & mu_0_Q 
  #k_0_P & k_0_Q
  #lambda_0_P & lambda_0_Q 
  
  # Variables initialization
  S_old <- S_init
  beta_old <- beta_init
  theta_old <- theta_init
  sigma_old <- sigma_init
  k_old <- length(unique(S_init))
  
  xi <- init_xi_star(S_init, xi_mu, xi_cov)
  xi_mu_star <- xi$xi_mu_star
  xi_cov_star <- xi$xi_cov_star
  
  
  for (r in 1:n_iter){
    # Updating the clusters
    clusters <- update_clusters(Y, xi_mu, xi_cov, xi_mu_star, xi_cov_star,
                                      S_old, beta_old, theta_old, sigma_old, k_old,...)
    
    # Updating the variables for next steps
    S_old = clusters$S_new
    k_old <- length(unique(S_old))
    
    xi_mu_star <- clusters$xi_mu_star
    xi_cov_star <- clusters$xi_cov_star
    
    # Updating the parameters of the groups’ distribution xi_star
    
    #' Input
    #' Y data
    #' S updated at previous step
    #' k obtained from new S 
    #' Q_param list of parameters for the Q0 distribution
    
    for (j in 1:k_old){
      xi <- update_xi(Y, S_old, k_old, Q_param)
      
      
      
      
    }
    
    
    
    
  }

list("xi_mean" = xi_mu, "xi_sigma" = xi_sigma)
  
update_cluster <- function(Y, xi_mu, xi_cov, xi_mu_star, xi_cov_star, S_old,
                           beta_old, theta_old, sigma_old, k_old,
                           nu_0_P, k_0_P, mu_0_P, lambda_0_P,
                           nu_0_Q, k_0_Q, mu_0_Q, lambda_0_Q)
  
 




my_list <-list("S_new"=S_new, "xi_mu"=xi_mu,"xi_cov"=xi_cov,
               "xi_mu_star"=xi_mu_star,"xi_cov_star"=xi_cov_star)



}
