

algorithm <- function(Y, S_init, sigma_init, theta_init, beta_init, xi_mu, xi_cov, n_iter){
  
  source("Update_clusters/Update_cluster.R")
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
  
  xi <- Initialization()
  xi_mu_star <- xi$xi_mu_star
  xi_cov_star <- xi$xi_cov_star
  
  
  for (r in 1:n_iter){
    # Updating the clusters
    clusters <- Update_cluster(Y, xi_mu, xi_cov, xi_mu_star, xi_cov_star,
                                      S_old, beta_old, theta_old, sigma_old, k_old,...)
    S_old = clusters$S_new
    
    
    k_old <- length(unique(clusters$S_new))
    
    
  }

update_cluster <- function(Y, xi_mu, xi_cov, xi_mu_star, xi_cov_star, S_old,
                           beta_old, theta_old, sigma_old, k_old,
                           nu_0_P, k_0_P, mu_0_P, lambda_0_P,
                           nu_0_Q, k_0_Q, mu_0_Q, lambda_0_Q)
  
 


source("Update_xi/Update_xi.R")

my_list <-list("S_new"=S_new, "xi_mu"=xi_mu,"xi_cov"=xi_cov,
               "xi_mu_star"=xi_mu_star,"xi_cov_star"=xi_cov_star)



}
