

algorithm <- function(Y, S_init, sigma_init, theta_init, beta_init, xi_mu, xi_cov, Q_param, P_param, n_iter ){
  
  source("update_clusters/update_clusters.R")
  source("update_xi/update_xi.R")
  source("update_sigma/update_sigma.R")
  source("auxiliary_functions/auxiliary_functions_xi.R")
  
  
  # CAMBIARE NOME FILES DI R CHE INIZINO TUTTI CON LA MINUSCOLA PER DIO
  
  #PASSO A
  
  
  # Variables initialization
  S_old <- S_init
  beta_old <- beta_init
  theta_old <- theta_init
  sigma_old <- sigma_init
  k_old <- length(unique(S_init))
  
  # Acceptance rates initialization
  acc_beta = 0
  acc_theta = 0
  acc_sigma = 0
  
  xi <- init_xi_star(S_init, xi_mu, xi_cov)
  xi_mu_star <- xi$xi_mu_star
  xi_cov_star <- xi$xi_cov_star
  
  
  for (r in 1:n_iter){
    # Step 2a: Updating the clusters
    #Input:
    # Y <- data.frame
    # xi_mu_star <- list of n vectors containing group means
    # xi_cov_star <- list of n matrices containing group covs
    # beta_old <- numero
    # theta_old <- numero
    # sigma_old <- numero
    # S_old <- clusters vector at previous iteration
    # k_old <- number of clusters at previous iteration
    # P_param <- list
    # Q_param <- list
    
    clusters <- update_clusters(Y, xi_mu, xi_cov, xi_mu_star, xi_cov_star,
                                      S_old, beta_old, theta_old, sigma_old, k_old, P_param, Q_param)
    
    # Updating the variables for next steps
    S_old = clusters$S_new
    k_old <- length(unique(S_old))
    
    for (j in 1:k_old){
      #' Step 2b: Updating the parameters of the groups’ distribution xi_star
      #' Input variables:
      #' Y data
      #' S updated at previous step
      #' k obtained from new S 
      #' Q_param list of parameters for the Q0 distribution
      xi <- update_xi(Y, S_old, k_old, Q_param)
    
      xi_mu_star <- xi$xi_mean
      xi_cov_star <- xi$xi_sigma
      
      # m1 and m1_bar computation
      m1_bar <- m1_bar(S_old)
      m1 <- m1(S_old)
      
      # Step 2c: Updating the parameters of the discrete component sigma
      #' Input variables:
      #' 
      sigma_list <- update_sigma(m1, m1_bar, k_old, sigma_old)
      
      
      
      
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
