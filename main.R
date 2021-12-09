source("update_clusters/update_clusters.R")
source("update_xi/update_xi.R")
source("update_sigma/update_sigma.R")
source("update_theta/update_theta.R")
source("update_beta/update_beta.R")
source("auxiliary_functions/auxiliary_functions_xi.R")

algorithm <- function(Y, S_init, sigma_init, theta_init, beta_init, xi_mu, xi_cov, Q_param, P_param, n_iter ){
 
  # Variables initialization
  n <- dim(Y)[1]
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
    
    clusters <- update_clusters(Y, xi_mu_star, xi_cov_star,
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
      
      # Freq computation
      freq <- as.integer(table(S))
      
      # Step 2c: Updating the parameters of the discrete component
      
      # Update sigma
      sigma_list <- update_sigma(m1, m1_bar, k_old, sigma_old, theta_old, freq, acc_sigma)
      sigma_old <- sigma_list$sigma
      acc_sigma <- sigma_list$acc
      
      # Update theta
      theta_list <- update_theta(n, m1_bar, k_old, theta_old, sigma_old, acc_theta)
      theta_old <- theta_list$theta
      acc_theta <- theta_list$acc
      
      # Step 2d: Updating the weight parameter
      beta_list <- update_beta(n, m1_bar, beta_old, acc_beta)
      beta_old <- beta_list$beta
      acc_beta <- beta_list$acc
      
    }

  }
  
  acc_beta <- acc_beta/n_iter
  acc_sigma <- acc_sigma/n_iter
  acc_theta <- acc_theta/n_iter
  
}
