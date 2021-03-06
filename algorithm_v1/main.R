##### ALGORITHM IMPLEMENTATION

library(progress)

source("algorithm_v1/update_clusters.R")
source("algorithm_v1/update_xi.R")
source("algorithm_v1/Update_sigma.R")
source("algorithm_v1/update_theta.R")
source("algorithm_v1/update_beta.R")
source("algorithm_v1/auxiliary_functions.R")

# INPUT: Y -> data 
#        S_init -> initialization for the vector of labels Si i=1,...,n
#        sigma_init -> initialization for sigma 
#        theta_init -> initialization for theta 
#        beta_init -> initialization for beta 
#        xi_mu -> list xi_mui i=1,...,n which contains for each data point the value of mu initialized 
#        xi_cov -> list xi_covi i=1,...,n which contains for each data point the value of cov initialized 
#        Q_param -> parameters of the prior for Q
#        P_param -> parameters of the prior for P
#        n_iter -> number of iterations to run the algorithm

# OUTPUT: S -> matrix of cluster labels for all the iterations by row
#         xi_star -> list containing xi_mu_star and xi_cov_star related to the last iteration
#         sigma -> vector of sigma values for all the iterations
#         theta -> vector of theta values for all the iterations
#         beta -> vector of beta values for all the iterations
#         acc_sigma -> acceptance rate of sigma
#         acc_theta -> acceptance rate of theta
#         acc_beta -> acceptance rate of beta
algorithm <- function(Y, S_init, sigma_init, theta_init, beta_init, beta_param, sigma_param, theta_param, xi_mu, xi_cov, Q_param, P_param, n_iter, burnin, thinning ){
 
  # Variables initialization:
  # initialize the variables for the iterations from the input
  n <- dim(Y)[1]
  S_old <- S_init
  beta_old <- beta_init
  theta_old <- theta_init
  sigma_old <- sigma_init
  k_old <- length(unique(S_init))
  
  # Data storage initialization:
  # create the vectors that will store the values of beta, theta and sigma
  # generated at each iteration, plus the matrix storing the cluster labels S
  # at each iteration by row
  beta_vec = {}
  theta_vec = {}
  sigma_vec = {}
  S_matrix = matrix(nrow = 0, ncol = n)
  
  # Acceptance rates initialization:
  # initialize the counter variables for the number of accepted values over MH iterations
  # they will be returned divided by n_iter
  #acc_beta = 0
  acc_theta = 0
  acc_sigma = 0
  
  # Initialize the xi_star objects from the input variables xi_mu, xi_cov and S_init
  # further informations about the function init_xi_star are available
  # in the script "auxiliary_functions.R"
  xi <- init_xi_star(S_init, xi_mu, xi_cov)
  xi_mu_star <- xi$xi_mu_star
  xi_cov_star <- xi$xi_cov_star
  
  # Progress bar
  pb = progress_bar$new(total=n_iter)
  pb$tick(0)
  # Start to iterate
  for (r in 1:n_iter){
    
    # Step 2a: Updating the clusters
    S_old <- update_clusters(Y, xi_mu_star, xi_cov_star,
                                      S_old, beta_old, theta_old, sigma_old, P_param, Q_param)
    
   
    
    # Find the actual number of groups k
    k_old = max(S_old)
    
    # Store the current cluster labels by rows in the matrix
    if (r >= burnin + 1 && r %% thinning == 0) S_matrix = rbind(S_matrix, S_old)
    
    # If there are some groups, I need to update their parameters
    if(k_old > 0){
      # Step 2b: Updating the parameters of the groups’ distribution xi_star
      xi <- update_xi(Y, S_old, k_old, Q_param)
      
      # Updating the variables for next steps
      xi_mu_star <- xi$xi_mu_star
      xi_cov_star <- xi$xi_cov_star
    }
    
    # m1 and m1_bar computation
    # further informations about the functions are available in the script
    # "auxiliary_functions.R"
    m1_bar <- m1_bar(S_old)
    m1 <- m1(S_old)
    
    # Computation of the frequency vector for the cluster labels
    freq = rep(0,k_old)
    for (i in 1:k_old){
      freq[i] <- sum(S_old == i)
    }
    
    # Step 2c: Updating the parameters of the discrete component
    
    # Updating sigma
    sigma_list <- update_sigma(m1, m1_bar, k_old, sigma_old, theta_old, freq, acc_sigma, sigma_param)
    # Updating the variables
    sigma_old <- sigma_list$sigma
    if (r >= burnin + 1){
      if (r %% thinning == 0){
        sigma_vec <- append(sigma_vec, sigma_old)
      }
      acc_sigma <- sigma_list$acc
    }
    
    # Updating theta
    theta_list <- update_theta(n, m1_bar, k_old, theta_old, sigma_old, acc_theta, theta_param)
    # Updating the variables
    theta_old <- theta_list$theta
    if (r >= burnin + 1){
      if (r %% thinning == 0){
        theta_vec <- append(theta_vec, theta_old)
      }
      acc_theta <- theta_list$acc
    }
    
    # Step 2d: Updating the weight parameter
    beta_new <- update_beta(n, m1_bar, beta_param)
    # Updating the variables
    beta_old <- beta_new
    if (r >= burnin + 1 && r %% thinning == 0) beta_vec <- append(beta_vec, beta_old)
    
    # Progress bar
    pb$tick()
    
  }

  # Return all final values in a list
  return (list("S"= S_matrix, "xi_star" = xi, "sigma" = sigma_vec, "theta" = theta_vec,
               "beta" = beta_vec, "acc_sigma" = acc_sigma/n_iter, 
               "acc_theta" = acc_theta/n_iter))
  
}
