##### UPDATE BETA

# Function to update the beta parameter
# INPUT: n -> number of data points
#        m1_bar -> number of points in group 0 (to be computed using the auxiliary function m1_bar)
#        param$a -> first parameter for the beta prior for beta
#        param$b -> second parameter for the beta prior for beta


# OUTPUT: beta -> value of the parameter beta at the r iteration

update_beta <- function(n, m1_bar, beta_param) { 
  
  # Assuming a conjugate beta prior on beta for the bernoulli likelihood, we obtain a 
  # posterior beta distribution from which sampling the new value for beta
  
  beta_new <- rbeta(1, beta_param$a+n-m1_bar, beta_param$b+m1_bar)
  
  return(beta_new)

  }