##### UPDATE THETA

# Function to update the theta parameter
# INPUT: n -> number of data points
#        m1_bar -> number of points in group 0 (to be computed using the auxiliary function m1_bar)
#        k -> number of distinct groups (excluding the groups from the contaminated part)
#        theta_old -> theta at the previous iteration
#        sigma -> sigma at the current iteration
#        n_acc -> number of accepted proposals in the previous iterations
#        theta_param -> list containing the prior parameter of theta (which is a gamma)
#                       a -> shape
#                       b -> rate
#        sd -> standard deviation for the proposal density in the MHRW (default value is 2)

# OUTPUT: theta -> value of the parameter theta at the r iteration
#         acc -> number of accepted proposals at the current iteration

update_theta <- function(n, m1_bar, k, theta_old, sigma, n_acc, theta_param, sd = 1){ 
  # METROPOLIS HASTINGS RANDOM WALK 
  
  # Extraction of a new value from the proposal distribution, doing an appropriate transformation 
  # to correct the fact that theta is in (0,+inf)
  y <- inv_theta( change_theta(theta_old) + rnorm(1,0,sd))
  
  # Computation the alpha of the new proposal wrt the old one
  aprob <- compute_alpha_theta(theta_old, y, k, m1_bar, sigma, n, theta_param)
  
  
  # Sampling from a U(0,1)
  u <- runif(1) 
  
  # If u < aprob, the proposed value is accepted and the number of accepted values is increased 
  if (u < aprob){
    n_acc = n_acc+1
    
    # Returning new value of sigma (equal to the proposed value y) and the accuracy
    return(list("theta" = y, "acc" = n_acc))
  } 
  else {
    
    # Returning new value of sigma (equal to the previous value sigma_old) and the accuracy
    return(list("theta" = theta_old, "acc" = n_acc))
    
  }
}

# Function to change the value of theta so theta_star is in (-inf, +inf) to do MH
# INPUT: x -> value of theta, x is in (0,+inf)
# OUTPUT: x_star -> value of transformed theta, x_star is in (-inf, +inf)
change_theta <- function(x) {
  return (log(x))
}

# Function to change the value of theta_star so theta is in (0, +inf) after MH
# Function to change the value of theta_star so theta is in (0, +inf) after MH
# inv_theta is the inverse function of change_theta
# INPUT: x -> value of theta_star, x is in (-inf, +inf)
# OUTPUT: x_star -> value of theta, x_star is in (0, +inf)
inv_theta <- function(x) {
  return (exp(x))
}

# Function to compute the alpha(x,y) for the proposed value y and the old value x
# INPUT: x -> value of theta at the previous step
#        y -> proposed value for theta
#        k -> number of distinct groups (excluding the groups from the contaminated part)
#        m1_bar -> number of points in group 0
#        sigma -> sigma at the current iteration
#        n -> number of data points
# OUTPUT: alpha -> the alpha needed to perform the acceptance/rejection in MH
compute_alpha_theta <- function(x, y, k, m1_bar, sigma, n, theta_param) {
  
  #We are using the log-form 
  #The prior chosen is a Gamma
  log_dens_y <- dgamma(y, theta_param$a, rate=theta_param$b, log = TRUE) + lgamma(y) + lgamma(y/sigma + k) - lgamma(y/sigma) - lgamma(y + n - m1_bar) + log(y)
  log_dens_x <- dgamma(x, theta_param$a, rate=theta_param$b, log = TRUE) + lgamma(x) + lgamma(x/sigma + k) - lgamma(x/sigma) - lgamma(x + n - m1_bar) + log(x)
  
  
  # Computation of alpha 
  return (min(1, exp(log_dens_y - log_dens_x)))
}
