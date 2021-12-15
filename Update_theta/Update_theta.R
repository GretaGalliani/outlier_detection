##### UPDATE THETA

# Function to update the theta parameter
# INPUT: n -> number of data points
#        m1_bar -> number of points in group 0 (to be computed using the auxiliary function m1_bar)
#        k -> number of distinct groups
#        theta_old -> theta at the previous iteration
#        sigma -> sigma at the current iteration
#        sd -> standard deviation for the proposal density in the MHRW (default value is 2)
#        n_acc -> number of accepted proposals in the previous iterations

# OUTPUT: theta -> value of the parameter theta at the r iteration
#         acc -> number of accepted proposals at the current iteration
update_theta <- function(n, m1_bar, k, theta_old, sigma, n_acc, sd = 2) { 
  
  # METROPOLIS HASTINGS RANDOM WALK 
  
  # Extraction of a new value from the proposal distribution, doing an appropriate transformation 
  # to correct the fact that theta is in (0,+inf)
  y <- inv_theta( change_theta(theta_old) + rnorm(1,0,sd))
  
  # Computation the alpha of the new proposal wrt the old one
  aprob <- compute_alpha_theta(theta_old, y, k, m1_bar, sigma, n)
  
  # Sampling from a U(0,1)
  u <- runif(1) 
  
  # If u < aprob, the proposed value is accepted and the number of accepted values is increased 
  if (u < aprob){
    n_acc = n_acc+1
    
    # Return of the new value of sigma (equal to the proposed value y) and the accuracy
    return(list("theta" = y, "acc" = n_acc))
  } 
  else {
    # Return of the new value of sigma (equal to the previous value sigma_old) and the accuracy
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
#        k -> number of distinct groups
#        m1_bar -> number of points in group 0
#        sigma -> sigma at the current iteration
#        n -> number of data points
# OUTPUT: alpha -> the alpha needed to perform the acceptance/rejection in MH
compute_alpha_theta <- function(x, y, k, m1_bar, sigma, n) {
  # Computation of the partial posterior density for y and x
  pi_y <- dens_theta(y, k, m1_bar, sigma, n)
  pi_x <- dens_theta(x, k, m1_bar, sigma, n)
  
  rapp <- pi_y/pi_x
  
  # Computation of alpha 
  return (min(1, rapp))
}


# Function to compute the partial posterior density (up to the normalizing constant) for theta
# INPUT: x -> point where the partial density is evaluated
#        k -> number of distinct groups
#        m1_bar -> number of points in group 0
#        sigma -> sigma at the current iteration
#        n -> number of data points
# OUTPUT: f -> the evaluation of the partial posterior density at point x
dens_theta <- function(x, k, m1_bar, sigma, n) {
  return ( gamma(x) * gamma(x/sigma + k) / (gamma(x/sigma) * gamma(x + n - m1_bar)) * (1/x) )
}




