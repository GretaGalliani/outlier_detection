##### UPDATE SIGMA

# Function to update the beta parameter
# INPUT: m1 -> number of singletons (number of points in group 0 + number of groups with one point)
#                                   (to be computed using the auxiliary function m1)
#        m1_bar -> number of points in group 0 (to be computed using the auxiliary function m1_bar)
#        k -> number of distinct groups (excluding the groups from the contaminated part)
#        sigma_old -> sigma at the previous iteration
#        theta_old -> theta at the current iteration
#        freq -> vector containing the number of points for each group j=1,...k
#        n_acc -> number of accepted proposals before the current iterations
#        sigma_param -> list containing the parameters of the prior for sigma (which is a beta)
#                       a -> shape1
#                       b -> shape2
#        sd -> standard deviation for the proposal density in the MHRW (default value is 2)

# OUTPUT: sigma -> value of the parameter sigma at the current iteration
#         acc -> number of accepted proposals at the current iteration

update_sigma <- function(m1, m1_bar, k, sigma_old, theta, freq, n_acc, sigma_param, sd = 1) { 

  # METROPOLIS HASTINGS RANDOM WALK
  
  # Extraction of a new value from the proposal distribution, doing an appropriate transformation 
  # to correct the fact that sigma is in (0,1)
  y <- inv_sigma(change_sigma(sigma_old) + rnorm(1,0,sd))
  
  # Computation of alpha of the new proposal wrt the old one 
  aprob <- compute_alpha_sigma(sigma_old, y, k, m1, m1_bar, theta, freq, sigma_param)
  
  # Sampling from a U(0,1)
  u <- runif(1) 
  
  # If u < aprob, the proposed value is accepted and the number of accepted values is increased 
  if (u < aprob){
    n_acc = n_acc+1
    
    # Return of the new value of sigma (equal to the proposed value y) and the accuracy 
    return(list("sigma" = y, "acc" = n_acc))
  } 
  else {
    
    # Return of the new value of sigma (equal to the previous value sigma_old) and the accuracy
    return(list("sigma" = sigma_old, "acc" = n_acc))
  }
}


# Function to change the value of sigma so sigma_star is in (-inf, +inf) to do MH
# INPUT: x -> value of sigma, x is in (0,1)
# OUTPUT: x_star -> value of transformed sigma, x_star is in (-inf, +inf)
change_sigma <- function(x) {
  return (log(x) - log(1-x))
}

# Function to change the value of sigma_star so sigma is in (0, 1) after MH
# inv_sigma is the inverse function of change_sigma
# INPUT: x -> value of sigma_star, x is in (-inf, +inf)
# OUTPUT: x_star -> value of sigma, x_star is in (0, 1)
inv_sigma <- function(x) {
  return (exp(x) / (exp(x) + 1))
}

# Function to compute the alpha(x,y) for the proposed value y and the old value x
# INPUT: x -> value of sigma at the previous step
#        y -> proposed value for sigma
#        k -> number of distinct groups (excluding the groups from the contaminated part)
#        m1 -> number of singletons (number of points in group 0 + number of groups with one point)
#        m1_bar -> number of points in group 0
#        theta -> theta at the current iteration
#        freq -> vector containing the number of points for each group j=1,...k
#        sigma_param -> list containing the prior parameter of sigma (which is a beta)
#                       a -> shape1
#                       b -> shape2

# OUTPUT: alpha -> the alpha needed to perform the acceptance/rejection in MH

compute_alpha_sigma <- function(x, y, k, m1, m1_bar, theta, freq, sigma_param) {
  freq_m1 = freq[freq>1]
  
  #We are using the log-form 
  #Pay attention the prior chosen is beta on (1,1) but it is egual to a uniform on (0,1)
  log_dens_y = dbeta(y, sigma_param$a, sigma_param$b, log = TRUE) + k*log(y) + lgamma(theta/y + k) - lgamma(theta/y) + log(y) + log(1-y)
  log_dens_x = dbeta(x, sigma_param$a, sigma_param$b, log = TRUE) + k*log(x) + lgamma(theta/x + k) - lgamma(theta/x) + log(x) + log(1-x)
  
  #computation the product of a sequence's part of the formula
  for (elem in freq_m1){
    log_dens_y = log_dens_y + lgamma(elem - y) - lgamma(1-y)
    log_dens_x = log_dens_x + lgamma(elem - x) - lgamma(1-x)
  }
  
  # Computation of alpha 
  return (min(1, exp(log_dens_y - log_dens_x)))
}
