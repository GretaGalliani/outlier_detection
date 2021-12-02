##### UPDATE BETA


# Function to update the beta parameter
# INPUT: n -> number of data points
#        m1_bar -> number of points in group 0 (to be computed using the auxiliary function m1_bar)
#        beta_old -> beta at the previous iteration
#        sd -> standard deviation for the proposal density in the MHRW 
#        n_acc -> number of accepted proposals in the previous iterations

# OUTPUT: beta -> value of the parameter beta at the r iteration
#         acc -> number of accepted proposals at the current iteration

update_beta <- function(n, m1_bar, beta_old, sd = 2, n_acc) { 
  # Extraction of a new value from the proposal distribution, doing an appropriate transformation 
  # to correct the fact that beta is in (0,1)
  y <- inv_beta( change_beta(beta_old) + rnorm(1,0,sd))
  
  # Computation the alpha of the new proposal wrt the old one 
  aprob <- calcolo_alpha_beta(beta_old, y, n, m1_bar)
  
  # Sampling from a U(0,1)
  u <- runif(1) 
  
  # If u < aprob, the proposed value is accepted and the number of accepted values is increased 
  if (u < aprob){
    n_acc = n_acc+1
    
    # Return of the new value of beta (equal to the proposed value y) and the accuracy 
    return(list("beta" = y, "acc" = n_acc))
  } 
  else {
    # Return of the new value of beta (equal to the previous value beta_old) and the accuracy
    return(list("beta" = beta_old, "acc" = n_acc))
  }
}




# Function to change the value of beta so beta_star is in (-inf, +inf) to do MH
# INPUT: x -> value of beta, x is in (0,1)
# OUTPUT: x_star -> value of transformed beta, x_star is in (-inf, +inf)

change_beta <- function(x) {
  return (log(x) - log(1-x))
}

# Function to change the value of beta_star so beta is in (0, 1) after MH
# inv_beta is the inverse function of change_beta
# INPUT: x -> value of beta_star, x is in (-inf, +inf)
# OUTPUT: x_star -> value of beta, x_star is in (0, 1)

inv_beta <- function(x) {
  return (exp(x) / (exp(x) + 1))
}

# Function to compute the alpha(x,y) for the proposed value y and the old value x
# INPUT: x -> value of beta at the previous step
#        y -> proposed value for beta
#        n -> number of data points
#        m1_bar -> number of points in group 0
# OUTPUT: alpha -> the alpha needed to perform the acceptance/rejection in MH

calcolo_alpha_beta <- function(x, y, n, m1_bar) {
  # Computation of the partial posterior density for y and x
  pi_y <- dens_beta(y, n, m1_bar)
  pi_x <- dens_beta(x, n, m1_bar)
  
  rapp <- pi_y/pi_x
  
  # Computation of alpha 
  return (min(1, rapp))
}



# Function to compute the partial posterior density (up to the normalizing constant) for beta
# INPUT: x -> point where the partial density is evaluated
#        n -> number of data points 
#        m1_bar -> number of points in group 0
# OUTPUT: alpha -> the alpha needed to perform the acceptance/rejection in MH

dens_beta <- function(x, n, m1_bar) {
  return ( x^(n - m1_bar) * (1 - x)^(m1_bar) * 1/(x*(1-x)))
}


