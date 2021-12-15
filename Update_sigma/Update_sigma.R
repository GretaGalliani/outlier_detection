##### UPDATE SIGMA

# Function to update the beta parameter
# INPUT: m1 -> number of singletons (number of points in group 0 + number of groups with one point)
#                                   (to be computed using the auxiliary function m1)
#        m1_bar -> number of points in group 0 (to be computed using the auxiliary function m1_bar)
#        k -> number of distinct groups
#        sigma_old -> beta at the previous iteration
#        theta -> theta at the current iteration
#        freq -> vector containing the number of points for each group j=1,...k
#        sd -> standard deviation for the proposal density in the MHRW (default value is 2)
#        n_acc -> number of accepted proposals until the previous iterations
#        a -> first shape parameter for the beta prior for sigma
#        b -> second shape parameter for the beta prior for sigma
#        sd -> standard deviation for the proposal density in the MHRW (default value is 2)

# OUTPUT: sigma -> value of the parameter sigma at the r iteration
#         acc -> number of accepted proposals at the current iteration
update_sigma <- function(m1, m1_bar, k, sigma_old, theta, freq, n_acc, a, b, sd = 2) { 
  # METROPOLIS HASTINGS RANDOM WALK 
  
  # Extraction of a new value from the proposal distribution, doing an appropriate transformation 
  # to correct the fact that sigma is in (0,1)
  y <- inv_sigma( change_sigma(sigma_old) + rnorm(1,0,sd))
  
  # Computation the alpha of the new proposal wrt the old one 
  aprob <- compute_alpha_sigma(sigma_old, y, k, m1, m1_bar, theta, freq, a, b)
  #print("Aprob sigma")
  #print(aprob)
  
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
#        k -> number of distinct groups
#        m1 -> number of singletons (number of points in group 0 + number of groups with one point)
#        m1_bar -> number of points in group 0
#        theta -> theta at the current iteration
#        freq -> vector containing the number of points for each group j=1,...k
# OUTPUT: alpha -> the alpha needed to perform the acceptance/rejection in MH
compute_alpha_sigma <- function(x, y, k, m1, m1_bar, theta, freq, a, b) {
  # Computation of the partial posterior density for y and x
  pi_y <- dens_sigma(y, k, m1, m1_bar, theta, freq, a, b)
  pi_x <- dens_sigma(x, k, m1, m1_bar, theta, freq, a, b)
  
  rapp <- pi_y/pi_x
  
  # Computation of alpha 
  return (min(1, rapp))
}


# Function to compute the partial posterior density (up to the normalizing constant) for sigma
# INPUT: x -> point where the partial density is evaluated
#        k -> number of distinct groups
#        m1 -> number of singletons (number of points in group 0 + number of groups with one point)
#        m1_bar -> number of points in group 0
#        theta -> theta at the current iteration
#        freq -> vector containing the number of points for each group j=1,...k
#        m1_bar -> number of points in group 0
# OUTPUT: f -> the evaluation of the partial posterior density at point x
dens_sigma <- function(x, k, m1, m1_bar, theta, freq, a, b) {
  
  # I select the groups which are not singletons
  single = TRUE
  for (elem in freq){
    if (elem > 1){
      single = FALSE
      break
    }
  }
  
  print("freq")
  print(freq)
  
  print("k")
  print(k)
  
  print("m1")
  print(m1)
  
  print("m1_bar")
  print(m1_bar)
  
  if(single)
    return ( dbeta(x, a, b) * x^(k) * gamma(theta/x + k - m1_bar)/gamma(theta/x) * 1/(x*(1-x)))
  else{
    freq_m1 = freq[freq>1]
    
    print("vettore produttoria")
    print(gamma(freq_m1-x)/gamma(1-x))
    
    return ( dbeta(x, a, b) * x^(k) * (gamma(theta/x + k)/gamma(theta/x)) * prod(gamma(freq_m1-x)/gamma(1-x)) * 1/(x*(1-x)))
  }  
}

