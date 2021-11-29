##### UPDATE THETA

# Function to update the theta parameter
update_theta <- function(n, m1_bar, k, theta_old, sigma, sd = 2) { 
  y <- inv_theta( change_theta(theta_old) + rnorm(1,0,sd))
  print(y)
  aprob <- calcolo_alpha_theta(theta_old, y, k, m1_bar, sigma, n)
  print(aprob)
  u <- runif(1) # If condition "(u < aprob)" is NOT met, we'll skip command "x <- y", #  so that the MC does not move from x                 
  if (u < aprob){
    return(y)
  } 
  else {
    return(theta_old)
  }
}


# Function to change the value of theta so theta_star is in (-inf, +inf) to do MH
# FUNZIONA
change_theta <- function(x) {
  return (log(x))
}

# Function to change the value of theta_star so theta is in (0, +inf) after MH
# FUNZIONA
inv_theta <- function(x) {
  return (exp(x))
}

# Function to compute the alpha(x,y) for theta
# We assume theta ~ U(0,1)
# DA TESTARE
calcolo_alpha_theta <- function(x, y, k, m1_bar, sigma, n) {
  pi_y <- dens_theta(y, k, m1_bar, sigma, n)
  pi_x <- dens_theta(x, k, m1_bar, sigma, n)
  rapp <- pi_y/pi_x
  return (min(1, rapp))
}


# Function to compute the partial posterior density (up to the normalizing constant) for theta
# INPUT: x (which is already inverted so that x is in (0,+inf))
# DA TESTARE
dens_theta <- function(x, k, m1_bar, sigma, n) {
  return ( gamma(x) * gamma(x/sigma + k - m1_bar) / gamma(x/sigma) / gamma(x + n - m1_bar) * (1/x) )
}




