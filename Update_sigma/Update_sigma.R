##### UPDATE SIGMA

# MAIN IMPLEMENTATION: pass as input for sigma_old and n_acc the element of
# the list using $

# Function to update the sigma parameter
# NB input variable freq is a vector containing the frequencies of the unique
# values xi* for i=1,...,k
update_sigma <- function(m1, m1_bar, k, sigma_old, theta, freq, sd = 2, n_acc) { 
  y <- inv_sigma( change_sigma(sigma_old) + rnorm(1,0,sd))
  print(y)
  aprob <- calcolo_alpha_sigma(sigma_old, y, k, m1, m1_bar, theta, freq)
  print(aprob)
  u <- runif(1) # If condition "(u < aprob)" is NOT met, we'll skip command "x <- y", #  so that the MC does not move from x                 
  if (u < aprob){
    n_acc = n_acc+1
    return(list("sigma" = y, "acc" = n_acc))
  } 
  else {
    return(list("sigma" = sigma_old, "acc" = n_acc))
  }
}


# Function to change the value of sigma so sigma_star is in (-inf, +inf) to do MH
# FUNZIONA
change_sigma <- function(x) {
  return (log(x) - log(1-x))
}

# Function to change the value of sigma_star so sigma is in (0, 1) after MH
# FUNZIONA
inv_sigma <- function(x) {
  return (exp(x) / (exp(x) + 1))
}

# Function to compute the alpha(x,y) for sigma
# We assume sigma ~ U(0,1)
# DA TESTARE
calcolo_alpha_sigma <- function(x, y, k, m1, m1_bar, theta, freq) {
  pi_y <- dens_sigma(y, k, m1, m1_bar, theta, freq)
  pi_x <- dens_sigma(x, k, m1, m1_bar, theta, freq)
  rapp <- pi_y/pi_x
  return (min(1, rapp))
}


# Function to compute the partial posterior density (up to the normalizing constant) for sigma
# INPUT: x (which is already inverted so that x is in (0,1))
# DA TESTARE
dens_sigma <- function(x, k, m1, m1_bar, theta, freq) {
  freq_m1 = freq[freq>1]
  return ( x^(k - m1_bar) * (gamma(theta/x + k - m1_bar)/gamma(theta/x))^(k-m1) * prod(gamma(freq_m1-x)/gamma(1-x)) * 1/(x*(1-x)))
}


