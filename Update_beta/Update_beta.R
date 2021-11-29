##### UPDATE BETA

# Function to update the beta parameter
# DA TESTARE
update_beta <- function(n, m1_bar, beta_old, sd = 2) { 
  y <- inv_beta( change_beta(beta_old) + rnorm(1,0,sd))
  print(y)
  aprob <- calcolo_alpha_beta(beta_old, y, n, m1_bar)
  print(aprob)
  u <- runif(1) # If condition "(u < aprob)" is NOT met, we'll skip command "x <- y", #  so that the MC does not move from x                 
  if (u < aprob){
    #acc <- acc + 1
    #my_list <- list("beta" = y, "acc" = acc)
    #return(my_list) 
    return(y)
  } 
  else {
    # my_list <- list("beta" = beta_old, "acc" = acc)
    # return(my_list)
    return(beta_old)
  }
}


# Function to change the value of beta so beta_star is in (-inf, +inf) to do MH
# FUNZIONA
change_beta <- function(x) {
  return (log(x) - log(1-x))
}

# Function to change the value of beta_star so beta is in (0, 1) after MH
# FUNZIONA
inv_beta <- function(x) {
  return (exp(x) / (exp(x) + 1))
}

# Function to compute the alpha(x,y) for beta
# We assume beta ~ U(0,1)
# DA TESTARE
calcolo_alpha_beta <- function(x, y, n, m1_bar) {
  pi_y <- dens_beta(y, n, m1_bar)
  pi_x <- dens_beta(x, n, m1_bar)
  rapp <- pi_y/pi_x
  return (min(1, rapp))
}



# Function to compute the partial posterior density (up to the normalizing constant) for beta
# INPUT: x (which is already inverted so that x is in (0,1))
# DA TESTARE
dens_beta <- function(x, n, m1_bar) {
  return ( x^(n - m1_bar) * (1 - x)^(m1_bar) * 1/(x*(1-x)))
}


