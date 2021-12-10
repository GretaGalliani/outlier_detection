### UPDATE CLUSTERS
library(MASS)
library(mvtnorm)

#NOTA: cosa succede se due dati sono uguali? Devo esplicitamente dirgli di non
#samplarlo in j=0?

#Importing some auxiliary functions
source("auxiliary_functions/auxiliary_functions.R")

#r=1,...,R

# Function to update the clusters
# INPUT: Y <- data.frame

#Input:
# Y <- data.frame
# xi_mu_star <- list of n vectors --> contains means of the groups
# xi_cov_star <- list of n matrices --> contains cov. matrices of the groups
# beta_old <- number
# theta_old <- number
# sigma_old <- number
# S_old <- vettore cluster passo precedente
# k_old <- numero
#nu_0_P & nu_0_Q
#mu_0_P & mu_0_Q 
#k_0_P & k_0_Q
#lambda_0_P & lambda_0_Q

update_clusters <- function(Y, xi_mu_star, xi_cov_star, S_old,
                              beta_old, theta_old, sigma_old, k_old,
                              P_param, Q_param)
{
  N <- dim(Y)[1]
  curr <- S_old 

  k_new <- k_old
  
  for (i in 1:N)
  {
    #We compute the probability that the point 
    # - comes from the contaminated distribution --> j=0
    # - comes from an already existing group --> j=1:K
    # - comes from a new group --> j=K+1
    
    prob <- rep(0,k_new+2) #0,1,...k_new, k_new+1
    
    m1_bar <- m1_bar(curr[-i])
    
    #j=0
    p = dim(Y)[2]
    
    prob[1] <- dens_contaminated(Y[i,], p,beta_old, P_param) 
    
    print("curr ")
    print(curr)
    #j=1:K
    for (t in (2:(k_new+1))) 
    {
      n_j <- as.integer(table(curr)[t-1])
      print("n_j")
      print(n_j)
      print("curr_i ")
      print (curr[i])
      if (n_j == 1 & curr[i] == t-1){
        prob[t-1]=0
      }
      else{
        n_j <- as.integer(table(which(curr[-i]==t-1)))
        n <- dim(data)[1]
        
        prob[t] <- dens_cluster_old(Y[i,],n_j, n, m1_bar, sigma_old, theta_old, beta_old, xi_mu_star[[t-1]],
                                   xi_cov_star[[t-1]])
      }
    }
    
    #j=K+1
    prob[k_new+2]<- dens_cluster_new(Y[i,], p, N, beta_old, sigma_old, theta_old, m1_bar, k_new,
                                     Q_param)
    
    #I sample the new assignment
    print(prob)
    j <- sample(0:(k_new+1),size=1,prob=prob)
    
    curr[i] <- j
    
    
      #If I sampled a new group, I need to sample also the new parameters
    if(j==k_new+1)
      {
          k_new=k_new+1
          
          curr[i] <- j
          
          #I sample from the predictive distribution considering only the current point 
          #The predictive distribution is known in closed form
          xi_mu_star <- construct_mu_new(Y[i],xi_mu_star, Q_param)
          
          xi_cov_star <- construct_cov_new(Y[i],xi_mu_star, Q_param)
          
      }
    
      
   } 
  
  # I delete and shift all the empty groups
  delete_list <- delete_and_shift(curr, xi_mu_star, xi_cov_star)
  
  curr <- delete_list$curr
  xi_mu_star <- delete_list$xi_mu_star
  xi_cov_star <- delete_list$xi_cov_star
  
  S_new <- curr
  my_list <-list("S_new"=S_new, "xi_mu_star"=xi_mu_star,"xi_cov_star"=xi_cov_star)
  return (my_list) 
}



dens_contaminated <- function(data, p, beta_old, P_param)
{
  weight <- (1-beta_old) *LaplacesDemon::dmvt(data,mu= P_param$mu_0, S = ((P_param$lambda_0 * (P_param$k_0+1))/(P_param$k_0*(P_param$nu_0-p+1))), 
                                  df = P_param$nu_0-p+1)  
  return (weight)
}


dens_cluster_old <- function(data, n_j, n, m1_bar, sigma_old, theta_old, beta_old, xi_mu_act, xi_cov_act)
{
  coeff <- beta_old * (n_j-sigma_old)/(theta_old+n-m1_bar-1)
  
  
  weight <- coeff * dmvnorm(data, mean=xi_mu_act, sigma=xi_cov_act)
  return (weight)
}


dens_cluster_new <- function(data, p, n, beta_old, sigma_old, theta_old, m1_bar, k_old, Q_param)
{
  coeff <- beta_old* (theta_old+(k_old-m1_bar)*sigma_old)/(theta_old+n-m1_bar-1)
  
  weight <- coeff *LaplacesDemon::dmvt(data, mu = Q_param$mu_0, S = (Q_param$lambda_0 * (Q_param$k_0+1))/(Q_param$k_0*(Q_param$nu_0-p+1)), 
                                  df = Q_param$nu_0-p+1) 
  return (weight)
}


delete_and_shift <- function(curr, xi_mu_star, xi_cov_star)
{
  end <- max(curr)
  # Iteration from 1 to the max of curr 
  j = 1
  while (j <= end)
  {
    # If the j-th group is empty
    if (j<=max(curr) & !(j %in% curr)){
      xi_mu_star <- xi_mu_star[-j]
      xi_cov_star <- xi_cov_star[-j]
      for (i in 1:length(curr)){
        if (curr[i]>j){
          curr[i] = curr[i]-1
        }
        
      }
    }
    else
    {
      j = j + 1
    }
  }
  }
  
  return(list("curr" = curr, "xi_mu_star" = xi_mu_star, "xi_cov_star" = xi_cov_star))
}

construct_mu_new <- function(data,xi_mu_star, Q_param)
{
  p <- dim(data)[2]
  nu_n <- Q_param$nu_0 + 1
  k_n <- Q_param$k_0 + 1
  mu_n <- Q_param$k_0/(k_n) * Q_param$mu_0 + 1/k_n * data 
  lambda_n <- Q_param$lambda_0 + (Q_param$k_0)/(k_n)*(data-Q_param$mu_0)*t(data-Q_param$mu_0)
  sigma = (lambda_n * 1)/(k_n*(nu_n-p+1))
                          
  mu_new <- rmvt(1, mu = mu_n, sigma = sigma, df = nu_n-p+1)
  xi_mu_star<-append(xi_mu_star,mu_new)
  return (xi_mu_star)
}

construct_cov_new <- function(data,xi_mu_star, Q_param)
{
  p <- dim(data)[2]
  nu_n <- Q_param$nu_0 + 1
  k_n <- Q_param$k_0 + 1
  
  lambda_n <- Q_param$lambda_0 + (Q_param$k_0)/(k_n)*(data-Q_param$mu_0)*t(data-Q_param$mu_0)
  sigma = (lambda_n * 1)/(k_n*(nu_n-p+1))
  
  lambda_n_chol <- chol(lambda_n)
  inv_lambda_n_chol <- chol2inv(lambda_n_chol)
  
  cov_inv <- inv(lambda_n)                        
  cov_new <-  rinvwishartc(nu_n, chol(inv_lambda_n_chol))
  
  xi_cov_star<-append(xi_cov_star,cov_new)
  return (xi_cov_star)
}


