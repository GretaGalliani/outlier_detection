### UPDATE CLUSTER
library(MASS)

#NOTA: cosa succede se due dati sono uguali? Devo esplicitamente dirgli di non
#samplarlo in j=0

#NEL MAIN: 
#Devo iniziare anche Xi_star (entrambi)

#Importing some auxiliary functions
source("Auxiliary_functions/Auxiliary_functions.R")

#r=1,...,R


#Input:
# Y <- data.frame
# Xi_mu_star <- lista di n vettori --> struttura che contiene medi gruppi
# Xi_cov_star <- lista di n matrici --> struttura che contiene matrici cov dei gruppi
# beta_old <- numero
# theta_old <- numero
# sigma_old <- numero
# S_old <- vettore cluster passo precedente
# k_old <- numero
#nu_0_P & nu_0_Q
#mu_0_P & mu_0_Q 
#k_0_P & k_0_Q
#lambda_0_P & lambda_0_Q

update_cluster <- function(Y, Xi_mu, Xi_cov, Xi_mu_star, Xi_cov_star, S_old,
                                     beta_old, theta_old, sigma_old, k_old,
                                    nu_0_P, k_0_P, mu_0_P, lambda_0_P,
                                    nu_0_Q, k_0_Q, mu_0_Q, lambda_0_Q)
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
    
    prob <- rep(0:k_new+2,0) #0,1,...k_new, k_new+1
    
    m1_bar <- m1_bar(curr[-i])
    #j=0
    prob[1] <- dens_contaminated(Y[i],beta_old, nu_0_P, k_0_P, mu_0_P, lambda_0_P) 
    
    #j=1:K
    for (t in (2:k_new+1)) 
    {
      n_j <- table(curr[-i])
      
      prob[t] <- dens_cluster_old(Y[i],n_j, m1_bar, sigma_old, theta_old, beta_old, Xi_mu_star[[t]],
                                 Xi_cov_star[[t]])
      
      }
    
    #j=K+1
    prob[k_new+2]<- dens_cluster_new(Y[i], beta_old, sigma_old, theta_old, m1_bar, k_new,
                                     nu_0_Q,mu_0_Q,k_0_Q,lambda_0_Q)
    
    #I sample the new assignment
    j <- sample(0:k_new+1,size=1,prob=prob)
    
    #If I sampled the same group, I can skip all the passages below
    if(j != curr[i])
    {
    
      #If i sampled an already existent group
      if(j %in% 1:k_new )
      {
        #I memorize the old group of the point
        old_group <- curr[i]
        #I assign the point to the j group
        curr[i] <- j
        
        #I assign the correct group parameters
        Xi_mu[[i]] <- Xi_mu_star[[j]]
        Xi_cov[[i]] <- Xi_cov_star[[j]]
        
        #I check if the old group is now empty
        if(!(old_group %in% curr))
          {
            #If it is, I need to shift all the groups' label and to cancel the old group's parameters
            curr <- shift(curr, old_group)
            Xi_mu_star <- Xi_mu_star[-old_group]
            Xi_cov_star <- Xi_cov_star[-old_group]
            k_new=k_new-1
          }
       }
     }
    
      #If I sampled a new group, I need to sample also the new parameters
      else if(j==k_new+1)
      {
          k_new=k_new+1
          
          curr[i] <- j
          
          #I sample from the predictive distribution considering only the current point 
          #The predictive distribution is known in closed form
          Xi_mu_star <- construct_mu_new(Y[i],Xi_mu_star, nu_0_Q,mu_0_Q,k_0_Q,lambda_0_Q)
          Xi_mu[[i]] <- Xi_mu_star[[j]]
          
          Xi_cov_star <- construct_cov_new(Y[i],Xi_mu_star, nu_0_Q,mu_0_Q,k_0_Q,lambda_0_Q)
          Xi_cov[[i]] <- Xi_cov_star[[j]]
          
      }
    
      else
      #If no one of the cases above is met, then j=0 and I do the assignment
      {
        curr[i] <- j
      }
   } 
  
  S_new <- curr
  my_list <-list("S_new"=S_new, "Xi_mu"=Xi_mu,"Xi_cov"=Xi_cov,
                            "Xi_mu_star"=Xi_mu_star,"Xi_cov_star"=Xi_cov_star)
  return (my_list) 
}



dens_contaminated <- function(data, beta_old, nu_0_P, k_0_P, mu_0_P, lambda_0_P)
{
  p <- dim(data)[1]
  x_medio <- mean(data)
  S = sum((data[i]-x_medio)*t((data[i]-x_medio)))
  
  weight <- (1-beta_old) *rmvt(1,mu= mu_0_P, sigma = (lambda_0_P * (k_0_P+1))/(k_0*(nu_0_P-p+1)), 
                                  df = nu_0_P-p+1)  
  
  return (weight)
}


dens_cluster_old <- function(data, n_j, m1_bar, sigma_old, theta_old, beta_old, Xi_mu_star, Xi_cov_star)
{
  
  coeff <- beta_old * (n_j-sigma_old)/(theta_old+N-m1_bar-1)
  
  weight <- coeff * mvrnorm(1, mean=Xi_mu_star, varcov=Xi_cov_star)
  return (weight)
}


dens_cluster_new <- function(data, beta_old, sigma_old, theta_old, m1_bar, k_old,
                 nu_0_Q,mu_0_Q,k_0_Q,lambda_0_Q)
{
  coeff <- beta_old* (theta_old+(k_old-m1_bar)*sigma_old)/(theta_old+n-m1_bar-1)
  
  p <- dim(data)[1]
  
  weight <- coeff *rmvt(1, mu = mu_0, sigma = (lambda_0_Q * (k_0_Q+1))/(k_0_Q*(nu_0_Q-p+1)), 
                                  df = nu_n-p+1) 
  return (weight)
}


shift <- function(curr, old_group)
{
  for(i in (1:dim(curr)[1]))
  {
    if(curr[i]>old_group)
    {
      curr[i]=curr[i]-1
    }
  }
  
  return(curr)
}

construct_mu_new <- function(data,Xi_mu_star, nu_0_Q,mu_0_Q,k_0_Q,lambda_0_Q)
{
  p <- dim(data)[2]
  nu_n <- nu_0_Q + 1
  k_n <- k_0_Q + 1
  mu_n <- k_0_Q/(k_n) * mu_0_Q + 1/k_n * data 
  lambda_n <- lambda_0_Q + (k_0_Q)/(k_n)*(data-mu_0_Q)*t(data-mu_0_Q)
  sigma = (lambda_n * 1)/(k_n*(nu_n-p+1))
                          
  mu_new <- rmvt(1, mu = mu_n, sigma = sigma, df = nu_n-p+1)
  Xi_mu_star<-append(Xi_mu_star,mu_new)
  return (Xi_mu_star)
}

construct_cov_new <- function(data,Xi_mu_star, nu_0_Q,mu_0_Q,k_0_Q,lambda_0_Q)
{
  p <- dim(data)[2]
  nu_n <- nu_0_Q + 1
  k_n <- k_0_Q + 1
  
  lambda_n <- lambda_0_Q + (k_0_Q)/(k_n)*(data-mu_0_Q)*t(data-mu_0_Q)
  sigma = (lambda_n * 1)/(k_n*(nu_n-p+1))
  
  cov_inv <- inv(lambda_n)                        
  cov_new <-  rinvwishart(nu_n, cov_inv)
  
  Xi_cov_star<-append(Xi_cov_star,cov_new)
  return (Xi_cov_star)
}


