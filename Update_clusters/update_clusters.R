### UPDATE CLUSTERS
library(MASS)
library(mvtnorm)
library(LaplacesDemon)

#NOTA: cosa succede se due dati sono uguali? Devo esplicitamente dirgli di non
#samplarlo in j=0?

#Importing some auxiliary functions
source("auxiliary_functions/auxiliary_functions.R")

# Function to update the clusters
# INPUT: Y -> data 
#        xi_mu_star -> list of k vectors, contains means of the groups
#        xi_cov_star -> list of k matrices, contains cov. matrices of the groups
#        beta_old -> beta at previous iteration
#        theta_old -> theta at previous iteration
#        sigma_old -> sigma at previous iteration
#        S_old -> n vector containing the clusters at previous iteration
#        k_old -> number of distinct groups
#        P_param -> list of parameters nu_0_P, mu_0_P, k_0_P, lambda_0_P
#        Q_param -> list of parameters nu_0_Q, mu_0_Q, k_0_Q, lambda_0_Q


# OUTPUT: list of the following:
#         S_new -> n vector containing the new clusters
#         xi_mu_star -> list of k vectors, contains means of the new groups
#         xi_cov_star -> list of k matrices, contains cov. matrices of the new groups

update_clusters <- function(Y, xi_mu_star, xi_cov_star, S_old,
                              beta_old, theta_old, sigma_old, k_old,
                              P_param, Q_param)
{
  n <- dim(Y)[1]
  
  # Creation of the new S and k initialized like the ones at previous iteration
  curr <- S_old 
  k_new <- k_old
  
  # Cycling over all the elements
  for (i in 1:n)
  {
    # Computation of the probability that the point 
    # - comes from the contaminated distribution --> j=0
    # - comes from an already existing group --> j=1:K
    # - comes from a new group --> j=K+1
    
    # Initialization of all probabilities to 0
    prob <- rep(0,k_new+2) #0,1,...k_new, k_new+1
    
    # Computation of m_bar discarding the i-th element
    m1_bar <- m1_bar(curr[-i])
    
    #j=0
    p = dim(Y)[2]
    
    # Probability that the i-th element is in the contaminated part
    prob[1] <- dens_contaminated(Y[i,], beta_old, P_param) 
    
    # print("prob contaminated")
    # print(prob[1])
    
    
    #j=1:K
    #note: t goes from 2 to k_new+1 because group j will be in position j+1 of the prob 
    #vector as the position 1 is occupied by the probability of group j being 0
    for (t in (2:(k_new+1))) 
    {
      # If the current group is empty the probability of sampling that group is 0
      if (!((t-1) %in% curr)){
        prob[t]=0
      }
      
      else
      {
        # Frequency of the current group
        n_j <- sum(curr == t-1)
      
        # print("n_j con i")
        # print(n_j)
        #print("gruppo")
        #print(t-1)
        # 
        # print("curr[i] è")
        # print(curr[i])
        # 
        # print("tutto curr è")
        # print(curr)
        # 
        
        # If data point i is the only point in the current group, then the probability 
        # of sampling that group is 0
        if (n_j == 1 & curr[i] == t-1){
          prob[t]=0
        }
        
        else{
          
          # Frequency of the current group without point i 
          n_j <- sum(curr[-i] == t-1)
          
          # Computation of the probability that data i is sampled in group j
          prob[t] <- dens_cluster_old(Y[i,],n_j, n, m1_bar, sigma_old, theta_old, beta_old, xi_mu_star[[t-1]],
                                      xi_cov_star[[t-1]])
          
          
        }
      }
      
      # print("gruppo")
      # print(t-1)
      # print( "prob gruppo ")
      # print(prob[t])
    }
    
    #j=K+1
    
    # Computation of the probability that data i is sampled in a new group 
    prob[k_new+2]<- dens_cluster_new(Y[i,], n, beta_old, sigma_old, theta_old, m1_bar, k_new,
                                     Q_param)
    
    # print("prob gruppo nuovo")
    # print(prob[k_new+2])
    
    
    print("vector of prob")
    print(prob)
    
    # Sampling of the new assignment
    j <- sample(1:(k_new+2),size=1,prob=prob) - 1
    
    
    # print("sampling")
    # print(j)
    
    # print("beta_old")
    # print(beta_old)
    
    curr[i] <- j
    
    
    # If a new group is sampled, then the new parameters also need to be sampled
    if(j==k_new+1)
      {
          # A new group has been formed 
          k_new=k_new+1
          
          # Sampling from the predictive distribution considering only the current point 
          # The predictive distribution is known in closed form
          xi_mu_star <- construct_mu_new(Y[i,],xi_mu_star, Q_param)
          
          xi_cov_star <- construct_cov_new(Y[i,],xi_cov_star, Q_param)
          
      }
    
      
   } 

  
  # At the end all the empty groups are deleted and the others are shifted
  delete_list <- delete_and_shift(curr, xi_mu_star, xi_cov_star)
  curr <- delete_list$curr
  xi_mu_star <- delete_list$xi_mu_star
  xi_cov_star <- delete_list$xi_cov_star
  
  
  print("curr")
  print(curr)
  
  print("beta")
  print(beta_old)
  
  print("sigma")
  print(sigma_old)
  
  print("theta")
  print(theta_old)
  
  S_new <- curr
  my_list <-list("S_new"=S_new, "xi_mu_star"=xi_mu_star,"xi_cov_star"=xi_cov_star)
  return (my_list) 
}


# Function to compute the probability that a point belongs to the contaminated part
# INPUT: data -> datapoint 
#        p -> dimension of the multivariate normal 
#        beta_old -> beta at previous iteration
#        P_param -> list of parameters nu_0_P, mu_0_P, k_0_P, lambda_0_P

# OUTPUT: weight -> the computed probability

dens_contaminated <- function(data, beta_old, P_param)
{
  data <- as.matrix(data)
  data <- t(data)
  p = dim(data)[2]
  
  mu_n = as.matrix(P_param$k_0/(P_param$k_0+1)*P_param$mu_0 + (1/P_param$k_0 + 1)*data)
  
  k_n = P_param$k_0 + 1
  nu_n = P_param$nu_0 + 1
  
  r <- as.matrix(data-P_param$mu_0)
  
  lambda_n = P_param$lambda_0 + (P_param$k_0/(P_param$k_0 + 1))*t(r)%*%r
  
  # Evaluation of a multivariate t-Student
  weight <- (1-beta_old) *LaplacesDemon::dmvt(data,mu= mu_n, S = lambda_n*(k_n+1)/(k_n*(nu_n-p+1)), 
                                              df = nu_n-p+1)
  
  return (weight)
}

# Function to compute the probability that a point belongs to the current group
# INPUT: data -> datapoint 
#        n_j -> frequency of the current group without data i
#        n -> number of data points
#        m1_bar -> number of points in group j=0
#        sigma_old -> sigma at previous iteration
#        theta_old -> theta at previous iteration
#        beta_old -> beta at previous iteration
#        xi_mu_act -> mean vector of the current group 
#        xi_cov_act -> cov. matrix of the current group

# OUTPUT: weight -> the computed probability

dens_cluster_old <- function(data, n_j, n, m1_bar, sigma_old, theta_old, beta_old, xi_mu_act, xi_cov_act)
{

  coeff <- beta_old * (n_j-sigma_old)/(theta_old+n-m1_bar-1)
  
  # Computation of the density of a multivariate normal 
  weight <- coeff * dmvnorm(data, mean=xi_mu_act, sigma=xi_cov_act)
  
  return (weight)
}

# Function to compute the probability that a point belongs to a new group 
# INPUT: data -> datapoint 
#        p -> dimension of the multivariate normal 
#        n -> number of data points
#        beta_old -> beta at previous iteration
#        sigma_old -> sigma at previous iteration
#        theta_old -> theta at previous iteration
#        m1_bar -> number of points in group j=0
#        k_old -> number of distinct groups
#        Q_param -> list of parameters nu_0_Q, mu_0_Q, k_0_Q, lambda_0_Q

# OUTPUT: weight -> the computed probability

dens_cluster_new <- function(data, n, beta_old, sigma_old, theta_old, m1_bar, k_old, Q_param)
{
  data <- as.matrix(data)
  p = dim(data)[2]
  
  data <- t(data)
  
  mu_n = as.matrix(Q_param$k_0/(Q_param$k_0+1)*Q_param$mu_0 + (1/Q_param$k_0 + 1)*data)
  k_n = Q_param$k_0 + 1
  nu_n = Q_param$nu_0 + 1
  
  r <- as.matrix(data-Q_param$mu_0)
  
  lambda_n = Q_param$lambda_0 + (Q_param$k_0/(Q_param$k_0 + 1))*t(r)%*%r
  
  coeff <- beta_old * (theta_old+(k_old-m1_bar)*sigma_old)/(theta_old+n-m1_bar-1)
  
  # Evaluation of a multivariate t-Student
  weight <- coeff *LaplacesDemon::dmvt(data,mu= mu_n, S = lambda_n*(k_n+1)/(k_n*(nu_n-p+1)), 
                                                        df = nu_n-p+1)  
  return (weight)
}

# Function to delete the empty groups and shift the other groups when needed
# INPUT: curr -> n vector containing the clusters at current iteration 
#        xi_mu_star -> list of k vectors, contains means of the current groups
#        xi_cov_star -> list of k matrices, contains cov. matrices of the current groups

# OUTPUT: curr -> n vector containing the updated clusters
#         xi_mu_star -> list of k vectors, contains means of the updated groups
#         xi_cov_star -> list of k matrices, contains cov. matrices of the updated groups

delete_and_shift <- function(curr, xi_mu_star, xi_cov_star)
{
  # end is the highest possible group index
  end <- max(curr)
  
  # Iteration from 1 to the max of curr 
  j = 1
  while (j <= end)
  {
    # If the j-th group is empty
    if (j<=max(curr) & !(j %in% curr)){
      
      # Elimination of group parameters
      xi_mu_star <- xi_mu_star[-j]
      xi_cov_star <- xi_cov_star[-j]
      
      # Shifting of groups with higher index
      for (i in 1:length(curr)){
        if (curr[i]>j){
          curr[i] = curr[i]-1
        }
        
      }
    }
    
    # If the current group is deleted, j is not augmented since the group needs to checked
    # again (after the shifting the current group corresponds to the next one)
    # If the is no deletion j is augmented to check the next group
    else
    {
      j = j + 1
    }
  }
  return(list("curr" = curr, "xi_mu_star" = xi_mu_star, "xi_cov_star" = xi_cov_star))
}

# Function to sample the parameter mu for a new group 
construct_mu_new <- function(data,xi_mu_star, Q_param)
{
  data <- as.matrix(data)
  p <- dim(data)[2]
  nu_n <- Q_param$nu_0 + 1
  k_n <- Q_param$k_0 + 1
  mu_n <- as.vector(Q_param$k_0/(k_n) * Q_param$mu_0 + 1/k_n * data)
  
  r <- as.matrix(data-Q_param$mu_0)

  lambda_n <- Q_param$lambda_0 + (Q_param$k_0)/(k_n)*r%*%t(r)
  
  sigma = (1/(k_n*(nu_n-p+1)))*(lambda_n)
                    
  mu_new <- LaplacesDemon::rmvt(mu = mu_n, S = sigma, df = nu_n-p+1)
  
  mu_new <- as.vector(mu_new)
  
  xi_mu_star<-append(xi_mu_star,list(mu_new))
  return (xi_mu_star)
}

construct_cov_new <- function(data,xi_cov_star, Q_param)
{
  data <- as.matrix(data)
  p <- dim(data)[2]
  
  nu_n <- Q_param$nu_0 + 1
  k_n <- Q_param$k_0 + 1
  
  r <- as.matrix(data-Q_param$mu_0)
  
  lambda_n <- Q_param$lambda_0 + (Q_param$k_0)/(k_n)*r%*%t(r)
  
  sigma = (1/(k_n*(nu_n-p+1)))*(lambda_n)
  
  # cov_new <-  LaplacesDemon::rinvwishartc(nu_n, chol2inv(chol(sigma)))
  
  cov_new <- LaplacesDemon::rinvwishart(nu_n, as.inverse(sigma))
  
  xi_cov_star<-append(xi_cov_star,list(cov_new))
  return (xi_cov_star)
}


