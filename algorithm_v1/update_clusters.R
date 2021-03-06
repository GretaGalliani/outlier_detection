### UPDATE CLUSTERS
library(MASS)
library(mvtnorm)
library(LaplacesDemon)



#Importing some auxiliary functions
source("algorithm_v1/auxiliary_functions.R")

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
                            beta_old, theta_old, sigma_old,
                            P_param, Q_param)
{
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  
  # Creation of the new S and k initialized like the ones at previous iteration
  curr <- S_old
  
  # Cycling over all the elements
  for (i in 1:n)
  {
    # Computation of the probability that the point 
    # - comes from the contaminated distribution --> j=0
    # - comes from an already existing group --> j=1:K
    # - comes from a new group --> j=K+1
    
    #print(paste0("I'm considering sample ", i))
    
    
    # Initialization of all probabilities to 0
    prob <- rep(0,max(curr)+2) #0,1,...k_new, k_new+1
    
    # Computation of m_bar discarding the i-th element
    m1_bar <- m1_bar(curr[-i])
    
    # Dimension of the multivariate data
    p = dim(Y)[2]
    
    #j=0
    # Probability that the i-th element is in the contaminated part
    prob[1] <- dens_contaminated(Y[i,], beta_old, P_param, p) 
    
    
    if(max(curr) > 0){
      for (t in 2:(max(curr)+1))
      {
        # Frequency of the current group
        n_j <- sum(curr == t-1)
        
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
    }
    #   
    
    
    
    #j=K+1
    # Computation of the probability that data i is sampled in a new group 
    prob[max(curr)+2]<- dens_cluster_new(Y[i,], n, beta_old, sigma_old, theta_old, m1_bar, max(curr[-i]),
                                         Q_param, p)
    

    
    # I save the old group
    old_group <- curr[i]
    
    
    # Sampling of the new assignment
    j <- sample(0:(max(curr)+1),size=1,prob=prob)
    
    old_k <- max(curr)
    
    # Assignment to the new group
    curr[i] <- j
    
    
    
    # If a new group is sampled, then the new parameters also need to be sampled
    if(j==old_k+1)
    {
      
      # Sampling from the predictive distribution considering only the current point 
      # The predictive distribution is known in closed form
      xi_mu_star <- construct_mu_new(Y[i,],xi_mu_star, Q_param)
      xi_cov_star <- construct_cov_new(Y[i,],xi_cov_star, Q_param)
      
    }
    
    # If the old group is now empty
    if (sum(curr==old_group)==0 & old_group!=0){
      # I call the function which delete the groups' parameters and shift the groups higher than the old one
      #print(paste0("I'm deleting the group ", old_group))
      delete_list <- delete_and_shift(curr, xi_mu_star, xi_cov_star, old_group)
      curr <- delete_list$curr
      xi_mu_star <- delete_list$xi_mu_star
      xi_cov_star <- delete_list$xi_cov_star
      
      #print(paste0("Ho cancellato il gruppo ", old_group, " a cui apparteneva il sample ", i))
    }
    
    
  }
  

  # At the end all the empty groups are deleted and the others are shifted
  # delete_list <- delete_and_shift(curr, xi_mu_star, xi_cov_star)
  # curr <- delete_list$curr
  # xi_mu_star <- delete_list$xi_mu_star
  # xi_cov_star <- delete_list$xi_cov_star
  
   
 
   
  S_new <- curr
  return (S_new) 
}



# Function to compute the probability that a point belongs to the contaminated part
# INPUT: data -> datapoint 
#        p -> dimension of the multivariate normal 
#        beta_old -> beta at previous iteration
#        P_param -> list of parameters nu_0_P, mu_0_P, k_0_P, lambda_0_P

# OUTPUT: weight -> the computed probability

dens_contaminated <- function(data, beta_old, P_param, p)
{
  
  # Evaluation of a multivariate t-Student
  weight <- log(1-beta_old) + LaplacesDemon::dmvt(data,mu = P_param$mu_0, 
                                              S =(P_param$k_0 + 1)/(P_param$k_0*(P_param$nu_0-p+1))*P_param$lambda_0, 
                                              df = P_param$nu_0-p+1,
                                              log = TRUE)
  
  
  return (exp(weight))
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
  
  coeff <- log(beta_old) + log(n_j-sigma_old) - log(theta_old+n-m1_bar-1)
  
  # Computation of the density of a multivariate normal 
  weight <- coeff + dmvnorm(data, mean=xi_mu_act, sigma=xi_cov_act, log = TRUE)
  
  return (exp(weight))
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

dens_cluster_new <- function(data, n, beta_old, sigma_old, theta_old, m1_bar, k_old, Q_param, p)
{
  # Transformation of data into a matrix

  # Compute mu_n as matrix
  # mu_n = as.matrix(Q_param$k_0/(Q_param$k_0+1)*Q_param$mu_0 + (1/(Q_param$k_0 + 1))*data)
  
  # Compute k_n and nu_n
  # k_n = Q_param$k_0 + 1
  # nu_n = Q_param$nu_0 + 1
  
  # r <- as.matrix(data-Q_param$mu_0)
  
  # Compute lambda_n
  # lambda_n = Q_param$lambda_0 + (Q_param$k_0/(Q_param$k_0 + 1))*t(r)%*%r
  
  coeff <- log(beta_old) + log(theta_old+k_old*sigma_old) - log(theta_old+n-m1_bar-1)
  

  # Evaluation of a multivariate t-Student
  weight <- coeff + LaplacesDemon::dmvt(data,mu = Q_param$mu_0, 
                                       S = (Q_param$k_0 + 1)/(Q_param$k_0*(Q_param$nu_0-p+1))*Q_param$lambda_0, 
                                       df = Q_param$nu_0-p+1,
                                       log = TRUE)  
  return (exp(weight))
}

# Function to delete the empty groups and shift the other groups when needed
# INPUT: curr -> n vector containing the clusters at current iteration 
#        xi_mu_star -> list of k vectors, contains means of the current groups
#        xi_cov_star -> list of k matrices, contains cov. matrices of the current groups
#        old_group -> index of the group to be cancelled

# OUTPUT: curr -> n vector containing the updated clusters
#         xi_mu_star -> list of k vectors, contains means of the updated groups
#         xi_cov_star -> list of k matrices, contains cov. matrices of the updated groups

delete_and_shift <- function(curr, xi_mu_star, xi_cov_star, old_group)
{
  # Elimination of groups' parameters
  xi_mu_star_upd <- xi_mu_star[-old_group]
  xi_cov_star_upd <- xi_cov_star[-old_group]
  
  # Shifting of groups with higher index
  for (i in 1:length(curr)){
    if (curr[i]>old_group){
      curr[i] = curr[i]-1
    }
  }
  
  
  return(list("curr" = curr, "xi_mu_star" = xi_mu_star_upd, "xi_cov_star" = xi_cov_star_upd))
}

# Function to sample the parameter mu for a new group 
# INPUT: data -> datapoint 
#        xi_mu_star -> list of k vectors, contains means of the current groups
#        Q_param -> list of parameters nu_0_Q, mu_0_Q, k_0_Q, lambda_0_Q

# OUTPUT: xi_mu_star -> list of k vectors, contains means of the updated groups
construct_mu_new <- function(data,xi_mu_star, Q_param)
{
  # Transformation of data into a matrix
  data <- as.matrix(data)
  
  # Dimension of data
  p <- dim(data)[2]
  
  # Compute nu_n, k_n and mu_n
  nu_n <- Q_param$nu_0 + 1
  k_n <- Q_param$k_0 + 1
  mu_n <- as.vector(Q_param$k_0/(k_n) * Q_param$mu_0 + 1/k_n * data)
  
  r <- as.matrix(data-Q_param$mu_0)
  
  # Compute lambda_n
  lambda_n <- Q_param$lambda_0 + (Q_param$k_0)/(k_n)*r%*%t(r)
  
  # Compute sigma
  sigma = (1/(k_n*(nu_n-p+1)))*(lambda_n)
  
  # Sampling from a multivariate t-Student                  
  mu_new <- LaplacesDemon::rmvt(mu = mu_n, S = sigma, df = nu_n-p+1)
  mu_new <- as.vector(mu_new)
  
  # Updating the list containing the mean parameters
  xi_mu_star <- append(xi_mu_star,list(mu_new))
  
  return (xi_mu_star)
}


# Function to sample the parameter cov for a new group 
# INPUT: data -> datapoint 
#        xi_cov_star -> list of k vectors, contains cov. matrices of the current groups
#        Q_param -> list of parameters nu_0_Q, mu_0_Q, k_0_Q, lambda_0_Q

# OUTPUT: xi_cov_star -> list of k vectors, contains cov. matrices of the updated groups
construct_cov_new <- function(data,xi_cov_star, Q_param)
{
  # Transformation of data into a matrix
  data <- as.matrix(data)
  
  p <- dim(data)[2]
  
  # Compute nu_n, k_n
  nu_n <- Q_param$nu_0 + 1
  k_n <- Q_param$k_0 + 1
  
  r <- as.matrix(data-Q_param$mu_0)
  
  # Compute lambda_n
  lambda_n <- Q_param$lambda_0 + (Q_param$k_0)/(k_n)*r%*%t(r)
  
  # Compute sigma
  sigma = (1/(k_n*(nu_n-p+1)))*(lambda_n)
  
  #Sampling from an IW distribution
  # cov_new <- LaplacesDemon::rinvwishart(nu_n, as.inverse(sigma))
  cov_new <- LaplacesDemon::rinvwishart(nu_n, sigma)
  
  # Updating the list containing the cov. matrices parameters
  xi_cov_star<-append(xi_cov_star,list(cov_new))
  
  return (xi_cov_star)
}

