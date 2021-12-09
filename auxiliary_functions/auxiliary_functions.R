##### AUXILIARY FUNCTIONS

# Function to compute m1_bar (number of points in group 0)
# INPUT: S -> vector of labels Si i=1,...,n
# OUTPUT: m1_bar -> number of points in group 0
m1_bar <- function(S){
  return (sum(S==0))
}

# Function to compute m1 (number of singletons (number of points in group 0 + number of groups with one point))
# INPUT: S -> vector of labels Si i=1,...,n
# OUTPUT: m1 -> number of singletons
m1 <- function(S){
  # Computation of the number of points in group 0, which are all singletons
  m1_bar = m1_bar(S)
  
  # Computation of the frequencies of all groups, excluding group 0
  x = S[S!=0]
  
  # I return the number of singletons
  return (m1_bar + sum(table(x)==1))
}


# Function to construct xi_mu_star and xi_cov_star from xi_mu and xi_cov at the beginning of the algorithm 
# INPUT: S -> vector of labels Si i=1,...,n
#        xi_mu -> list xi_mui i=1,...,n which contains for each data point the value of mu initialized 
#        xi_cov -> list xi_covi i=1,...,n which contains for each data point the value of cov initialized 
# WE ASSUME COHERENCE IN THE INITIALIZATION OF S AND XI
# OUTPUT: xi_mu_star -> list xi_mu_starj j=1,...,k which contains for each group the common value of mu initialized 
#         xi_cov_star -> list xi_cov_starj j=1,...,k which contains for each group the common value of cov initialized 
init_xi_star <-function(S, xi_mu, xi_cov)
{
  # Initialization of the lists
  xi_mu_star=list()
  xi_cov_star=list()
  
  # for each group found from the structure of S
  for (i in unique(S))
  {
    # Identification of one data point which belongs to group i
    l=which(S==i)
    
    # Appending of the mu and cov of the corresponding data point
    xi_mu_star=append(xi_mu_star,which(xi_mu[[l[1]]]))
    xi_cov_star=append(xi_cov_star,which(xi_cov[[l[1]]]))
  }
  
  
  my_list = list("xi_mu_star"=xi_mu_star, "xi_cov_star"=xi_cov_star)
  return (my_list)                              
}