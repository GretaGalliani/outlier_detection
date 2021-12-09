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
