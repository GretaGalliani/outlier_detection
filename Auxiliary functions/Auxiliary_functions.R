##### AUXILIARY FUNCTIONS

# Function computing m1_bar
# i.e. computing the number of samples classified as outliers (labeled with 0)
# Input: vector of labels Si i=1,...,n
# ALREADY TESTED
m1_bar <- function(S){
  return (sum(S==0))
}

# Function computing m1
# i.e. computing the number of singletons: the outliers and the groups
# containing only one sample
# ALREADY TESTED
m1 <- function(S){
  m1_bar = m1_bar(S)
  x = S[S!=0]
  return (m1_bar + sum(table(x)==1))
}
