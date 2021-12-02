#initialization

Initialization <-function(S, Xi_mu, Xi_cov)
  {
  k=lenght(unique(S))
  for (i in unique(S))
    {
    l=which(S==i)
    Xi_mu_star=append(Xi_mu_star,which(Xi_mu[[l[1]]]))
    Xi_cov_star=append(Xi_cov_star,which(Xi_cov[[l[1]]]))
}
  my_list = list("Xi_mu_star"=Xi_mu_star, "Xi_cov_star"=Xi_cov_star)
  return (my_list)                              
}

