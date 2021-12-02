

source("Update_clusters/Update_cluster.R")

#PASSO A
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
  
 
  