### UPDATE CLUSTER

#NEL MAIN: 
#Creo una struttura che è lista di liste --> support <-- inizializzazione (ogni dato è un gruupo)
#--> lista(liste(0:n))
#Devo inizializzare anche le frequenze iniziali
#Devo iniziare anche Xi_star (entrambi)

support <-list()
# K --> numero gruppi iniziali
k=3
empty_list <- vector(mode="list",length=k)

for(i in 1:N)
{
  for(j in 1:k)
    if(S[i]==j){
    empty_list[j] <- append(empty_list[j],i)
    }
}

#r=1,...,R


#Input:
# support <- lista di liste che contiene i cluster e li gestisce ad ogni passo
# Y <- data.frame
# Xi_mu_star <- lista di n vettori --> struttura che contiene medi gruppi
# Xi_cov_star <- lista di n matrici --> struttura che contiene matrici cov dei gruppi
# beta_old <- numero
# theta_old <- numero
# sigma_old <- numero
# k_old <- numero
# freq <- vettore contenente le frequenze di ogni cluster
# support_0 <- supporto per j=0



update_cluster <- function(support, support_0, Y, Xi_mu, Xi_cov, beta_old, Xi_mu_star, 
                                    Xi_cov_star, theta_old, sigma_old, k_old, freq)
{
  N <- dim(Y)[1]
  S_new <-rep(0,n)
  Xi_mu_new <- list(rep(0,n)) #lista( di n vettori)
  Xi_cov_new <- #list( di n matrici)
  k_new <- k_old
  
  for (i in 1:N)
  {
    prob <- rep(0:k_new+2,0) #0,1,...k_new, k_new+1
    
    prob[1] <- dens_contaminated(data[i],beta_old, nu_0, k_0, mu_0, delta_0, n) #j=0
    
    for (t in (2:k_new+1))
    {
      prob[t] <- dens_cluster_old(data[i], n_j, sigma_old, theta_old, Xi_mu_old_star[t-1],
                                 Xi_cov_old_star[t-1],n,beta_old)
    }
    
    prob[k_new+2]<- dens_cluster_new(data[i], n_j, sigma_old, theta_old, beta_old,n, k_old)
    
    j <- sample(0:k_new+1,size=1,prob=prob)
    
    S_new[i] <- j
    
    if (support == 1)
    {
      #Se frequenza si S_old[i] (o support) = 1 --> chiamo cancella (sistema support eliminando il gruppo j
      #e spostando indietro gli altri, aggiornando k e xi & traslo Xi_mu_star, Xi_cov_star)

    }
    else if ()
    Se nuovo gruppo --> creo nuova lista nel support
    Altrimenti --> sposto dato nel support se già esiste il gruppo
    
    
    if(j==k_new+1)
    {
      k_new=k_new+1
      
      mu_new <- construct_mu_new(??)
      Xi_mu_new[i] <- mu_new
      Xi_mu_star <- aggiungo mu_new
      
      cov_new <- construct_cov_new(??)
      Xi_cov_new[i] <- cov_new
      cov_new_star <- aggiungo cov_new
    }
    if(j in 1:k_new)
    {
      Xi_mu_new[i] <- Xi_mu_star[j]
      Xi_cov_new[i] <- Xi_cov_star[j]
    }
  } 
  
  #Calcolo nuove frequenze --> dopo lo spostamento   

}