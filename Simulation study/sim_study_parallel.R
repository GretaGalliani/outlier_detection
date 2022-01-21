library(parallel)
library(MASS)

starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)

numCores

system.time(
  results <- mclapply(starts, fx, mc.cores = numCores)
)

sim_study <- function(lista_di_merda)
{
  c <- lista_di_merda$c
  d <- lista_di_merda$d
  prova <- mvrnorm(1,mu= rep(0,d),Sigma= 3^2*diag(d))
  output <- list(prova)
  return (output)
  
}

list_core_1 <- list("d"=,"c"=,"k0_P"=,"k0_Q"=,"m"=) 
list_core_2 <- list("d"=,"c"=,"k0_P"=,"k0_Q"=,"m"=) 
list_core_3 <- list("d"=,"c"=,"k0_P"=,"k0_Q"=,"m"=) 
list_core_4 <- list("d"=,"c"=,"k0_P"=,"k0_Q"=,"m"=) 

input <- list("1" = list_core_1,
              "2" = list_core_2,
              "3" = list_core_3,
              "4" = list_core_4)

#Versione 1
numCores <- detectCores()
results <- mclapply(input, sim_study, mc.cores = numCores,
                    mc.set.seed=TRUE)

#Versione 2
cl=makeCluster(parallel::detectCores())
prova <- parSapply(cl=cl,X=input,FUN=sim_study,USE.NAMES=TRUE)


