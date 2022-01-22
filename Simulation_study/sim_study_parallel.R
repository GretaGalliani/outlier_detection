library(parallel)
library(MASS)

# sim_study <- function(lista_di_merda)
# {
#   c <- lista_di_merda$c
#   d <- lista_di_merda$d
#   prova <- mvrnorm(1,mu= rep(0,d),Sigma= 3^2*diag(d))
#   output <- list(prova)
#   return (output)
#   
# }

source("Simulation_study/function_simulation_study.R")

input_cores <- list("d"=2,"c"=1,"k0_Q"=1,"k0_P"=0.25,"m"=90) 
input <- list("1" = input_cores,
              "2" = input_cores,
              "3" = input_cores,
              "4" = input_cores)

#Versione 1
numCores <- detectCores()
output_model <- mclapply(input, sim_study, mc.cores = numCores,
                    mc.set.seed=TRUE)

#Versione 2
cl=makeCluster(parallel::detectCores())
prova <- parSapply(cl=cl,X=input,FUN=sim_study,USE.NAMES=TRUE)


