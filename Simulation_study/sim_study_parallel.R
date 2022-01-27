library(parallel)
library(MASS)
library(pbapply)
library(pbmcapply)
library(future)
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

input_cores <- list("d"=2,"c"=1.25,"k0_Q"=1,"k0_P"=0.25,"m"=90) 
input <- list("1" = input_cores,
              "2" = input_cores,
              "3" = input_cores,
              "4" = input_cores,
              "5" = input_cores,
              "6" = input_cores,
              "7" = input_cores,
              "8" = input_cores,
              "9" = input_cores,
              "10" = input_cores,
              "11" = input_cores,
              "12" = input_cores,
              "13" = input_cores,
              "14" = input_cores,
              "15" = input_cores,
              "16" = input_cores,
              "17" = input_cores,
              "18" = input_cores,
              "19" = input_cores,
              "20" = input_cores,
              "21" = input_cores,
              "22" = input_cores,
              "23" = input_cores,
              "24" = input_cores)



#Versione 1
numCores <- detectCores()
# tempo <- system.time({
# output_model <- mclapply(input, sim_study, mc.cores = numCores,mc.set.seed=TRUE)
# })
output_model <- pbmclapply(input, FUN=sim_study, mc.cores = numCores,mc.set.seed=TRUE)

vettore_beta <- NULL
vettore_k <- NULL
vettore_singletons <- NULL
vettore_BL <- NULL
vettore_VI <- NULL
vettore_BL_k <- NULL
vettore_VI_k <- NULL
for(i in (1:24))
{
  vettore_beta <- c(vettore_beta,output_model[[i]]$beta_mean)
  vettore_k <- c(vettore_k,output_model[[i]]$k_mean)
  vettore_singletons <- c(vettore_singletons,output_model[[i]]$singletons_mean)
  vettore_BL <- c(vettore_BL,output_model[[i]]$BL)
  vettore_VI <- c(vettore_VI,output_model[[i]]$VI)
  vettore_BL_k <- c(vettore_BL_k,output_model[[i]]$BL_k)
  vettore_VI_k <- c(vettore_VI_k,output_model[[i]]$VI_k)
  
}

data <- data.frame("beta"=vettore_beta, "k"=vettore_k, "singletons"=vettore_singletons,
                   "BL_k"=vettore_BL_k,"VI_k"=vettore_VI_k,
                   "BL"=vettore_BL,"VI"=vettore_VI)
#install.packages("xlsx")
library(xlsx)

percorso <- "Simulation_study/valori_modelli_sim_study.xlsx"
write.xlsx(data, file=percorso, sheetName = "Sheet3", 
           col.names = TRUE, row.names = FALSE, append = TRUE)



#Versione 2
# cl=makeCluster(parallel::detectCores())
# prova <- parSapply(cl=cl,X=input,FUN=sim_study,USE.NAMES=TRUE)


