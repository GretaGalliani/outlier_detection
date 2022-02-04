library(MASS)
library(RColorBrewer)
library(robustbase)
library(TSstudio)
library(parallel)
library(pbapply)
library(pbmcapply)
library(future)
library(doParallel) 

source("real_dataset_doge/doge-usd_function.R")
#SETTA I PARAMETRI CHE VUOI DARE IN INPUT HAI DIVERSI MODELLI
input_cores_1 <- list("costante"=6,"k0_Q"=1,"k0_P"=0.1) 
input_cores_2 <- list("costante"=6,"k0_Q"=1,"k0_P"=0.01) 
input_cores_3 <- list("costante"=6,"k0_Q"=0.5,"k0_P"=0.5) 
input_cores_4 <- list("costante"=6,"k0_Q"=0.5,"k0_P"=0.1) 
input <- list("1" = input_cores_1,
              "2" = input_cores_2,
              "3" = input_cores_3,
              "4" = input_cores_4)


############################################################
# WINDOWS

cl <- makeCluster(detectCores(), type='PSOCK')
result <-parLapply (cl, input, fun = doge_function )

result_1 = result[[1]]$Result
save(result_1, file='postcor6Q1P025.RData')
result_2 = result[[2]]$Result
save(result_2, file='postcor6Q05P05.RData')
result_3 = result[[3]]$Result
save(result_3, file='postcor6Q1P05.RData')
result_4 = result[[4]]$Result
save(result_4, file='postcor6Q05P025.RData')

stopCluster (cl)
#############################################################

#numCores <- detectCores()
numCores=4
output_model <- pbmclapply(input, FUN=doge_function, mc.cores = numCores,mc.set.seed=TRUE)


result_1 <- output_model[[1]]$Result
result_2 <- output_model[[1]]$Result
result_3 <- output_model[[3]]$Result
result_4 <- output_model[[4]]$Result


#SE VUOI SALVARE SU FILE EXCELL
#Costruisci il dataframe inserendo ogni pezzetto di ogni lista
data <- data.frame()
#install.packages("xlsx")
library(xlsx)
#Crei prima un file excell vuoto lo metti da qualche parte e ci metti il percorso
percorso <- ""
write.xlsx(data, file=percorso, sheetName = "Sheet3", 
           col.names = TRUE, row.names = FALSE, append = TRUE)

#Stai attento ogni volta che runni tutto a cambiare lo sheetName perchè senno sovrascrivi

