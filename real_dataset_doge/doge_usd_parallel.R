library(MASS)
library(RColorBrewer)
library(robustbase)
library(TSstudio)
library(parallel)
library(pbapply)
library(pbmcapply)
library(future)

source("real_dataset_doge/doge-usd_function.R")
#SETTA I PARAMETRI CHE VUOI DARE IN INPUT HAI DIVERSI MODELLI
input_cores_1 <- list("costante"=1,"k0_Q"=1,"k0_P"=0.25) 
input_cores_2 <- list("costante"=1,"k0_Q"=1,"k0_P"=0.25) 
input_cores_3 <- list("costante"=1,"k0_Q"=1,"k0_P"=0.25) 
input_cores_4 <- list("costante"=1,"k0_Q"=1,"k0_P"=0.25) 
input <- list("1" = input_cores_1,
              "2" = input_cores_2,
              "3" = input_cores_3,
              "4" = input_cores_4)


numCores <- detectCores()

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

