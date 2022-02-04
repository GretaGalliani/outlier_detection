library(MASS)
library(RColorBrewer)
library(robustbase)
library(TSstudio)
library(parallel)
library(pbapply)
library(pbmcapply)
library(future)

source("country_parallel/country_function.R")
#SETTA I PARAMETRI CHE VUOI DARE IN INPUT HAI DIVERSI MODELLI
input_cores_1 <- list("costante"=1,"k0_Q"=0.8,"k0_P"=0.05) 
input_cores_2 <- list("costante"=0.9,"k0_Q"=0.95,"k0_P"=0.05) 
input_cores_3 <- list("costante"=0.8,"k0_Q"=0.95,"k0_P"=0.05) 
input_cores_4 <- list("costante"=0.9,"k0_Q"=0.7,"k0_P"=0.05) 
input <- list("1" = input_cores_1,
              "2" = input_cores_2,
              "3" = input_cores_3,
              "4" = input_cores_4)


numCores <- detectCores()

output_model <- pbmclapply(input, FUN=country_function, mc.cores = numCores,mc.set.seed=TRUE)


result_1 <- output_model[[1]]$Result
result_2 <- output_model[[1]]$Result
result_3 <- output_model[[3]]$Result
result_4 <- output_model[[4]]$Result

save(result_1, file='country_Q08_P005.RData')
save(result_2, file='country_Q095_P005_Cov_09.RData')
save(result_3, file='country_Q095_P005_Cov_08.RData')
save(result_4, file='country_Q07_P005_Cov_09.RData')

#SE VUOI SALVARE SU FILE EXCELL
#Costruisci il dataframe inserendo ogni pezzetto di ogni lista
data <- data.frame()
#install.packages("xlsx")
library(xlsx)
#Crei prima un file excel vuoto lo metti da qualche parte e ci metti il percorso
percorso <- ""
write.xlsx(data, file=percorso, sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE, append = TRUE)

#Stai attento ogni volta che runni tutto a cambiare lo sheetName perchè senno sovrascrivi

