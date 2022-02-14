library(MASS)
library(RColorBrewer)
library(robustbase)
library(TSstudio)
library(parallel)
library(pbapply)
library(pbmcapply)
library(future)
library(doParallel) 

source("Dogecoin_parallel/Dogecoin_function.R")

# Set the input parameters for the models
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
result_2 = result[[2]]$Result
result_3 = result[[3]]$Result
result_4 = result[[4]]$Result

stopCluster (cl)
#############################################################
#MAC

numCores <- detectCores()
output_model <- pbmclapply(input, FUN=doge_function, mc.cores = numCores,mc.set.seed=TRUE)


result_1 <- output_model[[1]]$Result
result_2 <- output_model[[1]]$Result
result_3 <- output_model[[3]]$Result
result_4 <- output_model[[4]]$Result