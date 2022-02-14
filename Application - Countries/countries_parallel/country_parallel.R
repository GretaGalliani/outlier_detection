library(MASS)
library(RColorBrewer)
library(robustbase)
library(TSstudio)
library(parallel)
library(pbapply)
library(pbmcapply)
library(future)

source("country_parallel/country_function.R")

# Set the input parameters for the models
input_cores_1 <- list("const"=1,"k0_Q"=0.8,"k0_P"=0.05) 
input_cores_2 <- list("const"=0.9,"k0_Q"=0.95,"k0_P"=0.05) 
input_cores_3 <- list("const"=0.8,"k0_Q"=0.95,"k0_P"=0.05) 
input_cores_4 <- list("const"=0.9,"k0_Q"=0.7,"k0_P"=0.05) 
input <- list("1" = input_cores_1,
              "2" = input_cores_2,
              "3" = input_cores_3,
              "4" = input_cores_4)


numCores <- detectCores()

output_model <- pbmclapply(input, FUN=country_function, mc.cores = numCores)

result_1 <- output_model[[1]]$result
result_2 <- output_model[[1]]$result
result_3 <- output_model[[3]]$result
result_4 <- output_model[[4]]$result
