#This function implement the algorithm for the simulation study:
#Generate a set of n data from a mixture of two Gaussian distribution 
#Generate a set of S outliers from an over-disperse truncated gaussian

#Input = list(d,c,K0_Q,k0_P), we pass this parameters since we want to parallelized this function
sim_study <- function(input) 
{
  d <- input$d #Dimension: {2,4}
  c <- input$c # {1,1.25,1.5}
  k0_P <- input$k0_P #{0.5,0.25}
  k0_Q <- input$k0_Q #{0.5,1}
  m <- input$m #Number of samples {90,240}
  
  ## SAMPLING FROM THE DENSITY####
  # BASE MEASURE
  mean_a = rep(-3,d) #Vector of means
  sigma_b = diag(d) #Sigma
  m1 = rbinom(1, size=m, prob = 0.5) #Number of samples coming from the first gaussian distribution 
  m2 = m-m1 #Number of samples coming from the second gaussian distribution 
  
  val1 = mvrnorm (m1,mu = mean_a, Sigma = sigma_b) #Samples from first multivariate function
  val2 = mvrnorm (m2,mu = -mean_a, Sigma = sigma_b) #Samples from second multivariate function
  allval = rbind(val1,val2) #Combine
  
  # CONTAMINATED MEASURE ####
  s=10 #Number of outliers
  i=0 
  while(i<s){ #cycle to find s outliers
    value=mvrnorm(1,mu= rep(0,d),Sigma= 3^2*diag(d)) #sampling from a multivariate normal distribution
    module = norm(as.matrix(value), type="2")
    chi=qchisq(0.9, df = d)
    if(module^2>9*chi) #If we are sampling from the over-disperse truncated Gaussian distribution
    {
      i=i+1
      value = c*value #C has the role to shrink or expand the nuisance observations towards the origin
      allval = rbind(allval,value)
    }
  }
  # Constructing the dataset
  rownames(allval)=NULL
  data <- allval
  
  # Plotting the data with the original groups
  pal = brewer.pal(n = 9, name = "Set1")
  col_real = c(rep(pal[1], m1), rep(pal[3], m-m1), rep(pal[2],s))
  pc = c(rep(16,m), rep(17,s))
  
  pairs(data, col = col_real, pch = pc, main = "Real data")
  
  
  ##INITIALIZATION for Q0 and P0 ####
  Q_param = list()
  P_param = list()
  
  Q_param$k_0 = k0_Q 
  Q_param$mu_0 = c(0,0)
  Q_param$nu_0 = d + 3 # it must be > (p-1)
  Q_param$lambda_0 = diag(diag(cov(data)))
  
  P_param$k_0 = k0_P #{0.25,0.5}
  P_param$mu_0 = c(0,0)
  P_param$nu_0 = d + 3 # it must be > (p-1)
  P_param$lambda_0 = diag(diag(cov(data)))
  
  n = dim(data)[1]
  S_init = rep(1,n)
  
  beta_init <- 0.5
  sigma_init <- 0.5
  theta_init <- 1
  
  beta_param = list()
  sigma_param = list()
  theta_param = list()
  beta_param$a = 1
  beta_param$b = 1
  sigma_param$a = 1
  sigma_param$b = 1
  #Gamma con picco in 1
  theta_param$a = 2 
  theta_param$b = 1
  
  #Create two list empty for xi_mu and xi_cov
  xi_mu <- list()
  xi_cov <- list()
  
  init_mu <- colMeans(data)
  init_var <- cov(data)
  
  for (i in 1:n){
    xi_mu <- append(xi_mu, list(init_mu))
    xi_cov <- append(xi_cov, list(init_var))
  }
  
  ## RUNNING THE ALGORITHM ####
  source("algorithm_v1/main.R")
  result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, beta_param, sigma_param, theta_param, xi_mu, xi_cov, Q_param, P_param, 6000, 1000, 1)
  
  #### CLUSTER ANALYSIS ####
  max <- c()
  for (i in 1:dim(result$S)[1]){
    max <- c(max, max(result$S[i,]))
  }
  
  # Ww count the numbers of singletons
  source("algorithm_v1/auxiliary_functions.R")
  
  n_singletons <- c()
  for (i in 1:dim(result$S)[1]){
    n_singletons <- c(n_singletons, m1(result$S[i,]))
  }
  
  ## LOSS FUNCTION ####
  #Binder loss function
  library(mcclust)
  
  aux = result$S
  
  # These functions needs to have indexes of the groups >=1
  for (i in 1:dim(aux)[1]){
    for (j in 1:dim(aux)[2]){
      if (aux[i,j]==0){
        aux[i,j] = max(aux[i,]) + 1
      }
    }
  }
  
  # For a sample of clusterings of the same objects the proportion of clusterings 
  # in which observation $i$ and $j$ are together in a cluster is computed and 
  # a matrix containing all proportions is given out. 
  
  psm <- comp.psm(aux)
  
  # finds the clustering that minimizes the posterior expectation of Binders loss function
  min_bind <-  minbinder(psm, cls.draw = NULL, method = c("avg", "comp", "draws", 
                                                          "laugreen","all"), max.k = NULL, include.lg = FALSE, 
                                                          start.cl = NULL, tol = 0.001)
  
  # best cluster according binder loss 
  tab_bind <- table(min_bind$cl) 

  
  # Variation of Information
  #devtools::install_github("sarawade/mcclust.ext")
  library(mcclust.ext)
  
  # finds the clustering that minimizes  the lower bound to the posterior expected Variation of Information from Jensen's Inequality
  min_vi <- minVI(psm, cls.draw=NULL, method=c("avg","comp","draws","greedy","all"), 
                  max.k=NULL, include.greedy=FALSE, start.cl=NULL, maxiter=NULL,
                  l=NULL, suppress.comment=TRUE)
  
  
  # best cluster according to binder loss 
  tab_vi <- table(min_vi$cl)
  
  #We don't report the diagnostic since when we have the result from the parallelized function
  #we upload the ouput of the alghoritm in the file simulation_study_d=4 or simulation_study_d=2
  
  output_model<-list("beta_mean"=mean(result$beta), "k_mean"=mean(max),"singletons_mean"=mean(n_singletons),
                    "BL_k"=length(which(tab_bind>1)),"VI_k"=length(which(tab_vi>1)),
                   "BL"=length(which(tab_bind==1)),"VI"=length(which(tab_vi==1)))
  
return (output_model)
}