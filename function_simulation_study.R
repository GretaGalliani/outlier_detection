sim_study <- function(d,c,k0_P,k0_Q,m)
  mean_a = rep(-3,d)#vector of means
  sigma_b = diag(d) #sigma
  m1 = rbinom(1, size=m, prob = 0.5) #number of samples coming from the first gaussian 
  m2 = m-m1 # from the second one
  
  val1 = mvrnorm (m1,mu = mean_a, Sigma = sigma_b) #samples from first multivariate function
  val2 = mvrnorm (m2,mu = -mean_a, Sigma = sigma_b) #samples from second multivariate function
  allval = rbind(val1,val2) #combine
  
  # allval <- NULL
  
  s=10 #number of outliers
  i=0 
  while(i<s){ #cycle to find s outliers
    value=mvrnorm(1,mu= rep(0,d),Sigma= 3^2*diag(d)) #sampling from a multivariate normal distribution
    module = norm(as.matrix(value), type="2")
    chi=qchisq(0.9, df = d)
    # if(module^2>3*sqrt(chi))
    if(module^2>9*chi) #If we are sampling from the over-disperse truncated Gaussian distribution
    {
      i=i+1
      value = c*value
      allval = rbind(allval,value)
    }
  }
  # allval = allval[sample(m+s,m+s),] #randomizing rows
  
  rownames(allval)=NULL
  data <- allval
  
  #### INIZIALIZZAZIONE - P0 DIVERSO DA Q0 ####
  
  Q_param = list()
  P_param = list()
  
  d = 2
  
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
  
  xi_mu <- list()
  xi_cov <- list()
  
  init_mu <- colMeans(data)
  init_var <- cov(data)
  
  for (i in 1:n){
    xi_mu <- append(xi_mu, list(init_mu))
    xi_cov <- append(xi_cov, list(init_var))
  }
  
  #### RUNNING THE ALGORITHM ####
  source("main.R")
  result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, beta_param, sigma_param, theta_param, xi_mu, xi_cov, Q_param, P_param, 15000, 1000, 10)
  #SE VOGLIAMO SALVARE I FILE! BISOGNA PASSARE UN ITERATORE i
  nome_file <- paste('output_salvati_prova',as.character(i),sep="_")
  nome_file <-paste(nome_file,'dat',sep=".")
  #save(result,file=nome_file)
  
  #### CLUSTER ANALYSIS ####
  max <- c()
  for (i in 1:1000){
    max <- c(max, max(result$S[i,]))
  }
  
  # WE COUNT THE NUMBER OF SINGLETONS
  source("auxiliary_functions/auxiliary_functions.R")
  
  n_singletons <- c()
  for (i in 1:dim(result$S)[1]){
    n_singletons <- c(n_singletons, m1(result$S[i,]))
  }
  
  # IMPLEMENTING MIN BINDER LOSS
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

  
  # IMPLEMENTING MIN VARIATION OF INFORMATION
  #devtools::install_github("sarawade/mcclust.ext")
  library(mcclust.ext)
  
  # finds the clustering that minimizes  the lower bound to the posterior expected Variation of Information from Jensen's Inequality
  min_vi <- minVI(psm, cls.draw=NULL, method=c("avg","comp","draws","greedy","all"), 
                  max.k=NULL, include.greedy=FALSE, start.cl=NULL, maxiter=NULL,
                  l=NULL, suppress.comment=TRUE)
  
  
  # best cluster according to binder loss 
  tab_vi <- table(min_vi$cl)

 output_model<-list("beta_mean"=result$beta, "beta_sd"=sd(result$beta), "k_mean"=mean(max),
                   "k_sd"=sd(max),"singletons_mean"=n_singletons,"singletons_sd"=sd(n_singletons),
                   "BL"=length(which(tab_bind==1)),"VI"=length(which(tab_vi==1)))
  
return (output_model)
}