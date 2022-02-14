# Define the function to be called by the different cores

doge_function <- function(input)
{
  
  doge = read.csv('Dogecoin_dataset.csv')
  date = as.Date.character(doge$Date)
  doge$LogReturn = rep(0, dim(doge)[1])
  for (i in 2:dim(doge)[1]){doge$LogReturn[i] = 
    log(doge$Adj.Close[i]/doge$Adj.Close[i-1])*100}
  data = as.matrix(doge$LogReturn)
  data=data[1:365,]
  data= scale(data)
  n = dim(data)[1]
  d = dim(data)[2]

  # Initialization of the parameters for the priors
  Q_param = list()
  P_param = list()

  # Contaminated component
  Q_param$k_0 = input$k0_Q ###VARIABILE CHE GLI PASSI IN INPUT
  Q_param$mu_0 = mean(data)
  Q_param$nu_0 = d+3
  Q_param$lambda_0 = var(data)/input$costante ###VARIABILE CHE GLI PASSI IN INPUT
  
  # Contaminant diffuse component
  P_param$k_0 = input$k0_P ###VARIABILE CHE GLI PASSI IN INPUT
  P_param$mu_0 = mean(data)
  P_param$nu_0 = d+3
  P_param$lambda_0 = var(data)/input$costante ###VARIABILE CHE GLI PASSI IN INPUT
  
  # Initialization of the parameters for the Pitman-Yor and initial partition
  S_init = rep(1, n)
  beta_init <- 0.5
  sigma_init <- 0.2
  theta_init <- 0.2
  
  beta_param = list()
  sigma_param = list()
  theta_param = list()
  beta_param$a = 1
  beta_param$b = 1
  sigma_param$a = 1
  sigma_param$b = 1
  theta_param$a = 2
  theta_param$b = 1
  
  xi_mu <- list()
  xi_cov <- list()
  
  init_mu <- mean(data)
  init_var <- var(data)

  
  for (i in 1:dim(data)[1]){
    xi_mu <- append(xi_mu, list(init_mu))
    xi_cov <- append(xi_cov, list(init_var))
  }
  
  source("main.R")
  result <- algorithm(data, S_init, sigma_init, theta_init, beta_init, 
                      beta_param, sigma_param, theta_param, xi_mu, xi_cov,
                      Q_param, P_param, 2000, 100, 1)
  
    output_model<-list("Result"=result)
  
  return (output_model)
}
