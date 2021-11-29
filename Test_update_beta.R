n = 1000
m1_bar = 500

#x <- runif(n)

beta <- rep(0.1,1000)
for (r in 2:1000)
{
  beta[r] <- update_beta(n, m1_bar, beta[r-1])
}

