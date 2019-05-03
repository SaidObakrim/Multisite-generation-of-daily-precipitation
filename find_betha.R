find_betha = function(x,sim,sigma){
  # This function specifies the scale parameter of the mixed exponential
  # distribution at every location on everyday
  library("MixtureInf")
  n = dim(x)[1]
  m = dim(x)[2]
  beta = matrix(nrow = m ,ncol=2)
  alpha = c()
  for (j in 1:m) {
    b = c()
    for (i in 1:n) {
      if(x[i,j]!=0){
        b = c(b, x[i,j])
      }
    }
    beta[j,1] = emtest.exp(b, m0= 2)$`MLE of parameters under null hypothesis (order = m0)`[2,1]
    beta[j,2] = emtest.exp(b, m0= 2)$`MLE of parameters under null hypothesis (order = m0)`[2,2]
    alpha[j] = emtest.exp(b, m0= 2)$`MLE of parameters under null hypothesis (order = m0)`[1,2]
  }
  betha = matrix(nrow = n-1,ncol= m)
  mu = rep(0,m)
  set.seed(1000)
  mv = pnorm(mvrnorm(n-1,mu, sigma)) 
  for (j in 1:m) {
    for (i in 2:n) {
      fit = markovchainFit(data = x[,j], method = "mle", confidencelevel = 0.95)
      tran = fit$estimate
      u = mv[i-1,j]
      if(sim[i-1,j]==0){
        if((u/tran[1,2])<= alpha[j]){
          betha[i-1,j] = beta[j,2]
        } else{
          betha[i-1,j] = beta[j,1]
        }
      }else{
        if((u/tran[2,2])<= alpha[j]){
          betha[i-1,j] = beta[j,2]
        } else{
          betha[i-1,j] = beta[j,1]
        }
        
      }
    }
  }
  return(betha)
}