occurrence  = function(x,sigma){
  #This function uses us inputs the observed climate data and the variance
  #covariance matrix, sigma, and outputs the synthetic occurrences sequences
  n = dim(x)[1]# number of days
  m = dim(x)[2]# number of sites
  # creating dummy variables 
  for (i in 1:n) {
    for (j in 1:m) {
      if(x[i,j]<= 1.23){
        x[i,j] = 0
      } else{
        x[i,j] = 1
      }
    }
  }
  sim = matrix(nrow = n,ncol = m)#the matrix in which we will put daily occurrences simulations
  for (i in 1:m) {
    sim[1,i] = x[1,i]
  }
   mu = rep(0,m)# the mean vector of the multivariate normal distribution
  set.seed(1000)# the random seed
  mv = pnorm(mvrnorm(n-1,mu, sigma)) 
  for (j in 1:m) {
    for (i in 2:n) {
      fit = markovchainFit(data = x[,j], method = "mle", confidencelevel = 0.95)
      tran = fit$estimate# the estimation of transition probabilities
      u = mv[i-1,j]
      if(sim[i-1,j]==0){
        if(u>tran[1,2]){
          sim[i,j] = 0
        } else {
          sim[i,j] = 1
        } 
      } else {
        if(u>tran[2,2]){
          sim[i,j] = 0
        } else {
          sim[i,j] = 1
        }
      }
    }
  }
  return(sim)
}