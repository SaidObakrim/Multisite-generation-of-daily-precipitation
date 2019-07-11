amounts = function(x,sim,sigm){
  #This function simulates the nonzero precipitation amounts.
  # x: the observed precipitation amounts
  # sim : the simulated occurrence process
  # sigm : The output of find_sig function
  n = dim(sim)[1]
  m = dim(sim)[2]
  for (j in 1:m) {
    for (i in 1:n) {
      if(x[i,j]<=0.25){
        x[i,j]=0
      }
    }
  }
  mu = rep(0,m)
  set.seed(100)
  mv = pnorm(mvrnorm(n,mu, sigm))
  for (j in 1:m) {
    mod1= fitdistr(x[,j][x[,j]!=0],"gamma")
    for (i in 1:n) {
      if (sim[i,j]!=0){
        sim[i,j] = qgamma(mv[i,j],shape = mod1$estimate[1], rate =mod1$estimate[2])
      }
    }
  }
  return(sim)
}
