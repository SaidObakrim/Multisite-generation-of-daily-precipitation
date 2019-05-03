find_d = function(d,sim, betha){
  #This function illustrates the difference between correlations of 
  # simulated and observed daily rainfall data.
  library("MASS")
  d1= d[1]; d2 = d[2]
  n = dim(sim)[1]
  m = dim(sim)[2]
  cor = matrix(nrow = m, ncol= m)
  for (i in 1:m) {
    for(j in 1:m){
      if(i!=j){
        cor[i,j] = exp(-(d1*((t(x[,i]- x[,j]) %*% (x[,i]-x[,j]))^(1/2))^d2))
      } else{
        cor[i,j] = 1
      }
    }
  }
  mu = rep(0,m)
  set.seed(1000)
  mv = pnorm(mvrnorm(n,mu, cor))
  for (j in 1:m) {
    for (i in 2:n) {
      if (sim[i,j]!=0){
        sim[i,j] = 1.23 - betha[i-1,j]* log(mv[i,j])
      }
    }
  }
  v = cor(x);
  vec = c() 
  for (i in 2:m) {
    vec = c(vec, v[i-1,i:m])
  }
  vv = cor(sim)
  vvec = c()  
  for (i in 2:m) {
    vvec = c(vvec, vv[i-1,i:m])
  }
  return(sum(abs(vec-vvec)))
}
