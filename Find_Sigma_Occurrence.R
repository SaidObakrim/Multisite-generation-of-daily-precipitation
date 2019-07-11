find_sigma = function(x){
  #This function is used to find the correlation structure of the multivariate
  #normal distribution that gives the disered correlation for the occurrence process. 
  # x: the observed precipitation amounts
  library("MASS")
  n = dim(x)[1] # number of days
  m = dim(x)[2] # number of sites
  p  = m*(m-1)*0.5 # number of s
  # creating dummy variables 
  for (i in 1:n) {
    for (j in 1:m) {
      if(x[i,j]<= 0.25){
        x[i,j] = 0
      } else{
        x[i,j] = 1
      }
    }
  }
  un = c()
  v = cor(x);
  vec = c() # a vector contains correlations between sites 
  for (i in 2:m) {
    vec = c(vec, v[i-1,i:m])
  }
  seqq = seq(0,1, 0.01)
  ma = matrix(nrow = n, ncol = m)# the matrix in which we will put daily simulations
  for (i in 1:m) { 
    ma[1,i] = x[1,i]
  }
  cor1 = matrix(nrow = 100, ncol = p)#the matrix in which we will store the correlatins between sites
  cor2 = c()
  bet = matrix(nrow = m, ncol = m)
  for (j in 1:100) {
    for (v in 1:m) {
      set.seed(100)
      sigma = matrix(c(1, seqq[j+1],seqq[j+1],1), ncol = 2, nrow = 2)
      mu = c(0,0)
      set.seed(100)
      z = mvrnorm(n-1,mu, sigma)
      mv = pnorm(z)
      for (c in 1:m) {
        if (c != v){
          library("markovchain")
          fit1 = markovchainFit(data = x[,v], method = "mle", confidencelevel = 0.95)
          tran1 = fit1$estimate # estimate the transitions probabilities
          fit2 = markovchainFit(data = x[,c], method = "mle", confidencelevel = 0.95)
          tran2 = fit2$estimate
          for (i in 2:n) {
            u1 = mv[i-1,1]
            u2 = mv[i-1,2]
            if(ma[i-1,v]==0){
              if(u1>tran1[1,2]){
                ma[i,v] = 0
              } else {
                ma[i,v] = 1
              } 
            } else {
              if(u1>tran1[2,2]){
                ma[i,v] = 0
              } else {
                ma[i,v] = 1
              }
            }
            if(ma[i-1,c]==0){
              if(u2>tran2[1,2]){
                ma[i,c] = 0
              } else {
                ma[i,c] = 1
              } 
            } else {
              if(u2>tran2[2,2]){
                ma[i,c] = 0
              } else {
                ma[i,c] = 1
              }
            }
          }
          bet[v,c] = cor(ma[,v], ma[,c])# store the correlations in a trigonal matrix
        }
      }
    }
    vect = c()
    for (i in 2:m) {
      vect = c(vect, bet[i-1,i:m])# restore the previous correltions in a vector
    }
    for (i in 1:p) {
      cor1[j,i] = vect[i]# restore the correlations in the j'st cor1's row
    }
    cor2[j] = cor(z[,1],z[,2])
  }
  for (i in 1:p) {
    # estimating the relationship between correlations of 2 normal variates and 
    # the corresponding simulated daily occurences using linear regression
    y = cor1[,i]; x1 = cor2^1;x2 = cor2^16
    mod  = lm(y~ x1+x2)
    fu = function(x){
      mod$coefficients[1]+ mod$coefficients[2]*x^1+ mod$coefficients[3]*x^16- vec[i]
    }
    un[i] = bisection(fu, 0,1) # the value of cor(z[,1],z[,2]) which 
    #gives the desired correlation between the corresponding sites
  }
  sigma = matrix(nrow=m,ncol=m); j=1; v=m-2;# the desired covariance matrix
  for (i in 2:m) {
    k = j+v; sigma[i-1,i:m] = un[j:k]; j=k+1 ; v=v-1
  }
  for (i in 1:m) {
    for (j in 1:m) {
      sigma[j,i] = sigma[i,j]
      sigma[i,i] = 1;
    }
  }
  return(sigma)
}