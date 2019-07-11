find_sg = function(x,sim){
  #This function is used to find the correlation structure of the multivariate
  #normal distribution that gives the disered correlation for the amounts process.
  # x: the observed precipitation amounts
  # sim : the simulated occurrence process
  library("fitdistrplus")
  library("MASS")
  n = dim(x)[1] # number of days
  m = dim(x)[2] # number of sites
  p  = m*(m-1)*0.5 # number of pairwise correlations
  for (j in 1:m) {
    for (i in 1:n) {
      if(x[i,j]<=0.25){
        x[i,j]=0
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
  cor1 = matrix(nrow = 100, ncol = p)#the matrix in which we will store the correlatins between sites
  cor2 = c()
  bet = matrix(nrow = m, ncol = m)
  for (j in 1:100) {
    for (v in 1:m) {
      set.seed(100)
      sigma = matrix(c(1, seqq[j+1],seqq[j+1],1), ncol = 2, nrow = 2)
      mu = c(0,0)
      set.seed(100)
      z = mvrnorm(n,mu, sigma)
      mv = pnorm(z)
      for (c in 1:m) {
        if (c != v){
          mod1= fitdistr(x[,c][x[,c]!=0],"gamma")
          mod2= fitdistr(x[,v][x[,v]!=0],"gamma")
          for (i in 1:n) {
            u1 = mv[i,1]
            u2 = mv[i,2]
            sim[i,c] = qgamma(u1,shape = mod1$estimate[1], rate =mod1$estimate[2]) 
            sim[i,v] = qgamma(u2,shape = mod2$estimate[1], rate =mod2$estimate[2]) 
          }
          bet[v,c] = cor(sim[,v], sim[,c])# store the correlations in a trigonal matrix
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
    if(vec[i]<=0.8){
      # estimating the relationship between correlations of 2 normal variates and 
      # the corresponding simulated daily occurences using linear regression
      y = cor1[,i]; x1 = cor2
      mod  = lm(y~ x1)
      fu = function(x){
        mod$coefficients[1]+ mod$coefficients[2]*x- vec[i]-0.07
      }
      un[i] = bisection(fu, 0,1) # the value of cor(z[,1],z[,2]) which 
      #gives the desired correlation between the corresponding sites
    } else {
      # estimating the relationship between correlations of 2 normal variates and 
      # the corresponding simulated daily occurences using linear regression
      y = cor1[,i]; x1 = cor2
      mod  = lm(y~ x1)
      fu = function(x){
        mod$coefficients[1]+ mod$coefficients[2]*x- vec[i]
      }
      un[i] = bisection(fu, 0,1) # the value of cor(z[,1],z[,2]) which 
      #gives the desired correlation between the corresponding sites
    }
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