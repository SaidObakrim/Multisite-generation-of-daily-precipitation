mn = function(x,sm,m1){
  #This function simulates the monthly/annual data
  #x: the observed monthly/annual data
  #sm: the simulated data before adjustement 
  #m1: the lag one cross correlation matrix of the observed monthly/annual data
  mt = matrix(nrow = dim(sm)[1], ncol= dim(sm)[2])
  sm= scale(sm)
  x= scale(x)
  mt[1,]= sm[1,]
  m0 = cor(x)
  vc = c()
  for(i in 1:6){
    vc = c(vc, acf(x[,i],lag=1,plot=FALSE)$acf[2])
  }
  a = m1%*%solve(m0)
  a = diag(vc)
  ff = m0 - (a%*%t(m1))
  sv = svd(ff)
  f = (sv$u)%*%((diag(sv$d))^(1/2))
  dd = cor(sm)
  svv = svd(dd)
  d = (svv$u)%*%((diag(svv$d))^(1/2))
  b = f%*% solve(d)
  for (k in 2:dim(sm)[1]) {
    mt[k,]= a%*%mt[(k-1),]+ b%*%sm[k,]
  }
  return(mt)
}