Error.G <- function(n,q,rho){
  if (rho==0){
    sigma <- diag(q)
  }else{
    sigma <- matrix(0,q,q)
    for(i in 1:q){for(j in 1:q){sigma[i,j] <- rho^(abs(i-j))}}
  }
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2
  Y <- mvrnorm(n,rep(0,q),sigma)
  ref1 <- dcov.test(G,Y,R=199)$p.value<0.05
  ref2 <- wdcor.test(G,Y,R=199,G.list=seq(1,25,2))$p.value<0.05
  return(c(ref1=ref1,ref2=ref2))
}
