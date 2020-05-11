Wei_G_case1 <- function(n,q,q1){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2
  
  sigma <- diag(q)
  Y <- mvrnorm(n, mu=rep(0,q), Sigma = sigma)
  
  wei<- autoweight(G,Y)
  return(wei)
}
Wei_G_case2 <- function(n,q,q1){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2
  
  sigma <- diag(q)
  Beta <- c(rep(0.3,q1),rep(0,q-q1))  
  E <- mvrnorm(n,rep(0,q),sigma)  
  Y <- G%*%Beta + E 
  
  wei<- autoweight(G,Y)
  return(wei)
}
Wei_G_case3 <- function(n,q,q1){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2
  
  sigma <- diag(q)
  E <- mvrnorm(n,rep(0,q),sigma)
  
  Beta <- c(rep(0.3,q1-1),rep(0,q-q1))  
  
  Y1 <- 0.3*G^2 
  Y2 <- G%*%Beta 
  Y <- cbind(Y1,Y2)+E
  
  wei<- autoweight(G,Y)
  return(wei)
}
Wei_G_case4 <- function(n,q,q1){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2
  
  # gengrate Y
  sigma <- diag(q)
  E <- mvrnorm(n,rep(0,q),sigma)
  
  Beta <- c(rep(0.3,q1-1),rep(0,q-q1))  
  
  Y1 <- sin(1/6*pi*G) 
  Y2 <- G%*%Beta 
  Y <- cbind(Y1,Y2)+E
  
  wei<- autoweight(G,Y)
  return(wei)
}
Wei_G_case5 <- function(n,q,q1){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2
  
  # gengrate Y
  sigma <- diag(q)
  E <- mvrnorm(n,rep(0,q),sigma)
  
  Beta <- c(rep(0.3,q1-3),rep(0,q-q1))  
  
  Y1 <- 0.3*G^2
  Y2 <- sin(1/6*pi*G)
  Y3 <-  G%*%Beta 
  Y <- cbind(Y1,Y2,Y2,Y3) + E
  
  wei<- autoweight(G,Y)
  return(wei)
}
Wei_G_case6 <- function(n,q,q1){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2
  
  # gengrate Y
  sigma <- diag(q)
  E <- mvrnorm(n,rep(0,q),sigma)
  
  Beta <- c(rep(0.3,q1-4),rep(0,q-q1))
  
  Y1 <- tan(G/2)
  Y2 <- log(G+1) 
  Y3 <- 0.2^G
  Y4 <- 0.1+0.1*G+0.1*G^2
  Y5 <- G%*%Beta 
  Y <- cbind(Y1,Y2,Y3,Y4,Y5) + E
  
  wei<- autoweight(G,Y)
  return(wei)
}