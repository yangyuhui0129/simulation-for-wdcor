# GX1_power : beta = 0.3 
# GX2_power : beta ~ N(0,0.4)
G21_power <- function(n,q,q1,rho){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2  
  
  if(rho==0){
    sigma <- diag(q)
  }else{
    sigma<- matrix(0,q,q)
    for(i in 1:q){
      for(j in 1:q){
        sigma[i,j] <- rho^(abs(i-j))
      }
    }
  }
  
  Beta <- c(rep(0.3,q1),rep(0,q-q1))  # rnorm(q1,0,2) # ind.nz <- c(1:6,9:10)
  
  E <- mvrnorm(n,rep(0,q),sigma)  
  Y <- G%*%Beta + E
  ref1 <- dcov.test(G,Y,R=199)$p.value<0.05
  ref2 <- wdcor.test(G,Y,R=199,G.list=seq(1,25,2))$p.value<0.05
  
  return(c(ref1=ref1,ref2=ref2))
}
G31_power <- function(n,q,q1,rho){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2  
  
  if(rho==0){
    sigma <- diag(q)
  }else{
    sigma<- matrix(0,q,q)
    for(i in 1:q){
      for(j in 1:q){
        sigma[i,j] <- rho^(abs(i-j))
      }
    }
  }
  
  E <- mvrnorm(n,rep(0,q),sigma)
  
  Beta <- c(rep(0.3,q1-1),rep(0,q-q1))  # q1-1
  
  Y1 <- 0.3*G^2 
  Y2 <- G%*%Beta 
  Y <- cbind(Y1,Y2)+E
  
  ref1 <- dcov.test(G,Y,R=199)$p.value<0.05
  ref2 <- wdcor.test(G,Y,R=199,G.list=seq(1,25,2))$p.value<0.05
  
  return(c(ref1=ref1,ref2=ref2))
}
G41_power <- function(n,q,q1,rho){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2  
  
  if(rho==0){
    sigma <- diag(q)
  }else{
    sigma<- matrix(0,q,q)
    for(i in 1:q){
      for(j in 1:q){
        sigma[i,j] <- rho^(abs(i-j))
      }
    }
  }
  
  E <- mvrnorm(n,rep(0,q),sigma)
  
  Beta <- c(rep(0.3,q1-1),rep(0,q-q1))  # q1-1
  
  Y1 <- sin(1/6*pi*G)
  Y2 <- G%*%Beta 
  Y <- cbind(Y1,Y2)+E
  
  ref1 <- dcov.test(G,Y,R=199)$p.value<0.05
  ref2 <- wdcor.test(G,Y,R=199,G.list=seq(1,25,2))$p.value<0.05
  
  return(c(ref1=ref1,ref2=ref2))
}
G51_power <- function(n,q,q1,rho){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2  
  
  if(rho==0){
    sigma <- diag(q)
  }else{
    sigma<- matrix(0,q,q)
    for(i in 1:q){
      for(j in 1:q){
        sigma[i,j] <- rho^(abs(i-j))
      }
    }
  }
  
  E <- mvrnorm(n,rep(0,q),sigma)
  
  Beta <- c(rep(0.3,q1-3),rep(0,q-q1))  # q1-6
  
  Y1 <- 0.3*G^2
  Y2 <- sin(1/6*pi*G)
  Y3 <- G%*%Beta 
  Y <- cbind(Y1,Y2,Y2,Y3) + E
  
  ref1 <- dcov.test(G,Y,R=199)$p.value<0.05
  ref2 <- wdcor.test(G,Y,R=199,G.list=seq(1,25,2))$p.value<0.05
  
  return(c(ref1=ref1,ref2=ref2))
}
G61_power <- function(n,q,q1,rho){
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
  
  #Beta <- c(rep(0.3,2),rep(0,q-q1))
  Beta <- c(rep(0.3,q1-4),rep(0,q-q1))
  
  Y1 <- tan(G/2)
  Y2 <- log(G+1) 
  Y3 <- 0.2^G
  Y4 <- 0.1+0.1*G+0.1*G^2
  Y5 <- G%*%Beta 
  #Y <- cbind(Y1,Y1,Y2,Y2,Y3,Y3,Y4,Y4,Y5) + E
  Y <- cbind(Y1,Y2,Y3,Y4,Y5) + E
  
  ref1 <- dcov.test(G,Y,R=199)$p.value<0.05
  ref2 <- wdcor.test(G,Y,R=199,G.list=seq(1,25,2))$p.value<0.05
  
  return(c(ref1=ref1,ref2=ref2))
}
G22_power <- function(n,q,q1,rho){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2  
  
  if(rho==0){
    sigma <- diag(q)
  }else{
    sigma<- matrix(0,q,q)
    for(i in 1:q){
      for(j in 1:q){
        sigma[i,j] <- rho^(abs(i-j))
      }
    }
  }
  
  Beta <- c(rnorm( q1 , 0 , 0.4 ) ,rep(0,q-q1))
  
  E <- mvrnorm(n,rep(0,q),sigma)  
  Y <- G%*%Beta + E
  ref1 <- dcov.test(G,Y,R=199)$p.value<0.05
  ref2 <- wdcor.test(G,Y,R=199,G.list=seq(1,25,2))$p.value<0.05
  
  return(c(ref1=ref1,ref2=ref2))
}
G32_power <- function(n,q,q1,rho){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2  
  
  if(rho==0){
    sigma <- diag(q)
  }else{
    sigma<- matrix(0,q,q)
    for(i in 1:q){
      for(j in 1:q){
        sigma[i,j] <- rho^(abs(i-j))
      }
    }
  }
  
  E <- mvrnorm(n,rep(0,q),sigma)
  
  Beta <- c(rnorm( q1-1 , 0 , 0.4) ,rep(0,q-q1))  # q1-1
  
  Y1 <- (rnorm(1,0,0.4)) * G^2 
  Y2 <- G%*%Beta 
  Y <- cbind(Y1,Y2)+E
  
  ref1 <- dcov.test(G,Y,R=199)$p.value<0.05
  ref2 <- wdcor.test(G,Y,R=199,G.list=seq(1,25,2))$p.value<0.05
  
  return(c(ref1=ref1,ref2=ref2))
}
G42_power <- function(n,q,q1,rho){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2  
  
  if(rho==0){
    sigma <- diag(q)
  }else{
    sigma<- matrix(0,q,q)
    for(i in 1:q){
      for(j in 1:q){
        sigma[i,j] <- rho^(abs(i-j))
      }
    }
  }
  
  E <- mvrnorm(n,rep(0,q),sigma)
  
  Beta <- c(rnorm( q1-1 , 0 , 0.4 ) ,rep(0,q-q1))  # q1-1
  
  Y1 <- sin(1/6*pi*G)
  Y2 <- G%*%Beta 
  Y <- cbind(Y1,Y2) + E
  
  ref1 <- dcov.test(G,Y,R=199)$p.value<0.05
  ref2 <- wdcor.test(G,Y,R=199,G.list=seq(1,25,2))$p.value<0.05
  
  return(c(ref1=ref1,ref2=ref2))
}
G52_power <- function(n,q,q1,rho){
  X <- matrix(runif(n),n,1)
  G <- matrix(1,n,1)
  maf <- runif(1,0.3,0.5)
  Q1 <- maf^2
  Q2 <- 1-(1-maf)^2
  G[X<Q1] <-0
  G[X>Q2] <-2  
  
  if(rho==0){
    sigma <- diag(q)
  }else{
    sigma<- matrix(0,q,q)
    for(i in 1:q){
      for(j in 1:q){
        sigma[i,j] <- rho^(abs(i-j))
      }
    }
  }
  
  E <- mvrnorm(n,rep(0,q),sigma)
  
  Beta <- c(rnorm(q1-3,0,0.4) ,rep(0,q-q1))  # q1-6
  
  Y1 <- G^2 *rnorm(1,0,0.4)
  Y2 <- sin(1/6*pi*G)
  Y3 <- G%*%Beta 
  Y <- cbind(Y1,Y2,Y2,Y3) + E
  
  ref1 <- dcov.test(G,Y,R=199)$p.value<0.05
  ref2 <- wdcor.test(G,Y,R=199,G.list=seq(1,25,2))$p.value<0.05
  
  return(c(ref1=ref1,ref2=ref2))
}
G62_power <- function(n,q,q1,rho){
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
  
  #Beta <- c(rep(0.3,2),rep(0,q-q1))
  Beta <- c(rep(0.3,q1-4),rep(0,q-q1))
  
  Y1 <- tan(G/2)
  Y2 <- log(G+1) 
  Y3 <- 0.2^G
  Y4 <- 0.1+0.1*G+0.1*G^2
  Y5 <- G%*%Beta 
  #Y <- cbind(Y1,Y1,Y2,Y2,Y3,Y3,Y4,Y4,Y5) + E
  Y <- cbind(Y1,Y2,Y3,Y4,Y5) + E
  
  ref1 <- dcov.test(G,Y,R=199)$p.value<0.05
  #ref2 <- wdcor.test(G,Y,R=199,G.list=seq(1,25,2))$p.value<0.05
  
  return(c(ref1=ref1,ref2=0))
}