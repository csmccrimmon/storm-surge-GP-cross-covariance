#code to generate data for the simulation study
#script generates data that is already saved in the data folder

library(clhs)

#level 2 function
environ <- function(xx, s=c(0.5, 1, 1.5, 2, 2.5), t=seq(from=0.3, to=60, by=0.3), t_within_s=T)
{
  
  M   <- xx[1]
  D   <- xx[2]
  L   <- xx[3]
  tau <- xx[4]
  
  ds <- length(s)
  dt <- length(t)
  dY <- ds * dt
  Y <- matrix(0, ds, dt)
  
  # Create matrix Y, where each row corresponds to si and each column
  # corresponds to tj.
  for (ii in 1:ds) {
    si <- s[ii]
    for (jj in 1:dt) {
      tj <- t[jj]
      
      term1a <- M / sqrt(4*pi*D*tj)
      term1b <- exp(-si^2 / (4*D*tj))
      term1 <- term1a * term1b
      
      term2 <- 0
      if (tau < tj) {
        term2a <- M / sqrt(4*pi*D*(tj-tau))
        term2b <- exp(-(si-L)^2 / (4*D*(tj-tau)))
        term2 <- term2a * term2b
      }
      
      C <- term1 + term2
      Y[ii, jj] <- sqrt(4*pi) * C
    }
  }
  
  # Convert the matrix into a vector (by rows).
  #if sorting based on timewithin location
  if(t_within_s==T){
    Yrow <- t(Y)
    y <- t(as.vector(Yrow))
  }
  #sorting based on location within time
  else{
    y <- t(as.vector(Y))
  }
  return(y)
}


#function returns a one term taylor series expansion of first exponential term
#xx is ordered as M, D, L, tau
#xx.0, s.0 and t.0 are the points around which we apply the taylor series expansion
exp.tay.1 <- function(xx, s, t, xx.0=NA, s.0=NA, t.0=NA, method=1){
  
  #taylor series expansion around some point near the true point
  if(method==1){
    
    #in the exponential term for each part
    p1 <- -(s.0^2)/(4*xx.0[2]*t.0)
    #part of the partial derivative with respect to s
    p2 <- -s.0/(2*xx.0[2]*t.0)
    #part of the partial derivative with respect to D
    p3 <- (s.0^2)/(4*(xx.0[2]^2)*t.0)
    #part of the partial derivative with respect to t
    p4 <- (s.0^2)/(4*xx.0[2]*(t.0^2))
    
    final <- exp(p1) +(s-s.0)*exp(p1)*p2 + (xx[2]-xx.0[2])*exp(p1)*p3 + (t - t.0)*exp(p1)*p4
  }
  
  #use the maclaurin series representation of exp(x)
  else{
    p1 <- -(s^2)/(4*xx[2]*t)
    final <- 1+p1+((p1^2)/2) + ((p1^3)/6) 
  }
  
  return(final)
}

#function returns a one term taylor series expansion of first exponential term
#xx is ordered as M, D, L, tau
#xx.0, s.0 and t.0 are the points around which we apply the taylor series expansion
exp.tay.2 <- function(xx, s, t, xx.0=NA, s.0=NA, t.0=NA, method=1){
  
  #use multivariate taylor series expansion around a point near the actual point
  if(method==1){
    #in the exponential term for each part
    p1 <- -((s.0-xx.0[3])^2)/(4*xx.0[2]*(t.0-xx.0[4]))
    #derivative with respect to s (extra part)
    p2 <- -(s.0-xx.0[3])/(2*xx.0[2]*(t.0-xx.0[4]))
    #derivative with respect to L(extra part)
    p3 <- (s.0-xx.0[3])/(2*xx.0[2]*(t.0-xx.0[4]))
    #derivative with respect to D (extra part)
    p4 <- ((s.0-xx.0[3])^2)/(4*(xx.0[2]^2)*(t.0-xx.0[4]))
    #derivative with respect to t (exptra part)
    p5 <- ((s.0-xx.0[3])^2)/(4*xx.0[2]*(t.0-xx.0[4])^2)
    #derivative with respect to tau (extra part)
    p6 <- -((s.0-xx.0[3])^2)/(4*xx.0[2]*(t.0-xx.0[4])^2)
    
    final <- exp(p1) + (s-s.0)*exp(p1)*p2 + (xx[3] - xx.0[3])*exp(p1)*p3 + (xx[2] - xx.0[2])*exp(p1)*p4 + (t- t.0)*exp(p1)*p5 + (xx[4]-xx.0[4])*exp(p1)*p6
  }
  
  #use maclaurin series for exp(x)
  else{
    p1 <- -((s-xx[3])^2)/(4*xx[2]*(t-xx[4]))
    
    final <- 1+ p1 + ((p1^2)/2)+ ((p1^3)/6) 
  }
  return(final)
}


#function is based on a one term taylor series expansion for the two exponential terms
#same as first function other than that
#dif.s, dif.t control how far away from the true point we are when we use the taylor series expansion
#xx.0 is the point around which we apply the taylor series for xx
environ.lf <- function(xx, s=c(0.5, 1, 1.5, 2, 2.5), t=seq(from=0.3, to=60, by=0.3), dif.s=NA, dif.t=NA, xx.0=NA, method=1, t_within_s=T)
{
  M   <- xx[1]
  D   <- xx[2]
  L   <- xx[3]
  tau <- xx[4]
  
  ds <- length(s)
  dt <- length(t)
  dY <- ds * dt
  Y <- matrix(0, ds, dt)
  
  # Create matrix Y, where each row corresponds to si and each column
  # corresponds to tj.
  for (ii in 1:ds) {
    si <- s[ii]
    for (jj in 1:dt) {
      tj <- t[jj]
      
      if(method==1){
        s.0.cur <- si + dif.s
        t.0.cur <- tj + dif.t
      }
      else{
        s.0.cur <-NA
        t.0.cur <-NA
      }
      term1a <- M / sqrt(4*pi*D*tj)
      term1b <- exp.tay.1(xx=xx, s=si, t=tj, xx.0=xx.0, s.0=s.0.cur, t.0=t.0.cur, method = method)
      term1 <- term1a * term1b
      
      term2 <- 0
      if (tau < tj) {
        term2a <- M / sqrt(4*pi*D*(tj-tau))
        term2b <- exp.tay.2(xx=xx, s=si, t=tj, xx.0=xx.0, s.0=s.0.cur, t.0=t.0.cur, method=method)
        term2 <- term2a * term2b
      }
      
      C <- term1 + term2
      Y[ii, jj] <- sqrt(4*pi) * C
    }
  }
  
  # Convert the matrix into a vector (by rows).
  
  #if sorting based on timewithin location
  if(t_within_s==T){
    Yrow <- t(Y)
    y <- t(as.vector(Yrow))
  }
  #sorting based on location within time
  else{
    y <- t(as.vector(Y))
  }
  
  return(y)
}


#function that generates data using the low and high fidelity levels
dat.gen <- function(input.mat, s, t, dif.s=NA, dif.t=NA, dif.xx=NA, method=1, t_within_s=T){
  col.num <- length(s)*length(t)
  mat.lf <- matrix(NA, nrow=nrow(input.mat), ncol=col.num)
  mat.hf <- matrix(nrow=nrow(input.mat), ncol=col.num)
  
  for(x in 1:nrow(input.mat)){
    #determine which points we're using to calculate the taylor series approximation
    xx.0.cur <- input.mat[x,] + dif.xx
    mat.lf[x,] <- environ.lf(xx=input.mat[x,], s=s, t=t, dif.s=dif.s, dif.t=dif.t, xx.0=xx.0.cur, method=method, t_within_s=t_within_s)
    mat.hf[x,] <- environ(xx=input.mat[x,], s=s, t=t, t_within_s=t_within_s)
  }
  return(list(lf=mat.lf, hf=mat.hf))
}


#latin hypercube design from surrogates book
#code modified to handle different input space than [0,1] x [0,1] x ... x [0,1]
#n is number of model runs while 
#m is the size of the input space
#low.vec is a vector of lower bounds for each input
#high.vec is a vector of upper bounds for each input
mylhs <- function(n, m, low.vec, high.vec)
{
  ## generate the Latin hypercube 
  l <- (-(n - 1)/2):((n - 1)/2)
  L <- matrix(NA, nrow=n, ncol=m)
  for(j in 1:m){
    L[,j] <- sample(l, n)
  } 
  mat.min <- matrix(low.vec,nrow=n, ncol=m, byrow=T)
  mat.maxmin <- matrix((high.vec-low.vec),nrow=n, ncol=m, byrow=T)
  
  ## draw the random uniforms and turn the hypercube into a sample
  U <- matrix(runif(n*m), ncol=m)
  X <- (mat.maxmin*((L + (n - 1)/2 + U)/n)) +mat.min
  
  ## return the design and the grid it lives on for visualization
  return(list(X=X, g=c((l + (n - 1)/2)/n,1)))
}


#function to generate new inputs uniformly over the input space
in.uni.gen <- function(n,m, low.vec = c(7, 0.02, 0.01, 30.01), high.vec=c(13, 0.12, 3, 30.295)){
  
  mat <- matrix(nrow=n, ncol=m)
  for(i in 1:n){
    mat[i,] <- runif(n=m, low.vec, high.vec)
  }
  return(mat)
}


#generate the inputs
#use latin hypercube sample to generate the inputs that will be used at level 1
#we'll have 60 inputs at this level
set.seed(28925)
in1 <- mylhs(n=60,m=3, low.vec = c(7, 0.02, 0.01), high.vec=c(13, 0.12, 3))



#write.csv(in1$X, "C:/Users/smccr/Desktop/research_bayes/data/environ_sim_inputs/in1.csv", row.names=F)

#now use a conditional latin hypercube sample to generate the inputs that will be used at the second level
set.seed(1408726)
clhs.l2 <- clhs(as.data.frame(in1$X), size=30)


#generate the inputs for testing uniformly over the entire space
#generate 30 new inputs
set.seed(79425)
in.test <- in.uni.gen(30, m=3, low.vec = c(7, 0.02, 0.01), high.vec=c(13, 0.12, 3))

#spatial coordinates
s <- seq(from=0.5, to=5, length.out=20)
t <- seq(from=35, to=60, length.out=50)

dat1 <- dat.gen(cbind(in1$X, matrix(rep(30, times=60), nrow=60, ncol=1)), s=s, t=t, dif.s=0.075, dif.t=0.075, dif.xx=rep(.075, times=4), method=1, t_within_s=F)

#generate the data for making predictions using the randomly generated inputs
dat.pred.env <- dat.gen(cbind(in.test, matrix(rep(30, times=30), nrow=30, ncol=1)), s=s, t=t, dif.s=0.075, dif.t=0.075, dif.xx=rep(.075, times=4), method=1, t_within_s = F)

#60 inputs at level 1 and 30 inputs at level 2 using conditional Latin hypercube sampling
input.sim.env <- list(in1$X, in2)

#now look at the outputs
output.sim.env <- list(dat1$lf, dat1$hf[clhs.l2,])
#get the outputs at the lower fidelity level that we condition on for the second level 
output.sim.env.W <- dat1$lf[clhs.l2,]

#data used for simulation study that are already saved in the data folder
in.out.1k <- list(input=input.sim.env, output=output.sim.env, output.W=output.sim.env.W)
in.out.pred.1k <- list(input=in.test, output=dat.pred.env)

