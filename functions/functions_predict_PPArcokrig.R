#functions to sample from PP-ARCokrig and PP model

library(extraDistr)
library(Matrix)
library(Rfast)
library(evaluate)
library(farver)
library(backports)
library(scoringutils)


#separable correlation function with matern 
#for nu=2.5
R.mat.matern <- function(X, phi, nu=2.5){
  
  mat.d <- matrix(1, nrow(X), nrow(X))
  
  for(i in 1:(nrow(X)-1)){
    for(j in (i+1):nrow(X)){
      
      dist <- abs(X[i,]-X[j,])
      
      z <- sqrt(5) * dist / phi
      
      prod.ij <- prod((1 + z + z^2/3) * exp(-z))
      
      mat.d[i,j] <- prod.ij
      mat.d[j,i] <- prod.ij
    }
  }
  
  return(mat.d)
}

#y is a vector
b.hat <- function(H, R.inv, y){
  b <- chol2inv(chol(t(H)%*%R.inv%*%H))%*%t(H)%*%R.inv%*%y
  return(matrix(b, nrow=length(b), ncol=1))
}


pred.mom.t <- function(H, H.star, R.inv, R.0, R.00, y){
  if(is.null(ncol(H))==T){
    num.par <- 1
  }
  else{
    num.par <- ncol(H)
  }
  #compute the predictive mean
  coef.mat <- b.hat(H=H, R.inv=R.inv, y)
  m.p1 <- t(R.0)%*%R.inv%*%(y-H%*%coef.mat)
  pred.mean <-H.star%*%coef.mat + m.p1
  #compute the predictive variance
  deg.f <- length(y)-num.par
  sig2.hat <- t(y-H%*%coef.mat)%*%R.inv%*%(y-H%*%coef.mat)/deg.f
  c.p1 <- R.00 -t(R.0)%*%R.inv%*%R.0
  HRH.inv <- chol2inv(chol(t(H)%*%R.inv%*%H))
  c.p2 <- t(H.star)-t(H)%*%R.inv%*%R.0
  c.p3 <- t(c.p2)%*%HRH.inv%*%c.p2
  c.star <- diag(c.p1 + c.p3)
  pred.var <- as.numeric(sig2.hat)*c.star
  return(list(mean=as.vector(pred.mean), var=as.vector(pred.var), df=deg.f))
}

pred.mom.all.loc.1 <- function(H.mat,H.star.mat, R.inv, R.0, R.00, Y){
  mean.mat <- matrix(nrow=nrow(H.star.mat), ncol=ncol(Y))
  sd.mat <- matrix(nrow=nrow(H.star.mat), ncol=ncol(Y))
  #collect samples
  for(j in 1:ncol(Y)){
    mom.cur <- pred.mom.t(H=H.mat[,j,drop=F], H.star = H.star.mat[,j,drop=F], R.inv=R.inv, R.0=R.0, R.00=R.00, y=Y[,j,drop=F])
    mean.mat[,j] <- mom.cur$mean
    sd.mat[,j] <- sqrt(mom.cur$var)
  }
  return(list(mean=mean.mat, sd=sd.mat))
}



pred.samp.lev <- function(H.mat,H.star.mat, R.inv, R.0, R.00, Y, U=NULL,Y1.star=NULL, t){
  samp.mat <- matrix(nrow=nrow(H.star.mat), ncol=ncol(Y))
  if(t==1){
    H.ar <- array(dim=c(nrow(Y), ncol(Y), 1))
    H.ar[,,1] <- H.mat
    H.star.ar <- array(dim=c(nrow(H.star.mat), ncol(Y),1))
    H.star.ar[,,1] <- H.star.mat
  }
  else{
    H.ar <- array(dim=c(nrow(Y), ncol(Y), 2))
    H.ar[,,1] <- H.mat
    H.ar[,,2] <- U
    H.star.ar <- array(dim=c(nrow(H.star.mat), ncol(Y),2))
    H.star.ar[,,1] <- H.star.mat
    H.star.ar[,,2] <- Y1.star
  }
  #collect samples
  for(j in 1:ncol(Y)){
    mom.cur <- pred.mom.t(H=H.ar[,j,], H.star = H.star.ar[,j,], R.inv=R.inv, R.0=R.0, R.00=R.00, y=Y[,j,drop=F])
    samp.mat[,j] <- extraDistr::rlst(n=nrow(R.00), df=mom.cur$df, mu=mom.cur$mean, sigma=sqrt(mom.cur$var))
  }
  return(samp.mat)
}

pred.samp.2L <- function(H.list, H.star.list, R.inv.list, R.0.list, R.00.list, Y.list, U.list){
  samp.1L <- pred.samp.lev(H.mat=H.list[[1]], H.star.mat =H.star.list[[1]], R.inv=R.inv.list[[1]], 
                           R.0=R.0.list[[1]], R.00=R.00.list[[1]], Y=Y.list[[1]], U=NULL, 
                           Y1.star = NULL, t=1)
  samp.2L <- pred.samp.lev(H.mat=H.list[[2]], H.star.mat = H.star.list[[2]], R.inv=R.inv.list[[2]],
                           R.0=R.0.list[[2]], R.00=R.00.list[[2]], Y=Y.list[[2]], U=U.list[[2]], 
                           Y1.star=samp.1L, t=2)
  return(samp.2L)
}

#note that the prediction inputs are the first ones added to the matrix 
#R.11 is for X.pred (predictions)
#R.22 is for X (data)
R.partition <- function(X, X.pred, phi, nu=1.5,eta=0, nug_inc=F){
  if(nug_inc==T){
    R.all <- R.mat.matern(X=rbind(X.pred, X), phi=phi, nu=nu)+ diag(eta, nrow=(nrow(X.pred) + nrow(X)))
  }
  else{
    R.all <- R.mat.matern(X=rbind(X.pred, X), phi=phi, nu=nu) 
  }
  #handle the case where X.pred has only one row
  R.11<-R.all[1:nrow(X.pred), 1:nrow(X.pred)]
  R.12 <-R.all[1:nrow(X.pred), (nrow(X.pred)+1):(nrow(X)+nrow(X.pred))]
  R.22<- R.all[(nrow(X.pred)+1):(nrow(X)+nrow(X.pred)), (nrow(X.pred)+1):(nrow(X)+nrow(X.pred))]
  
  if(nrow(X.pred)==1){
    R.11 <- matrix(R.11, nrow=1,ncol=1)
    R.12 <- matrix(R.12, nrow=1, ncol=nrow(X))
  }
  
  return(list(R.11=R.11, R.12=R.12, R.22=R.22))
  # return(list(R.11=R.all[1:nrow(X.pred), 1:nrow(X.pred)], R.12 =R.all[1:nrow(X.pred), (nrow(X.pred)+1):(nrow(X)+nrow(X.pred))],
  #             R.22=R.all[(nrow(X.pred)+1):(nrow(X)+nrow(X.pred)), (nrow(X.pred)+1):(nrow(X)+nrow(X.pred))]))
}

pred.samp.2L.mult <- function(par.L1, par.L2, H.list, H.star.list, Y.list, U.list,X.list, X.pred,nu=2.5, times){
  Rp.1 <- R.partition(X=X.list[[1]], X.pred=X.pred, phi=par.L1,nu=nu, eta=0, nug_inc = F)
  Rp.2 <- R.partition(X=X.list[[2]], X.pred=X.pred, phi=par.L2,nu=nu, eta=0, nug_inc=F)
  R.inv.list <- list(chol2inv(chol(Rp.1$R.22)), chol2inv(chol(Rp.2$R.22)))
  R.00.list <- list(Rp.1$R.11, Rp.2$R.11)
  R.0.list <- list(t(Rp.1$R.12), t(Rp.2$R.12))
  # print(dim(R.00.list[[1]]))
  # print(dim(R.00.list[[2]]))
  # print(dim(R.0.list[[1]]))
  # print(dim(R.0.list[[2]]))
  # print(dim(R.inv.list[[1]]))
  # print(dim(R.inv.list[[2]]))
  # print(t(R.0.list))
  # 
  
  samp.ar <- array(dim=c(nrow(X.pred), ncol(Y.list[[2]]), times))
  for(i in 1:times){
    samp.ar[,,i] <- pred.samp.2L(H.list=H.list, H.star.list=H.star.list, R.inv.list=R.inv.list, 
                                 R.0.list=R.0.list, R.00.list=R.00.list, Y.list=Y.list, U.list=U.list)
    print(paste0("on iteration ",i))
  }
  
  return(samp.ar)
}


#compute energy score
fun.es.fast <- function(samp, Y.true){
  d1 <- dim(samp)[1]
  d2 <- dim(samp)[2]
  es.f <- rep(NA, times=d1)
  for(i in 1:d1){
    es.f[i] <- scoringutils::energy_score_multivariate(observed=Y.true[i,], predicted = samp[i,,], mv_group_id = rep(1, times=d2))
  }
  return(es.f)
}


pred.es.2L <- function(par.L1, par.L2, H.list, H.star.list, Y.list, U.list,X.list, X.pred,nu=2.5, times, Y.true){
  Rp.1 <- R.partition(X=X.list[[1]], X.pred=X.pred, phi=par.L1,nu=nu, eta=0, nug_inc = F)
  Rp.2 <- R.partition(X=X.list[[2]], X.pred=X.pred, phi=par.L2,nu=nu, eta=0, nug_inc=F)
  R.inv.list <- list(chol2inv(chol(Rp.1$R.22)), chol2inv(chol(Rp.2$R.22)))
  R.00.list <- list(Rp.1$R.11, Rp.2$R.11)
  R.0.list <- list(t(Rp.1$R.12), t(Rp.2$R.12))
  
  samp.ar <- array(dim=c(nrow(X.pred), ncol(Y.list[[2]]), times))
  for(i in 1:times){
    samp.ar[,,i] <- pred.samp.2L(H.list=H.list, H.star.list=H.star.list, R.inv.list=R.inv.list, 
                                 R.0.list=R.0.list, R.00.list=R.00.list, Y.list=Y.list, U.list=U.list)
    print(paste0("on iteration ",i))
  }
  
  es <- fun.es.fast(samp=samp.ar, Y.true=Y.true)
  pred.mean <- apply(samp.ar, MARGIN=c(1,2), mean)
  pred.sd <- apply(samp.ar, MARGIN = c(1,2), sd)
  rm(samp.ar)
  
  return(list(es=es, mean=pred.mean, sd=pred.sd))
}



#function to compute energy scores for the parallel partial model
pred.es.PP <- function(par.L1, H.mat, H.star.mat, Y, X, X.pred,nu=2.5, times, Y.true){
  Rp.1 <- R.partition(X=X, X.pred=X.pred, phi=par.L1,nu=nu, eta=0, nug_inc = F)
  
  R.inv <- chol2inv(chol(Rp.1$R.22))
  R.00 <- Rp.1$R.11
  R.0 <- t(Rp.1$R.12)
  
  #compute the matrix with the predictive mean and scale parameter for each location
  
  mom.L1 <- pred.mom.all.loc.1(H.mat=H.mat, H.star.mat=H.star.mat, R.inv=R.inv, R.0=R.0, R.00=R.00, Y=Y)
  
  samp.ar <- array(dim=c(nrow(X.pred), ncol(Y), times))
  for(i in 1:times){
    #extraDistr::rlst(n=nrow(R.00), df=mom.cur$df, mu=mom.cur$mean, sigma=sqrt(mom.cur$var))
    samp.cur <- matrix(extraDistr::rlst(n=(nrow(Y.true)*ncol(Y.true)), df=(nrow(Y)-1), mu=0, sigma=1), nrow=nrow(Y.true), ncol=ncol(Y.true))
    samp.cur <- samp.cur*mom.L1$sd + mom.L1$mean
    samp.ar[,,i] <- samp.cur
    print(paste0("on iteration ",i))
  }
  
  es <- fun.es.fast(samp=samp.ar, Y.true=Y.true)
  pred.mean <- apply(samp.ar, MARGIN=c(1,2), mean)
  pred.sd <- apply(samp.ar, MARGIN = c(1,2), sd)
  rm(samp.ar)
  
  return(list(es=es, mean=pred.mean, sd=pred.sd))
}


