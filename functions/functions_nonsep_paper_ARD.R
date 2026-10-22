#functions for the NONSEP emulator with ARD form for matern kernel

library(extraDistr)
library(Matrix)
library(MCMCpack)
library(fields)
library(Rfast)
library(evaluate)
library(farver)
library(backports)
library(scoringutils)



#function to setup everything to be used for the pca approach
#same input means model is trained using the same inputs for both low and high fidelity levels
fun.pca.setup <- function(Y.list, var.cutoff.list){
  pca.L1 <- prcomp(Y.list[[1]], center=T, scale=T)
  pca.L2 <- prcomp(Y.list[[2]], center=T, scale=T)
  pca.L1.summary <- summary(pca.L1)
  pca.L2.summary <- summary(pca.L2)
  print(pca.L1.summary)
  print(pca.L2.summary)
  dim.red.L1 <- min(which(pca.L1.summary$importance["Cumulative Proportion",]>var.cutoff.list[[1]]))
  dim.red.L2 <- min(which(pca.L2.summary$importance["Cumulative Proportion", ]>var.cutoff.list[[2]]))
  #get matrices with centered and scaled values of the original data
  L1.center <- matrix(pca.L1$center, nrow=nrow(Y.list[[1]]), ncol=ncol(Y.list[[1]]), byrow=T)
  L1.scale <- matrix(pca.L1$scale, nrow=nrow(Y.list[[1]]), ncol=ncol(Y.list[[1]]), byrow=T)
  #level 2 centered and scaled matrices
  L2.center <- matrix(pca.L2$center, nrow=nrow(Y.list[[2]]), ncol=ncol(Y.list[[2]]), byrow=T)
  L2.scale <- matrix(pca.L2$scale, nrow=nrow(Y.list[[2]]), ncol=ncol(Y.list[[2]]), byrow=T)
  #get eigenvectors that will be used for the reduced space
  L1.eig.red <- pca.L1$rotation[,1:dim.red.L1]
  L2.eig.red <- pca.L2$rotation[,1:dim.red.L2]
  #get the data on the reduced space
  L1.red <-((Y.list[[1]] -L1.center)/L1.scale)%*%L1.eig.red
  L2.red <- ((Y.list[[2]] -L2.center)/L2.scale)%*%L2.eig.red
  #get the pca scores from the two levels
  Z.red <- list(L1.red, L2.red)
  colnames(Z.red) <- NULL
  return(list(Y.red.list=list(L1.red, L2.red), Z.red=Z.red, eig.red.list=list(L1.eig.red, L2.eig.red), 
              dim.list=list(dim.red.L1, dim.red.L2), center.list=list(L1.center, L2.center), 
              scale.list=list(L1.scale, L2.scale), eig.all.list=list(pca.L1$rotation, pca.L2$rotation), 
              eig.val.list=list(pca.L1$sdev, pca.L2$sdev), pca.L1=pca.L1, pca.L2=pca.L2))
  
}

#calculate the input correlation function using the Matern function
#with ARD form
R.mat.matern <- function(X, phi, nu=2.5){
  X.scale <- sweep(X, MARGIN = 2, phi, FUN="/")
  dist.X <- fields::rdist(X.scale)
  mat.d <- fields::Matern(d=dist.X, smoothness=nu, phi=1)
  return(mat.d)
}

beta.gls <- function(H, R.inv, Y, W=NA, gam=NA,t=1 ){
  #Y.tilde is just Y for fidelity level 1
  if(t==1){
    Y.tilde <- Y
  }
  #if the fidelity level is greater than 1, we need to incorporate gam W into Y.tilde
  else{
    Y.tilde <- Y-gam*W
  }
  #find the inverse of H^T R.inv H
  p1.inv <- chol2inv(chol(t(H) %*% R.inv %*% H))
  #calculate beta.gls
  beta.gls <- p1.inv %*% t(H) %*% R.inv %*% Y.tilde
  return(beta.gls)
}


#function to give the log posterior of the range parameters
#w is a column vector of the weights (pca scores)
#W is the "design matrix" with the weights (pca scores) from level 1 if applicable
#W is set NULL if it's not used (level 1)
theta.post.pca <- function(phi, w, W=NULL, X, nu, nug_fix=0, phi.hp){
  
  #compute R^-1
  R.chol <- chol(R.mat.matern(X=X, phi=phi, nu=nu) + diag(nug_fix, nrow=nrow(X)))
  R.inv.cur <- chol2inv(R.chol)
  
  n <- nrow(X)
  
  #compute the determinant term
  R.det <- -1*sum(log(diag(R.chol)))
  
  if(is.null(W)==T){
    #compute the quadratic term
    log.q <- -((n)/2)*log(t(w)%*%R.inv.cur%*%w)
    WRW.det <- 0
  }
  else{
    #note that function is based on earlier notation, so W doesn't mean the same thing it does
    #as in this application
    beta.hat <- beta.gls(H=W, R.inv=R.inv.cur, Y=w, W=NA, gam=NA, t=1)
    p1 <- w - W%*%beta.hat
    log.q <- -((n-length(beta.hat))/2)*log(t(p1)%*%R.inv.cur%*%p1)
    WRW.chol <- chol(t(W)%*%R.inv.cur%*%W)
    WRW.det <- -1*sum(log(diag(WRW.chol)))
  }
  
  pl <- R.det + WRW.det + log.q + sum(dhcauchy(phi, sigma=phi.hp, log=T))
  return(list(post=pl, R.inv=R.inv.cur))
}

#function to calculate the map estimator of the log posterior to
#get starting values for mcmc procedure
#function to give the log posterior of the range parameters
#w is a column vector of the weights (pca scores)
#W is the "design matrix" with the weights (pca scores) from level 1 if applicable
#W is set NULL if it's not used (level 1)
#function returns -1*log posterior
theta.post.pca.map <- function(phi, w, W=NULL, X, nu, nug_fix=0, phi.hp){
  
  #compute R^-1
  R.chol <- chol(R.mat.matern(X=X, phi=phi, nu=nu) + diag(nug_fix, nrow=nrow(X)))
  R.inv.cur <- chol2inv(R.chol)
  
  n <- nrow(X)
  
  #compute the determinant term
  R.det <- -1*sum(log(diag(R.chol)))
  
  if(is.null(W)==T){
    #compute the quadratic term
    log.q <- -((n)/2)*log(t(w)%*%R.inv.cur%*%w)
    WRW.det <- 0
  }
  else{
    #note that function is based on earlier notation, so W doesn't mean the same thing it does
    #as in this application
    beta.hat <- beta.gls(H=W, R.inv=R.inv.cur, Y=w, W=NA, gam=NA, t=1)
    p1 <- w - W%*%beta.hat
    log.q <- -((n-length(beta.hat))/2)*log(t(p1)%*%R.inv.cur%*%p1)
    WRW.chol <- chol(t(W)%*%R.inv.cur%*%W)
    WRW.det <- -1*sum(log(diag(WRW.chol)))
  }
  
  pl <- R.det + WRW.det + log.q + sum(dhcauchy(phi, sigma=phi.hp, log=T))
  return(-pl)
}


#function to sample theta for the pca model
mh.step.pca  <- function(post.prev, phi.prev, R.inv.prev, w,W=NULL, X, nu, nug_fix=0, phi.hp, c.0.phi){
  
  #next, update the range parameters phi
  #initialize quantities that will be saved to speed up computations
  phi.up <- phi.prev
  post.up <- post.prev
  a.phi <- 0
  R.inv.up <- R.inv.prev
  
  #generate a proposal for phi
  #values are generated independently of one another
  par.prop <- rnorm(length(phi.prev), mean=phi.prev, sd=c.0.phi)
  par.gz <- which(par.prop<=0)
  
  #make sure the proposal falls in the domain for phi (>0)
  if(length(par.gz)==0){
    phi.prop <- par.prop
    
    #now compute the value of the integrated posterior at the proposal
    post.cur <- theta.post.pca(phi=phi.prop, w=w,W=W, X=X, nu=nu, nug_fix=nug_fix, phi.hp=phi.hp)
    
    post.prop <- post.cur$post
    
    #entire probability ratio for determining acceptance
    pr.phi.cur <- post.prop-post.prev
    
    u <- runif(1, 0, 1)
    if(exp(pr.phi.cur) >u){
      #update the phi parameter
      phi.up <- phi.prop
      post.up <- post.prop
      a.phi <- 1
      R.inv.up <- post.cur$R.inv
    }
    else{
      #don't need to update phi or 
      a.phi <- 0
    }
  }
  #phi proposal isn't in the support of phi
  else{
    #phi, R and S don't need to be updated
    a.phi <- 0
  }
  
  #can add a nugget effect update later if it would be helpful
  return(list(phi=phi.up, post=post.up, a.phi=a.phi, R.inv=R.inv.up))
}


mh.pca <- function(post.st, phi.st, R.inv.st, w,W=NULL, X, nu, nug_fix=0, phi.hp, c.0.phi, times){
  
  phi.mat <- matrix(ncol=ncol(X), nrow=times)
  samp.1 <- mh.step.pca(post.prev=post.st, phi.prev = phi.st, R.inv.prev = R.inv.st, w=w,W=W, X=X, nu=nu, nug_fix=nug_fix, phi.hp=phi.hp, c.0.phi=c.0.phi)
  #get first sample
  phi.mat[1,] <- samp.1$phi
  post.prev <- samp.1$post
  R.inv.prev <- samp.1$R.inv
  for(i in 2:times){
    samp.cur <- mh.step.pca(post.prev=post.prev, phi.prev = phi.mat[i-1,], R.inv.prev = R.inv.prev, w=w,W=W, X=X, nu=nu, nug_fix=nug_fix, phi.hp=phi.hp, c.0.phi=c.0.phi)
    phi.mat[i,] <- samp.cur$phi
    post.prev <- samp.cur$post
    R.inv.prev <- samp.cur$R.inv
    
    if(i%%100 ==0){
      
      print(paste0("On iteration ", i))
    }
  }
  
  return(phi.mat)
}


#functions for a half-normal prior on the range parameters

#function to give the log posterior of the range parameters
#w is a column vector of the weights (pca scores)
#W is the "design matrix" with the weights (pca scores) from level 1 if applicable
#W is set NULL if it's not used (level 1)
#norm.ind gives the indices of phi that should follow a half normal prior
theta.post.pca.npr <- function(phi, w, W=NULL, X, nu, nug_fix=0, phi.hp, norm.ind=c(2)){
  
  #compute R^-1
  R.chol <- chol(R.mat.matern(X=X, phi=phi, nu=nu) + diag(nug_fix, nrow=nrow(X)))
  R.inv.cur <- chol2inv(R.chol)
  
  n <- nrow(X)
  
  #compute the determinant term
  R.det <- -1*sum(log(diag(R.chol)))
  
  if(is.null(W)==T){
    #compute the quadratic term
    log.q <- -((n)/2)*log(t(w)%*%R.inv.cur%*%w)
    WRW.det <- 0
  }
  else{
    #note that function is based on earlier notation, so W doesn't mean the same thing it does
    #as in this application
    beta.hat <- beta.gls(H=W, R.inv=R.inv.cur, Y=w, W=NA, gam=NA, t=1)
    p1 <- w - W%*%beta.hat
    log.q <- -((n-length(beta.hat))/2)*log(t(p1)%*%R.inv.cur%*%p1)
    WRW.chol <- chol(t(W)%*%R.inv.cur%*%W)
    WRW.det <- -1*sum(log(diag(WRW.chol)))
  }
  
  pl <- R.det + WRW.det + log.q + sum(dhcauchy(phi[-norm.ind], sigma=phi.hp[-norm.ind], log=T)) + sum(dnorm(phi[norm.ind], sd=phi.hp[norm.ind], log=T))
  return(list(post=pl, R.inv=R.inv.cur))
}


#function that give negative log posterior to find the map 
#estimator to start the mcmc algorithm
#w is a column vector of the weights (pca scores)
#W is the "design matrix" with the weights (pca scores) from level 1 if applicable
#W is set NULL if it's not used (level 1)
#norm.ind gives the indices of phi that should follow a half normal prior
theta.post.pca.npr.map <- function(par, w, W=NULL, X, nu, nug_fix=0, phi.hp, norm.ind=c(2)){
  phi <- par
  #compute R^-1
  R.chol <- chol(R.mat.matern(X=X, phi=phi, nu=nu) + diag(nug_fix, nrow=nrow(X)))
  R.inv.cur <- chol2inv(R.chol)
  
  n <- nrow(X)
  
  #compute the determinant term
  R.det <- -1*sum(log(diag(R.chol)))
  
  if(is.null(W)==T){
    #compute the quadratic term
    log.q <- -((n)/2)*log(t(w)%*%R.inv.cur%*%w)
    WRW.det <- 0
  }
  else{
    #note that function is based on earlier notation, so W doesn't mean the same thing it does
    #as in this application
    beta.hat <- beta.gls(H=W, R.inv=R.inv.cur, Y=w, W=NA, gam=NA, t=1)
    p1 <- w - W%*%beta.hat
    log.q <- -((n-length(beta.hat))/2)*log(t(p1)%*%R.inv.cur%*%p1)
    WRW.chol <- chol(t(W)%*%R.inv.cur%*%W)
    WRW.det <- -1*sum(log(diag(WRW.chol)))
  }
  
  pl <- R.det + WRW.det + log.q + sum(dhcauchy(phi[-norm.ind], sigma=phi.hp[-norm.ind], log=T)) + sum(dnorm(phi[norm.ind], sd=phi.hp[norm.ind], log=T))
  return(-pl)
}

#function to sample theta for the pca model
mh.step.pca.npr  <- function(post.prev, phi.prev, R.inv.prev, w,W, X, nu, nug_fix=0, phi.hp, c.0.phi, norm.ind=c(2)){
  
  #next, update the range parameters phi
  #initialize quantities that will be saved to speed up computations
  phi.up <- phi.prev
  post.up <- post.prev
  a.phi <- 0
  R.inv.up <- R.inv.prev
  
  #generate a proposal for phi
  #values are generated independently of one another
  par.prop <- rnorm(length(phi.prev), mean=phi.prev, sd=c.0.phi)
  par.gz <- which(par.prop<=0)
  
  #make sure the proposal falls in the domain for phi (>0)
  if(length(par.gz)==0){
    phi.prop <- par.prop
    
    #now compute the value of the integrated posterior at the proposal
    post.cur <- theta.post.pca.npr(phi=phi.prop, w=w,W=W, X=X, nu=nu, nug_fix=nug_fix, phi.hp=phi.hp, norm.ind=norm.ind)
    
    post.prop <- post.cur$post
    
    #entire probability ratio for determining acceptance
    pr.phi.cur <- post.prop-post.prev
    
    u <- runif(1, 0, 1)
    if(exp(pr.phi.cur) >u){
      #update the phi parameter
      phi.up <- phi.prop
      post.up <- post.prop
      a.phi <- 1
      R.inv.up <- post.cur$R.inv
    }
    else{
      #don't need to update phi or 
      a.phi <- 0
    }
  }
  #phi proposal isn't in the support of phi
  else{
    #phi, R and S don't need to be updated
    a.phi <- 0
  }
  
  #can add a nugget effect update later if it would be helpful
  return(list(phi=phi.up, post=post.up, a.phi=a.phi, R.inv=R.inv.up))
}


mh.pca.npr <- function(post.st, phi.st, R.inv.st, w,W, X, nu, nug_fix=0, phi.hp, c.0.phi, norm.ind=c(2), times){
  
  phi.mat <- matrix(ncol=ncol(X), nrow=times)
  samp.1 <- mh.step.pca.npr(post.prev=post.st, phi.prev = phi.st, R.inv.prev = R.inv.st, w=w,W=W, X=X, nu=nu, nug_fix=nug_fix, phi.hp=phi.hp, c.0.phi=c.0.phi, norm.ind=norm.ind)
  #get first sample
  phi.mat[1,] <- samp.1$phi
  post.prev <- samp.1$post
  R.inv.prev <- samp.1$R.inv
  for(i in 2:times){
    samp.cur <- mh.step.pca.npr(post.prev=post.prev, phi.prev = phi.mat[i-1,], R.inv.prev = R.inv.prev, w=w,W=W, X=X, nu=nu, nug_fix=nug_fix, phi.hp=phi.hp, c.0.phi=c.0.phi, norm.ind=norm.ind)
    phi.mat[i,] <- samp.cur$phi
    post.prev <- samp.cur$post
    R.inv.prev <- samp.cur$R.inv
    
    if(i%%100 ==0){
      
      print(paste0("On iteration ", i))
    }
  }
  
  return(phi.mat)
}


#X.pred are the inputs for the predictions
#note that the function doesn't include R.21 since it can be computed from R.12
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

#mean and variance for predictive distributions of weights at a single level
pred.pca.uni <- function(phi, w, X,X.pred, nu, nug_fix=0){
  
  R.all <- R.partition(X=X, X.pred=X.pred, phi=phi, nu=nu, eta=nug_fix, nug_inc=T)
  R.inv <- chol2inv(chol(R.all$R.22))
  pred.m <- R.all$R.12%*%R.inv%*%w
  
  df.pred <- nrow(X)
  pred.var <- R.all$R.11 - R.all$R.12%*%R.inv%*%t(R.all$R.12)
  sig.est <- t(w)%*%R.inv%*%w
  var.final <- (as.numeric(sig.est)/(df.pred))*pred.var
  return(list(mean=pred.m, var=var.final))
  
}

#function to compute the mean and variance of the predictive distribution for level 2
pred.pca.mom.l2 <- function(phi, w, W, W.star, X,X.pred, nu, nug_fix=0){
  
  R.all <- R.partition(X=X, X.pred=X.pred, phi=phi, nu=nu, eta=nug_fix, nug_inc=T)
  R.inv <- chol2inv(chol(R.all$R.22))
  beta.hat <- beta.gls(H=W, R.inv=R.inv, Y=w, W=NA, gam=NA, t=1)
  pred.m <- W.star%*%matrix(beta.hat, nrow=ncol(W.star), ncol=1) + R.all$R.12%*%R.inv%*%(w-W%*%beta.hat)
  
  df.pred <- nrow(X) -ncol(W)
  pred.var1 <- R.all$R.11 - R.all$R.12%*%R.inv%*%t(R.all$R.12)
  pred.var2 <- (W.star -R.all$R.12%*%R.inv%*%W)%*%chol2inv(chol(t(W)%*%R.inv%*%W))%*%t((W.star -R.all$R.12%*%R.inv%*%W))
  sig.est <- t(w- W%*%beta.hat)%*%R.inv%*%(w- W%*%beta.hat)
  var.final <- (as.numeric(sig.est)/(df.pred))*(pred.var1 +pred.var2)
  return(list(mean=pred.m, var=var.final))
}

#function to predict the weights or pca scores
#function is set up to make predictions at each weight at levels 1 and 2
pred.wt.pca <- function(phi1.mat, phi2.mat,  w1, w2, W.mat, X1, X2, X.pred, nu, nug_fix=0){
  mat.1 <- matrix(ncol=nrow(phi1.mat), nrow=nrow(X.pred))
  mat.2 <- matrix(ncol=nrow(phi2.mat), nrow=nrow(X.pred))
  for(i in 1:nrow(phi1.mat)){
    mom.1 <- pred.pca.uni(phi=phi1.mat[i,], w=w1, X=X1, X.pred, nu=nu, nug_fix=nug_fix)
    mat.1[,i] <- extraDistr::rlst(n=nrow(X.pred), df=(nrow(X1)), mu=mom.1$mean, sigma=sqrt(diag(mom.1$var)))
    # mat.1[,i] <- mvtnorm::rmvt(n=1, delta=mom.1$mean, sigma=diag(diag(mom.1$var)), df=(nrow(X1)-1), type="shifted")
    mom.2 <- pred.pca.mom.l2(phi=phi2.mat[i,], w=w2, W=W.mat, W.star=mat.1[,i,drop=F], X=X2, X.pred=X.pred, nu=nu, nug_fix=nug_fix)
    #degrees of freedom are n-1-1 (since we're only using one prediction at each level 2 weight)
    # mat.2[,i] <- mvtnorm::rmvt(n=1, delta=mom.2$mean, sigma=diag(diag(mom.2$var)), df=(nrow(X2)-2), type="shifted")
    mat.2[,i] <- extraDistr::rlst(n=nrow(X.pred), df=(nrow(X2)-1), mu=mom.2$mean, sigma=sqrt(diag(mom.2$var)))
    if(i %% 100 == 0){
      print(paste0("On iteration ",i))
    }
  }
  return(list(L1=mat.1, L2=mat.2))
}

#function to make predictions at the weights in level 1
pred.wt.L1.pca <- function(phi1.mat,  w1, X1, X.pred, nu, nug_fix=0){
  mat.1 <- matrix(ncol=nrow(phi1.mat), nrow=nrow(X.pred))
  for(i in 1:nrow(phi1.mat)){
    mom.1 <- pred.pca.uni(phi=phi1.mat[i,], w=w1, X=X1, X.pred, nu=nu, nug_fix=nug_fix)
    mat.1[,i] <- extraDistr::rlst(n=nrow(X.pred), df=(nrow(X1)), mu=mom.1$mean, sigma=sqrt(diag(mom.1$var)))
    if(i %% 100 == 0){
      print(paste0("On iteration ",i))
    }
  }
  return(mat.1)
}


pred.wt.L2.pca <- function(phi2.mat, w2, W.mat,W.star.ar, X2, X.pred, nu, nug_fix=0){
  mat.2 <- matrix(ncol=nrow(phi2.mat), nrow=nrow(X.pred))
  for(i in 1:nrow(phi2.mat)){
    mom.2 <- pred.pca.mom.l2(phi=phi2.mat[i,], w=w2, W=W.mat, W.star=W.star.ar[,,i], X=X2, X.pred=X.pred, nu=nu, nug_fix=nug_fix)
    #degrees of freedom are n-1-1 (since we're only using one prediction at each level 2 weight)
    # mat.2[,i] <- mvtnorm::rmvt(n=1, delta=mom.2$mean, sigma=diag(diag(mom.2$var)), df=(nrow(X2)-2), type="shifted")
    mat.2[,i] <- extraDistr::rlst(n=nrow(X.pred), df=(nrow(X2)-1), mu=mom.2$mean, sigma=sqrt(diag(mom.2$var)))
    if(i %% 100 == 0){
      print(paste0("On iteration ",i))
    }
  }
  return(mat.2)
}



#function setup so we don't need to save every prediction
#set up for making predictions at level 2 only
#if predictive intervals are going to be used, create a pool of elements for lower and upper
pred.large.pca <-  function(wt.list,eig.red,scale.vec, center.vec, Y.true, pred_int=F, sing_val_add=F, eig.val=NULL, eig.vec.dis=NULL){
  ar.sum <- matrix(0, nrow=nrow(wt.list[[1]]), ncol=ncol(Y.true))
  ar.sum.sq <-matrix(0, nrow=nrow(wt.list[[1]]), ncol=ncol(Y.true))
  ar.dif.sq <- matrix(0, nrow=nrow(wt.list[[1]]), ncol=ncol(Y.true))
  #initialize values so they will be udpated
  ar.quant <- matrix(0, nrow=nrow(wt.list[[1]]), ncol=ncol(Y.true))
  mat.spatial.sum <- matrix(0, nrow=nrow(wt.list[[1]]), ncol=ncol(wt.list[[1]]))
  #ar.L2 <- array(NA,dim= c(nrow(X.pred.list[[2]]), ncol(Y.true), nrow(phi.list[[2]])))
  
  if(pred_int==T){
    #pool for both levels is determined by the total number of samples in samp.ar
    #based on 95% prediction intervals
    pool.size <- round(.025*ncol(wt.list[[1]]))
    mat.pi.low <- matrix(nrow=pool.size, ncol=(nrow(Y.true)*ncol(Y.true)))
    mat.pi.high <- matrix(nrow=pool.size, ncol=(nrow(Y.true)*ncol(Y.true)))
  }
  else{
    mat.pi.low <- NA
    mat.pi.high <-NA
  }
  
  for(i in 1:ncol(wt.list[[1]])){
    
    wt.mat <- matrix(wt.list[[1]][,i], nrow=1, ncol=nrow(wt.list[[1]]))
    for(p in 2:length(wt.list)){
      wt.mat <- rbind(wt.mat, matrix(wt.list[[p]][,i], nrow=1, ncol=nrow(wt.list[[p]])))
    }
    
    if(sing_val_add==F){
      pred_1 <- eig.red%*%wt.mat
      pred.cur <- pred_1*matrix(scale.vec, nrow=ncol(Y.true), ncol=nrow(Y.true), byrow=F) +matrix(center.vec, nrow=ncol(Y.true), ncol=nrow(Y.true), byrow=F)
      pred.cur <- t(pred.cur)
    }
    #add a reconstruction error based on the eigenvalues from the original 
    else{
      pred_1 <- eig.red%*%wt.mat
      #get eigen values for the columns that were discarded for dimension reducution
      eig.val.d <- eig.val[(ncol(eig.red) + 1):length(eig.val)]
      #sample from normal distribution for reconstruction
      eig.norm <- rnorm(n=length(eig.val.d)*ncol(wt.mat), sd=eig.val.d)
      #reconstruct the data using samples of principal components for GP and random noise with sd proportional to eigenvalues
      dat.recon <- t(pred_1) +matrix(eig.norm, nrow=ncol(wt.mat), ncol=length(eig.val.d), byrow=T)%*%t(eig.vec.dis)
      pred.cur <-  t(dat.recon)*matrix(scale.vec, nrow=ncol(Y.true), ncol=nrow(Y.true), byrow=F) +matrix(center.vec, nrow=ncol(Y.true), ncol=nrow(Y.true), byrow=F)
      pred.cur <- t(pred.cur)
    }
    
    
    if(pred_int==T){
      if(i < pool.size){
        mat.pi.low[i,] <- as.vector(pred.cur)
        mat.pi.high[i,] <-as.vector(pred.cur)
      }
      else if(i==pool.size){
        mat.pi.low[i,] <- as.vector(pred.cur)
        mat.pi.high[i,] <-as.vector(pred.cur)
        mat.pi.low <- colSort(mat.pi.low)
        mat.pi.high <- colSort(mat.pi.high)
      }
      else{
        ind.low <- which(mat.pi.low[nrow(mat.pi.low),]>as.vector(pred.cur))
        ind.high <- which(mat.pi.high[1,]<as.vector(pred.cur))
        if(length(ind.low)>0){
          mat.pi.low[nrow(mat.pi.low), ind.low] <- as.vector(pred.cur)[ind.low]
          #sort the columns that were updated
          mat.pi.low[,ind.low] <- colSort(as.matrix(mat.pi.low[,ind.low]))
          
          # for(j in ind.low){
          #   mat.pi.low[,j] <- insert_sorted(vec=mat.pi.low[,j], x=as.vector(pred.cur)[ind.low[j]])
          # }
        }
        if(length(ind.high)>0){
          #row 1 as they are sorted low to high
          mat.pi.high[1, ind.high] <- as.vector(pred.cur)[ind.high]
          #sort the columns that were updated
          mat.pi.high[,ind.high] <- colSort(as.matrix(mat.pi.high[,ind.high]))
          # for(j in ind.high){
          #   mat.pi.high[,j] <- insert_sorted(vec=mat.pi.high[,j], x=as.vector(pred.cur)[ind.high[j]])
          # }
        }
      }
    }
    
    
    
    #quantities needed for mse
    ar.sum <- ar.sum + pred.cur
    ar.sum.sq <- ar.sum.sq + pred.cur^2
    ar.dif.sq <- ar.dif.sq + (pred.cur - Y.true)^2
    
    #helper matrix to keep track of quantiles
    mat.quant.mid <- matrix(0, nrow=nrow(wt.list[[1]]), ncol=ncol(Y.true))
    ind.larger <- which((Y.true-pred.cur)>0, arr.ind = T)
    #anytime Y.true is larger than ar.pred, add a 1
    mat.quant.mid[ind.larger] <- 1
    #keep track of sum of indicators for Y.star.true being bigger
    ar.quant <- ar.quant + mat.quant.mid
    
    #take mean over all spatial locations
    mat.spatial.sum[,i] <- rowmeans(pred.cur)
    
    print(paste0("on iteration ", i))
  }
  pred.mean <- ar.sum/(ncol(wt.list[[1]]))
  pred.var <- ar.sum.sq/(ncol(wt.list[[1]])-1) - ((ncol(wt.list[[1]]))/(ncol(wt.list[[1]])-1))*(ar.sum/ncol(wt.list[[1]])^2)
  pred.mse <- ar.dif.sq/ncol(wt.list[[1]])
  
  if(pred_int==T){
    pred.high <- mat.pi.high[1,]
    pred.low <- mat.pi.low[nrow(mat.pi.low),]
    pred.high.mat <- matrix(pred.high, nrow=nrow(Y.true),ncol=ncol(Y.true),byrow=F)
    pred.low.mat <- matrix(pred.low, nrow=nrow(Y.true), ncol=ncol(Y.true))
  }
  else{
    pred.high <- NA
    pred.low <- NA
    pred.high.mat <- NA
    pred.low.mat <-NA
  }
  rm(mat.pi.low)
  rm(mat.pi.high)
  
  return(list(mean=pred.mean, var=pred.var, mse=pred.mse,quant.ind=(ar.quant/ncol(wt.list[[1]])),mean.spatial=mat.spatial.sum, low=pred.low, high=pred.high))
}


rmspe.by.pcomp <- function(eig.red, wt.mat, center.vec, scale.vec, y.true){
  
  center.mat <- matrix(center.vec, nrow=nrow(y.true), ncol=ncol(y.true), byrow = T)
  scale.mat <- matrix(scale.vec, nrow=nrow(y.true), ncol=ncol(y.true), byrow = T)
  
  rmspe.final <- rep(NA, times=ncol(eig.red))
  for(i in 1:length(rmspe.final)){
    proj1 <- eig.red[,1:i,drop=F]%*%t(wt.mat[,1:i, drop=F])
    proj2 <- t(proj1)*scale.mat +center.mat
    rmspe.final[i] <- sqrt(mean((proj2-y.true)^2))
  }
  return(rmspe.final)
}


mean.cvg.large <- function(mat.mean, mat.true){
  sd.mean <- apply(mat.mean, 1, sd)
  mean.mean <- apply(mat.mean, 1, mean)
  mean.true <- apply(mat.true, 1, mean)
  mse <- (mean.mean-mean.true)^2 +sd.mean^2
  rmse.new <- sqrt(mean((mean.mean-mean.true)^2))
  mean.quant <- apply(mat.mean, 1, quantile, probs=c(.025, .975))
  cvg.ind <- which(mean.true>mean.quant[1,] & mean.true<mean.quant[2,])
  cvg.len <- mean.quant[2,] -mean.quant[1,]
  return(list(cvg=length(cvg.ind)/length(mean.true), cvg.ind=cvg.ind, mse=mse, cvg.len=cvg.len, rmse=rmse.new, sd=sd.mean))
}

#samp is organized with dim (# of test inputs, # spatial locations, # samples)
fun.es <- function(samp, Y.true){
  d1 <- dim(samp)[1]
  d2 <- dim(samp)[2]
  d3 <- dim(samp)[3]
  sum.es.p1 <- rep(0, times=d1)
  sum.es.p2 <-rep(0, times=d1)
  for(i in 1:d3){
    print(paste0("On sample",i))
    s.p1 <- sqrt(rowSums((Y.true-samp[,,i])^2))
    sum.es.p1 <- sum.es.p1 + s.p1
    for(j in i:d3){
      s.p2 <- sqrt(rowSums((samp[,,i] - samp[,,j])^2))
      if(i==j){
        sum.es.p2 <- sum.es.p2 + s.p2
      }
      else{
        sum.es.p2 <- sum.es.p2 + 2*s.p2
      }
    }
  }
  es.f <- sum.es.p1/d3 -.5*sum.es.p2/(d3^2)
  return(es.f)
}

fun.es.fast <- function(samp, Y.true){
  d1 <- dim(samp)[1]
  d2 <- dim(samp)[2]
  es.f <- rep(NA, times=d1)
  for(i in 1:d1){
    es.f[i] <- scoringutils::energy_score_multivariate(observed=Y.true[i,], predicted = samp[i,,], mv_group_id = rep(1, times=d2))
  }
  return(es.f)
}



#function to compute the energy score at level two
pred.energy.score.pca <-  function(wt.list,eig.red,scale.vec, center.vec, Y.true, sing_val_add=F, eig.val=NULL, eig.vec.dis=NULL){
  
  ar.samp <- array(NA, dim=c(nrow(Y.true), ncol(Y.true), ncol(wt.list[[1]])))
  
  for(i in 1:ncol(wt.list[[1]])){
    
    wt.mat <- matrix(wt.list[[1]][,i], nrow=1, ncol=nrow(wt.list[[1]]))
    for(p in 2:length(wt.list)){
      wt.mat <- rbind(wt.mat, matrix(wt.list[[p]][,i], nrow=1, ncol=nrow(wt.list[[p]])))
    }
    
    if(sing_val_add==F){
      pred_1 <- eig.red%*%wt.mat
      pred.cur <- pred_1*matrix(scale.vec, nrow=ncol(Y.true), ncol=nrow(Y.true), byrow=F) +matrix(center.vec, nrow=ncol(Y.true), ncol=nrow(Y.true), byrow=F)
      pred.cur <- t(pred.cur)
    }
    #add a reconstruction error based on the eigenvalues from the original 
    else{
      pred_1 <- eig.red%*%wt.mat
      #get eigen values for the columns that were discarded for dimension reducution
      eig.val.d <- eig.val[(ncol(eig.red) + 1):length(eig.val)]
      #sample from normal distribution for reconstruction
      eig.norm <- rnorm(n=length(eig.val.d)*ncol(wt.mat), sd=eig.val.d)
      #reconstruct the data using samples of principal components for GP and random noise with sd proportional to eigenvalues
      dat.recon <- t(pred_1) +matrix(eig.norm, nrow=ncol(wt.mat), ncol=length(eig.val.d), byrow=T)%*%t(eig.vec.dis)
      pred.cur <-  t(dat.recon)*matrix(scale.vec, nrow=ncol(Y.true), ncol=nrow(Y.true), byrow=F) +matrix(center.vec, nrow=ncol(Y.true), ncol=nrow(Y.true), byrow=F)
      pred.cur <- t(pred.cur)
    }
    
    ar.samp[,,i] <- pred.cur
    print(paste0("on iteration ", i))
  }
  
  #compute the energy score
  time.st <- Sys.time()
  es <- fun.es.fast(samp=ar.samp, Y.true=Y.true)
  time.end <- Sys.time()
  print(time.end-time.st)
  #get rid of large array of samples
  rm(ar.samp)
  #return summary statistics
  return(es)
}


#samp.L1 contains indices of samples for L1 only and
#samp.L2 contains indices of samples for L2 and L1
#hold out one sample from L2 when fitting the model
pca.cv <- function(Y.list, X, comp.num, samp.L1, samp.L2, par.st.L1, par.st.L2, nu=2.5, nug_fix=.00001, phi.hp.L1, phi.hp.L2, norm.ind=1:6){
  
  res.rmspe <- matrix(nrow=length(samp.L2), ncol=comp.num)
  map.L1.ar <- array(dim=c(comp.num, ncol(X), length(samp.L2)))
  map.L2.ar <- array(dim=c(comp.num, ncol(X), length(samp.L2)))
  for(i in 1:length(samp.L2)){
    ic.L2 <- samp.L2[-i]
    ic.L1 <- c(ic.L2, samp.L1)
    
    #perform PCA
    pca.L1 <- prcomp(Y.list[[1]][ic.L1,], center=T, scale=T)
    pca.L2 <- prcomp(Y.list[[2]][ic.L2, ], center=T, scale=T)
    
    #get the weights for the two levels
    w.L1 <- pca.L1$x[,1:comp.num]
    w.L2 <- pca.L2$x[,1:comp.num]
    
    #compute the MAP estimators for each level-1 weight
    pwt.L1 <- rep(NA, times=comp.num)
    for(j in 1:comp.num){
      m1.cur <- optim(par=par.st.L1[j,], fn=theta.post.pca.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                      w=w.L1[,j, drop=F], X=X[ic.L1,], nu=nu, nug_fix=nug_fix, phi.hp=phi.hp.L1)
      map.L1.ar[j,,i] <- m1.cur$par
      #compute the expected value of the L1 weight for the held-out sample
      p1.cur <- pred.pca.uni(phi=m1.cur$par, w=w.L1[,j, drop=F], X=X[ic.L1,], 
                             X.pred= X[samp.L2[i], ,drop=F],nu=nu, nug_fix=nug_fix)
      pwt.L1[j] <- as.numeric(p1.cur$mean)
    }
    #compute the MAP estimators for each level-2 weight
    pwt.L2 <- rep(NA, times=comp.num)
    for(j in 1:comp.num){
      # print(par.start.L2[j,])
      # print(w.L2[,j,drop=F])
      # print(X[ic.L2,])
      # print(w.L1[1:length(ic.L2),j,drop=F])
      m2.cur <- optim(par=par.st.L2[j,], fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                      w=w.L2[,j, drop=F], W=w.L1[1:length(ic.L2),j,drop=F], X=X[ic.L2,], nu=nu, 
                      nug_fix=nug_fix, phi.hp=phi.hp.L2, norm.ind=norm.ind)
      map.L2.ar[j,,i] <- m2.cur$par
      #compute the expected value of the L2 weight for the held-out sample
      # print(dim(w.L2[,j,drop=F]))
      # print(dim(w.L1[1:length(ic.L2),,drop=F]))
      # print(pwt.L1[j])
      p2.cur <- pred.pca.mom.l2(phi=m2.cur$par, w=w.L2[,j, drop=F], W=w.L1[1:length(ic.L2),j,drop=F],  
                                W.star=matrix(pwt.L1[j], nrow=1, ncol=1), X=X[ic.L2,], X.pred=X[samp.L2[i],,drop=F], 
                                nu=nu, nug_fix=nug_fix)
      pwt.L2[j] <- as.numeric(p2.cur$mean)
    }
    wt.mat.L2 <- matrix(pwt.L2, nrow=1, ncol=length(pwt.L2))
    rmspe.prcomp.cur <- rmspe.by.pcomp(eig.red = pca.L2$rotation[,1:comp.num], wt.mat=wt.mat.L2, 
                                       center.vec=pca.L2$center, scale.vec=pca.L2$scale, y.true = Y.list[[2]][samp.L2[i],,drop=F])
    
    res.rmspe[i,] <- rmspe.prcomp.cur
    print(paste0("on input ", i))
  }
  return(list(rmspe=res.rmspe, map.L1=map.L1.ar, map.L2=map.L2.ar))
}


#X are the inputs for level 2
#W has the weights for all nested inputs in both fidelity levels
#while w.L1 has the weights for all inputs at level 1
mean.var.pca.L2 <- function(phi, w, W, W.star, X,X.pred, nu, nug_fix=0,var.L1, X.L1, phi.L1, w.L1){
  
  R.all <- R.partition(X=X, X.pred=X.pred, phi=phi, nu=nu, eta=nug_fix, nug_inc=T)
  R.inv <- chol2inv(chol(R.all$R.22))
  beta.hat <- beta.gls(H=W, R.inv=R.inv, Y=w, W=NA, gam=NA, t=1)
  pred.m <- W.star%*%matrix(beta.hat, nrow=ncol(W.star), ncol=1) + R.all$R.12%*%R.inv%*%(w-W%*%beta.hat)
  
  df.pred.L1 <- nrow(X)
  df.pred <- nrow(X) -ncol(W)
  pred.var1 <- R.all$R.11 - R.all$R.12%*%R.inv%*%t(R.all$R.12)
  
  R.all.L1 <- R.partition(X=X.L1, X.pred=X.pred, phi=phi.L1, nu=nu, eta=nug_fix, nug_inc=T)
  R.inv.L1 <- chol2inv(chol(R.all.L1$R.22))
  
  pred.var2 <- ((R.all.L1$R.12 %*%R.inv.L1%*%w.L1 - R.all$R.12%*%R.inv%*%W))^2
  
  #code only works right now for a single wt in W
  pred.var3 <-  (var.L1 + pred.var2)/as.numeric(t(W)%*%R.inv%*%W)
  
  pred.var4 <- diag(pred.var1) + pred.var3
  
  sig.est <- t(w- W%*%beta.hat)%*%R.inv%*%(w- W%*%beta.hat)
  pred.var5 <- ((as.numeric(sig.est)/(df.pred))*(pred.var4))*(df.pred/(df.pred-2))
  
  pred.var6 <- as.numeric(beta.hat^2)*var.L1*(df.pred.L1/(df.pred.L1-2))
  var.final <- pred.var6 + pred.var5
  
  return(list(mean=pred.m, var=var.final))
}


#can be used with emulators for either one or two levels
#set all L1 inputs to NULL if there's only one level

pred.mean.var.L2 <- function(phi.L1=NULL, phi.L2, w1=NULL, w2, W=NULL, X.L1=NULL, X.L2, X.pred, nu, nug_fix=0, one_L =F){
  
  if(one_L ==T){
    mom.L2 <- pred.pca.uni(phi=phi.L2, w=w2, X=X.L2,X.pred=X.pred, nu=nu, nug_fix=nug_fix)
    mean.wt <- mom.L2$mean
    #degrees of freedom adjustment
    #subtract by 1 because of integrating out variance parameter
    var.wt <- ((nrow(X.L1))/(nrow(X.L1)-2))*diag(mom.L2$var)
    mean.wt.L1 <- mom.L1$mean
    
  }
  else{
    mom.L1 <- pred.pca.uni(phi=phi.L1, w=w1, X=X.L1,X.pred=X.pred, nu=nu, nug_fix=nug_fix)
    #degrees of freedom adjustment
    #subtract by 1 because of integrating out variance parameter
    var.wt.L1 <- ((nrow(X.L1))/(nrow(X.L1)-2))*diag(mom.L1$var)
    mean.wt.L1 <- mom.L1$mean
    mom.L2 <- mean.var.pca.L2(phi=phi.L2, w=w2, W=W, W.star=matrix(mom.L1$mean, nrow=nrow(X.pred),ncol=1), 
                              X=X.L2,X.pred=X.pred, nu=nu, nug_fix=nug_fix, var.L1=var.wt.L1, X.L1=X.L1, phi.L1=phi.L1,w.L1=w1)
    mean.wt <- mom.L2$mean
    var.wt <- mom.L2$var
  }
  return(list(mean=mean.wt, var=var.wt, mean.L1=mean.wt.L1))
  
}


mean.var.proj.L2 <- function(phi.L1.list, phi.L2.list, w1.list, w2.list, W.list, X.list, X.pred, nu, nug_fix=0,one_L, eig.vec, eig.val, center.mat, scale.mat){
  
  mean.ar <- array(dim=c(nrow(X.pred),nrow(eig.vec), length(w2.list)))
  var.ar <- array(dim=c(nrow(X.pred),nrow(eig.vec), length(w2.list)))
  var.recon <- array(dim=c(nrow(X.pred),nrow(eig.vec), length(w2.list)))
  w.mean <- matrix(NA,nrow=nrow(X.pred),ncol=length(phi.L1.list))
  w.L1.mean <- matrix(NA,nrow=nrow(X.pred),ncol=length(phi.L1.list))
  
  for(i in 1:length(w2.list)){
    mom.cur <- pred.mean.var.L2(phi.L1=phi.L1.list[[i]], phi.L2=phi.L2.list[[i]], w1=w1.list[[i]], 
                                w2=w2.list[[i]], W=W.list[[i]], X.L1=X.list[[1]], X.L2=X.list[[2]], 
                                X.pred=X.pred, nu=nu, nug_fix=nug_fix, one_L =one_L)
    w.mean[,i] <- mom.cur$mean
    w.L1.mean[,i] <- mom.cur$mean.L1
    
    if(i==1){
      #project back for a rank 1 approximation
      mean.ar[,,i] <- matrix(mom.cur$mean,nrow=length(mom.cur$mean),ncol=1)%*%t(eig.vec[,1,drop=F])
      #compute the variance (multiply the weight variance by the squared elements of k)
      var.ar[,,i] <- matrix(mom.cur$var,nrow=length(mom.cur$mean),ncol=1)%*%t(eig.vec[,1,drop=F]^2)
    }
    else{
      #project back for a rank i approximation
      mean.ar[,,i] <- mean.ar[,,(i-1)] + matrix(mom.cur$mean,nrow=length(mom.cur$mean),ncol=1)%*%t(eig.vec[,i,drop=F])
      #compute the variance (multiply the weight variance by the squared elements of k)
      var.ar[,,i] <- var.ar[,,(i-1)] + matrix(mom.cur$var,nrow=length(mom.cur$mean),ncol=1)%*%t(eig.vec[,i,drop=F]^2)
    }
    
    
    var.recon[,,i] <- matrix(eig.val[(i+1):length(eig.val)]^2, nrow=1)%*%t(eig.vec[,(i+1):length(eig.val)]^2)
    
  }
  
  #rescale and center the outputs
  mean.ar <- mean.ar*array(scale.mat, c(nrow(scale.mat),ncol(scale.mat),dim(mean.ar)[3])) + array(center.mat, c(nrow(center.mat),ncol(center.mat),dim(mean.ar)[3]))
  var.ar <- var.ar*array(scale.mat, c(nrow(scale.mat),ncol(scale.mat),dim(mean.ar)[3]))^2
  
  var.recon <- var.recon*array(scale.mat, c(nrow(scale.mat),ncol(scale.mat),dim(mean.ar)[3]))^2
  
  var.tot <- var.ar + var.recon
  
  return(list(mean=mean.ar, var=var.tot, var.nr=var.ar, var.recon=var.recon, wt=w.mean, wt.L1=w.L1.mean))
  
}


#function to save samples used to generate plots
pred.pca.all <-  function(wt.list,eig.red,scale.vec, center.vec, Y.true, sing_val_add=F, eig.val=NULL, eig.vec.dis=NULL){
  
  ar.samp <- array(NA, dim=c(nrow(Y.true), ncol(Y.true), ncol(wt.list[[1]])))
  
  for(i in 1:ncol(wt.list[[1]])){
    
    wt.mat <- matrix(wt.list[[1]][,i], nrow=1, ncol=nrow(wt.list[[1]]))
    for(p in 2:length(wt.list)){
      wt.mat <- rbind(wt.mat, matrix(wt.list[[p]][,i], nrow=1, ncol=nrow(wt.list[[p]])))
    }
    
    if(sing_val_add==F){
      pred_1 <- eig.red%*%wt.mat
      pred.cur <- pred_1*matrix(scale.vec, nrow=ncol(Y.true), ncol=nrow(Y.true), byrow=F) +matrix(center.vec, nrow=ncol(Y.true), ncol=nrow(Y.true), byrow=F)
      pred.cur <- t(pred.cur)
    }
    #add a reconstruction error based on the eigenvalues from the original 
    else{
      pred_1 <- eig.red%*%wt.mat
      #get eigen values for the columns that were discarded for dimension reducution
      eig.val.d <- eig.val[(ncol(eig.red) + 1):length(eig.val)]
      #sample from normal distribution for reconstruction
      eig.norm <- rnorm(n=length(eig.val.d)*ncol(wt.mat), sd=eig.val.d)
      #reconstruct the data using samples of principal components for GP and random noise with sd proportional to eigenvalues
      dat.recon <- t(pred_1) +matrix(eig.norm, nrow=ncol(wt.mat), ncol=length(eig.val.d), byrow=T)%*%t(eig.vec.dis)
      pred.cur <-  t(dat.recon)*matrix(scale.vec, nrow=ncol(Y.true), ncol=nrow(Y.true), byrow=F) +matrix(center.vec, nrow=ncol(Y.true), ncol=nrow(Y.true), byrow=F)
      pred.cur <- t(pred.cur)
    }
    
    ar.samp[,,i] <- pred.cur
    print(paste0("on iteration ", i))
  }
  #return samples
  return(ar.samp)
}




