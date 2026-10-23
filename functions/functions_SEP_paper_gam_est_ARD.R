#script with functions to use ARD kernel for the input correlation function
#estimate gamma as well

library(extraDistr)
library(Matrix)
library(Rfast)
library(fields)
library(evaluate)
library(farver)
library(backports)
library(scoringutils)


#Matern correlation function with ARD form
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




#function returns numeric value for dhat and det of V hat that show up in the posterior of theta|y^obs
d.hat.V.det <- function(S, c.j,j, tau_sq, lam_t){
  if(length(c.j)==0){
    dh <- as.numeric(t(S[,j])%*%S[,j] + lam_t)
    det.log <- 0
  }
  else{
    if(tau_sq==0){
      S.sub <- chol(t(S[,c.j])%*%S[,c.j])
    }
    else{
      S.sub <- chol(t(S[,c.j])%*%S[,c.j] + diag(1/tau_sq, nrow=length(c.j)))
    }
    
    dh.1 <- t(S[,j]%*%S[,j])
    #uses faster solving procedure to get second sum of squares
    dh.2 <- backsolve(t(S.sub), x= t(S[,c.j])%*%S[,j], upper.tri=F)
    dh.3 <- t(dh.2)%*%dh.2
    dh <- as.numeric(dh.1 - dh.3)
    dh <- dh + lam_t
    
    #now compute log(det(V.hat))*(-0.5)
    det.log <- -1*sum(log(diag(S.sub)))
  }
  
  return(list(dhat=dh, det=det.log))
}

#compute d_hat and det(V.hat) for every location
#q_t is number of basis functions used to estimate the mean
dV.for <- function(S, c.j.list, q_t, tau_sq, lam_t, eta_t, df_cj=F){
  dh.vec <- rep(NA, times=ncol(S))
  det.vec <- rep(NA, times=ncol(S))
  #compute the degrees of freedom
  df.dh <- -.5*(nrow(S) +eta_t-q_t)
  for(j in 1:ncol(S)){
    dv.cur <- d.hat.V.det(S=S, c.j=c.j.list[[j]], j=j, tau_sq=tau_sq, lam_t=lam_t)
    if(df_cj==T){
      df.cur <- df.dh -.5*length(c.j.list[[j]])
    }
    else{
      df.cur <- df.dh
    }
    dh.vec[j] <- log(dv.cur$dhat)*df.cur
    det.vec[j] <- dv.cur$det
  }
  return(list(dhat=dh.vec, det=det.vec))
}

#compute the posterior distribution of theta and gamma 
#(function used specifically for computing the MAP estimator)
#function returns -1* posterior density
#run at level 1 with U=NULL
#run at level 2 with U= previous outputs of Y_t-1 aligned with inputs from Y
theta.gam.int.post <- function(par, Y, U=NULL,  H, X, c.j.list, nu, tau_sq, lam_t, eta_t, nug_fix=0, df_cj=F, phi.hp){
  if(is.null(U)==T){
    phi <- par
    t <- 1
  }
  else{
    # print("here")
    phi <- par[1:ncol(X)]
    gam <- par[length(par)]
    t <- 2
    # print(phi)
    # print(gam)
  }
  
  
  
  #compute R^-1
  R.chol <- chol(R.mat.matern(X=X, phi=phi, nu=nu) + diag(nug_fix, nrow=nrow(X)))
  R.inv.cur <- chol2inv(R.chol)
  R.inv.chol <- t(backsolve(R.chol, diag(1, nrow=nrow(R.chol))))
  
  #setup if we want to fix the beta parameters
  if(is.null(H)==F){
    Beta.cur <- beta.gls(H=H, R.inv=R.inv.cur, Y=Y, W=U,gam=gam, t=t)
    if(t==1){
      S.cur <- R.inv.chol%*%(Y - H%*%Beta.cur)  
    }
    else{
      S.cur <- R.inv.chol%*%(Y -gam*U - H%*%Beta.cur) 
    }
    
  }
  else{
    if(t==1){
      S.cur <- R.inv.chol%*%Y   
    }
    else{
      S.cur <- R.inv.chol%*%(Y-gam*U)
    }
    
  }
  
  N <- ncol(S.cur)
  if(is.null(H)==FALSE){
    HRH.chol <- chol(t(H) %*% R.inv.cur%*% H)
    HRH.det <- -(N)*sum(log(diag(HRH.chol)))
    q_t <- ncol(H)
  }
  else{
    HRH.det <- 0
    q_t <- 0
  }
  
  R.det <- -N*sum(log(diag(R.chol)))
  
  #compute the log determinant etc
  logdet <- dV.for(S=S.cur, c.j.list=c.j.list, q_t=q_t,tau_sq=tau_sq,lam_t=lam_t, eta_t=eta_t, df_cj=df_cj)
  
  pl <- sum(logdet$dhat) + sum(logdet$det)+HRH.det +R.det +sum(dhcauchy(phi, sigma=phi.hp, log=T)) 
  return(-pl)
}


#computes the integrated value of the posterior for theta and gamma for the Metropolis hastings
#algorithm
theta.gam.int.post.mhg <- function(par, Y, U=NULL, H, X, c.j.list, nu, tau_sq, lam_t, eta_t, nug_fix=0, df_cj=F, phi.hp){
  if(is.null(U)==T){
    phi <- par
    t <- 1
  }
  else{
    phi <- par[1:ncol(X)]
    gam <- par[length(par)]
    t <- 2
  }
  
  #compute R^-1
  R.chol <- chol(R.mat.matern(X=X, phi=phi, nu=nu) + diag(nug_fix, nrow=nrow(X)))
  R.inv.cur <- chol2inv(R.chol)
  R.inv.chol <- t(backsolve(R.chol, diag(1, nrow=nrow(R.chol))))
  
  #setup if we want to fix the beta parameters
  if(is.null(H)==F){
    Beta.cur <- beta.gls(H=H, R.inv=R.inv.cur, Y=Y, W=U,gam=gam, t=t)
    if(t==1){
      S.cur <- R.inv.chol%*%(Y - H%*%Beta.cur)  
    }
    else{
      S.cur <- R.inv.chol%*%(Y -gam*U - H%*%Beta.cur) 
    }
    
  }
  else{
    if(t==1){
      S.cur <- R.inv.chol%*%Y   
    }
    else{
      S.cur <- R.inv.chol%*%(Y-gam*U)
    }
    
  }
  
  N <- ncol(S.cur)
  if(is.null(H)==FALSE){
    HRH.chol <- chol(t(H) %*% R.inv.cur%*% H)
    HRH.det <- -(N)*sum(log(diag(HRH.chol)))
    q_t <- ncol(H)
  }
  else{
    HRH.det <- 0
    q_t <- 0
  }
  
  R.det <- -N*sum(log(diag(R.chol)))
  
  #compute the log determinant etc
  logdet <- dV.for(S=S.cur, c.j.list=c.j.list, q_t=q_t,tau_sq=tau_sq,lam_t=lam_t, eta_t=eta_t, df_cj=df_cj)
  
  pl <- sum(logdet$dhat) + sum(logdet$det)+HRH.det +R.det +sum(dhcauchy(phi, sigma=phi.hp, log=T)) 
  
  return(list(prob=pl, S=S.cur, R.inv=R.inv.cur))
}


#function to update all of the parameters at once
#phi.hp are hyperparameters for the half-cauchy priors on phi
#post.prev gives the value of the integrated posterior from the previous sample
#par.prev is a vector of the previous values of theta and gamma if t>1

mh.gam.gibbs.tInt  <- function(post.prev, par.prev, S.prev, R.inv.prev, Y, U=NULL, H, X, c.j.list, nu, tau_sq, lam_t, eta_t, nug_fix=0, df_cj=F, phi.hp, c.0.phi, c.0.gam=NULL){
  #initialize quantities that will be saved to speed up computation
  if(is.null(U)==T){
    par.up <- par.prev
    t<- 1
    #generate a proposal for phi
    #values are generated independently of one another
    phi.prop <- rnorm(length(par.prev), mean=par.prev, sd=c.0.phi)
    par.gz <- which(phi.prop<=0)
    par.prop <- phi.prop
  }
  else{
    phi.prev <- par.prev[1:ncol(X)]
    gam.prev <- par.prev[length(par.prev)]
    par.up <- par.prev
    t <- 2
    #generate a proposal for phi
    #values are generated independently of one another
    #phi.up is the previous values of phi
    phi.prop <- rnorm(length(phi.prev), mean=phi.prev, sd=c.0.phi)
    par.gz <- which(phi.prop<=0)
    gam.prop <- rnorm(n=1, mean=gam.prev, sd=c.0.gam)
    par.prop <- c(phi.prop, gam.prop)
  }
  post.up <- post.prev
  a.phi <- 0
  S.up <- S.prev
  R.inv.up <- R.inv.prev
  
  
  
  #make sure the proposal falls in the domain for phi (>0)
  if(length(par.gz)==0){
    #now compute the value of the integrated posterior at the proposal
    post.cur <- theta.gam.int.post.mhg(par=par.prop, Y=Y, U=U, H=H, X=X, c.j.list = c.j.list, nu=nu, 
                                       tau_sq = tau_sq, lam_t=lam_t, eta_t=eta_t, nug_fix=nug_fix, df_cj=df_cj, phi.hp=phi.hp)
    post.prop <- post.cur$prob
    
    #entire probability ratio for determining acceptance
    pr.phi.cur <- post.prop-post.prev
    
    u <- runif(1, 0, 1)
    if(exp(pr.phi.cur) >u){
      #update the parameters
      par.up <- par.prop
      post.up <- post.prop
      a.phi <- 1
      S.up <- post.cur$S
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
  return(list(par=par.up, post=post.up, a.phi=a.phi, S=S.up, R.inv=R.inv.up))
}



#function to loop through all fidelity levels indexed using t
#post.list, phi.list are the values of the posterior density and the parameters from the previous iteration
#the dimension of the list for each parameter type corresponds to the number of levels for t
#s is the number of levels
t.gam.loop.mhg.tInt <- function(post.list, par.list,S.prev.list,R.inv.prev.list, Y.list, U.list, H.list, X.list, c.j.list, nu, tau_sq.list, lam_t.list, eta_t.list, nug_fix=0, df_cj=F, phi.hp.list, c.0.phi.list, c.0.gam.list, s){
  #list to store results
  res.t <- list()
  
  for(t in 1:s){
    res.cur <- mh.gam.gibbs.tInt(post.prev=post.list[[t]], par.prev=par.list[[t]],S.prev= S.prev.list[[t]],R.inv.prev=R.inv.prev.list[[t]], Y=Y.list[[t]], U=U.list[[t]], H=H.list[[t]], 
                                 X=X.list[[t]], c.j.list = c.j.list[[t]], nu=nu, tau_sq=tau_sq.list[[t]], lam_t=lam_t.list[[t]], 
                                 eta_t=eta_t.list[[t]], nug_fix=nug_fix, df_cj = df_cj, phi.hp=phi.hp.list[[t]], c.0.phi = c.0.phi.list[[t]], c.0.gam=c.0.gam.list[[t]])
    res.t <- append(res.t, list(res.cur))
  }
  return(res.t)
}

#parameter lists are the starting values for each of the parameters
#H.list, Y.list, X.list, n.list, q.list are lists with the values for each set of parameters at t=1, ..., s
#times is the number of samples we need to collect
#predictions are commented out
#use separate prediction function that doesn't save every predictive sample
algo1.gam.mhg.tInt <- function(post.list, par.list, S.prev.list,R.inv.prev.list, Y.list, U.list, H.list, X.list, c.j.list, nu, tau_sq.list, lam_t.list, eta_t.list, nug_fix=0, df_cj=F, phi.hp.list, c.0.phi.list, c.0.gam.list,s, times, X.pred.list, H.star.list){
  #initialize everything that will change after each iteration
  par.cur <- par.list
  post.cur <- post.list
  S.cur <- S.prev.list
  R.inv.cur <- R.inv.prev.list
  
  #save phi and gamma parameters for one level only
  num.par <- length(par.cur[[1]]) 
  mat.par.1 <- matrix(nrow=times, ncol=(num.par))
  #save predictions using the following array
  #ar.pred.1 <- array(dim=c(nrow(X.pred.list[[1]]), ncol(Y.list[[1]]), times))
  
  #setup so the function can be used with only a single fidelity level
  #save parameters from level 2 (theta and gamma)
  if(s>1){
    #ar.pred.2 <- array(dim=c(nrow(X.pred.list[[2]]), ncol(Y.list[[2]]),times))
    num.par2 <- length(par.cur[[2]])
    mat.par.2 <- matrix(nrow=times, ncol=num.par2)
  }
  else{
    ar.pred.2 <- NA
    mat.par.2 <- NA
  }
  
  
  for(l in 1:times){
    #collect samples from current iteration
    iter.cur <- t.gam.loop.mhg.tInt(post.list=post.cur, par.list=par.cur,S.prev.list=S.cur, R.inv.prev=R.inv.cur, 
                                    Y.list=Y.list, U.list=U.list, H.list=H.list, X.list=X.list, c.j.list=c.j.list, 
                                    nu=nu, tau_sq.list=tau_sq.list, lam_t.list=lam_t.list, eta_t.list=eta_t.list, 
                                    nug_fix=nug_fix, df_cj=df_cj, phi.hp.list=phi.hp.list, c.0.phi.list=c.0.phi.list,
                                    c.0.gam.list=c.0.gam.list, s=s)
    
    mat.par.1[l,] <- c(iter.cur[[1]]$par)
    
    if(s>1){
      mat.par.2[l,] <- c(iter.cur[[2]]$par)
    }
    
    #update the current parameters
    par.cur <- list(iter.cur[[1]]$par)
    post.cur <- list(iter.cur[[1]]$post)
    S.cur <- list(iter.cur[[1]]$S)
    R.inv.cur <- list(iter.cur[[1]]$R.inv)
    
    #update the parameters at the remaining fidelity levels
    if(s > 1){
      for(x in 2:s){
        par.cur <- append(par.cur, list(iter.cur[[x]]$par))
        post.cur <- append(post.cur, list(iter.cur[[x]]$post))
        S.cur <- append(S.cur, list(iter.cur[[x]]$S))
        R.inv.cur <- append(R.inv.cur, list(iter.cur[[x]]$R.inv))
      }
    }
    
    #make predictions using samples from current iteration
    # pred.cur <- t.pred.tInt(X.list=X.list, X.pred.list=X.pred.list, H.list=H.list, H.star.list=H.star.list,
    #                         S.list=S.cur,R.inv.list = R.inv.cur, Y.list=Y.list, phi.list=phi.cur, s=s, c.j.list=c.j.list, nu=nu, tau_sq.list=tau_sq.list,
    #                         lam_t.list=lam_t.list, eta_t.list=eta_t.list, nug_fix=nug_fix, df_cj=df_cj)
    # 
    # #save predictions
    # ar.pred.1[,,l] <- as.matrix(pred.cur$res.list[[1]])
    # 
    # if(s>1){
    #   ar.pred.2[,,l] <- as.matrix(pred.cur$res.list[[2]])
    # }
    
    print(paste0("on iteration", l))
  }
  
  ar.pred.1 <- NA
  ar.pred.2 <- NA
  
  return(list(par.l1=mat.par.1, par.l2=mat.par.2, pred1 = ar.pred.1, pred2=ar.pred.2))
}


#function to compute ahat, dhat, V.hat
ad.hat <-  function(S, c.j,j, tau_sq, lam_t){
  if(length(c.j)==0){
    dh <- as.numeric(t(S[,j])%*%S[,j] + lam_t)
    a.hat <- NA
    V.hat <-NA
  }
  else{
    if(tau_sq==0){
      S.sub <- chol(t(S[,c.j])%*%S[,c.j])
    }
    else{
      S.sub <- chol(t(S[,c.j])%*%S[,c.j] + diag(1/tau_sq, nrow=length(c.j)))
      #print("here")
    }
    dh.1 <- t(S[,j]%*%S[,j])
    #uses faster solving procedure to get second sum of squares
    dh.2 <- backsolve(t(S.sub), x= t(S[,c.j])%*%S[,j], upper.tri=F)
    dh.3 <- t(dh.2)%*%dh.2
    dh <- as.numeric(dh.1 - dh.3)
    dh <- dh + lam_t
    
    #now compute a_hat
    S.sub.inv <- backsolve(S.sub, x=diag(1, nrow=nrow(S.sub)), upper.tri = T)
    #get the mean of each a_t,j
    #note that S.sub.inv %*% dh.2 will give a_hat mean
    a.hat <- S.sub.inv%*%dh.2
    
    #compute V.hat
    V.hat <- chol2inv(S.sub)
  }
  
  return(list(dhat=dh,ahat=a.hat, V.hat=V.hat))
}


#pred_num gives the number of realizations to generate
pred.S.star.j <- function(S.star, S, c.j, j,q_t=1, tau_sq, lam_t, eta_t, df_cj=F, pred_num=1){
  if(df_cj==T){
    df.cur <- nrow(S) -q_t + eta_t + length(c.j)
  }
  else{
    df.cur <- nrow(S) - q_t + eta_t 
  }
  
  if(length(c.j)>0){
    if(nrow(S.star)>1 & length(c.j)>1){
      S.star.cj <- S.star[,c.j]
    }
    else if(nrow(S.star)>1 & length(c.j)==1){
      S.star.cj <- matrix(S.star[,c.j], nrow=nrow(S.star), ncol=1)
    }
    #fixes problem with subsetting matrices with only one row
    else{
      S.star.cj <- matrix(S.star[,c.j], nrow=1, ncol=length(c.j))
    }
    if(length(c.j)==1){
      S.cj <- matrix(S[,c.j], nrow=nrow(S), ncol=1)
      S.j <- matrix(S[,j], nrow=nrow(S), ncol=1)
    }
    else{
      S.cj <- S[,c.j]
      S.j <- S[,j]
    }
    
    ad.cur <- ad.hat(S=S, c.j=c.j, j=j, tau_sq=tau_sq, lam_t=lam_t)
    
    #mu.star.j will be a vector with mean at pred_num points
    mu.star.j <- S.star.cj%*%matrix(ad.cur$ahat, nrow=length(ad.cur$ahat), ncol=1)
    
    #vector of variance terms for each pred_num point
    Sigma.star.j <- (diag(S.star.cj%*% ad.cur$V.hat%*%t(S.star.cj)) + 1)*ad.cur$dhat/df.cur
    
    
  }
  else{
    ad.cur <- ad.hat(S=S, c.j=c.j, j=j, tau_sq=tau_sq, lam_t=lam_t)
    #print(ad.cur)
    #print(paste0("dhat is ", ad.cur$dhat/df.cur))
    mu.star.j <- 0
    Sigma.star.j <- ad.cur$dhat/df.cur
    #print(Sigma.star.j)
  }
  
  #generate new realizations for pred_num points
  #note that the realizations are generated independently
  S.star.samp <- extraDistr::rlst(n=pred_num, df=df.cur, mu=mu.star.j, sigma=sqrt(Sigma.star.j))
  
  return(matrix(S.star.samp, nrow=length(S.star.samp), ncol=1))
  
}


loop.samp.S.star <- function(S, c.j.list, q_t=1, tau_sq, lam_t, eta_t, df_cj=F, pred_num=1){
  
  S.star <- matrix(nrow=pred_num, ncol=ncol(S))
  S.star.cur <- pred.S.star.j(S.star=NA, S=S, c.j=c(), j=1, q_t=q_t, tau_sq=tau_sq, lam_t=lam_t, eta_t=eta_t, df_cj=df_cj, pred_num=pred_num)
  S.star[,1] <- S.star.cur
  
  #update samples for the remaining spatial locations
  for(j in 2:ncol(S)){
    
    #sample from a new location
    S.star.cur <- pred.S.star.j(S.star=S.star, S=S, c.j=c.j.list[[j]], j=j,q_t=q_t, tau_sq=tau_sq, lam_t=lam_t, eta_t=eta_t, pred_num=pred_num)
    
    #update S.star so it's ready for the next sample
    #S.star <- cbind(S.star, S.star.cur)
    S.star[,j] <- S.star.cur
    
  }
  return(S.star)
}

#testing that the function works
# S.star.num10 <- loop.samp.S.star(S=sim2$dat, c.j.list = cj.sim2, q_t=1, tau_sq=.0001, lam_t=.01, eta_t=2, df_cj=F, pred_num=2)

S.star.samp.mult <- function(S, c.j.list, q_t=1,tau_sq, lam_t, eta_t,times){
  
  dat.final <- matrix(nrow=times, ncol=ncol(S))
  for(i in 1:times){
    dat.final[i,] <- loop.samp.S.star(S=S, c.j.list=c.j.list,q_t=q_t, tau_sq=tau_sq, lam_t=lam_t, eta_t=eta_t)
    #print(paste0("on iteration ",i))
  }
  return(dat.final)
}




#function that computes the correlation function for both the original inputs and the inputs for the predictions
#function returns a partitioned version of the correlation matrix
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
}

R.star.star.H <- function(H, H.star, R.inv, R.11, R.12){
  R.star <- R.11 - R.12 %*% R.inv %*% t(R.12)
  if(is.null(H)==F){
    p1 <- t(H) %*%R.inv%*%H
    chol.p1 <- chol(p1)
    p2 <- H.star - R.12%*%R.inv%*%H
    cor.final <- R.star + p2%*%chol2inv(chol.p1)%*%t(p2)
  }
  else{
    cor.final <- R.star
  }
  return(cor.final)
}



#make predictions
#predictions are based on using the product of t distributions
#par includes both the range parameters and scale discrepancy parameter at level-2
pred.gam.lev.tInt <- function(par, X, X.pred,H, H.star, S, R.inv, Y, U=NULL, Y.pred.L1=NULL, eta,t, c.j.list, nu, tau_sq, lam_t, eta_t, nug_fix=0, df_cj=F){
  
  #compute R using old+pred inputs
  #R.11 give correlation matrix at given range parameters for predictions
  #R.22 gives correlation matrix of test inputs at given range parameters
  #R.inv is passed in so inverse doesn't need to be computed again
  if(nug_fix==0){
    nug_inc <- F
    
  }
  else{
    nug_inc <- T
  }
  
  if(is.null(U)==T){
    phi <- par
    t <- 1
  }
  else{
    phi <- par[1:ncol(X)]
    gam <- par[length(par)]
    t <- 2
  }
  
  R.all <- R.partition(X=X, X.pred=X.pred, phi=phi,eta=eta, nu=nu, nug_inc=nug_inc)
  
  if(is.null(H)==F){
    beta.hat <- beta.gls(H=H, R.inv=R.inv, Y=Y,W=U, gam=gam, t=t)
    p1 <- H.star%*%beta.hat
    if(t>1){
      p2 <- Y-gam*U -H%*%beta.hat
      M.ss <- gam*Y.pred.L1 + p1 +R.all$R.12%*%R.inv%*%p2
    }
    else{
      p2 <- Y-H%*%beta.hat
      M.ss <- p1 +R.all$R.12%*%R.inv%*%p2
    }
    q_t <- ncol(H)
  }
  else{
    if(t>1){
      M.ss <- R.all$R.12%*%R.inv%*%(Y-gam*U)  
    }
    else{
      M.ss <- R.all$R.12%*%R.inv%*%Y
    }
    q_t=0
  }
  #compute the conditional correlation matrix
  R.ss <- R.star.star.H(H, H.star, R.inv=R.inv, R.11 = R.all$R.11, R.12=R.all$R.12)
  
  #use values of S to generate a new realization of Y.star
  #note that predictions are independent
  samp.cur <- loop.samp.S.star(S=S, c.j.list=c.j.list, q_t=q_t, tau_sq=tau_sq, lam_t=lam_t, eta_t=eta_t, df_cj=df_cj, pred_num=nrow(X.pred))
  #now multiply the sqrt of the diagonal elements along R.ss and add the sample mean
  Y.star <- diag(sqrt(diag(R.ss)))%*%samp.cur + M.ss
  
  return(list(Y.star=Y.star))
}



#make predictions at multiple fidelity levels
t.gam.pred.tInt <- function(par.list, X.list, X.pred.list, H.list, H.star.list,S.list,R.inv.list, Y.list, U.list, s, c.j.list, nu, tau_sq.list, lam_t.list, eta_t.list, nug_fix=0, df_cj=F){
  res.list <- vector("list", s)
  Y.st.prev <- NULL
  for(t in 1:s){
    pred.cur <- pred.gam.lev.tInt(par=par.list[[t]], X=X.list[[t]], X.pred=X.pred.list[[t]],H=H.list[[t]], H.star=H.star.list[[t]], 
                                  S=S.list[[t]], R.inv=R.inv.list[[t]], Y=Y.list[[t]],Y.pred.L1=Y.st.prev, U=U.list[[t]], eta=nug_fix,t=t, c.j.list=c.j.list[[t]], 
                                  nu=nu, tau_sq=tau_sq.list[[t]], lam_t=lam_t.list[[t]], eta_t=eta_t.list[[t]], nug_fix=nug_fix, df_cj=df_cj)
    #get the Y star values from the last level
    Y.st.prev <- pred.cur$Y.star
    #assumes that scale discrepancy parameter is fixed
    res.list[[t]] <- pred.cur$Y.star
  }
  return(list(res.list=res.list))
}


#function setup so we don't need to save every prediction
#set up for making predictions at level 2 only
#if predictive intervals are going to be used, create a pool of elements for lower and upper
pred.gam.large.tInt <-  function(par.list, X.list, X.pred.list, H.list, H.star.list, Y.list, U.list, s, c.j.list, nu, tau_sq.list,
                                 lam_t.list, eta_t.list, nug_fix, df_cj,Y.true, pred_int=F){
  ar.sum <- matrix(0, nrow=nrow(X.pred.list[[1]]), ncol=ncol(Y.list[[1]]))
  ar.sum.sq <-matrix(0, nrow=nrow(X.pred.list[[1]]), ncol=ncol(Y.list[[1]]))
  ar.dif.sq <- matrix(0, nrow=nrow(X.pred.list[[1]]), ncol=ncol(Y.list[[1]]))
  #initialize values so they will be udpated
  ar.quant <- matrix(0, nrow=nrow(X.pred.list[[1]]), ncol=ncol(Y.list[[1]]))
  mat.spatial.sum <- matrix(0, nrow=nrow(X.pred.list[[1]]), ncol=nrow(par.list[[1]]))
  #ar.L2 <- array(NA,dim= c(nrow(X.pred.list[[2]]), ncol(Y.true), nrow(phi.list[[2]])))
  
  if(pred_int==T){
    #pool for both levels is determined by the total number of samples in samp.ar
    #based on 95% prediction intervals
    pool.size <- round(.025*nrow(par.list[[1]]))
    mat.pi.low <- matrix(nrow=pool.size, ncol=(nrow(X.pred.list[[1]])*ncol(Y.list[[1]])))
    mat.pi.high <- matrix(nrow=pool.size, ncol=(nrow(X.pred.list[[1]])*ncol(Y.list[[1]])))
  }
  else{
    mat.pi.low <- NA
    mat.pi.high <-NA
  }
  
  for(i in 1:nrow(par.list[[1]])){
    
    #compute R and S for the current parameters
    R.inv.list <- vector(mode = "list", length = 2)
    S.list <- vector(mode="list", length=2)
    for(t in 1:2){
      if(t==1){
        phi.cur <- par.list[[t]][i,]
      }
      else{
        phi.cur <- par.list[[t]][i, 1:ncol(X.list[[t]])]
        gam.cur <- par.list[[t]][i,ncol(par.list[[t]])]
      }
      R.chol <- chol(R.mat.matern(X=X.list[[t]], phi=phi.cur, nu=nu) + diag(nug_fix, nrow=nrow(X.list[[t]])))
      R.inv.cur <- chol2inv(R.chol)
      R.inv.chol <- t(backsolve(R.chol, diag(1, nrow=nrow(R.chol))))
      
      R.inv.list[[t]] <- R.inv.cur
      
      #setup if we want to fix the beta parameters (H is NULL)
      if(is.null(H.list[[t]])==F){
        if(t==1){
          Beta.cur <- beta.gls(H=H.list[[t]], R.inv=R.inv.cur, Y=Y.list[[t]],t=t)
          S.cur <- R.inv.chol%*%(Y.list[[t]] - H.list[[t]]%*%Beta.cur)
        }
        else{
          Beta.cur <- beta.gls(H=H.list[[t]], R.inv=R.inv.cur, Y=Y.list[[t]],W=U.list[[t]], gam=gam.cur, t=t)
          S.cur <- R.inv.chol%*%(Y.list[[t]] -gam.cur*U.list[[t]] - H.list[[t]]%*%Beta.cur)
        }
      }
      else{
        if(t==1){
          S.cur <- R.inv.chol%*%Y.list[[t]]
        }
        else{
          S.cur <- R.inv.chol%*%(Y.list[[t]]-gam.cur*U.list[[t]])
        }
      }
      S.list[[t]] <- S.cur
    }
    
    
    #make a prediction using the current values of the range parameter
    #make predictions using samples from current iteration
    #predictions commented out to make sure the sampling is working properly
    pc <- t.gam.pred.tInt(par.list=list(par.list[[1]][i,], par.list[[2]][i,]), X.list=X.list, X.pred.list=X.pred.list, H.list=H.list, H.star.list=H.star.list,
                          S.list=S.list,R.inv.list = R.inv.list, Y.list=Y.list, U.list=U.list, s=s, c.j.list=c.j.list, nu=nu, tau_sq.list=tau_sq.list,
                          lam_t.list=lam_t.list, eta_t.list=eta_t.list, nug_fix=nug_fix, df_cj=df_cj)
    
    pred.cur <- pc$res.list[[2]]
    
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
        }
        if(length(ind.high)>0){
          #row 1 as they are sorted low to high
          mat.pi.high[1, ind.high] <- as.vector(pred.cur)[ind.high]
          #sort the columns that were updated
          mat.pi.high[,ind.high] <- colSort(as.matrix(mat.pi.high[,ind.high]))
        }
      }
    }
    
    
    
    #quantities needed for mse
    ar.sum <- ar.sum + pred.cur
    ar.sum.sq <- ar.sum.sq + pred.cur^2
    ar.dif.sq <- ar.dif.sq + (pred.cur - Y.true)^2
    
    #helper matrix to keep track of quantiles
    mat.quant.mid <- matrix(0, nrow=nrow(X.pred.list[[1]]), ncol=ncol(Y.list[[1]]))
    ind.larger <- which((Y.true-pred.cur)>0, arr.ind = T)
    #anytime Y.true is larger than ar.pred, add a 1
    mat.quant.mid[ind.larger] <- 1
    #keep track of sum of indicators for Y.star.true being bigger
    ar.quant <- ar.quant + mat.quant.mid
    
    #take mean over all spatial locations
    mat.spatial.sum[,i] <- rowmeans(pred.cur)
    
    print(paste0("on iteration ", i))
  }
  pred.mean <- ar.sum/(nrow(par.list[[1]]))
  pred.var <- ar.sum.sq/(nrow(par.list[[1]])-1) - ((nrow(par.list[[1]])/(nrow(par.list[[1]])-1))*(ar.sum/nrow(par.list[[1]]))^2)
  pred.mse <- ar.dif.sq/nrow(par.list[[1]])
  
  if(pred_int==T){
    pred.high <- mat.pi.high[1,]
    pred.low <- mat.pi.low[nrow(mat.pi.low),]
  }
  else{
    pred.high <- NA
    pred.low <- NA
  }
  rm(mat.pi.low)
  rm(mat.pi.high)
  
  return(list(mean=pred.mean, var=pred.var, mse=pred.mse,quant.ind=(ar.quant/nrow(par.list[[1]])),mean.spatial=mat.spatial.sum, low=pred.low, high=pred.high))
}

#mean.cvg with samples from large prediction function
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


fun.es.fast <- function(samp, Y.true){
  d1 <- dim(samp)[1]
  d2 <- dim(samp)[2]
  es.f <- rep(NA, times=d1)
  for(i in 1:d1){
    es.f[i] <- scoringutils::energy_score_multivariate(observed=Y.true[i,], predicted = samp[i,,], mv_group_id = rep(1, times=d2))
  }
  return(es.f)
}


#function setup to compute the energy score
#set up for making predictions at level 2 only
#saves predicted samples to compute the energy score
pred.gam.energy.score.tInt <-  function(par.list, X.list, X.pred.list, H.list, H.star.list, Y.list, U.list, s, c.j.list, nu, tau_sq.list,
                                        lam_t.list, eta_t.list, nug_fix, df_cj,Y.true){
  ar.samp <- array(NA, dim=c(nrow(X.pred.list[[1]]), ncol(Y.list[[1]]), nrow(par.list[[1]])))
  
  time.st <- Sys.time()
  
  for(i in 1:nrow(par.list[[1]])){
    
    #compute R and S for the current parameters
    R.inv.list <- vector(mode = "list", length = 2)
    S.list <- vector(mode="list", length=2)
    for(t in 1:2){
      if(t==1){
        phi.cur <- par.list[[t]][i,]
      }
      else{
        phi.cur <- par.list[[t]][i, 1:ncol(X.list[[t]])]
        gam.cur <- par.list[[t]][i,ncol(par.list[[t]])]
      }
      R.chol <- chol(R.mat.matern(X=X.list[[t]], phi=phi.cur, nu=nu) + diag(nug_fix, nrow=nrow(X.list[[t]])))
      R.inv.cur <- chol2inv(R.chol)
      R.inv.chol <- t(backsolve(R.chol, diag(1, nrow=nrow(R.chol))))
      
      R.inv.list[[t]] <- R.inv.cur
      
      #setup if we want to fix the beta parameters (H is NULL)
      if(is.null(H.list[[t]])==F){
        if(t==1){
          Beta.cur <- beta.gls(H=H.list[[t]], R.inv=R.inv.cur, Y=Y.list[[t]],t=t)
          S.cur <- R.inv.chol%*%(Y.list[[t]] - H.list[[t]]%*%Beta.cur)
        }
        else{
          Beta.cur <- beta.gls(H=H.list[[t]], R.inv=R.inv.cur, Y=Y.list[[t]],W=U.list[[t]], gam=gam.cur, t=t)
          S.cur <- R.inv.chol%*%(Y.list[[t]] -gam.cur*U.list[[t]] - H.list[[t]]%*%Beta.cur)
        }
      }
      else{
        if(t==1){
          S.cur <- R.inv.chol%*%Y.list[[t]]
        }
        else{
          S.cur <- R.inv.chol%*%(Y.list[[t]]-gam.cur*U.list[[t]])
        }
      }
      S.list[[t]] <- S.cur
    }
    
    time.end <- Sys.time()
    print(time.end-time.st)
    #make a prediction using the current values of the range parameter
    #make predictions using samples from current iteration
    #predictions commented out to make sure the sampling is working properly
    pc <- t.gam.pred.tInt(par.list=list(par.list[[1]][i,], par.list[[2]][i,]), X.list=X.list, X.pred.list=X.pred.list, H.list=H.list, H.star.list=H.star.list,
                          S.list=S.list,R.inv.list = R.inv.list, Y.list=Y.list, U.list=U.list, s=s, c.j.list=c.j.list, nu=nu, tau_sq.list=tau_sq.list,
                          lam_t.list=lam_t.list, eta_t.list=eta_t.list, nug_fix=nug_fix, df_cj=df_cj)
    
    ar.samp[,,i] <- pc$res.list[[2]]
    
    print(paste0("on iteration ",i))
  }
  
  time.end <- Sys.time()
  print(time.end-time.st)
  print("done with samples")
  time.st <- Sys.time()
  #use another function to calculate the energy score
  # es <- fun.es(samp=ar.samp, Y.true=Y.true)
  es <- fun.es.fast(samp=ar.samp, Y.true=Y.true)
  time.end <- Sys.time()
  print(time.end-time.st)
  #get rid of large array of samples
  rm(ar.samp)
  #return summary statistics
  return(es)
}

