#code to obtain results for the SEP emulator in the simulation
#results are based on estimating gamma and using the product form for the correlation function
#predictions are made using the map estimators

rm(list=ls())
#set working directory to folder with code
setwd("/supplementary_material_code_R1")

library(GpGp)
library(FNN)
library(fields)

#load the data 
load("data/in.out.1k.rda")
load("data/in.out.pred.1k.rda")

#load helper functions
source("functions/functions_SEP_paper_gam_est.R")

#spatial coordinates
s <- seq(from=0.5, to=5, length.out=20)
t <- seq(from=35, to=60, length.out=50)

#60 inputs at level 1 and 30 inputs at level 2 using conditional Latin hypercube sampling
input.sim.env <- in.out.1k$input

#now look at the outputs
output.sim.env <- in.out.1k$output
#get the outputs at the lower fidelity level that we condition on for the second level 
output.sim.env.W <- in.out.1k$output.W

#inputs for predictions
in.test <- in.out.pred.1k$input

#outputs for predictions
dat.pred.env <- in.out.pred.1k$output


#ordering functions from GpGp package

#takes in vector of spatial locations that was randomly chosen
#use set.seed first to make sure we get the same initial point
sl.order.fun.environ <- function(mat.loc, neigh){
  #use maxmin ordering for the reduced set of locations
  st.order <- order_maxmin(mat.loc)
  mat.red <- mat.loc[st.order, ]
  #find the neigh nearest neighbors
  mat.nn <- find_ordered_nn(locs=mat.red, m=neigh)
  #return a list with the samples ordered based on maxmin ordering and the neigh nearest neighbors
  return(list(order=st.order, nn=mat.nn))
}


#function to get initial values to run algorithm for sampling
#function is modified for the environ function simulation
#in.list is list of inputs at each fidelity level
#out.list is list of outputs at each fidelity level
fun.setup.env <- function(in.list, out.list,out.W,mat.loc=NA, neigh=NULL){
  #use max-min distance ordering for the outputs
  sortnn <- sl.order.fun.environ(mat.loc=mat.loc, neigh=neigh)
  sl <- sortnn$order
  
  #outputs after sorting the locations
  out.list.up <- list(out.list[[1]][,sl], out.list[[2]][,sl])
  out.W.up <- out.W[,sl]
  
  #set up to use only constant basis functions for the mean
  H.list <- list(matrix(1, nrow=nrow(in.list[[1]]), ncol=1), matrix(1,nrow=nrow(in.list[[2]]), ncol=1))
  
  #number of neighbors for each locations
  nj.p1 <- seq(from=0, to=neigh, by=1)
  nj.p2 <- rep(neigh, times=(length(sl)-length(nj.p1)))
  
  #get the neighbors for each locations
  cj <- vector("list",ncol(out.list[[2]]))
  cj[[1]] <- integer(0)
  for(i in 2:ncol(out.list[[2]])){
    ind.col.i <- 1+c(nj.p1, nj.p2)[i]
    loc.cur <- sort(sortnn$nn[i,2:ind.col.i])
    if(is.null(loc.cur)==T){
      cj[[i]] <- integer(0)
    }
    else{
      cj[[i]] <- loc.cur
    }
  }
  
  cj.list <- list(cj, cj)
  
  return(list(H.list=H.list, out.list=out.list.up, out.W=out.W.up,cj.list=cj.list, sl=sl))
}


phi.1.sim <- c(1.81, 0.0135, 0.1767)
phi.2.sim <- c(3.25, 0.004, 0.09)
#now create a list of locations to use for sorting the samples
mat.loc <- matrix(c(rep(s, times=50), rep(t, each=20)), nrow=1000, ncol=2, byrow=F)

#use a single nearest neighbor
#get initial values to run the sampling algorithm
set.seed(376282)
env.loc1 <- fun.setup.env(in.list=input.sim.env, out.list=output.sim.env, out.W = output.sim.env.W, mat.loc=mat.loc, neigh=1)

#calculate MAP estimators
st.e <- Sys.time()
map.L1 <- optim(par=c(11.45542135, 0.03793705, 1.05211330), fn=theta.gam.int.post, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=3),
                Y=env.loc1$out.list[[1]],U=NULL, H=env.loc1$H.list[[1]],X=input.sim.env[[1]],
                c.j.list=env.loc1$cj.list[[1]], nu=2.5, tau_sq=1, lam_t=2, eta_t=4, nug_fix = .000001,
                df_cj=F, phi.hp=c(3,.05,1.5))

map.L2 <- optim(par=c(12.73331421, 0.02252476, 0.75618546,1), fn=theta.gam.int.post, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=3),
                Y=env.loc1$out.list[[2]], U=env.loc1$out.W, H=env.loc1$H.list[[2]],X=input.sim.env[[2]],
                c.j.list=env.loc1$cj.list[[2]], nu=2.5, tau_sq=1, lam_t=2, eta_t=4, nug_fix = .000001,
                df_cj=F, phi.hp=c(3,.05,1.5))

et.e <- Sys.time()
print("Time for estimation: ")
print(et.e-st.e)


#make predictions
pred.num <- 1000
mat.par.L1 <- matrix(map.L1$par, nrow=pred.num, ncol=3, byrow=T)
mat.par.L2 <- matrix(map.L2$par, nrow=pred.num, ncol=4, byrow=T)

st.p <- Sys.time()
set.seed(739203)
pred.a1 <- pred.gam.large.tInt(par.list=list(mat.par.L1, mat.par.L2), X.list=input.sim.env, X.pred.list=list(in.test, in.test),
                               H.list =env.loc1$H.list, H.star.list=list(matrix(1, nrow=30,ncol=1), matrix(1,nrow=30,ncol=1)),
                               Y.list=list(env.loc1$out.list[[1]], env.loc1$out.list[[2]]), U.list =list(NULL, env.loc1$out.W), s=2,
                               c.j.list=env.loc1$cj.list, nu=2.5, tau_sq.list = list(1,1), lam_t.list = list(2, 2), eta_t=list(4,4),
                               nug_fix=.000001, df_cj=F, Y.true=dat.pred.env$hf[,env.loc1$sl], pred_int=T)

et.p <- Sys.time()
print("Time for prediction: ")
print(et.p-st.p)


#summary statistics for predictions
rmspe.m <- round(sqrt(mean((pred.a1$mean-dat.pred.env[[2]][,env.loc1$sl])^2)), digits=3)
cvg.m <- 100*round(length(which(as.vector(dat.pred.env[[2]][,env.loc1$sl])>pred.a1$low & as.vector(dat.pred.env[[2]][,env.loc1$sl])<pred.a1$high))/(nrow(pred.a1$mean)*ncol(pred.a1$mean)),digits=3)
alci.m <- round(mean(pred.a1$high - pred.a1$low), digits=3)

mean.sp <- mean.cvg.large(mat.mean = pred.a1$mean.spatial, mat.true=dat.pred.env[[2]][,env.loc1$sl])
rmspe.sp <- round(mean.sp$rmse, digits=3)
cvg.sp <- 100*round(mean.sp$cvg, digits=3)
alci.sp <- round(mean(mean.sp$cvg.len), digits=3)



#compute the energy score for the test inputs
set.seed(739203)
pred.es <- pred.gam.energy.score.tInt(par.list=list(mat.par.L1, mat.par.L2), X.list=input.sim.env, X.pred.list=list(in.test, in.test),
                               H.list =env.loc1$H.list, H.star.list=list(matrix(1, nrow=30,ncol=1), matrix(1,nrow=30,ncol=1)),
                               Y.list=list(env.loc1$out.list[[1]], env.loc1$out.list[[2]]), U.list =list(NULL, env.loc1$out.W), s=2,
                               c.j.list=env.loc1$cj.list, nu=2.5, tau_sq.list = list(1,1), lam_t.list = list(2, 2), eta_t=list(4,4),
                               nug_fix=.000001, df_cj=F, Y.true=dat.pred.env$hf[,env.loc1$sl])

es <- round(mean(pred.es), digits=2)


res.sep <- matrix(c(rmspe.m, cvg.m, alci.m,es, rmspe.sp, cvg.sp, alci.sp, round((et.e-st.e)/60, digits=2), round((et.p-st.p), digits=2)), nrow=1)
colnames(res.sep) <- c("RMSPE marginal", "CVG (95%)", "ALCI (95%)", "ES", "RMSPE spatial average", "CVG (95%) spatial average", "ALCI (95%) spatial average", "Run time estimation", "Run time prediction")

#results for SEP emulator in Table 1 of the manuscript
print(res.sep)

