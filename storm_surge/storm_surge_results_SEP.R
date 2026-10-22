#script to obtain results for sep emulator in the paper
#uses ARD form for input correlation function and MAP estimator

rm(list=ls())
#set working directory to folder with code
setwd("/supplementary_material_code_R1")

#libraries
library(fields)
library(GpGp)
library(FNN)
library(xtable)
library(evaluate)
library(knitr)

#load the data
load("data/StormSurge.RData")

#load helper functions
source("functions/functions_SEP_paper_gam_est_ARD.R")

#setup for the data
lon0 = min(inputs[ ,6])
lat0 = max(inputs[, 7])
reference = matrix(c(lon0, lat0), nrow=1)
dist = c(rdist.earth(reference, inputs[ ,c(6,7)]))
inputs.6Dim = cbind(inputs[ ,1:5], dist)

#sample of 120 low-fidelity test inputs
set.seed(193723)
samp.120 <- sample(1:226,120)

#available inputs for L2
samp.subl2 <- which(seq(from=1, to=226, by=1) %in% samp.120)

#select 60 inputs to use as high fidelity test inputs
set.seed(673627)
samp.60 <- sample(samp.subl2, 60)

#get test inputs
samp.pred <- which(!(seq(from=1, to=226, by=1) %in% samp.120))

#number of locations used
samp.loc <- 1:9284

#test inputs for predictions
X.new.comp <- list(inputs.6Dim[sort(samp.pred),], inputs.6Dim[sort(samp.pred),])

#test outputs for predictions
output.pred <- list(PSE.lf[sort(samp.pred),samp.loc], PSE.hf[sort(samp.pred), samp.loc])

#maximum-minimum distance ordering and selection of nearest neighbors
#use set.seed first to make sure we get the same initial point
#ll.mat is a matrix where each row is the location (longitude and lattitude) of the inputs
sl.order.fun <- function(samp.l, neigh, ll.mat){
  #use maxmin ordering for the reduced set of locations
  lonlat.order <- order_maxmin(locs=ll.mat[samp.l,], lonlat=T)
  samp.order <- samp.l[lonlat.order]
  lonlat.red <-ll.mat[samp.order, ]
  #find the neigh nearest neighbors
  lonlat.nn <- find_ordered_nn(locs=lonlat.red, m=neigh, lonlat=T)
  #return a list with the samples ordered based on maxmin ordering and the neigh nearest neighbors
  return(list(order=samp.order, nn=lonlat.nn))
}

#function that organizes the data needed for the sampling algorithm
#function to get initial values to run algorithm for sampling
#fix out.list for uneven number of training data points
fun.setup <- function(samp.l, ll.mat, neigh, samp.i.list, input.mat, out.L1, out.L2){
  #sort distances and find nearest neighbors
  sortnn <- sl.order.fun(samp.l=samp.l, neigh=neigh, ll.mat=ll.mat)
  sl <- sortnn$order
  
  in.list <- list(input.mat[samp.i.list[[1]],], input.mat[samp.i.list[[2]],])
  out.list <- list(out.L1[samp.i.list[[1]],sl], out.L2[samp.i.list[[2]],sl])
  out.W <- out.L1[samp.i.list[[2]], sl]
  H.list <- list(matrix(1, nrow=length(samp.i.list[[1]]), ncol=1), matrix(1,nrow=length(samp.i.list[[2]]), ncol=1))
  
  #create neighbor sets at each fidelity level
  #assumes that neighbor sets are identical for each fidelity level
  
  #number of neighbors for each location
  nj.p1 <- seq(from=0, to=neigh, by=1)
  nj.p2 <- rep(neigh, times=(length(sl)-length(nj.p1)))
  
  cj <- vector("list",length(samp.l))
  cj[[1]] <- integer(0)
  for(i in 2:length(samp.l)){
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
  
  return(list(in.list=in.list, out.list=out.list, out.W=out.W, H.list=H.list, cj.list=cj.list, sl=sl))
}

set.seed(123)
set.loc1 <- fun.setup( samp.l=1:9284, ll.mat=lonlat, neigh=1, samp.i.list = list(samp.120, samp.60), input.mat=inputs.6Dim, out.L1=PSE.lf, out.L2=PSE.hf)

range.L1.inputs <- apply(inputs.6Dim[samp.120,], 2, range)[2,]-apply(inputs.6Dim[samp.120,], 2, range)[1,]
range.L2.inputs <- apply(inputs.6Dim[samp.60,], 2, range)[2,]-apply(inputs.6Dim[samp.60,], 2, range)[1,]


st.e <- Sys.time()
#starting values for level 1 and level 2 range parameters based on map estimator
phi.map.st.L1 <- c(33.98, 15.24, 4.58, 18.57, 0.87,9.62)
phi.map.st.L2 <- c(30.04, 12.86, 12.14,25.87, 0.84,13.48)

#find MAP estimator of theta and gamma
#level 1 no gamma
map.L1 <- optim(par=phi.map.st.L1, fn=theta.gam.int.post, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                Y=set.loc1$out.list[[1]], U=NULL, H=set.loc1$H.list[[1]],X=set.loc1$in.list[[1]],
                c.j.list=set.loc1$cj.list[[1]], nu=2.5, tau_sq=10000000, lam_t=.02, eta_t=.2, nug_fix = .000001,
                df_cj=F, phi.hp=range.L1.inputs/2)
# #level 2 estimate gamma as well as theta
map.L2 <- optim(par=c(phi.map.st.L2,1), fn=theta.gam.int.post, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                Y=set.loc1$out.list[[2]], U=set.loc1$out.W, H=set.loc1$H.list[[2]],X=set.loc1$in.list[[2]],
                c.j.list=set.loc1$cj.list[[2]], nu=2.5, tau_sq=10000000, lam_t=.002, eta_t=.2, nug_fix = .000001,
                df_cj=F, phi.hp=range.L1.inputs/2)

et.e <- Sys.time()
print("Time for estimation: ")
print(et.e-st.e)


#make predictions 
pred.num <- 1000
mat.par.L1 <- matrix(phi.map.st.L1, nrow=pred.num, ncol=6, byrow=T)
mat.par.L2 <- matrix(c(phi.map.st.L2,1), nrow=pred.num, ncol=7, byrow=T)
st.p <- Sys.time()
set.seed(39403)
pred.gam.a4 <- pred.gam.large.tInt(par.list=list(mat.par.L1, mat.par.L2), X.list=set.loc1$in.list, X.pred.list=X.new.comp,
                                   H.list=set.loc1$H.list, H.star.list=list(matrix(1,nrow=106,ncol=1),matrix(1,nrow=106,ncol=1) ),
                                   Y.list=list(set.loc1$out.list[[1]], set.loc1$out.list[[2]]), U.list=list(NULL, set.loc1$out.W), s=2, c.j.list=set.loc1$cj.list, nu=2.5,
                                   tau_sq.list = list(10000000,10000000), lam_t.list = list(.02, .002), eta_t=list(.2,.2), nug_fix=.00001, df_cj=F,Y.true=output.pred[[2]][,set.loc1$sl], pred_int = T)
et.p <- Sys.time()
print("Time for prediction: ")
print(et.p-st.p)


#summary statistics for predictive performance of the SEP emulator
rmspe.m <- round(sqrt(mean((pred.gam.a4$mean-output.pred[[2]][,set.loc1$sl])^2)), digits=3)
cvg.m <- 100*round(length(which(pred.gam.a4$low<output.pred[[2]][,set.loc1$sl] & pred.gam.a4$high > output.pred[[2]][,set.loc1$sl]))/(9284*106), digits=4)
alci.m <- round(mean(pred.gam.a4$high-pred.gam.a4$low), digits=3)

mean.sum.a4 <- mean.cvg.large(mat.mean = pred.gam.a4$mean.spatial, mat.true=output.pred[[2]])

rmspe.sp <- round(mean.sum.a4$rmse, digits=3)
cvg.sp <- 100*round(mean.sum.a4$cvg, digits=4)
alci.sp <- round(mean(mean.sum.a4$cvg.len), digits=3)

et.e-st.e

res.sep <- matrix(c(rmspe.m, cvg.m, alci.m, rmspe.sp, cvg.sp, alci.sp, round((et.e-st.e), digits=2), round(et.p-st.p, digits=2)), nrow=1)
colnames(res.sep) <- c("RMSPE marginal", "CVG (95%)", "ALCI (95%)", "RMSPE spatial average", "CVG (95%) spatial average", "ALCI (95%) spatial average", "Run time estimation", "Run time prediction")


#results for SEP model in table 2 of the manuscript 
print(res.sep)

#results for the energy score are provided at the end of the script and require additional 
#memory to compute


#compute the energy score
#each of the 1,000 samples are saved for every test input and spatial location
#running this code requires around 8 GB of memory
pred.num <- 1000
mat.par.L1 <- matrix(map.L1$par, nrow=pred.num, ncol=6, byrow=T)
mat.par.L2 <- matrix(map.L2$par, nrow=pred.num, ncol=7, byrow=T)


set.seed(39403)
pred.gam.es <- pred.gam.energy.score.tInt(par.list=list(mat.par.L1, mat.par.L2), X.list=set.loc1$in.list, X.pred.list=X.new.comp,
                                          H.list=set.loc1$H.list, H.star.list=list(matrix(1,nrow=106,ncol=1),matrix(1,nrow=106,ncol=1) ),
                                          Y.list=list(set.loc1$out.list[[1]], set.loc1$out.list[[2]]), U.list=list(NULL, set.loc1$out.W), s=2, c.j.list=set.loc1$cj.list, nu=2.5,
                                          tau_sq.list = list(10000000,10000000), lam_t.list = list(.02, .002), eta_t=list(.2,.2), nug_fix=.00001, df_cj=F,Y.true=output.pred[[2]][,set.loc1$sl])

#results for energy score for SEP model in Table 2 of the manuscript
round(mean(pred.gam.es), digits=2)





