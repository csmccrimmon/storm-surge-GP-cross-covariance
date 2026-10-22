#code to obtain results for the nonsep emulator in the simulation study
#uses ARD form for the input correlation function and predictions are made
#using the map estimator

rm(list=ls())
#set working directory to folder with code
#setwd("/supplementary_material_code_R1")

#libraries
library(xtable)

#load the data 
load("data/in.out.1k.rda")
load("data/in.out.pred.1k.rda")

#load helper functions
source("functions/functions_nonsep_paper_ARD.R")

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


#determine which rows of the matrix with the inputs for level 1 match the order of the inputs in level 2
in.hf <- match(
  apply(input.sim.env[[2]], 1, function(x) paste0(sprintf("%.17g", x), collapse = ",")),
  apply(input.sim.env[[1]], 1, function(x) paste0(sprintf("%.17g", x), collapse = ","))
)
#check that it's right
input.sim.env[[1]][in.hf,]-input.sim.env[[2]]

#determine the row indices for the remaining low-fidelity only inputs
in.lf <- which(!(1:60 %in% in.hf))

sL1.sort <- c(in.hf, in.lf)

#use pca setup function to get everthing ready to run
set.pca <- fun.pca.setup(Y.list = list(output.sim.env[[1]][sL1.sort,],output.sim.env[[2]]),  var.cutoff.list=list(0.99999, 0.99995))

phi.start <- c(11.475, .038, 1.059)
range.L1.inputs <- c(6, .1, 2.99)
range.L2.inputs <- c(6, .1, 2.99)

#calculate MAP estimators
st.e <- Sys.time()
map.w1 <- optim(par=phi.start, fn=theta.post.pca.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                w=matrix(set.pca$Z.red[[1]][,1], ncol=1, nrow=nrow(set.pca$Z.red[[1]])), W=NULL,
                X=input.sim.env[[1]][sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w2 <- optim(par=phi.start, fn=theta.post.pca.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                w=matrix(set.pca$Z.red[[1]][,2], ncol=1, nrow=nrow(set.pca$Z.red[[1]])), W=NULL,
                X=input.sim.env[[1]][sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w3 <- optim(par=phi.start, fn=theta.post.pca.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                w=matrix(set.pca$Z.red[[1]][,3], ncol=1, nrow=nrow(set.pca$Z.red[[1]])), W=NULL,
                X=input.sim.env[[1]][sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w4 <- optim(par=phi.start, fn=theta.post.pca.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                w=matrix(set.pca$Z.red[[1]][,4], ncol=1, nrow=nrow(set.pca$Z.red[[1]])), W=NULL,
                X=input.sim.env[[1]][sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w5 <- optim(par=phi.start, fn=theta.post.pca.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                w=matrix(set.pca$Z.red[[1]][,5], ncol=1, nrow=nrow(set.pca$Z.red[[1]])), W=NULL,
                X=input.sim.env[[1]][sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w6 <- optim(par=phi.start, fn=theta.post.pca.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                w=matrix(set.pca$Z.red[[1]][,6], ncol=1, nrow=nrow(set.pca$Z.red[[1]])),W=NULL,
                X=input.sim.env[[1]][sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w7 <- optim(par=phi.start, fn=theta.post.pca.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                w=matrix(set.pca$Z.red[[1]][,7], ncol=1, nrow=nrow(set.pca$Z.red[[1]])), W=NULL,
                X=input.sim.env[[1]][sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w8 <- optim(par=phi.start, fn=theta.post.pca.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                w=matrix(set.pca$Z.red[[1]][,8], ncol=1, nrow=nrow(set.pca$Z.red[[1]])), W=NULL,
                X=input.sim.env[[1]][sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

#map estimators for level two
map.l2.w1 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                   w=matrix(set.pca$Z.red[[2]][,1], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:30,1,drop=F],
                   X=input.sim.env[[2]], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:3)

map.l2.w2 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                   w=matrix(set.pca$Z.red[[2]][,2], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:30,2,drop=F],
                   X=input.sim.env[[2]], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:3)

map.l2.w3 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                   w=matrix(set.pca$Z.red[[2]][,3], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:30,3,drop=F],
                   X=input.sim.env[[2]], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2,norm.ind=1:3)

map.l2.w4 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                   w=matrix(set.pca$Z.red[[2]][,4], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:30,4,drop=F],
                   X=input.sim.env[[2]], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2,norm.ind=1:3)

map.l2.w5 <- optim(par=phi.start/10, fn=theta.post.pca.npr.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                   w=matrix(set.pca$Z.red[[2]][,5], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:30,5,drop=F],
                   X=input.sim.env[[2]], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2,norm.ind=1:3)

map.l2.w6 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                   w=matrix(set.pca$Z.red[[2]][,6], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:30,6,drop=F],
                   X=input.sim.env[[2]], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2,norm.ind=1:3)

map.l2.w7 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                   w=matrix(set.pca$Z.red[[2]][,7], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:30,7,drop=F],
                   X=input.sim.env[[2]], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2,norm.ind=1:3)

map.l2.w8 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B", lower=rep(.0001, times=3),
                   w=matrix(set.pca$Z.red[[2]][,8], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:30,8,drop=F],
                   X=input.sim.env[[2]], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2,norm.ind=1:3)

et.e <- Sys.time()
print("Time for estimation: ")
print(et.e-st.e)

#make predictions
pred.num <- 1000
mh.w1 <- matrix(map.w1$par, nrow=pred.num, ncol=3, byrow=T)
mh.w2 <- matrix(map.w2$par, nrow=pred.num, ncol=3, byrow=T)
mh.w3 <- matrix(map.w3$par, nrow=pred.num, ncol=3, byrow=T)
mh.w4 <- matrix(map.w4$par, nrow=pred.num, ncol=3, byrow=T)
mh.w5 <- matrix(map.w5$par, nrow=pred.num, ncol=3, byrow=T)
mh.w6 <- matrix(map.w6$par, nrow=pred.num, ncol=3, byrow=T)
mh.w7 <- matrix(map.w7$par, nrow=pred.num, ncol=3, byrow=T)
mh.w8 <- matrix(map.w8$par, nrow=pred.num, ncol=3, byrow=T)

mh.l2.w1 <- matrix(map.l2.w1$par, nrow=pred.num, ncol=3, byrow=T)
mh.l2.w2 <- matrix(map.l2.w2$par, nrow=pred.num, ncol=3, byrow=T)
mh.l2.w3 <- matrix(map.l2.w3$par, nrow=pred.num, ncol=3, byrow=T)
mh.l2.w4 <- matrix(map.l2.w4$par, nrow=pred.num, ncol=3, byrow=T)
mh.l2.w5 <- matrix(map.l2.w5$par, nrow=pred.num, ncol=3, byrow=T)
mh.l2.w6 <- matrix(map.l2.w6$par, nrow=pred.num, ncol=3, byrow=T)
mh.l2.w7 <- matrix(map.l2.w7$par, nrow=pred.num, ncol=3, byrow=T)
mh.l2.w8 <- matrix(map.l2.w8$par, nrow=pred.num, ncol=3, byrow=T)

st.p <- Sys.time()
set.seed(2937302)
pred.wt1 <- pred.wt.pca(phi1.mat = mh.w1, phi2.mat = mh.l2.w1, w1=set.pca$Z.red[[1]][,1,drop=F], w2=set.pca$Z.red[[2]][,1, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:30,1,drop=F], X1=input.sim.env[[1]][sL1.sort,], X2=input.sim.env[[2]], X.pred=in.test, nu=2.5, nug_fix=.00001)

set.seed(2979836)
pred.wt2 <- pred.wt.pca(phi1.mat = mh.w2, phi2.mat = mh.l2.w2, w1=set.pca$Z.red[[1]][,2,drop=F], w2=set.pca$Z.red[[2]][,2, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:30,2,drop=F], X1=input.sim.env[[1]][sL1.sort,], X2=input.sim.env[[2]], X.pred=in.test, nu=2.5, nug_fix=.00001)

set.seed(7799559)
pred.wt3 <- pred.wt.pca(phi1.mat = mh.w3, phi2.mat = mh.l2.w3, w1=set.pca$Z.red[[1]][,3,drop=F], w2=set.pca$Z.red[[2]][,3, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:30,3,drop=F], X1=input.sim.env[[1]][sL1.sort,], X2=input.sim.env[[2]], X.pred=in.test, nu=2.5, nug_fix=.00001)

set.seed(2113617)
pred.wt4 <- pred.wt.pca(phi1.mat = mh.w4, phi2.mat = mh.l2.w4, w1=set.pca$Z.red[[1]][,4,drop=F], w2=set.pca$Z.red[[2]][,4, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:30,4,drop=F], X1=input.sim.env[[1]][sL1.sort,], X2=input.sim.env[[2]], X.pred=in.test, nu=2.5, nug_fix=.00001)

set.seed(2621428)
pred.wt5 <- pred.wt.pca(phi1.mat = mh.w5, phi2.mat = mh.l2.w5, w1=set.pca$Z.red[[1]][,5,drop=F], w2=set.pca$Z.red[[2]][,5, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:30,5,drop=F], X1=input.sim.env[[1]][sL1.sort,], X2=input.sim.env[[2]], X.pred=in.test, nu=2.5, nug_fix=.00001)

set.seed(5834664)
pred.wt6 <- pred.wt.pca(phi1.mat = mh.w6, phi2.mat = mh.l2.w6, w1=set.pca$Z.red[[1]][,6,drop=F], w2=set.pca$Z.red[[2]][,6, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:30,6,drop=F], X1=input.sim.env[[1]][sL1.sort,], X2=input.sim.env[[2]], X.pred=in.test, nu=2.5, nug_fix=.00001)

set.seed(9895865)
pred.wt7 <- pred.wt.pca(phi1.mat = mh.w7, phi2.mat = mh.l2.w7, w1=set.pca$Z.red[[1]][,7,drop=F], w2=set.pca$Z.red[[2]][,7, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:30,7,drop=F], X1=input.sim.env[[1]][sL1.sort,], X2=input.sim.env[[2]], X.pred=in.test, nu=2.5, nug_fix=.00001)

set.seed(2133924)
pred.wt8 <- pred.wt.pca(phi1.mat = mh.w8, phi2.mat = mh.l2.w8, w1=set.pca$Z.red[[1]][,8,drop=F], w2=set.pca$Z.red[[2]][,8, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:30,8,drop=F], X1=input.sim.env[[1]][sL1.sort,], X2=input.sim.env[[2]], X.pred=in.test, nu=2.5, nug_fix=.00001)


#predictions based on eight weights

set.seed(826292)
pred.L2.pca.recon <- pred.large.pca(wt.list=list(pred.wt1$L2, pred.wt2$L2,pred.wt3$L2, pred.wt4$L2, pred.wt5$L2, pred.wt6$L2, pred.wt7$L2, pred.wt8$L2), eig.red=set.pca$eig.red.list[[2]][,c(1:8)], scale.vec=set.pca$scale.list[[2]][1,],
                                      center.vec=set.pca$center.list[[2]][1,], Y.true=dat.pred.env[[2]], pred_int=T, sing_val_add = T, eig.val =set.pca$eig.val.list[[2]], eig.vec.dis = set.pca$eig.all.list[[2]][,9:30] )

et.p <- Sys.time()
print("Time for prediction: ")
print(et.p-st.p)


rmspe.m <- round(sqrt(mean((pred.L2.pca.recon$mean-dat.pred.env[[2]])^2)), digits=3)
cvg.m <- 100*round(length(which(as.vector(dat.pred.env[[2]])>pred.L2.pca.recon$low & as.vector(dat.pred.env[[2]])<pred.L2.pca.recon$high))/(nrow(pred.L2.pca.recon$mean)*ncol(pred.L2.pca.recon$mean)), digits=3)
alci.m <- round(mean(pred.L2.pca.recon$high - pred.L2.pca.recon$low), digits=3)

mean.sp.pca.recon <- mean.cvg.large(mat.mean = pred.L2.pca.recon$mean.spatial, mat.true=dat.pred.env[[2]])
rmspe.sp <- round(mean.sp.pca.recon$rmse, digits=3)
cvg.sp <- 100*round(mean.sp.pca.recon$cvg, digits=3)
alci.sp <- round(mean(mean.sp.pca.recon$cvg.len), digits=3)

# #compute the energy score for the model
set.seed(826292)
pred.pca.es <- pred.energy.score.pca(wt.list=list(pred.wt1$L2, pred.wt2$L2,pred.wt3$L2, pred.wt4$L2, pred.wt5$L2, pred.wt6$L2, pred.wt7$L2, pred.wt8$L2), eig.red=set.pca$eig.red.list[[2]][,c(1:8)], scale.vec=set.pca$scale.list[[2]][1,],
                                    center.vec=set.pca$center.list[[2]][1,], Y.true=dat.pred.env[[2]], sing_val_add = T, eig.val =set.pca$eig.val.list[[2]], eig.vec.dis = set.pca$eig.all.list[[2]][,9:30] )

es <- round(mean(pred.pca.es), digits=2)

res.nonsep <- matrix(c(rmspe.m, cvg.m, alci.m,es, rmspe.sp, cvg.sp, alci.sp, round((et.e-st.e)/60, digits=2), round((et.p-st.p)/60, digits=2)), nrow=1)
colnames(res.nonsep) <- c("RMSPE marginal", "CVG (95%)", "ALCI (95%)", "ES", "RMSPE spatial average", "CVG (95%) spatial average", "ALCI (95%) spatial average", "Run time estimation", "Run time prediction")

#results for NONSEP emulator in Table 1
print(res.nonsep)
