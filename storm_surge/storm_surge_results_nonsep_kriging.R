#code to calculate the MSPE for the NONSEP kriging model in Table 3 of the manuscript
#models estimated using ARD form for the input correlation function
#with map estimator for correlation parameters

rm(list=ls())
#set working directory to folder with code
#setwd("/supplementary_material_code_R1")


#libraries
library(fields)

#load the data
load("data/StormSurge.RData")

#load helper functions
source("functions/functions_nonsep_paper_ARD.R")

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

#sort the inputs so first 60 from level 1 are the same as for level 2
sL1.sort <- c(sort(samp.60), sort(which(!(samp.120 %in% samp.60))))
sL2.sort <- sort(samp.60)

#use pca setup function to get everthing ready to run
set.pca <- fun.pca.setup(Y.list = list(PSE.lf[sL1.sort,], PSE.hf[sL2.sort,]), var.cutoff.list=list(.997, 0.999))

#range of training inputs
range.L1.inputs <- apply(inputs.6Dim[samp.120,], 2, range)[2,]-apply(inputs.6Dim[samp.120,], 2, range)[1,]
range.L2.inputs <- apply(inputs.6Dim[samp.60,], 2, range)[2,]-apply(inputs.6Dim[samp.60,], 2, range)[1,]

#starting values for estimation of the range parameters
phi.L1 <- c(33.63, 16.07,5.48, 24.48,1.14,9.73)
phi.start <- c(phi.L1[1]/2,9.5, phi.L1[3:5]/2, 10)

#calculate the map estimators for each weight to get starting values to run the sampling algorithm

map.w1 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[2]][,1], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=NULL,
                X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:6)

map.w2 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[2]][,2], ncol=1, nrow=nrow(set.pca$Z.red[[2]])),W=NULL,
                X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:6)

map.w3 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B", control=list(maxit=500),lower=rep(.000001, times=6), 
                w=matrix(set.pca$Z.red[[2]][,3], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=NULL,
                X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind = 1:6)

map.w4 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[2]][,4], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=NULL,
                X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2,norm.ind=1:6)

map.w5 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6), 
                w=matrix(set.pca$Z.red[[2]][,5], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=NULL,
                X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:6)

map.w6 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[2]][,6], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=NULL,
                X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2,norm.ind=1:6)

map.w7 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[2]][,7], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=NULL,
                X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:6)

map.w8 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[2]][,8], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=NULL,
                X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:6)

map.w9 <- optim(par=phi.start, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[2]][,9], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=NULL,
                X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:6)


#use 1,000 samples to generate predictions
pred.num <- 1000

mh.w1 <- matrix(map.w1$par, nrow=pred.num, ncol=6, byrow=T)
mh.w2 <- matrix(map.w2$par, nrow=pred.num, ncol=6, byrow=T)
mh.w3 <- matrix(map.w3$par, nrow=pred.num, ncol=6, byrow=T)
mh.w4 <- matrix(map.w4$par, nrow=pred.num, ncol=6, byrow=T)
mh.w5 <- matrix(map.w5$par, nrow=pred.num, ncol=6, byrow=T)
mh.w6 <- matrix(map.w6$par, nrow=pred.num, ncol=6, byrow=T)
mh.w7 <- matrix(map.w7$par, nrow=pred.num, ncol=6, byrow=T)
mh.w8 <- matrix(map.w8$par, nrow=pred.num, ncol=6, byrow=T)
mh.w9 <- matrix(map.w9$par, nrow=pred.num, ncol=6, byrow=T)


#compute predictions
set.seed(6267625)
pred.wt1 <- pred.wt.L1.pca(phi1.mat = mh.w1, w1=set.pca$Z.red[[2]][,1, drop=F],
                           X1=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(5842375)
pred.wt2 <- pred.wt.L1.pca(phi1.mat = mh.w2, w1=set.pca$Z.red[[2]][,2, drop=F],
                           X1=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(2533430)
pred.wt3 <- pred.wt.L1.pca(phi1.mat = mh.w3, w1=set.pca$Z.red[[2]][,3, drop=F],
                           X1=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(2516322)
pred.wt4 <- pred.wt.L1.pca(phi1.mat = mh.w4, w1=set.pca$Z.red[[2]][,4, drop=F],
                           X1=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(534194)
pred.wt5 <- pred.wt.L1.pca(phi1.mat = mh.w5, w1=set.pca$Z.red[[2]][,5, drop=F],
                           X1=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(9200475)
pred.wt6 <- pred.wt.L1.pca(phi1.mat = mh.w6, w1=set.pca$Z.red[[2]][,6, drop=F],
                           X1=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(6977180)
pred.wt7 <- pred.wt.L1.pca(phi1.mat = mh.w7, w1=set.pca$Z.red[[2]][,7, drop=F],
                           X1=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(849302)
pred.wt8 <- pred.wt.L1.pca(phi1.mat = mh.w8, w1=set.pca$Z.red[[2]][,8, drop=F],
                           X1=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(637201)
pred.wt9 <- pred.wt.L1.pca(phi1.mat = mh.w9, w1=set.pca$Z.red[[2]][,9, drop=F],
                           X1=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

#calculate the average predicted weight
wt1.mean.L2 <- apply(pred.wt1,1, mean)
wt2.mean.L2 <- apply(pred.wt2,1, mean)
wt3.mean.L2 <- apply(pred.wt3,1, mean)
wt4.mean.L2 <- apply(pred.wt4,1, mean)
wt5.mean.L2 <- apply(pred.wt5,1, mean)
wt6.mean.L2 <- apply(pred.wt6,1, mean)
wt7.mean.L2 <- apply(pred.wt7,1, mean)
wt8.mean.L2 <- apply(pred.wt8, 1, mean)
wt9.mean.L2 <- apply(pred.wt9, 1, mean)

w.L2.pred.mat <- cbind(wt1.mean.L2, wt2.mean.L2, wt3.mean.L2, wt4.mean.L2, wt5.mean.L2, wt6.mean.L2, wt7.mean.L2, wt8.mean.L2, wt9.mean.L2)

#calculate the rmspe by number of components used
error.prcomp.1L <- rmspe.by.pcomp(eig.red=set.pca$eig.all.list[[2]][,1:9], wt.mat=w.L2.pred.mat, center.vec=set.pca$center.list[[2]][1,], scale.vec=set.pca$scale.list[[2]][1,], y.true=output.pred[[2]])

#results for MSPE for the NONSEP kriging model in Table 3 of the manuscript
round(100*error.prcomp.1L^2, digits=2)



