#script to obtain results for the NONSEP emulator in the storm surge application
#uses ARD form for the input correlation function
#results based on the map estimator with 1,000 samples generated for predictions

rm(list=ls())
#set working directory to folder with code
#setwd("/supplementary_material_code_R1")


#libraries
library(fields)
library(GpGp)
library(FNN)
library(xtable)

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
sL1.sort <- c(sort(samp.60), sort(samp.120[which(!(samp.120 %in% samp.60))]))
sL2.sort <- sort(samp.60)


#use pca setup function to get everthing ready to run
set.pca <- fun.pca.setup(Y.list = list(PSE.lf[sL1.sort,], PSE.hf[sL2.sort,]), var.cutoff.list=list(.999, 0.999))

#range of training inputs
range.L1.inputs <- apply(inputs.6Dim[samp.120,], 2, range)[2,]-apply(inputs.6Dim[samp.120,], 2, range)[1,]
range.L2.inputs <- apply(inputs.6Dim[samp.60,], 2, range)[2,]-apply(inputs.6Dim[samp.60,], 2, range)[1,]

#starting values for estimation of the range parameters
phi.L1 <- c(33.63, 16.07,5.48, 24.48,1.14,9.73)
phi.start <- c(phi.L1[1]/2,9.5, phi.L1[3:5]/2, 10)

#collect map estimates of model parameters
st.e <- Sys.time()
#determine starting values for each set of weights using the map estimator
map.w1 <- optim(par=phi.L1, fn=theta.post.pca.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[1]][,1], ncol=1, nrow=nrow(set.pca$Z.red[[1]])),
                X=inputs.6Dim[sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w2 <- optim(par=phi.L1, fn=theta.post.pca.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[1]][,2], ncol=1, nrow=nrow(set.pca$Z.red[[1]])),
                X=inputs.6Dim[sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w3 <- optim(par=phi.L1, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[1]][,3], ncol=1, nrow=nrow(set.pca$Z.red[[1]])),
                X=inputs.6Dim[sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2, norm.ind=c(2))

map.w4 <- optim(par=phi.L1, fn=theta.post.pca.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[1]][,4], ncol=1, nrow=nrow(set.pca$Z.red[[1]])),
                X=inputs.6Dim[sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w5 <- optim(par=phi.L1, fn=theta.post.pca.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[1]][,5], ncol=1, nrow=nrow(set.pca$Z.red[[1]])),
                X=inputs.6Dim[sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w6 <- optim(par=phi.L1, fn=theta.post.pca.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[1]][,6], ncol=1, nrow=nrow(set.pca$Z.red[[1]])),
                X=inputs.6Dim[sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w7 <- optim(par=phi.L1, fn=theta.post.pca.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[1]][,7], ncol=1, nrow=nrow(set.pca$Z.red[[1]])),W=NULL,
                X=inputs.6Dim[sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w8 <- optim(par=phi.L1, fn=theta.post.pca.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[1]][,8], ncol=1, nrow=nrow(set.pca$Z.red[[1]])),W=NULL,
                X=inputs.6Dim[sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

map.w9 <- optim(par=phi.L1, fn=theta.post.pca.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                w=matrix(set.pca$Z.red[[1]][,9], ncol=1, nrow=nrow(set.pca$Z.red[[1]])),W=NULL,
                X=inputs.6Dim[sL1.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L1.inputs/2)

# #calculate the MAP estimators for the level two weights
map.w1.L2 <- optim(par=phi.L1, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                   w=matrix(set.pca$Z.red[[2]][,1], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:60,1,drop=F],
                   X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:6)

map.w2.L2 <- optim(par=phi.L1, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                   w=matrix(set.pca$Z.red[[2]][,2], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:60,2,drop=F],
                   X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:6)

map.w3.L2 <- optim(par=phi.L1, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                   w=matrix(set.pca$Z.red[[2]][,3], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:60,3,drop=F],
                   X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2,norm.ind=1:6)

map.w4.L2 <- optim(par=phi.L1, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                   w=matrix(set.pca$Z.red[[2]][,4], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:60,4,drop=F],
                   X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2,norm.ind=1:6)

map.w5.L2 <- optim(par=phi.L1, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                   w=matrix(set.pca$Z.red[[2]][,5], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:60,5,drop=F],
                   X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:6)

map.w6.L2 <- optim(par=phi.L1, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                   w=matrix(set.pca$Z.red[[2]][,6], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:60,6,drop=F],
                   X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2,norm.ind=1:6)

map.w7.L2 <- optim(par=phi.L1, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                   w=matrix(set.pca$Z.red[[2]][,7], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:60,7,drop=F],
                   X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:6)

map.w8.L2 <- optim(par=phi.L1, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                   w=matrix(set.pca$Z.red[[2]][,8], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:60,8,drop=F],
                   X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:6)
map.w9.L2 <- optim(par=phi.L1, fn=theta.post.pca.npr.map, method="L-BFGS-B",  control=list(maxit=500),lower=rep(.000001, times=6),
                   w=matrix(set.pca$Z.red[[2]][,9], ncol=1, nrow=nrow(set.pca$Z.red[[2]])), W=set.pca$Z.red[[1]][1:60,9,drop=F],
                   X=inputs.6Dim[sL2.sort,], nu=2.5, nug_fix=.00001, phi.hp=range.L2.inputs/2, norm.ind=1:6)

et.e <- Sys.time()
print("Time for estimation: ")
print(et.e-st.e)


#make predictions using 1,000 Monte Carlo samples
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

mh.l2.w1 <- matrix(map.w1.L2$par, nrow=pred.num, ncol=6, byrow=T)
mh.l2.w2 <- matrix(map.w2.L2$par, nrow=pred.num, ncol=6, byrow=T)
mh.l2.w3 <- matrix(map.w3.L2$par, nrow=pred.num, ncol=6, byrow=T)
mh.l2.w4 <- matrix(map.w4.L2$par, nrow=pred.num, ncol=6, byrow=T)
mh.l2.w5 <- matrix(map.w5.L2$par, nrow=pred.num, ncol=6, byrow=T)
mh.l2.w6 <- matrix(map.w6.L2$par, nrow=pred.num, ncol=6, byrow=T)
mh.l2.w7 <- matrix(map.w7.L2$par, nrow=pred.num, ncol=6, byrow=T)
mh.l2.w8 <- matrix(map.w8.L2$par, nrow=pred.num, ncol=6, byrow=T)
mh.l2.w9 <- matrix(map.w9.L2$par, nrow=pred.num, ncol=6, byrow=T)


st.p <- Sys.time()
set.seed(11039302)
pred.wt1 <- pred.wt.pca(phi1.mat = mh.w1, phi2.mat = mh.l2.w1, w1=set.pca$Z.red[[1]][,1,drop=F], w2=set.pca$Z.red[[2]][,1, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:60,1,drop=F], X1=inputs.6Dim[sL1.sort,], X2=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(4752608)
pred.wt2 <- pred.wt.pca(phi1.mat = mh.w2, phi2.mat = mh.l2.w2, w1=set.pca$Z.red[[1]][,2,drop=F], w2=set.pca$Z.red[[2]][,2, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:60,2,drop=F], X1=inputs.6Dim[sL1.sort,], X2=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(5865850)
pred.wt3 <- pred.wt.pca(phi1.mat = mh.w3, phi2.mat = mh.l2.w3, w1=set.pca$Z.red[[1]][,3,drop=F], w2=set.pca$Z.red[[2]][,3, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:60,3,drop=F], X1=inputs.6Dim[sL1.sort,], X2=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(3182612)
pred.wt4 <- pred.wt.pca(phi1.mat = mh.w4, phi2.mat = mh.l2.w4, w1=set.pca$Z.red[[1]][,4,drop=F], w2=set.pca$Z.red[[2]][,4, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:60,4,drop=F], X1=inputs.6Dim[sL1.sort,], X2=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(7255663)
pred.wt5 <- pred.wt.pca(phi1.mat = mh.w5, phi2.mat = mh.l2.w5, w1=set.pca$Z.red[[1]][,5,drop=F], w2=set.pca$Z.red[[2]][,5, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:60,5,drop=F], X1=inputs.6Dim[sL1.sort,], X2=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(9804269)
pred.wt6 <- pred.wt.pca(phi1.mat = mh.w6, phi2.mat = mh.l2.w6, w1=set.pca$Z.red[[1]][,6,drop=F], w2=set.pca$Z.red[[2]][,6, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:60,6,drop=F], X1=inputs.6Dim[sL1.sort,], X2=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(1596876)
pred.wt7 <- pred.wt.pca(phi1.mat = mh.w7, phi2.mat = mh.l2.w7, w1=set.pca$Z.red[[1]][,7,drop=F], w2=set.pca$Z.red[[2]][,7, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:60,7,drop=F], X1=inputs.6Dim[sL1.sort,], X2=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)

set.seed(65108362)
pred.wt8 <- pred.wt.pca(phi1.mat = mh.w8, phi2.mat = mh.l2.w8, w1=set.pca$Z.red[[1]][,8,drop=F], w2=set.pca$Z.red[[2]][,8, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:60,8,drop=F], X1=inputs.6Dim[sL1.sort,], X2=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)


set.seed(1762093)
pred.wt9 <- pred.wt.pca(phi1.mat = mh.w9, phi2.mat = mh.l2.w9, w1=set.pca$Z.red[[1]][,9,drop=F], w2=set.pca$Z.red[[2]][,9, drop=F],
                        W.mat=set.pca$Z.red[[1]][1:60,9,drop=F], X1=inputs.6Dim[sL1.sort,], X2=inputs.6Dim[sL2.sort,], X.pred=X.new.comp[[1]], nu=2.5, nug_fix=.00001)


set.seed(3719382)
pred.L2.pca.recon <- pred.large.pca(wt.list=list(pred.wt1$L2, pred.wt2$L2,pred.wt3$L2, pred.wt4$L2, pred.wt5$L2, pred.wt6$L2, pred.wt7$L2, pred.wt8$L2, pred.wt9$L2), eig.red=set.pca$eig.red.list[[2]][,c(1:9)], scale.vec=set.pca$scale.list[[2]][1,],
                                    center.vec=set.pca$center.list[[2]][1,], Y.true=output.pred[[2]], pred_int=T, sing_val_add = T, eig.val =set.pca$eig.val.list[[2]], eig.vec.dis = set.pca$eig.all.list[[2]][,10:60] )

et.p <- Sys.time()
print("Time for prediction: ")
print(et.p-st.p)

#summarize results
rmspe.m <- round(sqrt(mean((pred.L2.pca.recon$mean-output.pred[[2]])^2)), digits=3)
cvg.m <- 100*round(length(which(as.vector(output.pred[[2]])>pred.L2.pca.recon$low & as.vector(output.pred[[2]])<pred.L2.pca.recon$high))/(nrow(pred.L2.pca.recon$mean)*ncol(pred.L2.pca.recon$mean)), digits=4)
alci.m <- round(mean(pred.L2.pca.recon$high - pred.L2.pca.recon$low), digits=3)
#average peak surge over spatial locations for nonsep emulator in Table 2 of the manuscript
mean.sp.pca.recon <- mean.cvg.large(mat.mean = pred.L2.pca.recon$mean.spatial, mat.true=output.pred[[2]])
rmspe.sp <- round(mean.sp.pca.recon$rmse, digits=3)
cvg.sp <- 100*round(mean(mean.sp.pca.recon$cvg), digits=4)
alci.sp <- round(mean(mean.sp.pca.recon$cvg.len), digits=3)

res.nonsep <- matrix(c(rmspe.m, cvg.m, alci.m, rmspe.sp, cvg.sp, alci.sp, round((et.e-st.e)/60, digits=2), round(et.p-st.p, digits=2)), nrow=1)
colnames(res.nonsep) <- c("RMSPE marginal", "CVG (95%)", "ALCI (95%)", "RMSPE spatial average", "CVG (95%) spatial average", "ALCI (95%) spatial average", "Run time estimation", "Run time prediction")

#code to recreate results in Table 3 of the manuscript
#matrices of center and scale from level 2 training outputs
mat.L2.scale <- matrix(set.pca$scale.list[[2]][1,],nrow=nrow(X.new.comp[[1]]), ncol=ncol(output.pred[[1]]), byrow=T)
mat.L2.center <- matrix(set.pca$center.list[[2]][1,],nrow=nrow(X.new.comp[[1]]), ncol=ncol(output.pred[[1]]), byrow=T)

#calculate the average predictive weight for each PC
wt1.mean.L2 <- apply(pred.wt1$L2,1, mean)
wt2.mean.L2 <- apply(pred.wt2$L2,1, mean)
wt3.mean.L2 <- apply(pred.wt3$L2,1, mean)
wt4.mean.L2 <- apply(pred.wt4$L2,1, mean)
wt5.mean.L2 <- apply(pred.wt5$L2,1, mean)
wt6.mean.L2 <- apply(pred.wt6$L2,1, mean)
wt7.mean.L2 <- apply(pred.wt7$L2,1, mean)
wt8.mean.L2 <- apply(pred.wt8$L2,1, mean)
wt9.mean.L2 <- apply(pred.wt9$L2,1, mean)

w.L2.pred.mat <- cbind(wt1.mean.L2, wt2.mean.L2, wt3.mean.L2, wt4.mean.L2, wt5.mean.L2, wt6.mean.L2, wt7.mean.L2, wt8.mean.L2, wt9.mean.L2)

#calculate the rmspe by number of components used

rmspe.pred <- rmspe.by.pcomp(eig.red=set.pca$eig.all.list[[2]][,1:9], wt.mat=w.L2.pred.mat, center.vec=set.pca$center.list[[2]][1,], scale.vec=set.pca$scale.list[[2]][1,], y.true=output.pred[[2]])

# #PCA with both the training and test data
pca.L2.all <- prcomp(rbind(PSE.hf[sL2.sort,], output.pred[[2]]), center=T, scale=T)
#scores for the predictions from optimal basis representation
scores.pca.all.L2 <- pca.L2.all$x[61:nrow(pca.L2.all$x),1:9]

#compute rmspe when optimal basis representation is used
rmspe.best.bf <- rmspe.by.pcomp(eig.red=pca.L2.all$rotation[,1:9], wt.mat=scores.pca.all.L2, center.vec=pca.L2.all$center, scale.vec=pca.L2.all$scale, y.true=output.pred[[2]])

#what are the pca scores based on the original basis representation
scores.pred.L2 <- ((output.pred[[2]]-mat.L2.center)/mat.L2.scale)%*%set.pca$eig.red.list[[2]][,1:9]

#rmspe based on projecting basis representation of predictions back to the original space
rmspe.proj <- rmspe.by.pcomp(eig.red=set.pca$eig.all.list[[2]][,1:9], wt.mat=scores.pred.L2, center.vec=mat.L2.center[1,], scale.vec=mat.L2.scale[1,], y.true=output.pred[[2]])

#pca for level 1 and 2 to get explained variance
pca.L1.sum <- summary(prcomp(PSE.lf[sL1.sort,], center=T, scale=T))
pca.L2.sum <- summary(prcomp(PSE.hf[sL2.sort,], center=T, scale=T))

sum.rmspe <- data.frame(pc_num=1:9,var_L1 = pca.L1.sum$importance[3,1:9]*100, var_L2=pca.L2.sum$importance[3,1:9]*100,
                        MRE=rmspe.best.bf, TRE=rmspe.proj, rmspe_cokriging=rmspe.pred)
sum.mse <- sum.rmspe
sum.mse[,4:6] <- 100*sum.rmspe[,4:6]^2

#load the MSPE from NONSEP kriging emulator from the data folder
#the results for the MSPE of the NONSEP kriging emulator can be checked
#using the script storm_surge/storm_surge_results_nonsep_kriging.R
load(file="data/error.prcomp.1L.rda")
sum.mse$rmspe_kriging <- 100*error.prcomp.1L^2




#final results for the NONSEP emulator
#results for NONSEP emulator in Table 2 of the manuscript
print(res.nonsep)

#results for the energy score are provided at the end of the script and require additional 
#memory to compute

#results for Table 3 of the manuscript
print(round(sum.mse, digits=2))



#code to compute the energy score
#each of the 1,000 samples are saved for every test input and spatial location
#running this code requires around 8 GB of memory
set.seed(3719382)
pred.L2.pca.es <- pred.energy.score.pca(wt.list=list(pred.wt1$L2, pred.wt2$L2,pred.wt3$L2, pred.wt4$L2, pred.wt5$L2, pred.wt6$L2, pred.wt7$L2, pred.wt8$L2, pred.wt9$L2), eig.red=set.pca$eig.red.list[[2]][,c(1:9)], scale.vec=set.pca$scale.list[[2]][1,],
                                        center.vec=set.pca$center.list[[2]][1,], Y.true=output.pred[[2]], sing_val_add = T, eig.val =set.pca$eig.val.list[[2]], eig.vec.dis = set.pca$eig.all.list[[2]][,10:60])

#results for energy score for NONSEP emulator in Table 2 of the manuscript
round(mean(pred.L2.pca.es), digits=2)

