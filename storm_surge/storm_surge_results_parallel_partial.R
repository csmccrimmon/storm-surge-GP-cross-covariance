#code to calculate results for PP-ARCokrig and PP emulators in Table 2 of the manuscript

rm(list=ls())
#set working directory to folder with code
#setwd("/supplementary_material_code_R1")

library(fields)
library(farver)
library(RobustGaSP)

#load the data
load("data/StormSurge.RData")

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

#results for the PP emulator trained using only the ADCIRC + SWAN outputs
input.pp = inputs.6Dim[samp.60,]
p.x = dim(input.pp)[2]
formula.pp = ~1
colnames(input.pp) = paste0("x", 1:p.x)
trend = model.matrix(formula.pp, data.frame(input.pp))
input.new <- X.new.comp[[2]]
colnames(input.new) = paste0("x", 1:p.x)
trend.testing = model.matrix(formula.pp, data.frame(input.new))

### fit pp model
st.e <- Sys.time()
fit.PP = ppgasp(design=input.pp, response=PSE.hf[samp.60,], trend=trend, kernel="matern_5_2", isotropic = F, optimization = "nelder-mead")
et.e <- Sys.time()
print("Time for estimation: ")
print(et.e-st.e)

#make predictions
st.p <- Sys.time()
pred.pp = predict(fit.PP, X.new.comp[[2]], testing_trend=trend.testing)
et.p <- Sys.time()
print(et.p-st.p)

nump <- function( y.testing, y.pred, Lower95, Upper95, y.predSE){
  RMSPE=sqrt(mean((y.testing-y.pred)^2))
  PCI95 = mean( (Lower95<y.testing) & (Upper95>y.testing) )
  LCI95 = mean(Upper95 - Lower95)
  
  measure = c(RMSPE, PCI95, LCI95)
  names(measure) <- c("RMSPE", "P95CI", "L95CI")
  return(measure)
}

#marginal predictive performance of PP emulator
pp.res = nump(y.testing=output.pred[[2]], 
              y.pred=pred.pp$mean,
              Lower95=pred.pp$lower95, 
              Upper95=pred.pp$upper95,
              y.predSE=pred.pp$sd)

rmspe.m <- round(pp.res[1], digits=3)
cvg.m <- 100*round(pp.res[2], digits=4)
alci.m <- round(pp.res[3], digits=3)


#results for average output over spatial locations
pp.mean.sd<- sqrt(apply(pred.pp$sd^2, 1, sum))/ncol(output.pred[[2]])
pp.mean.upper <- apply(pred.pp$mean, 1,mean) + 2*pp.mean.sd
pp.mean.lower <- apply(pred.pp$mean, 1,mean) - 2*pp.mean.sd
rmse.pp.mean <- sqrt(mean((apply(pred.pp$mean, 1,mean) - apply(output.pred[[2]], 1, mean))^2))

#predictive performance at average output over spatial locations
#rmspe
rmspe.sp <- round(rmse.pp.mean, digits=3)
cvg.sp <- 100*round(length(which(apply(output.pred[[2]], 1,mean) > pp.mean.lower & apply(output.pred[[2]], 1,mean) < pp.mean.upper))/nrow(output.pred[[2]]), digits=4)
alci.sp <- round(mean(4*pp.mean.sd), digits=3)


res.pp <- matrix(c(rmspe.m, cvg.m, alci.m, rmspe.sp, cvg.sp, alci.sp, round((et.e-st.e)/60, digits=2), round((et.p-st.p)/60, digits=2)), nrow=1)
colnames(res.pp) <- c("RMSPE marginal", "CVG (95%)", "ALCI (95%)", "RMSPE spatial average", "CVG (95%) spatial average", "ALCI (95%) spatial average", "Run time estimation", "Run time prediction")

#results for PP emulator in Table 2 of the manuscript
print(res.pp)


#results for the energy score are provided at the end of the script and require additional 
#memory to compute


#compute the energy score
#each of the 1,000 samples are saved for every test input and spatial location
#running this code requires around 8 GB of memory
source("functions/functions_predict_PPArcokrig.R")

pp.map.L1 <- c(66.12665, 35.94269, 19.10454, 63.73145, 3.524803, 39.50297)

set.seed(2048402)
pp.es <-  pred.es.PP(par.L1 = pp.map.L1, H.mat=matrix(1, nrow=60, ncol=9284), H.star.mat = matrix(1, nrow=106, ncol=9284), 
                     Y=PSE.hf[samp.60,], X=inputs.6Dim[samp.60,], X.pred=X.new.comp[[1]],nu=2.5, times=1000, Y.true=output.pred[[2]])

#energy score for PP emulator in Table 2 of the manuscript
round(mean(pp.es$es), digits=2)




