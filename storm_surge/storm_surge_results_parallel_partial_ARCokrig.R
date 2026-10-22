#storm surge results for parallel partial autoregressive cokriging emulator

rm(list=ls())
#set working directory to folder with code
#setwd("/supplementary_material_code_R1")

#libraries
library(fields)
library(farver)
library(ARCokrig)

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


# #results for PP cokriging model
cov.model = "matern_5_2"
#constants for both levels
formula = list(~1, ~1)

opt = list(maxit=800)
#different prior on the range parameters
prior = list(name="JR")

NestDesign = TRUE

out.list <- list(PSE.lf[samp.120,], PSE.hf[samp.60,])
in.list <- list(inputs.6Dim[samp.120,], inputs.6Dim[samp.60,])

obj = mvcokm(formula=formula, output=out.list, input=in.list,
             cov.model=cov.model, opt=opt, NestDesign = NestDesign, prior=prior)

st.e <- Sys.time()
#fit the model
obj.fitted = mvcokm.fit(obj)
et.e <- Sys.time()
print("Time for estimation: ")
print(et.e-st.e)




# #get predictions
st.p <- Sys.time()
pred.obj <- mvcokm.predict(obj.fitted, X.new.comp[[1]])
et.p <- Sys.time()
print("Time for prediction: ")
print(et.p-st.p)

ppco.up <- pred.obj$upper95$Level2
ppco.low <- pred.obj$lower95$Level2


nump <- function(y.training.hf, y.testing, y.pred, Lower95, Upper95, y.predSE){
  RMSPE=sqrt(mean((y.testing-y.pred)^2))
  PCI95 = mean( (Lower95<y.testing) & (Upper95>y.testing) )
  LCI95 = mean(Upper95 - Lower95)

  measure = c(RMSPE, PCI95, LCI95)
  names(measure) <- c("RMSPE", "P95CI", "L95CI")
  return(measure)
}


pp.res <-  nump(y.training.hf = out.list[[2]],
                y.testing=output.pred[[2]],
                y.pred=pred.obj$mu$Level2,
                Lower95=pred.obj$lower95$Level2,
                Upper95=pred.obj$upper95$Level2,
                y.predSE=pred.obj$SE$Level2)

#marginal predictive performance
rmspe.m <- round(pp.res[1], digits=3)
cvg.m <- 100*round(pp.res[2], digits=4)
alci.m <- round(pp.res[3], digits=3)

#get summary statistics for the mean storm surge over all locations
#follows the code from the function which uses 2*SE as the margin of error
#look at the function predict.ND to get the t quantile used for making predictions
pred.obj.mean.sd<- sqrt(apply(pred.obj$SE$Level2^2, 1, sum))/ncol(output.pred[[2]])
po.mean.upper <- apply(pred.obj$mu$Level2, 1,mean) + 2*pred.obj.mean.sd
po.mean.lower <- apply(pred.obj$mu$Level2, 1,mean) - 2*pred.obj.mean.sd
rmspe.sp <- round(sqrt(mean((apply(pred.obj$mu$Level2, 1,mean) - apply(output.pred[[2]], 1, mean))^2)), digits=3)
cvg.sp <- 100*round(length(which(apply(output.pred[[2]], 1,mean) > po.mean.lower & apply(output.pred[[2]], 1,mean) < po.mean.upper))/106, digits=4)
alci.sp <- round(mean(po.mean.upper-po.mean.lower), digits=3)

res.pparcokrig <- matrix(c(rmspe.m, cvg.m, alci.m, rmspe.sp, cvg.sp, alci.sp, round((et.e-st.e)/60, digits=2), round((et.p-st.p)/60, digits=2)), nrow=1)
colnames(res.pparcokrig) <- c("RMSPE marginal", "CVG (95%)", "ALCI (95%)", "RMSPE spatial average", "CVG (95%) spatial average", "ALCI (95%) spatial average", "Run time estimation", "Run time prediction")


#results for PP-ARCokrig model in Table 2 of the manuscript
print(res.pparcokrig)


#results for the energy score are provided at the end of the script and require additional 
#memory to compute



#compute the energy score
#code takes over two hours to run
#each of the 1,000 samples are saved for every test input and spatial location
#running this code requires around 8 GB of memory
sL1.sort <- c(sort(samp.60), sort(samp.120[which(!(samp.120 %in% samp.60))]))
sL2.sort <- sort(samp.60)


map.L1 <- c(58.923516, 29.960394, 12.417932, 52.978936, 1.271851, 32.721979)
map.L2 <- c(36.344035, 22.153126, 5.007215, 23.616352, 1.356960, 20.253545)

source("functions/functions_predict_PPArcokrig.R")
set.seed(2984020)
es.ppar <- pred.es.2L(par.L1=map.L1, par.L2=map.L2, H.list=list(matrix(1, nrow=120, ncol=9284),matrix(1, nrow=60, ncol=9284)),
                             H.star.list=list(matrix(1, nrow=106, ncol=9284), matrix(1, nrow=106, ncol=9284)),
                             Y.list=list(PSE.lf[sL1.sort,], PSE.hf[sL2.sort,]), U.list=list(NULL, PSE.lf[sL2.sort,]),
                             X.list=list(inputs.6Dim[sL1.sort,], inputs.6Dim[sL2.sort,]), X.pred=X.new.comp[[1]], times=10, Y.true=output.pred[[2]])

#results for energy score for PP-ARCokrig model
round(mean(es.ppar$es), digits=2)

