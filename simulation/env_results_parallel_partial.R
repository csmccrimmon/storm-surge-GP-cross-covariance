#code to obtain results for parallel partial model in the simulation

rm(list=ls())
#set working directory to folder with code
#setwd("/supplementary_material_code_R1")

library(RobustGaSP)

#load the data 
load("data/in.out.1k.rda")
load("data/in.out.pred.1k.rda")

#60 inputs at level 1 and 30 inputs at level 2 using conditional Latin hypercube sampling
input.sim.env <- in.out.1k$input

#list of outputs for the two levels
output.sim.env <- in.out.1k$output

#inputs for predictions
in.test <- in.out.pred.1k$input

#outputs for predictions
dat.pred.env <- in.out.pred.1k$output

#PP emulator for only the high fidelity outputs
input.pp = input.sim.env[[2]]
p.x = dim(input.pp)[2]
formula.pp = ~1
colnames(input.pp) = paste0("x", 1:p.x)
trend = model.matrix(formula.pp, data.frame(input.pp))

input.new <- in.test
colnames(input.new) = paste0("x", 1:p.x)
trend.testing = model.matrix(formula.pp, data.frame(input.new))

#estimate model parameters
st.e <- Sys.time()
fit.PP = ppgasp(design=input.pp, response=output.sim.env[[2]], trend=trend, optimization = "nelder-mead")
et.e <- Sys.time()
print("Time for estimation: ")
print(et.e-st.e)

#make predictions
st.p <- Sys.time()
pred.pp = predict(fit.PP, input.new, testing_trend=trend.testing)
et.p <- Sys.time()
print("Time for prediction: ")
print(et.p-st.p)


#function to summarize marginal predictive performance
nump <- function(y.testing, y.pred, Lower95, Upper95, y.predSE){
  RMSPE=sqrt(mean((y.testing-y.pred)^2))
  PCI95 = mean( (Lower95<y.testing) & (Upper95>y.testing) )
  LCI95 = mean(Upper95 - Lower95)
  measure = c(RMSPE, PCI95, LCI95)
  names(measure) <- c("RMSPE", "P95CI", "L95CI")
  return(measure)
}

#marginal predictive performance of PP emulator 
pp.res = nump(y.testing=dat.pred.env$hf, 
              y.pred=pred.pp$mean,
              Lower95=pred.pp$lower95, 
              Upper95=pred.pp$upper95,
              y.predSE=pred.pp$sd)

rmspe.m <- round(pp.res[1], digits=3)
cvg.m <- 100*round(pp.res[2], digits=3)
alci.m <- round(pp.res[3], digits=3)

#predictive performance at average output over spatial locations 
pp.mean.sd<- sqrt(apply(pred.pp$sd^2, 1, sum))/ncol(dat.pred.env$hf)
pp.mean.upper <- apply(pred.pp$mean, 1,mean) + 2*pp.mean.sd
pp.mean.lower <- apply(pred.pp$mean, 1,mean) - 2*pp.mean.sd
rmse.pp.mean <- sqrt(mean((apply(pred.pp$mean, 1,mean) - apply(dat.pred.env$hf, 1, mean))^2))

#predictive performance at average output over spatial locations
rmspe.sp <- round(rmse.pp.mean, digits=3)
cvg.sp <- 100*round(length(which(apply(dat.pred.env$hf, 1,mean) > pp.mean.lower & apply(dat.pred.env$hf, 1,mean) < pp.mean.upper))/nrow(dat.pred.env$hf), digits=3)
alci.sp <- round(mean(pp.mean.upper-pp.mean.lower), digits=3)


#compute the energy score the PP emulator
source("functions/functions_predict_PPArcokrig.R")
pp.map.L1 <- c(17.03096, 0.1067997, 2.839486)
set.seed(748303)
pp.es <-  pred.es.PP(par.L1 = pp.map.L1, H.mat=matrix(1, nrow=30, ncol=1000), H.star.mat = matrix(1, nrow=30, ncol=1000), 
                     Y=output.sim.env[[2]], X=input.sim.env[[2]], X.pred=in.test, nu=2.5, times=1000, Y.true=dat.pred.env$hf)

es <- round(mean(pp.es$es), digits=2)

res.pp <- matrix(c(rmspe.m, cvg.m, alci.m,es, rmspe.sp, cvg.sp, alci.sp, round((et.e-st.e)/60, digits=2), round((et.p-st.p)/60, digits=2)), nrow=1)
colnames(res.pp) <- c("RMSPE marginal", "CVG (95%)", "ALCI (95%)", "ES", "RMSPE spatial average", "CVG (95%) spatial average", "ALCI (95%) spatial average", "Run time estimation", "Run time prediction")

#results for PP emulator in Table 1
print(res.pp)
