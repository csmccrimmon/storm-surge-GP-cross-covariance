#code to obtain results for PPArcokrig model in the simulation

rm(list=ls())
#set working directory to folder with code
#setwd("/supplementary_material_code_R1")

library(fields)
library(farver)
library(ARCokrig)

#load the data 
load("data/in.out.1k.rda")
load("data/in.out.pred.1k.rda")


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


#results for parallel partial autoregressive cokriging model
cov.model = "matern_5_2"
#constants for both levels
formula = list(~1, ~1)

opt = list(maxit=800)
#different prior on the range parameters
prior = list(name="JR")

NestDesign = TRUE

obj = mvcokm(formula=formula, output=output.sim.env, input=input.sim.env,
             cov.model="matern_5_2", opt=opt, NestDesign = NestDesign, prior=prior)

st.e <- Sys.time()
#fit the model
obj.fitted = mvcokm.fit(obj)
et.e <- Sys.time()
print("Time for estimation: ")
print(et.e-st.e)




#get predictions
st.p <- Sys.time()
pred.obj <- mvcokm.predict(obj.fitted, in.test)
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

# #marginal predictive performance of PP-ARCokrig emulator
pp.ar.res <-  nump(y.testing=dat.pred.env$hf,
                   y.pred=pred.obj$mu$Level2,
                   Lower95=pred.obj$lower95$Level2,
                   Upper95=pred.obj$upper95$Level2,
                   y.predSE=pred.obj$SE$Level2)

rmspe.m <- round(pp.ar.res[1], digits=3)
cvg.m <- 100*round(pp.ar.res[2], digits=3)
alci.m <- round(pp.ar.res[3], digits=3)

#get summary statistics for predictive performance at average output over the spatial locations
pred.obj.mean.sd<- sqrt(apply(pred.obj$SE$Level2^2, 1, sum))/ncol(dat.pred.env$hf)
po.mean.upper <- apply(pred.obj$mu$Level2, 1,mean) + 2*pred.obj.mean.sd
po.mean.lower <- apply(pred.obj$mu$Level2, 1,mean) - 2*pred.obj.mean.sd
rmse.po.mean <- sqrt(mean((apply(pred.obj$mu$Level2, 1,mean) - apply(dat.pred.env$hf, 1, mean))^2))

#predictive performance at average output over spatial locations
rmspe.sp <- round(rmse.po.mean, digits=3)
cvg.sp <- 100*round(length(which(apply(dat.pred.env$hf, 1,mean) > po.mean.lower & apply(dat.pred.env$hf, 1,mean) < po.mean.upper))/nrow(dat.pred.env$hf),digits=3)
alci.sp <- round(mean(po.mean.upper-po.mean.lower), digits=3)

#compute the energy score for the model
map.L1 <- c(37.7729944, 0.1322507, 3.9683379)
map.L2 <- c(19.93946773, 0.05525765, 2.14161426)
source("functions/functions_predict_PPArcokrig.R")
set.seed(738201)
es.ppar <- pred.es.2L(par.L1=map.L1, par.L2=map.L2, H.list=list(matrix(1, nrow=60, ncol=1000),matrix(1, nrow=30, ncol=1000)),
                             H.star.list=list(matrix(1, nrow=30, ncol=1000), matrix(1, nrow=30, ncol=1000)),
                             Y.list=list(output.sim.env[[1]], output.sim.env[[2]]), U.list=list(NULL, output.sim.env.W),
                             X.list=list(input.sim.env[[1]], input.sim.env[[2]]), X.pred=in.test, times=1000, Y.true=dat.pred.env$hf)

es <- round(mean(es.ppar$es), digits=2)

res.pparcokrig <- matrix(c(rmspe.m, cvg.m, alci.m,es, rmspe.sp, cvg.sp, alci.sp, round((et.e-st.e)/60, digits=2), round((et.p-st.p)/60, digits=2)), nrow=1)
colnames(res.pparcokrig) <- c("RMSPE marginal", "CVG (95%)", "ALCI (95%)", "ES", "RMSPE spatial average", "CVG (95%) spatial average", "ALCI (95%) spatial average", "Run time estimation", "Run time prediction")

#results for PPARCokrig emulator in Table 1
print(res.pparcokrig)






