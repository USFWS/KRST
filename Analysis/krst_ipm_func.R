krst.scenario <- function(Time, Time.fore,
                          krst_con, krst_data,dat,
                          hist.bioc.means,
                          clim.var, ssp,slrp.summary,
                          sl.summary,slr.scenario,
                          nest.mange.fore,krst.ipm,
                          params){
  

prob.fem = krst_data$prob.fem
prob.fem[(Time+Time.fore)] = NA

n.fem = krst_data$n.fem
n.fem[(Time+Time.fore)] <- NA

bio5 = subset(hist.bioc.means,bioc == clim.var)
bio5.out = rep(NA,length(krst_data$J.ind))
bio5.fore.245 = subset(bioc.forecast,bioc == clim.var & ssp == ssp)
bio5.tmp = data.frame(mean = c(rep(bio5.fore.245$mean[1],length.out=19),
                               rep(bio5.fore.245$mean[2],length.out=20),
                               rep(bio5.fore.245$mean[3],length.out=20),
                               rep(bio5.fore.245$mean[4],length.out=20)),
                      year = seq(2022,2100))
bio5.scale = as.numeric(scale(c(bio5$mean,bio5.tmp$mean)))
bio5$mean.scale = bio5.scale[1:31]

slr.cov = rep(NA,length(krst_data$J.ind))
slr.fore.med = subset(slrp.summary, Scenario == slr.scenario)
slr.tmp = data.frame(mean = c(rep(slr.fore.med$msl.mean[1],length.out=19),
                              rep(slr.fore.med$msl.mean[2],length.out=20),
                              rep(slr.fore.med$msl.mean[3],length.out=20),
                              rep(slr.fore.med$msl.mean[4],length.out=20)),
                     year = seq(2022,2100))
slr.scale = as.numeric(scale(c(sl.summary$msl.mean[2:8],slr.tmp$mean)))
sl.summary$mean.scale = c(NA,slr.scale[1:7])

#may.temp = subset(hist.temp.means,month == 5)
#may.tmax = rep(NA,length(nest.pro$Year)-1)
for(i in 1:1597){
  #may.tmax[i] = may.temp$tmax[which(may.temp$year==(year+1990)[i])]
  bio5.out[i] = bio5$mean.scale[which(bio5$year==(year+1990)[i])]
  slr.cov[i] = sl.summary$msl.mean[which(sl.summary==(year+1990)[i])]
}
#may.tmax.ann = as.numeric(scale(may.temp$tmax))
bio5.ann = bio5$mean.scale
slr.cov = as.numeric(scale(slr.cov))

bio5.fore.cov = bio5.tmp$mean[1:Time.fore]
slr.cov.fore = slr.tmp$mean[1:Time.fore]

#forecast count data
y.fem = krst_data$y.fem
y.fem[32:(Time+Time.fore)] = NA

#need to truncate N.b + N.i to be at least as 
# large as number of observed individuals. Just
# have to use the y.fem data twice, once as data
# once as truncation constraint

trun = c(apply(dat,2,sum))

#### Data for NIMBLE ####

krst_data <- list(y = krst_data$y, 
                  y.fem = y.fem,
                  psiNBB.mom=krst_data$psiNBB.mom,
                  psiBNB.mom=krst_data$psiBNB.mom,
                  p.mom=krst_data$p.mom,
                  s.p.mom=krst_data$s.p.mom,
                  in.s = krst_data$in.s,
                  J = krst_data$J, 
                  E.R = krst_data$E.R,
                  E = krst_data$E,
                  N.p1.ini = 10,
                  N.p2.ini = 3,
                  N.p3.ini = 3,
                  N.p4.ini = 3,
                  N.p5.ini = 3,
                  N.p6.ini = 3,
                  N.p7.ini = 3,
                  N.p8.ini = 3,
                  N.p9.ini = 3,
                  N.p10.ini = 3,
                  N.b.ini=1, 
                  N.nb.ini=2,
                  prob.fem = prob.fem,
                  n.fem = n.fem,
                  J.ind = krst_data$J.ind,
                  E.ind = krst_data$E.ind,
                  E.ind.fore = 96)

krst_con2 <- list(N = nrow(y), Time = Time, 
                 first=krst_con$first,
                 Time.fore = Time.fore,
                 nest.mange = krst_con$nest.mange,
                 year = krst_con$year,
                 year.mean.idx = krst_con$year.mean.idx,
                 N.nest = length(J.ind),
                 may.tmax.ann = bio5.ann,
                 clim.cov = bio5.out,
                 clim.cov.fore = bio5.fore.cov,
                 nest.mange.fore = nest.mange.fore,
                 slr.cov=slr.cov,
                 slr.cov.fore = slr.cov.fore)

zinit <- dat
y <- krst_data$y
first <- krst_con$first
for (i in 1:nrow(y)) {
  for (j in 1:ncol(y)) {
    if (j > first[i] & dat[i,j]==0) {zinit[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1)}
    if (j < first[i]) {zinit[i,j] <- 0}
  }
}

J.ini = J.ind
J.ini[which(is.na(J.ini))] = runif(length(which(is.na(J.ini))),0,80)
J.ini[which(!is.na(J.ind))] = NA

prob.fem.ini = prob.fem
prob.fem.ini[which(is.na(prob.fem.ini))] = runif(length(which(is.na(prob.fem.ini))),0.1,0.9)
prob.fem.ini[which(!is.na(prob.fem))] = NA

#krst_ini = function() list(
#  psiNBB=runif(1,0,.5),
#  psiBNB=runif(1,0.5,1),
#  pB=runif(1,0,1), 
#  z=zinit,
#  s.p1 = runif(1,.1,.4), 
#  N.nb = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
  #N.b = rpois((Time+Time.fore),c(trun[1:Time]/2,rep(50,Time.fore))),
#  N.b = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
#  N.p1 = rpois((Time+Time.fore),seq(10,200,length.out=(Time+Time.fore))),
#  N.p2 = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
#  N.p3 = rpois((Time+Time.fore),seq(2,500,length.out=(Time+Time.fore))),
#  N.p4 = rpois((Time+Time.fore),seq(2,500,length.out=(Time+Time.fore))),
#  N.p5 = rpois((Time+Time.fore),seq(2,500,length.out=(Time+Time.fore))),
#  N.p6 = rpois((Time+Time.fore),seq(2,500,length.out=(Time+Time.fore))),
#  N.p7 = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
#  N.p8 = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
#  N.p9 = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
#  N.p10 = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
#  #N.i = rpois((Time+Time.fore),c(trun[1:Time]/2,rep(50,Time.fore))),
#  N.i = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
#  r = runif(1,.8,1.2),
#  ln.mean.fem = runif(1,-1,1),
#  mu.n.fem=2,
#  sd.n.fem = runif(1,.1,1),
#  mu.phi = runif(1,1,2.2),
#  sd.phi2=runif(1,.1,1),
#  mu.egg = runif(1,-1,1),
#  nest.cov = c(runif(5,-1,1),NA),
#  insitu.surv = runif(1,.1,.8),
#  mu.temp = runif(1,-1,1),
#  var.temp = runif(1,.1,3),
#  beta.temp.all = runif(1,-1,1),
#  beta.temp = runif(6,-1,1),
#  mu.imm = runif(1,-1,1),
#  sd.imm = runif(1,.01,1),
#  ln.mu.phi = runif((Time+Time.fore-1),1,2.2),
#  eps.imm.tmp = runif((Time+Time.fore),-1,1),
#  J.ind = round(J.ini),
#  prob.fem = prob.fem.ini,
#  n.fem = c(rep(NA,Time),rep(2,Time.fore)),
#  y.fem = c(rep(NA,Time),rep(100,Time.fore)),
#  beta.slr = c(rep(NA,5),runif(1,-1,1)),
#  f = rep(76,(Time+Time.fore)))

#krst.model <- nimbleModel(code = krst.ipm,
#                          constants = krst_con,
#                          data = krst_data,
#                          inits = krst_ini())

#krst.config.mcmc <- configureMCMC(krst.model)

#krst.config.mcmc$addMonitors("f","prob.egg.fore")

#krst.mcmc <- buildMCMC(krst.config.mcmc)

#krst.cmodel <- compileNimble(krst.model)

#krst.cmcmc <- compileNimble(krst.mcmc, project = krst.model)

#krst.samples <- runMCMC(krst.cmcmc, niter= 110000, 
#                        nburnin = 100000,nchains = 3,
#                        samplesAsCodaMCMC = TRUE)


nc = 3 #num of chains
clus = makeCluster(nc)

#nsamp2 = nsamp2
#nburnin1 = nburnin1

clusterExport(clus, c("krst.ipm","krst_ini",
                      "krst_data",
                      "krst_con2","params"),
              envir = environment())
                      

for(j in seq_along(clus)){
  set.seed(j)
  init <- krst_ini()
  clusterExport(clus[j],"init")
}

out <- clusterEvalQ(clus, {
  library(nimble)
  library(coda)
  model <- nimbleModel(code = krst.ipm, 
                       name = "krst.ipm",
                       constants = krst_con2,
                       data = krst_data, 
                       inits = init)
  Cmodel <- compileNimble(model)
  modelConf <- configureMCMC(model)
  modelConf$addMonitors(params)
  modelMCMC <- buildMCMC(modelConf)
  CmodelMCMC <- compileNimble(modelMCMC, 
                              project = model)
  CmodelMCMC$run(2500,reset=TRUE, resetMV=TRUE)
  return(as.mcmc(as.matrix(CmodelMCMC$mvSamples)))
})

#out.mcmc <- as.mcmc.list(out) 

out2 <- clusterEvalQ(clus,{
  CmodelMCMC$run(5000, reset=FALSE, resetMV=TRUE)
  return(as.mcmc(as.matrix(CmodelMCMC$mvSamples)))
})

stopCluster(clus)
out.up <- as.mcmc.list(out2)

 return(out.up = out.up)

}
