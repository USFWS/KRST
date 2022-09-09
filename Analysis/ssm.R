ssm.ipm <- nimbleCode({
  
  s.p1 ~ dbeta(s.p.mom[1],s.p.mom[2])
  #s.p1 ~ dbeta(1,1)
  
  #priors initial abundance
  N.p1[1] ~ dpois(N.p1.ini)
  N.p2[1] ~ dpois(N.p2.ini)
  N.p3[1] ~ dpois(N.p3.ini)
  N.p4[1] ~ dpois(N.p4.ini)
  N.p5[1] ~ dpois(N.p5.ini)
  N.p6[1] ~ dpois(N.p6.ini)
  N.p7[1] ~ dpois(N.p7.ini)
  N.p8[1] ~ dpois(N.p8.ini)
  N.p9[1] ~ dpois(N.p9.ini)
  N.p10[1] ~ dpois(N.p10.ini)
  #N.b[1] ~ dpois(N.b.ini)
  N.b[1] <- 1
  #N.nb[1] ~ dpois(N.nb.ini)
  N.nb[1] <- 1
  
  pre.br1[1] <- N.p1[1]
  no.breed[1] <- N.nb[1]
  breed[1] <- N.b[1]
  imm[1] <- 1
  N.i[1] <- 1
  N.tot[1] <- 2
  mu.imm ~ dnorm(0, sd = 2)
  #sd.imm ~ dinvgamma(2,1)
  
  for(t in 1:(Time+Time.fore)){
    #eps.imm.tmp[t] ~ dnorm(mu.imm, var=sd.imm)
    eps.imm.tmp[t] ~ dnorm(mu.imm, sd = 1)
    #eps.imm.tmp[t] <- 2
    #eps.imm[t] <- min(eps.imm.tmp[t],10)
    eps.imm[t] <- eps.imm.tmp[t]
  }
  
  mean.fem <- .7
  mu.n.fem <- 2
  f[1:73] <- 85
  phi[1:72] <- .9
  psiNBB <- .3
  psiBNB <- .9
  #prob.trans ~ dbeta(1,1)
  
  for(t in 2:(Time+Time.fore)){
    
    N.tot[t] <- N.b[t] + N.i[t]
    
    #s.p1 is transition probability * survival
    pre.br1[t] <- f[t-1] * mean.fem * N.tot[t-1] * s.p1 * mu.n.fem #+ 
    
    
    no.breed[t] <- (N.nb[t-1] * phi[t-1] * (1-psiNBB)) + 
                         (N.b[t-1] * phi[t-1] * psiBNB)
    breed[t] <- (N.p7[t-1] * s.p1) +
                      (N.nb[t-1] * phi[t-1] * psiNBB) +
                      (N.b[t-1] * phi[t-1] * (1-psiBNB))
    
    N.p1[t] ~ dpois(pre.br1[t])
    #p2[t] <- (s.p1*N.p1[t-1])
    #N.p2[t] ~ dpois(p2[t])
    N.p2[t] ~ dbinom(s.p1,N.p1[t-1])
    #p3[t] <- (s.p1*prob.trans*N.p2[t-1]) + (1-prob.trans)*N.p3[t-1]
    #N.p3[t] ~ dpois(p3[t])
    N.p3[t] ~ dbinom(s.p1,N.p2[t-1])
    #p4[t] <- (s.p1*prob.trans*N.p3[t-1]) + (1-prob.trans)*N.p4[t-1]
    #N.p4[t] ~ dpois(p4[t])
    N.p4[t] ~ dbinom(s.p1,N.p3[t-1])
    #p5[t] <- (s.p1*prob.trans*N.p4[t-1]) + (1-prob.trans)*N.p5[t-1]
    #N.p5[t] ~ dpois(p5[t])
    N.p5[t] ~ dbinom(s.p1,N.p4[t-1])
    #p6[t] <- (s.p1*prob.trans*N.p5[t-1]) + (1-prob.trans)*N.p6[t-1]
    #N.p6[t] ~ dpois(p6[t])
    N.p6[t] ~ dbinom(s.p1,N.p5[t-1])
    #p7[t] <- (s.p1*prob.trans*N.p6[t-1]) + (1-prob.trans)*N.p7[t-1]
    #N.p7[t] ~ dpois(p7[t])
    N.p7[t] ~ dbinom(s.p1,N.p6[t-1])
    #N.p8[t] ~ dpois(N.p7[t-1]*s.p1)
    #N.p8[t] ~ dbinom(s.p1,N.p7[t-1])
    #N.p9[t] ~ dpois(N.p8[t-1]*s.p1)
    #N.p9[t] ~ dbinom(s.p1,N.p8[t-1])
    #N.p10[t] ~ dpois(N.p9[t-1]*s.p1)
    #N.p10[t] ~ dbinom(s.p1,N.p9[t-1])
    N.nb[t] ~ dpois(no.breed[t])
    N.b[t] ~ dpois(breed[t])
    log(imm[t]) <- eps.imm[t]
    N.i[t] ~ dpois(imm[t])
    #N.i[t] <- round(exp(eps.imm[t]))
  }
  
  #number of females seen as count data
  
  #for(t in 1:(Time+Time.fore)){
  for(t in 1:Time){
    #y.fem[t] ~ dpois(N.b[t] + N.i[t])
    y.fem[t] ~ dpois(N.tot[t])
  }
  
})

y.fem[1] = 0
ssm_data <- list(y.fem = y.fem,
                              s.p.mom=s.p.mom,
                              N.p1.ini = 200, #initial values taken from Cailouet et al 1995, Table 1
                              N.p2.ini = 130, #survival of 0.6 assumed each year
                              N.p3.ini = 86,
                              N.p4.ini = 50,
                              N.p5.ini = 33,
                              N.p6.ini = 20,
                              N.p7.ini = 13,
                              N.p8.ini = 10,
                              N.p9.ini = 5,
                              N.p10.ini = 3)#,
                              #N.b.ini=2, 
                              #N.nb.ini=2)

ssm_con <- list(Time = ncol(dat), 
                 Time.fore = 45)


ssm_ini = function() list(
  s.p1 = runif(1,.1,.5), 
  prob.trans = runif(1,.1,.9),
  N.nb = c(NA,rpois((Time+Time.fore-1),seq(2,200,length.out=(Time+Time.fore-1)))),
  #N.b = rpois((Time+Time.fore),c(trun[1:Time]/2,rep(50,Time.fore))),
  N.b = c(NA,rpois((Time+Time.fore-1),seq(2,100,length.out=(Time+Time.fore-1)))),
  N.p1 = rpois((Time+Time.fore),seq(200,5000,length.out=(Time+Time.fore))),
  N.p2 = rpois((Time+Time.fore),seq(130,2000,length.out=(Time+Time.fore))),
  N.p3 = rpois((Time+Time.fore),seq(85,1000,length.out=(Time+Time.fore))),
  #N.p2 = rpois((Time+Time.fore),10),
  #N.p3 = rpois((Time+Time.fore),10),
  N.p4 = rpois((Time+Time.fore),seq(50,800,length.out=(Time+Time.fore))),
  N.p5 = rpois((Time+Time.fore),seq(30,400,length.out=(Time+Time.fore))),
  N.p6 = rpois((Time+Time.fore),seq(20,200,length.out=(Time+Time.fore))),
  N.p7 = rpois((Time+Time.fore),seq(13,100,length.out=(Time+Time.fore))),
  #N.p8 = rpois((Time+Time.fore),seq(10,80,length.out=(Time+Time.fore))),
  #N.p9 = rpois((Time+Time.fore),seq(8,60,length.out=(Time+Time.fore))),
  #N.p10 = rpois((Time+Time.fore),seq(6,45,length.out=(Time+Time.fore))),
  #N.i = rpois((Time+Time.fore),c(trun[1:Time]/2,rep(50,Time.fore))),
  N.i = c(NA,rpois((Time+Time.fore-1),seq(1,100,length.out=(Time+Time.fore-1)))),
  mu.imm = runif(1,-1,1),
  #sd.imm = runif(1,.01,1),
  eps.imm.tmp = runif((Time+Time.fore),-1,1),
  y.fem = c(rep(NA,Time),rep(150,Time.fore)))

#### nimble run ####

ssm.model <- nimbleModel(code = ssm.ipm,
                          constants = ssm_con,
                          data = ssm_data,
                          inits = ssm_ini())

ssm.config.mcmc <- configureMCMC(ssm.model)

#ssm.config.mcmc$addMonitors("N.i")

ssm.mcmc <- buildMCMC(ssm.config.mcmc)

ssm.cmodel <- compileNimble(ssm.model)

ssm.cmcmc <- compileNimble(ssm.mcmc, project = ssm.model)

ssm.samples <- runMCMC(ssm.cmcmc, niter= 11000, 
                        nburnin = 10000,
                        nchains = 3,
                        samplesAsCodaMCMC = TRUE)

out = as.mcmc(rbind(ssm.samples[[1]], ssm.samples[[2]], 
                    ssm.samples[[3]]))

library(parallel)

#### nimble cluster ####

nc = 3 #num of chains
clus = makeCluster(nc)

clusterExport(clus, c("ssm.ipm","ssm_ini",
                      "ssm_data",
                      "ssm_con"))

for(j in seq_along(clus)){
  set.seed(j)
  init <- ssm_ini()
  clusterExport(clus[j],"init")
}

out <- clusterEvalQ(clus, {
  library(nimble)
  library(coda)
  model <- nimbleModel(code = ssm.ipm, 
                       name = "ssm.ipm",
                       constants = ssm_con,
                       data = ssm_data, 
                       inits = init)
  Cmodel <- compileNimble(model)
  modelConf <- configureMCMC(model)
  #modelConf$addMonitors(params)
  modelConf$addMonitors(c("N.i","N.b"))
  #configureRJ(modelConf,
  #            targetNodes = 'beta',
  #            indicatorNodes = 'indA',
  #            control = list(mean=0, scale = 2))
  modelMCMC <- buildMCMC(modelConf)
  CmodelMCMC <- compileNimble(modelMCMC, 
                              project = model)
  #out1 <- runMCMC(CmodelMCMC, niter=15000, nburnin=14000)
  #return(as.mcmc(out1))
  CmodelMCMC$run(50000,reset=TRUE, resetMV=TRUE)
  return(as.mcmc(as.matrix(CmodelMCMC$mvSamples)))
})

out.mcmc <- as.mcmc.list(out) 

#stopCluster(clus)

out2 <- clusterEvalQ(clus,{
  CmodelMCMC$run(10000,reset=FALSE, resetMV=TRUE)
  return(as.mcmc(as.matrix(CmodelMCMC$mvSamples)))
})

out.up <- as.mcmc.list(out2)
