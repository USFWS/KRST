
nest.code <- nimbleCode({

  mu.egg ~ dnorm(0, sd = 2)
  for(s in 1:5){
    nest.cov[s] ~ dnorm(0, sd = 2)
    beta.temp[s] ~ dnorm(mu.temp, var = var.temp)
  }
  nest.cov[6] <- log(insitu.surv/(1-insitu.surv))
  beta.temp[6] ~ dnorm(-1,sd=1)
  beta.temp.all ~ dnorm(mu.temp,var = var.temp)
  #beta.temp.all2 ~ dnorm(0,sd=2)
  mu.temp ~ dnorm(0,sd = 1)
  var.temp ~ dinvgamma(1,1)
  
  #first 24 years (until 2015) only have data on
  # total number of eggs and total hatchlings released
  for(t in 1:24){
    #number of successful hatchlings released
    # out of total number of eggs
    
    #J is total hatchlings released
    J[t] ~ dbinom(size = E[t], prob = prob.egg[t]) 
    
    #prob.egg is probability of egg surviving to hatchling
    logit(prob.egg[t]) <- mu.egg + beta.temp.all*may.tmax.ann[t]
    
    #f is hatchlings per nest per female per year
    f[t] ~ dpois(E.R[t]*prob.egg[t])
  }
  
  #last 7 years have data on number of eggs, nests,
  # and incubation type for each nest
  
  #insitu survival (from Shaver et al. 2020)
  insitu.surv ~ dbeta(in.s[1],in.s[2])
  
  for(i in 1:N.nest){
    J.ind[i] ~ dbinom(size = E.ind[i], prob = prob.egg.ind[i])
    #E.ind2[i] ~ dpois(E.ind[i])
    logit(prob.egg.ind[i]) <- nest.cov[nest.mange[i]] + 
      beta.temp[nest.mange[i]]*may.tmax[i]
  }
  
  f[25] ~ dpois(mean(E.ind[1:year.mean.idx[25]])*mean(prob.egg.ind[1:year.mean.idx[25]]))
  
  for(t in 26:Time){
    f[t] ~ dpois(mean(E.ind[(year.mean.idx[t-1]+1):year.mean.idx[t]])*mean(prob.egg.ind[(year.mean.idx[t-1]+1):year.mean.idx[t]]))
  }
  
  for(t in (Time+1):(Time+Time.fore)){
    logit(prob.egg.fore.tmp[1:6,(t-24)]) <- nest.cov[1:6] + beta.temp[1:6]*may.tmax[t]
    prob.egg.fore[(t-24)] <- sum(prob.egg.fore.tmp[1:6,(t-24)]*nest.mange.fore[1:6])/sum(nest.mange.fore[1:6])
    f[t] ~ dpois(E.ind.fore*prob.egg.fore[(t-24)])
  }
  
})
    
#j.e = J/E
#j.e[2:3] = 0.01
incu = nest.pro$incubation.type
incu = incu[-which(incu=="S")]
incu[incu == "I" | incu == "I "] = 1
incu[incu == "C"] = 0
incu = as.numeric(incu)
incu=incu[-1448]

#nest condition index. G = good, others are "poor"
cond = nest.pro$condition.code[-which(nest.pro$incubation.type=="S")]
cond = cond[-1448]
cond[cond == "G"] = 1
cond[cond == "N"] = 0
cond[cond == "O"] = 0
cond[cond == "L"] = 0
cond = as.numeric(cond)

good.incu = rep(0,1569)
good.corr = rep(0,1569)
ngood.incu = rep(0,1569)
ngood.corr = rep(0,1569)
good.incu[cond == 1 & incu == 1] = 1
good.corr[cond == 1 & incu == 0] = 1
ngood.incu[cond == 0 & incu == 1] = 1
ngood.corr[cond == 0 & incu == 0] = 1

j.e.tmp = nest.pro$emerge.success[-which(nest.pro$incubation.type=="S")]/100
j.e.tmp = j.e.tmp[-1448]

J = nest.pro$hatch.released[-which(nest.pro$incubation.type=="S")]
J = J[-1448]
E = nest.pro$total.eggs[-which(nest.pro$incubation.type=="S")]
E = E[-1448]

egg.clim = rnorm(7,0,1)

year = nest.pro$Year[-which(nest.pro$incubation.type == "S")]
year = year - 2014
year = year[-1448]
year.mean.idx = as.numeric()
for(t in 1:7){
  year.mean.idx[t] = min(which(year == t))
}

nest_data <- list(#J = tx_nest_dat$hatch_release[-c(1:4)], 
                  E.R = E.R,
  J = J,                
  #j.e = j.e.tmp)
                  E = E,in.s = in.s,
  J.ind = J.ind,
  E.ind = E.ind)

nest_con <- list(N.nest = length(J.ind), Time = 31, 
                 nest.mange = nest.mange,
                 year = year,
                 year.mean.idx = year.mean.idx,
                 N.nest = length(J.ind),
                 #may.tmax.ann = may.tmax.ann,
                 may.tmax.ann = bio1.ann,
                 may.tmax = bio1.out)

E.ini = E.ind
E.ini[which(is.na(E.ini))] = runif(length(which(is.na(E.ini))),81,100)
E.ini[which(!is.na(E.ind))] = NA

nest_ini = function() list(
  mu.egg = runif(1,-1,1),
  mu.temp = runif(1,-1,1),
  var.temp = runif(1,.1,2),
  nest.cov = c(runif(5,-1,1),NA),
  f = rpois(31,76))

nest.model <- nimbleModel(code = nest.code,
                          constants = nestsim_con,
                          data = nestsim_data,
                          inits = nestsim_ini())

nest.config.mcmc <- configureMCMC(nest.model)

nest.config.mcmc$addMonitors("f","beta.temp.all")

nest.mcmc <- buildMCMC(nest.config.mcmc)

nest.cmodel <- compileNimble(nest.model)

nest.cmcmc <- compileNimble(nest.mcmc, project = nest.model)

nest.samples <- runMCMC(nest.cmcmc, niter= 2000, 
                        nburnin = 1000,nchains = 3,
                        samplesAsCodaMCMC = TRUE)


#sim data for binomial nest data
#6 nest treatments
b.clim = rnorm(7,-1,sd=2)
x = rnorm(10,0,1)
x2 = rnorm(24,0,1)
b1 = -1
b2 = 1
b3 = 1.5
b4 = -.5
b5 = .5
b6 = -2

#p1 = 1/(1+exp(-(b1)))
#p2 = 1/(1+exp(-(b2)))
#p3 = 1/(1+exp(-(b3)))
#p4 = 1/(1+exp(-(b4)))
#p5 = 1/(1+exp(-(b5)))
#p6 = 1/(1+exp(-(b6)))
#p7 = 1/(1+exp(-(b1)))

p1 = 1/(1+exp(-(b1+(b.clim[1]*x))))
p2 = 1/(1+exp(-(b2+b.clim[2]*x)))
p3 = 1/(1+exp(-(b3+b.clim[3]*x)))
p4 = 1/(1+exp(-(b4+b.clim[4]*x)))
p5 = 1/(1+exp(-(b5+b.clim[5]*x)))
p6 = 1/(1+exp(-(b6+b.clim[6]*x)))
p7 = 1/(1+exp(-(b1+b.clim[7]*x2)))

j.sim = c(rbinom(100,90,rep(p1,10)),
          rbinom(100,90,rep(p2,10)),
          rbinom(100,90,rep(p3,10)),
          rbinom(100,90,rep(p4,10)),
          rbinom(100,90,rep(p5,10)),
          rbinom(100,90,rep(p6,10)),
          rbinom(24,900,p7))

time = rep(seq(1,10),60)
dat.sim = data.frame(j.sim = j.sim[1:600], time = time,
                     cov = rep(x,60), nest = rep(1:6,each=100))
dat.sim = dat.sim[order(dat.sim$time),]

year = dat.sim$time+24
year.mean.idx = as.numeric()
for(t in 25:34){
  year.mean.idx[t] = max(which(year == t))
}

nests.fore = 250*Time.fore #250 nests each year in future
n.fore.fems = (225/250)*nests.fore

year[(year.mean.idx[34]+1):(year.mean.idx[34]+nests.fore)] = rep(35:(34+Time.fore),each=250)

for(t in 35:(34+Time.fore)){
  year.mean.idx[t] = max(which(year == t))
}

nest.mange = dat.sim$nest

nest.mange.fore = c(10,20,10,0,5,20)

J.ind = dat.sim$j.sim
E.ind = rep(90,600)

E.ind[(year.mean.idx[34]+1):(year.mean.idx[34]+nests.fore)] = round(mean(E.ind))
J.ind[(year.mean.idx[34]+1):(year.mean.idx[34]+nests.fore)] = NA

#forecast nest management scenarios
#%50 in corrals, %50 in situ
mange.fore = c(rep(1,63),rep(2,62),rep(6,125))
nest.mange[(year.mean.idx[34]+1):(year.mean.idx[34]+nests.fore)] = rep(mange.fore,Time.fore)

clim.cov = c(dat.sim$cov,rep(.1,10),rep(-.1,10),rep(-.3,10))

nestsim_data <- list(E.R = rep(90,24), J = j.sim[601:624],                
  E = rep(900,24),in.s = in.s,
  J.ind = J.ind,
  E.ind = E.ind,
  E.ind.fore =90)

nestsim_con <- list(N.nest = 600, Time = 34, 
                    Time.fore = 30,
                 nest.mange = nest.mange,
                 year = time,
                 year.mean.idx = year.mean.idx,
                 may.tmax.ann = x2,
                 may.tmax = clim.cov,
                 nest.mange.fore = nest.mange.fore)

J.ini=J.ind
J.ini[is.na(J.ini)] = 70
prob.egg.ini = c(rep(NA,600),runif(12500,.1,.9))

nestsim_ini = function() list(
  mu.egg = runif(1,-1,1),
  mu.temp = runif(1,-1,1),
  var.temp = runif(1,.1,2),
  nest.cov = c(runif(5,-1,1),NA),
  beta.temp = runif(6,-1,1),
  beta.temp.all = runif(1,-1,1),
  f = rpois(64,76),
  insitu.surv = runif(1,.1,.9))
  #J.ind = J.ini,
  #prob.egg.ind = prob.egg.ini)

nestsim.model <- nimbleModel(code = nest.code,
                             constants = nestsim_con,
                             data = nestsim_data,
                             inits = nestsim_ini())

nestsim.config.mcmc <- configureMCMC(nestsim.model)

nestsim.config.mcmc$addMonitors("f","beta.temp.all")

nestsim.mcmc <- buildMCMC(nestsim.config.mcmc)

nestsim.cmodel <- compileNimble(nestsim.model)

nestsim.cmcmc <- compileNimble(nestsim.mcmc, project = nestsim.model)

nestsim.samples <- runMCMC(nestsim.cmcmc, niter= 2000, 
                        nburnin = 1000,nchains = 3,
                        samplesAsCodaMCMC = TRUE)

nest.f.code <- nimbleCode({
  
  mu.egg ~ dnorm(0, sd = 2)
  for(s in 1:5){
    nest.cov[s] ~ dnorm(0, sd = 2)
    beta.temp[s] ~ dnorm(mu.temp, var = var.temp)
  }
  nest.cov[6] <- log(insitu.surv/(1-insitu.surv))
  beta.temp[6] ~ dnorm(-1,sd=1)
  beta.temp.all ~ dnorm(mu.temp,var = var.temp)
  #beta.temp.all2 ~ dnorm(0,sd=2)
  mu.temp ~ dnorm(0,sd = 1)
  var.temp ~ dinvgamma(1,1)
  
  #first 24 years (until 2015) only have data on
  # total number of eggs and total hatchlings released
  for(t in 1:24){
    #number of successful hatchlings released
    # out of total number of eggs
    
    #J is total hatchlings released
    #J[t] ~ dbinom(size = E[t], prob = prob.egg[t]) 
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*f[t]
    log(f[t]) <- mu.egg + beta.temp.all*may.tmax.ann[t]
    
    #prob.egg is probability of egg surviving to hatchling
    #logit(prob.egg[t]) <- mu.egg + beta.temp.all*may.tmax.ann[t]
    
    #f is hatchlings per nest per female per year
    #f[t] ~ dpois(E.R[t]*prob.egg[t])
  }
  
  #last 7 years have data on number of eggs, nests,
  # and incubation type for each nest
  
  #insitu survival (from Shaver et al. 2020)
  insitu.surv ~ dbeta(in.s[1],in.s[2])
  
  for(i in 1:N.nest){
    J.ind[i] ~ dpois(f.ind[i])
    log(f.ind[i]) <- nest.cov[nest.mange[i]] + beta.temp[nest.mange[i]]*may.tmax[i]
  }
  
  f[25] ~ dpois(mean(E.ind[1:year.mean.idx[25]])*mean(prob.egg.ind[1:year.mean.idx[25]]))
  
  #for(t in 26:(Time+Time.fore)){
  #  f[t] ~ dpois(mean(E.ind[(year.mean.idx[t-1]+1):year.mean.idx[t]])*mean(prob.egg.ind[(year.mean.idx[t-1]+1):year.mean.idx[t]]))
  #}
  
  
  #for(i in 1:N.nest){
  #  J.ind[i] ~ dbinom(size = E.ind[i], prob = prob.egg.ind[i])
  #  #E.ind2[i] ~ dpois(E.ind[i])
  #  logit(prob.egg.ind[i]) <- nest.cov[nest.mange[i]] + 
  #    beta.temp[nest.mange[i]]*may.tmax[i]
  #}
  
  #f[25] ~ dpois(mean(E.ind[1:year.mean.idx[25]])*mean(prob.egg.ind[1:year.mean.idx[25]]))
  
  #for(t in 26:(Time+Time.fore)){
  #  f[t] ~ dpois(mean(E.ind[(year.mean.idx[t-1]+1):year.mean.idx[t]])*mean(prob.egg.ind[(year.mean.idx[t-1]+1):year.mean.idx[t]]))
  #}
  
})