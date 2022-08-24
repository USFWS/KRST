library(janitor)
library(readxl)
library(lubridate)
library(nimble)
library(nimbleEcology)
library(stringr)

cmr = read_xls("NPS CMR data updated.xls",sheet = 1, col_names = TRUE)
cmr.df = as.data.frame(cmr)

#need to clean up id numbers

#first remove "2 transmitters" and "3 transmitters"
cmr.df$id[cmr.df$id=="2 transmitters"] = NA
cmr.df$id[cmr.df$id=="3 transmitters"] = NA

#need to add in 114 which got renumbered
cmr.df$id[715] = 114

idx = str_sort(cmr.df$id,na_last=TRUE,decreasing=TRUE)
cmr.df$id[which(cmr.df$id==idx[10])]=813
cmr.df$id[which(cmr.df$id==idx[7])]=814
cmr.df$id[which(cmr.df$id==idx[6])]=815
cmr.df$id[which(cmr.df$id==idx[5])]=816
cmr.df$id[which(cmr.df$id==idx[4])]=817
cmr.df$id[which(cmr.df$id==idx[3])]=818
cmr.df$id[which(cmr.df$id==idx[2])]=819
cmr.df$id[which(cmr.df$id==idx[1])]=820
cmr.df$id[which(cmr.df$id==idx[9])]=821
cmr.df$id[which(cmr.df$id==idx[8])]=822

#remove comments/turtle names
for(i in 11:27){
cmr.df$id[which(cmr.df$id==idx[i])]=NA}

cmr.df$id=as.numeric(cmr.df$id)

#remove text from bottom of spreadsheet
cmr.df = cmr.df[1:3694,]

#add in id numbers for NAs. In spreadsheet,
# only the first capture has an id number assigned

id = cmr.df$id

#manually fix mislabeled ids (813:822)
id[9] = 4
id[231] = 47
id[233] = 814
id[280:281] = 57
id[283] = 815
id[388:403] = 71
id[405] = 816
id[442:443] = 81
id[445] = 817
id[758:759] = 118
id[761] = 818
id[763] = 819
id[1998] = 347
id[2000:2001] = 820
id[1534] = 254
id[1536:1538] = 821
id[3693] = 812

skip.id = c(4,47,57,71,81,118,347,254)

#(missing ids: 114, 227, 291, 506, 556, 588, 604)

for(i in (1:811)[-skip.id]){
  tmp1 = which(id == i)
  tmp2 = which(id == i+1)[1]
  if(tmp2-tmp1==1) next
  id[(tmp1+1):(tmp2-1)] = rep(i,(tmp2-tmp1-1))
}

cmr.df$id2 = id

#fix dates using stringr
idx2 = str_split(cmr.df$clutch,"-",simplify = TRUE)[,1]
#cmr.df$date_nest[which(idx2== "MEXICO")]

#made up date (kept year) to match other dates
cmr.df$date_nest[which(idx2== "MEXICO")][7] = "16 May 2015"

#extract dates from date_nest and move them to idx2
idx2[which(idx2=="MEXICO")] = year(as.Date(cmr.df$date_nest[which(idx2== "MEXICO")],
        format = "%d %B %Y"))

idx.tmp = str_sort(idx2,na_last=TRUE,decreasing=TRUE)

for(i in 1:56){
idx2[which(idx2 == idx.tmp[i])] = NA}

cmr.df$year = as.numeric(idx2)

#cmr.df$date_nest[11] = "13 June 1997"
#cmr.df$date_nest[21] = "19 April 2000"
#cmr.df$date_nest[22] = "26 April 2002"
#cmr.df$date_nest[23] = "16 May 2002"
#cmr.df$date_nest[24] = "26 April 2006"
#cmr.df$date_nest[25] = "11 May 2006"
#cmr.df$date_nest[26] = "26 April 2008"
#cmr.df$date_nest[27] = "30 April 2008"
#cmr.df$date_nest[28] = "16 May 2008"
#cmr.df$date_nest[29] = "1 June 2008"
#cmr.df$date_nest[30] = "18 June 2012"
#cmr.df$date_nest[137] = "3 June 2008" #was 1008
#cmr.df$date_nest[215] = "1 May 2015" #no day or month in excel file
#cmr.df$date_nest[468] = "30 April 2007"
#cmr.df$date_nest[1761] = "1 May 2012" #only 'in situ' recorded for 2012
#cmr.df$date_nest[1818] = "7 May 2007"
#cmr.df$date_nest[1819] = "22 May 2007"
#cmr.df$date_nest[1921] = "1 May 2013" #only 'in situ' recorded for 2013
#cmr.df$date_nest[1941] = "3 July 2018"
#cmr.df$date_nest[2278] = '1 May 2008' #'in situ' for 2008
#cmr.df$date_nest[2776] = '8 August 2018'

dat = matrix(0, nrow=822, ncol = (2023-1991))
#wild = rep(NA, 812)

for(i in 1:822){
  tmp = subset(cmr.df,id2==i)
  #tmp2 = as.Date(tmp$date_nest, format="%d %B %Y")
  #years = year(tmp2)
  dat[i,(tmp$year-1990)] = 1
  #wild[i] = max(ifelse(tmp$yr_hatch == "W",1,0),na.rm=TRUE)
}

#wild=wild[-c(114, 227, 291, 506, 556, 588, 604,801)] #renumbered individuals
#wild[!is.finite(wild)] = 1

dat = dat[-(which(apply(dat,1,sum)==0)),] #remove individuals that were renumbered in dataset (missing ids: 114, 227, 291, 506, 556, 588, 604,801)

tx_nest_dat = read.csv("tx_nest_hatch_data.csv")
colnames(tx_nest_dat)[1] = "year"

#nest searching effort
nest.effort = read_xlsx("patrol effort data_NPI_1986 to 2021.xlsx",sheet = 1, col_names = TRUE)
kilo.patrol = nest.effort$`Kilometers patrolled`[nest.effort$Year>1990]
#nest protection data
#hatch success = prob of hatching
#emerge success = prob of being released
nest.pro = read.csv("Nesting Data 2015-2021.csv")

#get values for priors
source("~/LPC/beta_mom.R",echo=TRUE)
phi.mom <- beta.mom(.9,.1) #ad survival is high
psiBNB.mom <- beta.mom(.8,.2) #prob of B to NB is low
psiNBB.mom <- beta.mom(.7,.2) #prob of NB to B is higher
p.mom <- beta.mom(.8,.2) #detection prob is high
s.p.mom <- beta.mom(.43,.1)
s.p.mom2 <- beta.mom(.5,.1)
in.s <- beta.mom(.63,.1)

krst.ipm <- nimbleCode({
  
  #### priors for adult survival ####
  
  for(t in 1:((Time+Time.fore)-1)){
  #phi[t] ~ dbeta(phi.mom[1],phi.mom[2]) #surv prob of breeders
  l.phi[t] <- ln.mu.phi[t]
  l.lim[t] <- min(99, max(-99, l.phi[t]))
  phi[t] <- 1/(1+exp(-l.lim[t]))
  }
  
  #for(j in 1:max(phi.ind)){
    #ln.mu.phi[j] ~ dnorm(mu.phi,var = sd.phi2)
  for(t in 1:((Time+Time.fore)-1)){
    ln.mu.phi[t] ~ dnorm(mu.phi,var=sd.phi2)
  }
  mu.phi ~ dnorm(2.2, sd = 1)
  sd.phi2 ~ dinvgamma(1,1)

  psiBNB ~ dbeta(psiBNB.mom[1], psiBNB.mom[2]) #transition of Breeders to NonBreeders
  psiNBB ~ dbeta(psiNBB.mom[1], psiNBB.mom[2]) #transition of NB to B
  #pB ~ dbeta(p.mom[1], p.mom[2]) #det of breeders
  a.0 ~ dnorm(0,sd = 1)
  a.k ~ dnorm(0,sd = 1)
  
  #model for detection probability by kilo searched
  for(t in 1:(Time+Time.fore)){
  logit(pB[t]) <- a.0 + a.k*kilo.patrol[t]
  }
  
  #priors pre breeding survival
  #s.p1 ~ dbeta(1,1) #same table 3, hatch & 1yr surv
  s.p1 ~ dbeta(s.p.mom[1],s.p.mom[2])
  
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
  N.b[1] ~ dpois(N.b.ini)
  #N.b[1] <- 1
  N.nb[1] ~ dpois(N.nb.ini)
  
  #probabilities of state z(t+1) given z(t)
  for(t in 1:((Time+Time.fore)-1)){
  gamma[1,1,t] <- phi[t] * (1-psiBNB)
  gamma[1,2,t] <- phi[t] * psiBNB
  gamma[1,3,t] <- 1-phi[t]
  
  gamma[2,1,t] <- phi[t] * psiNBB
  gamma[2,2,t] <- phi[t] * (1-psiNBB)
  gamma[2,3,t] <- 1-phi[t]
  
  gamma[3,1,t] <- 0
  gamma[3,2,t] <- 0
  gamma[3,3,t] <- 1
  #}
  
  #probabilities of y(t) given z(t)
  omega[1,1,t] <- 1-pB[t] #pr(alive B t -> non-detected)
  omega[1,2,t] <- pB[t] #pr(alive B -> detected B t)
  
  omega[2,1,t] <- 1 #pr(alive NB -> non-detected)
  omega[2,2,t] <- 0 #pr(alive NB -> detected)
  
  omega[3,1,t] <- 1 #pr(dead t -> nondetected)
  omega[3,2,t] <- 0 #pr(dead t -> detected)
  }
  
  ##
  #####HMM for adult survival and breeding probability####
  ##
  
  #Need to fix "dynamic index out of bounds" warning??
  
  for(i in 1:N){
    #latent state at first capture
    z[i,first[i]] <- 1
    
    for(t in (first[i]+1):Time){
      #z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:3,t-1])
      #y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:2,t])
    }
  }

  ##
  #####Nest success model####
  ##
  
  mu.egg ~ dnorm(0, sd = 2)
  for(s in 1:5){
    nest.cov[s] ~ dnorm(0, sd = 0.5)
    beta.temp[s] ~ dnorm(mu.temp, var = var.temp)
    beta.slr[s] <- 0
  }
  nest.cov[6] <- log(insitu.surv/(1-insitu.surv))
  beta.temp[6] ~ dnorm(-1,sd=1)
  beta.slr[6] ~ dnorm(-1,sd=1)
  beta.temp.all ~ dnorm(mu.temp, var = var.temp)
  #beta.temp.all2 ~ dnorm(0,sd=2)
  mu.temp ~ dnorm(0, sd = 1)
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
    logit(prob.egg.ind[i]) <- nest.cov[nest.mange[i]] + 
      beta.temp[nest.mange[i]]*clim.cov[i] +
      beta.slr[nest.mange[i]]*slr.cov[i]
  }
  
  f[25] ~ dpois(mean(E.ind[1:year.mean.idx[25]])*mean(prob.egg.ind[1:year.mean.idx[25]]))
  
  for(t in 26:Time){
    f[t] ~ dpois(mean(E.ind[(year.mean.idx[t-1]+1):year.mean.idx[t]])*mean(prob.egg.ind[(year.mean.idx[t-1]+1):year.mean.idx[t]]))
  }
  
  for(t in (Time+1):(Time+Time.fore)){
    logit(prob.egg.fore.tmp[1:6,(t-31)]) <- nest.cov[1:6] + beta.temp[1:6]*clim.cov.fore[(t-31)] +
      beta.slr[1:6]*slr.cov.fore[(t-31)]
    prob.egg.fore[(t-31)] <- sum(prob.egg.fore.tmp[1:6,(t-31)]*nest.mange.fore[1:6])/sum(nest.mange.fore[1:6])
    f[t] ~ dpois(E.ind.fore*prob.egg.fore[(t-31)])
  }
  
  ##
  #####number of nests per female per year####
  ##
  
  for(t in 1:(Time+Time.fore)){  
  n.fem[t] ~ dgamma(mean = mu.n.fem, sd = sd.n.fem)
  }
  
  mu.n.fem ~ dunif(1,4)
  sd.n.fem ~ dunif(0,10)
  
  #Probability of hatchling being female
  for(t in 1:Time+Time.fore){
    prob.fem[t] ~ dbeta(r*mu.fem[t], r*(1-mu.fem[t])) #Data from Shaver & Calliway
    logit(mu.fem[t]) <- ln.mean.fem
  }
  
  ln.mean.fem ~ dnorm(0,sd=10) #prior, mean prob
  mean.fem <- exp(ln.mean.fem)/(1+exp(ln.mean.fem))
  r ~ dgamma(.1,.1) #overdispersion param for beta dist'n
  
  #state space model for abundance/pop matrix
  
  pre.br1[1] <- N.p1[1]
  pre.br2[1] <- N.p2[1]
  pre.br3[1] <- N.p3[1]
  pre.br4[1] <- N.p4[1]
  pre.br5[1] <- N.p5[1]
  pre.br6[1] <- N.p6[1]
  pre.br7[1] <- N.p7[1]
  pre.br8[1] <- N.p8[1]
  pre.br9[1] <- N.p9[1]
  pre.br10[1] <- N.p10[1]
  no.breed[1] <- N.nb[1]
  breed[1] <- N.b[1]
  imm[1] <- 1
  N.i[1] <- 1
  mu.imm ~ dnorm(0, sd = 2)
  sd.imm ~ dinvgamma(1,1)
  
  for(t in 1:(Time+Time.fore)){
    eps.imm.tmp[t] ~ dnorm(mu.imm, var=sd.imm)
    eps.imm[t] <- min(eps.imm.tmp[t],8)
  }
  
  for(t in 2:(Time+Time.fore)){
    
    #s.p1 is transition probability * survival
    pre.br1[t] <- f[t-1] * mean.fem * N.b[t-1] * s.p1 * mu.n.fem + 
      (N.p1[t-1]*(1-s.p1))
    pre.br2[t] <- N.p1[t-1]*s.p1 +
      N.p2[t-1]*(1-s.p1) 
    pre.br3[t] <- N.p2[t-1]*s.p1 +
      N.p3[t-1]*(1-s.p1)
    pre.br4[t] <- N.p3[t-1]*s.p1 +
      N.p4[t-1]*(1-s.p1)
    pre.br5[t] <- N.p4[t-1]*s.p1 +
      N.p5[t-1]*(1-s.p1)
    pre.br6[t] <- N.p5[t-1]*s.p1 +
      N.p6[t-1]*(1-s.p1)
    pre.br7[t] <- N.p6[t-1]*s.p1 +
      N.p7[t-1]*(1-s.p1)
    pre.br8[t] <- N.p7[t-1]*s.p1 +
      N.p8[t-1]*(1-s.p1)
    pre.br9[t] <- N.p8[t-1]*s.p1 +
      N.p9[t-1]*(1-s.p1)
    pre.br10[t] <- N.p9[t-1]*s.p1 +
      N.p10[t-1]*(1-s.p1)
    
    no.breed[t] <- max(1,(N.nb[t-1] * phi[t-1] * (1-psiNBB)) + 
      (N.b[t-1] * phi[t-1] * psiBNB))
    breed[t] <- max(1,(N.p10[t-1] * s.p1) +
      (N.nb[t-1] * phi[t-1] * psiNBB) +
      (N.b[t-1] * phi[t-1] * (1-psiBNB)))
                       
    N.p1[t] ~ dpois(pre.br1[t])
    N.p2[t] ~ dpois(pre.br2[t])
    N.p3[t] ~ dpois(pre.br3[t])
    N.p4[t] ~ dpois(pre.br4[t])
    N.p5[t] ~ dpois(pre.br5[t])
    N.p6[t] ~ dpois(pre.br6[t])
    N.p7[t] ~ dpois(pre.br7[t])
    N.p8[t] ~ dpois(pre.br8[t])
    N.p9[t] ~ dpois(pre.br9[t])
    N.p10[t] ~ dpois(pre.br10[t])
    N.nb[t] ~ dpois(no.breed[t])
    N.b[t] ~ dpois(breed[t])
    log(imm[t]) <- eps.imm[t]
    N.i[t] ~ dpois(imm[t])
  }
  
  for(t in 1:(Time+Time.fore)){
    Ntot[t] <- N.b[t] + N.nb[t]
  }
  
  #number of females seen as count data

  for(t in 1:(Time+Time.fore)){
    #tried to constrain distribution in some way
    # to ensure y.fem was the minimum known alive
    # rather than a typical state space model.
    # Could not get reasonable results.
    
    #y.fem[t] ~ dnorm(N.b[t] + N.i[t],var=sd.y2)
    #y.fem[t] ~ T(dnorm(N.b[t] + N.i[t], sd=sd.y2),1,500)
    #y.fem[t] ~ dpois(theta*(N.b[t] + N.i[t]))
    y.fem[t] ~ dpois(N.b[t] + N.i[t])
    #y.fem[t] ~ T(dpois(N.b[t] + N.i[t]),trun[t],1000)
    #constraint_data[t] ~ dconstraint(N.b[t] + N.i[t] > trun[t])
    }
  
})

f_ini <- mean(tx_nest_dat$hatch_per_nest)

#CMR data starts in 1991 - 2020
#Nest data starts w/ four years prior to 1991 - 2014, minus total in last row

#remove last column (2022) from dat, not complete data
dat = dat[,1:31]

y=dat+1
#z.ini=zinits_dcat+1
first = apply(dat, 1, function(x) min(which(x !=0)))

fem.dat = read.csv("hatch_sexes.csv")

#forecasting
Time = dim(dat)[2]
Time.fore = 45

prob.fem = c(fem.dat$prob_fem[-c(1:4)]/100,rep(NA,8))
prob.fem[(Time+Time.fore)] = NA

#phi.cov = rnorm(30,0,1)
#fec.cov = rnorm(31,0,1)
#phi.ind = c(rep(1,18), rep(2,3),rep(3,9))
#phi.ind = c(rep(1,5),rep(2,5),rep(3,5),rep(4,5),rep(5,5),rep(6,5))

n.fem = tx_nest_dat$nests[-c(1:4)]/apply(dat,2,sum)
n.fem[2:4] <-1
n.fem[(Time+Time.fore)] <- NA

#Nest success data

E = (tx_nest_dat$eggs_intact[-c(1:4,25:31)] + tx_nest_dat$eggs_broken[-c(1:4,25:31)])
J = tx_nest_dat$hatch_release[-c(1:4,25:31)]
#j.e = J/E
#j.e[2:3] = NA #beta data can't include 0
E.R = E/tx_nest_dat$nests[-c(1:4,25:31)]
E.R[2:3] <- 0.01

incu = nest.pro$incubation.type
#incu = incu[-which(incu=="S")]
incu[incu == "I" | incu == "I "] = 1
incu[incu == "C"] = 0
incu[incu == "S"] = 2
incu = as.numeric(incu)
incu=incu[-1472]

#nest condition index. G = good, others are "poor"
cond = nest.pro$condition.code
cond = cond[-1472]
cond[cond == "G"] = 1
cond[cond == "N"] = 0
cond[cond == "O"] = 0
cond[cond == "L"] = 0
cond[cond == ""] = 2
cond = as.numeric(cond)

nest.loc = nest.pro$General.Location
nest.loc = nest.loc[-1472]
nest.loc[nest.loc != "SPI"] = 0
nest.loc[nest.loc == "SPI"] = 1
nest.loc = as.numeric(nest.loc)

good.incu = rep(0,1597)
good.corr.pais = rep(0,1597)
good.corr.spi = rep(0,1597)
ngood.incu = rep(0,1597)
ngood.corr = rep(0,1597)
insitu = rep(0,1597)
good.incu[cond == 1 & incu == 1] = 1
good.corr.pais[cond == 1 & incu == 0 & nest.loc == 0] = 1
good.corr.spi[cond == 1 & incu == 0 & nest.loc == 1] = 1
ngood.incu[cond == 0 & incu == 1] = 1 
ngood.corr[cond == 0 & incu == 0] = 1
insitu[incu == 2] = 1

nest.mange = rep(0,1597)
nest.mange[good.incu==1] = 1
nest.mange[good.corr.pais==1] = 2
nest.mange[good.corr.spi==1] = 3
nest.mange[ngood.incu == 1] = 4
nest.mange[ngood.corr == 1] = 5
nest.mange[insitu == 1] = 6

J.ind = nest.pro$hatch.released
J.ind = J.ind[-1472]
E.ind = nest.pro$total.eggs
E.ind = E.ind[-1472]
E.ind[is.na(E.ind)] = round(mean(E.ind,na.rm=TRUE))

year = nest.pro$Year
year = year - 1990
year = year[-1472]

year.mean.idx = as.numeric()
for(t in 25:31){
  year.mean.idx[t] = max(which(year == t))
}

#nests.fore = 250*Time.fore #250 nests each year in future
#n.fore.fems = (225/250)*nests.fore

#year[(year.mean.idx[31]+1):(year.mean.idx[31]+nests.fore)] = rep(32:(31+Time.fore),each=250)

#for(t in 32:(31+Time.fore)){
#  year.mean.idx[t] = max(which(year == t))
#}

#inport climate data
#historic
load("~/KRST/Analysis/KRST climate histical/krst.climate.variables.Rdata")
#start with bio5
bioc.forecast = read.csv("~/KRST/Analysis/CMIP6_BIOC/CMIP6_BIOC/bioc_means_all.csv")

bio5 = subset(hist.bioc.means,bioc == "bio5")
bio5.out = rep(NA,length(nest.pro$Year)-1)
bio5.fore.585 = subset(bioc.forecast,bioc == "bio5" & ssp == 585)
bio5.tmp = data.frame(mean = seq(bio5.fore.585$mean[1],bio5.fore.585$mean[2],length.out=19),
                      year = seq(2022,2040))
bio5.scale = as.numeric(scale(c(bio5$mean,bio5.tmp$mean)))
bio5$mean.scale = bio5.scale[1:31]

#temp during may (most common nesting month)
#setwd("~/KRST/Analysis/KRST climate histical")

load("C:/Users/beross/OneDrive - DOI/Documents/KRST/Analysis/SLR Data/slr.krst.data.Rdata")
slr.cov = rep(NA,length(nest.pro$Year)-1)
slr.fore.med = subset(slrp.summary, Scenario == "0.5 - MED")
slr.tmp = data.frame(mean = seq(sl.summary$msl.mean[8],slr.fore.med$msl.mean[1],length.out=19),
                     year = seq(2022,2040))
slr.scale = as.numeric(scale(c(sl.summary$msl.mean[2:8],slr.tmp$mean)))
sl.summary$mean.scale = c(NA,slr.scale[1:7])

may.temp = subset(hist.temp.means,month == 5)
may.tmax = rep(NA,length(nest.pro$Year)-1)
for(i in 1:1597){
  may.tmax[i] = may.temp$tmax[which(may.temp$year==(year+1990)[i])]
  bio5.out[i] = bio5$mean.scale[which(bio5$year==(year+1990)[i])]
  slr.cov[i] = sl.summary$msl.mean[which(sl.summary==(year+1990)[i])]
  }
may.tmax.ann = as.numeric(scale(may.temp$tmax))
bio5.ann = bio5$mean.scale
slr.cov = as.numeric(scale(slr.cov))

#need to assign each 20 year block (21-40, 41-60, etc) to different climate
#bio5.fore.cov = c(rep(bio5.scale[32],20),rep(bio5.scale[33],20),rep(bio5.scale[34],10))
bio5.fore.cov = bio5.tmp$mean[1:Time.fore]
slr.cov.fore = slr.tmp$mean[1:Time.fore]
#100 nests in year T+1
#E.ind[(year.mean.idx[31]+1):(year.mean.idx[31]+nests.fore)] = round(mean(E.ind))
#J.ind[(year.mean.idx[31]+1):(year.mean.idx[31]+nests.fore)] = NA

#forecast nest management scenarios
#%50 in corrals, %50 in situ
nest.mange.fore = c(63,62,0,0,0,125)
#nest.mange[(year.mean.idx[31]+1):(year.mean.idx[31]+nests.fore)] = rep(mange.fore,Time.fore)

#forecast CMR data - 225 new individuals for 250 new nests
#y.fore = rbind(y,matrix(NA,n.fore.fems,31))
#y.fore = cbind(y.fore,matrix(NA,(814+n.fore.fems),Time.fore))

#first[815:12064] = rep(32:81,each=225)

#forecast count data
y.fem = apply(dat,2,sum)
y.fem[32:(Time+Time.fore)] = NA

#dat.fore = rbind(dat,matrix(0,100,31))
#dat.fore = cbind(dat.fore,c(rep(0,795),rep(1,100)))

#need to truncate N.b + N.i to be at least as 
# large as number of observed individuals. Just
# have to use the y.fem data twice, once as data
# once as truncation constraint
#trun = c(apply(dat,2,sum),NA)
trun = c(apply(dat,2,sum))

#### Data for NIMBLE ####

krst_data <- list(y = y, y.fem = y.fem,
                  psiNBB.mom=psiNBB.mom,
                  psiBNB.mom=psiBNB.mom,
                  p.mom=p.mom,
                  s.p.mom=s.p.mom,
                  in.s = in.s,
                  J = J, 
                  E.R = E.R,
                  E = E,
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
                  J.ind = J.ind,
                  E.ind = E.ind,
                  E.ind.fore = 96)

source("~/KRST/Analysis/make_fore_dat.R", echo=TRUE)
slr.scenario = "2.0 - MED"
nest.mange.fore = c(40,19,24,15,1,1) #low,incu; low,cor,PAIS; low,cor,SPI; high,incu; high,cor; insitu
bio5.585.2 = make.fore.dat(hist.bioc.means = hist.bioc.means, 
                           clim.var = "bio5", 
                           ssp = 585, 
                           slrp.summary=slrp.summary, 
                           sl.summary = sl.summary, 
                           slr.scenario = slr.scenario)

#future patrol effort
k.scale = as.numeric(scale(kilo.patrol))
kilo.patrol.fore = c(k.scale,rep(k.scale[31],Time.fore))

krst_con <- list(N = nrow(y), Time = ncol(dat), 
                 first=first,
                 Time.fore = 45,
                 nest.mange = nest.mange,
                 year = year,
                 year.mean.idx = year.mean.idx,
                 N.nest = length(J.ind),
                 kilo.patrol = kilo.patrol.fore,
                 may.tmax.ann = bio5.585.2$bio5.ann,
                 clim.cov = bio5.585.2$bio5.out,
                 clim.cov.fore = bio5.585.2$bio5.fore.cov,
                 nest.mange.fore = nest.mange.fore,
                 slr.cov=bio5.585.2$slr.cov,
                 slr.cov.fore = bio5.585.2$slr.cov.fore)

zinit <- dat
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

krst_ini = function() list(
  psiNBB= runif(1,.1,.5),
  psiBNB= runif(1,.5,.9),
  pB=runif(1,0.1,0.9), 
  z=zinit,
  s.p1 = runif(1,.1,.5), 
  N.nb = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
  #N.b = rpois((Time+Time.fore),c(trun[1:Time]/2,rep(50,Time.fore))),
  N.b = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
  N.p1 = rpois((Time+Time.fore),seq(10,200,length.out=(Time+Time.fore))),
  N.p2 = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
  N.p3 = rpois((Time+Time.fore),seq(2,500,length.out=(Time+Time.fore))),
  N.p4 = rpois((Time+Time.fore),seq(2,500,length.out=(Time+Time.fore))),
  N.p5 = rpois((Time+Time.fore),seq(2,500,length.out=(Time+Time.fore))),
  N.p6 = rpois((Time+Time.fore),seq(2,500,length.out=(Time+Time.fore))),
  N.p7 = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
  N.p8 = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
  N.p9 = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
  N.p10 = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
  #N.i = rpois((Time+Time.fore),c(trun[1:Time]/2,rep(50,Time.fore))),
  N.i = rpois((Time+Time.fore),seq(2,200,length.out=(Time+Time.fore))),
  r = runif(1,.8,10),
  ln.mean.fem = runif(1,-1,1),
  mu.n.fem=2,
  sd.n.fem = runif(1,.1,1),
  mu.phi = runif(1,1,2.2),
  sd.phi2=runif(1,.1,1),
  mu.egg = runif(1,-1,1),
  nest.cov = c(runif(5,-1,1),NA),
  insitu.surv = runif(1,.1,.8),
  mu.temp = runif(1,-1,1),
  var.temp = runif(1,.1,3),
  beta.temp.all = runif(1,-1,1),
  beta.temp = runif(6,-1,1),
  mu.imm = runif(1,-1,1),
  sd.imm = runif(1,.01,1),
  ln.mu.phi = runif((Time+Time.fore-1),1,2.2),
  eps.imm.tmp = runif((Time+Time.fore),-1,1),
  J.ind = round(J.ini),
  prob.fem = prob.fem.ini,
  n.fem = c(rep(NA,Time),rep(2,Time.fore)),
  y.fem = c(rep(NA,Time),rep(150,Time.fore)),
  beta.slr = c(rep(NA,5),runif(1,-1,1)),
  a.0 = runif(1,-1,1),
  a.k = runif(1,-1,1),
  f = rep(round(f_ini),(Time+Time.fore)))

krst.model <- nimbleModel(code = krst.ipm,
                           constants = krst_con,
                           data = krst_data,
                           inits = krst_ini())

krst.config.mcmc <- configureMCMC(krst.model)

krst.config.mcmc$addMonitors("f")

krst.mcmc <- buildMCMC(krst.config.mcmc)

krst.cmodel <- compileNimble(krst.model)

krst.cmcmc <- compileNimble(krst.mcmc, project = krst.model)

krst.samples <- runMCMC(krst.cmcmc, niter= 2000, 
                         nburnin = 1000,nchains = 3,
                         samplesAsCodaMCMC = TRUE)

library(coda)
library(parallel)

nc = 3 #num of chains
clus = makeCluster(nc)

clusterExport(clus, c("krst.ipm","krst_ini",
                      "krst_data",
                      "krst_con"))

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
                       constants = krst_con,
                       data = krst_data, 
                       inits = init)
  Cmodel <- compileNimble(model)
  modelConf <- configureMCMC(model)
  #modelConf$addMonitors(params)
  modelConf$addMonitors(c("phi","N.i","N.b","beta.temp.all","f","prob.fem"))
  #configureRJ(modelConf,
  #            targetNodes = 'beta',
  #            indicatorNodes = 'indA',
  #            control = list(mean=0, scale = 2))
  modelMCMC <- buildMCMC(modelConf)
  CmodelMCMC <- compileNimble(modelMCMC, 
                              project = model)
  #out1 <- runMCMC(CmodelMCMC, niter=15000, nburnin=14000)
  #return(as.mcmc(out1))
  CmodelMCMC$run(10000,reset=TRUE, resetMV=TRUE)
  return(as.mcmc(as.matrix(CmodelMCMC$mvSamples)))
})

out.mcmc <- as.mcmc.list(out) 

out2 <- clusterEvalQ(clus,{
  CmodelMCMC$run(100000,thin=10,reset=FALSE, resetMV=TRUE)
  return(as.mcmc(as.matrix(CmodelMCMC$mvSamples)))
})

out.up <- as.mcmc.list(out2)

save(out.up, file="KRST_out.RData")
