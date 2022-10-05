library(viridis)
library(coda)
library(vioplot)
library(gplots)

load("~/KRST/Analysis/results_currentnm_5852.RData")
load("~/KRST/Analysis/results_currentnm_2455_2100.RData")
load("~/KRST/Analysis/results_lowis_2455.RData")
load("~/KRST/Analysis/results_noincu_2455.RData")
load("~/KRST/Analysis/results_currentnm_2455_is1.RData")
load("~/KRST/Analysis/results_currentnm_2455_is0.RData")

out.2455 = as.mcmc(rbind(out.current.2455[[1]], out.current.2455[[2]], out.current.2455[[3]]))
out.5852 = as.mcmc(rbind(out.current.5852[[1]], out.current.5852[[2]], out.current.5852[[3]]))

out.is0 = as.mcmc(rbind(out.current.2455.is0[[1]], out.current.2455.is0[[2]], out.current.2455.is0[[3]]))
out.is1 =as.mcmc(rbind(out.current.2455.is1[[1]], out.current.2455.is1[[2]], out.current.2455.is1[[3]]))

out.lis =as.mcmc(rbind(out.lowis.2455[[1]], out.lowis.2455[[2]], out.lowis.2455[[3]]))
out.noincu =as.mcmc(rbind(out.noincu.2455[[1]], out.noincu.2455[[2]], out.noincu.2455[[3]]))

col.plot=viridis(5,alpha=.6)

Time = 28
Time.fore = 79
plot.year = 1993+Time+Time.fore

## Figure for abundance
Nb.2455 = apply(out.2455[,1:(Time+Time.fore)],2,quantile,.5)
Nb.low.2455 = apply(out.2455[,1:(Time+Time.fore)],2,quantile,0.025)
Nb.high.2455 = apply(out.2455[,1:(Time+Time.fore)],2,quantile,0.975)
Nb.5852 = apply(out.5852[,1:(Time+Time.fore)],2,quantile,.5)
Nb.low.5852 = apply(out.5852[,1:(Time+Time.fore)],2,quantile,0.025)
Nb.high.5852 = apply(out.5852[,1:(Time+Time.fore)],2,quantile,0.975)

imm = apply(out.2455[,(Time+Time.fore+1):(2*(Time+Time.fore))],2,quantile,.5)
imm.l = apply(out.2455[,(Time+Time.fore+1):(2*(Time+Time.fore))],2,quantile,0.025)
imm.h = apply(out.2455[,(Time+Time.fore+1):(2*(Time+Time.fore))],2,quantile,0.975)
Nnb = apply(out[,(2*(Time+Time.fore)+1):(3*(Time+Time.fore))],2,mean)
Nnb.low = apply(out[,(2*(Time+Time.fore)+1):(3*(Time+Time.fore))],2,quantile,0.025)
Nnb.high = apply(out[,(2*(Time+Time.fore)+1):(3*(Time+Time.fore))],2,quantile,0.975)
#Np1 = apply(out[,63:93],2,mean)
#Np2 = apply(out[,94:124],2,mean)
#Np3 = apply(out[,125:155],2,mean)
#Np4 = apply(out[,156:186],2,mean)
#Np5 = apply(out[,187:217],2,mean)
#Np6 = apply(out[,218:248],2,mean)
#Np7 = apply(out[,249:279],2,mean)

#Np = rowSums(cbind(Np1,Np2,Np3,Np4,Np5,Np6,Np7))
Nad = rowSums(cbind(Nb,Nnb))

Nad.low <- Nad.high <- rep(NA,(Time+Time.fore))
for(i in 1:(Time+Time.fore)){
Nad.low[i]=quantile(rowSums(cbind(out[,i],out[,i+(2*(Time+Time.fore))])),0.025)
Nad.high[i] =quantile(rowSums(cbind(out[,i],out[,i+(2*(Time+Time.fore))])),0.975)}

par(mfrow=c(2,1))

plot(Nad,x=seq(1994,plot.year), type="l", ylab = "Abundance", 
     xlab= "Year",col=col.plot[1],ylim=c(0,max(Nad.high)))
polygon(c(seq(1994,plot.year),rev(seq(1994,plot.year))),
     c(Nad.low,rev(Nad.high)),col=col.plot[1])
#mtext("A) Total Abundance", side=3, line=.5, adj=0)
#axis(side=1, at=seq(1991,2021), labels=FALSE)

plot(Nb, x = seq(1994,plot.year), type="l", ylab="Abundance",
     xlab="Year",ylim=c(0,max(Nb.high)))
polygon(c(seq(1994,plot.year), rev(seq(1994,plot.year))),
        c(Nb.low,rev(Nb.high)),col=col.plot[2])


lines(Nnb, x= seq(1991,plot.year),type="l")
polygon(c(seq(1991,plot.year),rev(seq(1991,plot.year))),
        c(Nnb.low,rev(Nnb.high)),col=col.plot[3])

#mtext("B) Breeders and Non-Breeders", side=3, line=.5, adj=0)


legend("topleft", 
       legend=c("Breeders","Non-Breeders"), pch=c(16,16),
       bty="n",col=c(col.plot[2],col.plot[3]),
       x.intersp=0.7, y.intersp=1.0, cex=1.1, pt.cex=1.5)

#dev.print(tiff,"tx_ad.tiff",res=600,width=9,units="in")

N.t = rowSums(cbind(Nb.2455,imm))

Nt.low <- Nt.high <- rep(NA,(Time+Time.fore))
for(i in 1:(Time+Time.fore)){
  Nt.low[i]=quantile(rowSums(cbind(out.2455[,i],out.2455[,i+(Time+Time.fore)])),0.025)
  Nt.high[i] =quantile(rowSums(cbind(out.2455[,i],out.2455[,i+(Time+Time.fore)])),0.975)}

par(
  mfrow=c(1,2), # panels will plot in 1 row with 3 columns
  oma=c(3,3.25,3,3), # outer margins (bottom, left, top, right)
  mai=c(0.5,0.1,0.1,0.1), # inner margins (between panels)
  cex=1, 	# make sure the text and points are a normal size even if your plot window is an unusual size (R likes to resize spontaneously)
  mgp=c(3,0.8,0))

plot(N.t,x=seq(1994,plot.year), type="l", 
     xlab= "Year",col="black",
     ylab = "Abundance",
     ylim=c(0,max(Nt.high)))
polygon(c(seq(1994,plot.year),rev(seq(1994,plot.year))),
        c(Nt.low,rev(Nt.high)),col=col.plot[4])


#axis(side=1, at=seq(1991,2021), labels=FALSE)


plot(Nb, x= seq(1994,plot.year),type="l",
     ylim=c(0,max(Nb.high)),
     xlab = "Year", ylab = "Abundance")
polygon(c(seq(1994,plot.year),rev(seq(1994,plot.year))),
        c(Nb.low,rev(Nb.high)),col=col.plot[3])
lines(imm, x = seq(1994,plot.year), type="l", ylab="Abundance",
     xlab="Year",ylim=c(0,max(imm.h)))
polygon(c(seq(1994,plot.year), rev(seq(1994,plot.year))),
        c(imm.l,rev(imm.h)),col=col.plot[5])
#lines(Nb, x = seq(1991,2021), type="l", ylab="Abundance",
#     xlab="Year",ylim=c(0,max(Nnb.high)))
#polygon(c(seq(1991,2021), rev(seq(1991,2021))),
#        c(Nb.low,rev(Nb.high)),col=col.plot[2])

#mtext("C) Imm, Breeders, Non-Breeders", side=3, line=.5, adj=0)

legend("topleft", 
       legend=c("Non-TX Breeders", "TX Breeders"), pch=c(16,16),
       bty="n",col=c(col.plot[5],col.plot[3]),
       x.intersp=0.7, y.intersp=1.0, cex=1.1, pt.cex=1.5)

Nb.is0 = apply(out.is0[,1:(Time+Time.fore)],2,quantile,.5)
Nb.low.is0 = apply(out.is0[,1:(Time+Time.fore)],2,quantile,0.025)
Nb.high.is0 = apply(out.is0[,1:(Time+Time.fore)],2,quantile,0.975)
Nb.is1 = apply(out.is1[,1:(Time+Time.fore)],2,quantile,.5)
Nb.low.is1 = apply(out.is1[,1:(Time+Time.fore)],2,quantile,0.025)
Nb.high.is1 = apply(out.is1[,1:(Time+Time.fore)],2,quantile,0.975)
Nb.lis = apply(out.lis[,1:(Time+Time.fore)],2,quantile,.5)
Nb.low.lis = apply(out.lis[,1:(Time+Time.fore)],2,quantile,0.025)
Nb.high.lis = apply(out.lis[,1:(Time+Time.fore)],2,quantile,0.975)
Nb.noincu = apply(out.noincu[,1:(Time+Time.fore)],2,quantile,.5)
Nb.low.noincu = apply(out.noincu[,1:(Time+Time.fore)],2,quantile,0.025)
Nb.high.noincu = apply(out.noincu[,1:(Time+Time.fore)],2,quantile,0.975)


plot(Nb.2455, x= seq(1994,plot.year),type="l",
     ylim=c(0,max(Nb.high.2455)),
     xlab = "Year", ylab = "Abundance")
polygon(c(seq(1994,plot.year),rev(seq(1994,plot.year))),
        c(Nb.low.2455,rev(Nb.high.2455)),col=col.plot[3])
lines(Nb.5852, x = seq(1994,plot.year), lty=3, ylab="Abundance",
      xlab="Year",ylim=c(0,max(Nb.high.2455)),lwd=2)
polygon(c(seq(1994,plot.year), rev(seq(1994,plot.year))),
        c(Nb.low.5852,rev(Nb.high.5852)),col=col.plot[5])

legend("topleft", 
       legend=c("SSP2-4.5, 0.5 m SLR", "SSP5-8.5, 2 m SLR"), pch=c(16,16),
       bty="n",col=c(col.plot[3],col.plot[5]),
       x.intersp=0.7, y.intersp=1.0, cex=1.1, pt.cex=1.5)

#all in sito, effect = 0
lines(Nb.is0, x = seq(1994,plot.year), lty=3,lwd=2)
polygon(c(seq(1994,plot.year), rev(seq(1994,plot.year))),
        c(Nb.low.is0,rev(Nb.high.is0)),col=col.plot[5])

legend("topleft", 
       legend=c("Current Management", "All In Situ PAIS"), pch=c(16,16),
       bty="n",col = c(col.plot[3],col.plot[5]),
       x.intersp=0.7, y.intersp=1.0, cex=1.1, pt.cex=1.5)

#all in sito, effect = -1
lines(Nb.is1, x = seq(1994,plot.year), lty=3,lwd=2)
polygon(c(seq(1994,plot.year), rev(seq(1994,plot.year))),
        c(Nb.low.is1,rev(Nb.high.is1)),col=col.plot[5])

legend("topleft", 
       legend=c("Current Management", "All In Situ PAIS"), pch=c(16,16),
       bty="n",col = c(col.plot[3],col.plot[5]),
       x.intersp=0.7, y.intersp=1.0, cex=1.1, pt.cex=1.5)

#PAIS low in situ
lines(Nb.lis, x = seq(1994,plot.year), lty=3, ylab="Abundance",
      xlab="Year",ylim=c(0,max(Nb.2455)),lwd=2)
polygon(c(seq(1994,plot.year), rev(seq(1994,plot.year))),
        c(Nb.low.lis,rev(Nb.high.lis)),col=col.plot[5])

legend("topleft", 
       legend=c("Current Management", "8% In Situ PAIS"), pch=c(16,16),
       bty="n",col = c(col.plot[3],col.plot[5]),
       x.intersp=0.7, y.intersp=1.0, cex=1.1, pt.cex=1.5)

#no incubation at PAIS
lines(Nb.noincu, x = seq(1994,plot.year), lty=3, ylab="Abundance",
      xlab="Year",ylim=c(0,max(Nb.2455)),lwd=2)
polygon(c(seq(1994,plot.year), rev(seq(1994,plot.year))),
        c(Nb.low.noincu,rev(Nb.high.noincu)),col=col.plot[5])

legend("topleft", 
       legend=c("Current Management", "No Incubation"), pch=c(16,16),
       bty="n",col = c(col.plot[3],col.plot[5]),
       x.intersp=0.7, y.intersp=1.0, cex=1.1, pt.cex=1.5)


#### LAMBDA ####
lam = matrix(NA,30000,(Time+Time.fore-1))
n.tot = matrix(NA,30000,(Time+Time.fore))
lam.b.2455 = matrix(NA,30000,(Time+Time.fore-1))
lam.b.5852 = matrix(NA,30000,(Time+Time.fore-1))
lam.b.is0 = matrix(NA,30000, (Time+Time.fore-1))
lam.b.is1 = matrix(NA,30000, (Time+Time.fore-1))
lam.b.lis = matrix(NA,30000, (Time+Time.fore-1))
lam.b.noin = matrix(NA,30000, (Time+Time.fore-1))

#for(t in 1:(Time+Time.fore)){
#  n.tot[,t] = rowSums(cbind(out[,t], out[,(2*(Time+Time.fore))+t]))
#}

for(t in 1:(Time+Time.fore-1)){ 
  #lam[,t] = n.tot[,t+1]/n.tot[,t]
  lam.b.2455[,t] = out.2455[,t+1]/out.2455[,t]
  lam.b.5852[,t] = out.5852[,t+1]/out.5852[,t]
  #lam.b.is0[,t] = out.is0[,t+1]/out.is0[,t]
  #lam.b.is1[,t] = out.is1[,t+1]/out.is1[,t]
  #lam.b.lis[,t] = out.lis[,t+1]/out.lis[,t]
  lam.b.noin[,t] = out.noincu[,t+1]/out.noincu[,t]
  }
#remove Inf from lam
lam.b.2455[!is.finite(lam.b.2455)] = NA
lam.b.5852[!is.finite(lam.b.5852)] = NA
lam.b.is0[!is.finite(lam.b.is0)] = NA
lam.b.is1[!is.finite(lam.b.is1)]= NA
lam.b.lis[!is.finite(lam.b.lis)] = NA
lam.b.noin[!is.finite(lam.b.noin)] = NA

#geo mean of lambda through future
geo_mean <- function(data){
  gm = rep(NA,30000)
  for(i in 1:length(data[,1])){
  gm[i] <- exp(mean(log(data[i,])))
  }
  return(gm)
  }

gm.lam.2455 = geo_mean(lam.b.2455[,29:72])
gm.lam.5852 = geo_mean(lam.b.5852[,29:72])
gm.lam.is1 = geo_mean(lam.b.is1[,29:72])
gm.lam.is0 = geo_mean(lam.b.is0[,29:72])
gm.lam.lis = geo_mean(lam.b.lis[,29:72])
gm.lam.noin = geo_mean(lam.b.noin[,29:72])

boxplot(cbind(gm.lam.2455,gm.lam.5852,gm.lam.is1,
              gm.lam.is0,gm.lam.lis,gm.lam.noin),
        ylab = expression(lambda),
        names = c("Current,SSP2","Current,SSP5",
                  "All In Situ -1","All In Situ 0",
                  "Limited In Situ","No Incubation"),outline=FALSE)
abline(h=1)

zero.gm.2455 = length(which(gm.lam.2455==0))
zero.gm.5852 = length(which(gm.lam.5852==0))
zero.gm.is1 = length(which(gm.lam.is1==0))
zero.gm.is0 = which(gm.lam.is0==0)
zero.gm.lis = which(gm.lam.lis==0)
zero.gm.noin = which(gm.lam.noin==0)

lam.low = apply(lam[,1:(Time-1)],2,quantile,0.025,na.rm=TRUE)
lam.high = apply(lam[,1:(Time-1)],2,quantile,0.975,na.rm=TRUE)
lam.low.fore.2455 = apply(lam.b.2455[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.5852 = apply(lam.b.5852[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.is0 = apply(lam.b.is0[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.is1 = apply(lam.b.is1[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.lis = apply(lam.b.lis[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.noin = apply(lam.b.noin[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)

lam.high.fore.2455 = apply(lam.b.2455[,(Time):(Time+Time.fore-1)],2,quantile,0.975)
lam.high.fore.5852 = apply(lam.b.5852[,(Time):(Time+Time.fore-1)],2,quantile,0.975)
lam.high.fore.is0 = apply(lam.b.is0[,(Time):(Time+Time.fore-1)],2,quantile,0.975,na.rm=TRUE)
lam.high.fore.is1 = apply(lam.b.is1[,(Time):(Time+Time.fore-1)],2,quantile,0.975,na.rm=TRUE)
lam.high.fore.lis = apply(lam.b.lis[,(Time):(Time+Time.fore-1)],2,quantile,0.975,na.rm=TRUE)
lam.high.fore.noin = apply(lam.b.noin[,(Time):(Time+Time.fore-1)],2,quantile,0.975,na.rm=TRUE)

lam.med.all.2455 = apply(lam.b.2455[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.5852 = apply(lam.b.5852[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.is0 = apply(lam.b.is0[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.is1 = apply(lam.b.is1[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.lis = apply(lam.b.lis[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.noin = apply(lam.b.noin[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)

par(mfrow=c(1,1))
#vioplot(lam[,-1],ylim=c(0,3),names=seq(1993,(1991+Time+Time.fore)))
boxplot(cbind(lam.b.2455[,1:(Time-1)],matrix(NA,30000,(Time.fore-1))),
        ylab = expression(lambda),
        names = seq(1995,(1993+Time+Time.fore-1)),outline=FALSE)

#plot(y=lam.med.all, x=seq(1992,(1991+Time+Time.fore-1)),
#     type="l",ylim=c(min(lam.low),max(lam.high)),
#     lty=2,lwd=3,ylab=expression(lambda),
#     xlab="Year")
#polygon(c(seq(1995,(1994+Time-1)), 
#          rev(seq(1995,(1994+Time-1)))),
#        c(lam.low,rev(lam.high)),
#        col=col.plot[1])

polygon(c(seq(Time,(Time+Time.fore-1)), 
          rev(seq(Time,(Time+Time.fore-1)))),
        c(lam.low.fore.2455,rev(lam.high.fore.2455)),
          col=col.plot[3])
polygon(c(seq(Time,(Time+Time.fore-1)), 
          rev(seq(Time,(Time+Time.fore-1)))),
        c(lam.low.fore.5852,rev(lam.high.fore.5852)),
        col=col.plot[5])
#in situ effects = 0
polygon(c(seq(Time,(Time+Time.fore-1)), 
          rev(seq(Time,(Time+Time.fore-1)))),
        c(lam.low.fore.is0,rev(lam.high.fore.is0)),
        col=col.plot[2])
#in situ effects = -1
polygon(c(seq(Time,(Time+Time.fore-1)), 
          rev(seq(Time,(Time+Time.fore-1)))),
        c(lam.low.fore.is1,rev(lam.high.fore.is1)),
        col=col.plot[5])
#low in situ, 8% at PAIS
polygon(c(seq(Time,(Time+Time.fore-1)), 
          rev(seq(Time,(Time+Time.fore-1)))),
        c(lam.low.fore.lis,rev(lam.high.fore.lis)),
        col=col.plot[5])
#no incubation at PAIS
polygon(c(seq(Time,(Time+Time.fore-1)), 
          rev(seq(Time,(Time+Time.fore-1)))),
        c(lam.low.fore.noin,rev(lam.high.fore.noin)),
        col=col.plot[5])

lines(y=lam.med.all.5852[Time:(Time+Time.fore-1)],
      x=seq(Time,(Time+Time.fore-1)),lty=3,lwd=2)
lines(y=lam.med.all.2455[Time:(Time+Time.fore-1)],
      x=seq(Time,(Time+Time.fore-1)),lty=1,lwd=2)
lines(y=lam.med.all.is0[Time:(Time+Time.fore-1)],
      x=seq(Time,(Time+Time.fore-1)),lty=1,lwd=2)
lines(y=lam.med.all.is1[Time:(Time+Time.fore-1)],
      x=seq(Time,(Time+Time.fore-1)),lty=2,lwd=2)
lines(y=lam.med.all.lis[Time:(Time+Time.fore-1)],
      x=seq(Time,(Time+Time.fore-1)),lty=2,lwd=2)
lines(y=lam.med.all.noin[Time:(Time+Time.fore-1)],
      x=seq(Time,(Time+Time.fore-1)),lty=2,lwd=2)


abline(h=1)

legend("topright", 
       legend=c("SSP2-4.5, 0.5 m", "SSP5-8.5, 2 m"), pch=c(16,16),
       bty="n",col=c(col.plot[3],col.plot[5]),
       x.intersp=0.7, y.intersp=1.0, cex=1.1, pt.cex=1.5)

#Detection probability
det.prob = out[,665]
plot(density(det.prob))

boxplot(det.prob,outline=FALSE,
        names="",
        ylab = "Detection Probability",
        col="white",ylim = c(.5,1))

#transition probabilities
boxplot(out[,711],out[,712],
        names=c("Breeder to Non","Non to Breeder"),
        outline=FALSE)

#survival rates through time
#666-710
phi.plot = as.matrix(out[,666:710])
boxplot(phi.plot,names=seq(1992,(1990+Time+Time.fore)),
        ylab = "Adult Survival",outline=FALSE)
abline(h=(1/(1+exp(-(mean(out[,657]))))))

##### covs plots ####
#climate effects
clim.cov.plot = out[,311]

#yrs 1-24
seq.plot = seq(min(c(bio5.245.5$bio5.ann,bio5.245.5$bio5.fore.cov)),
               max(c(bio5.245.5$bio5.ann,bio5.245.5$bio5.fore.cov)),
               length.out=100)

#demonstration of effect of -1 on prob of nest success
demo.plot = 1/(1+exp(-(mean(out.is1[,973])+seq.plot*-1)))
clim = seq(min(subset(hist.bioc.means,bioc == "bio5")$mean),39,length.out=100)
plot(y = demo.plot, x = clim,
     type= "l", ylab = "Prob. of Egg Success", xlab = "Max Temp")

plot.low.incu <- plot.low.cor.pais <-
  plot.low.cor.spi <- plot.high.incu <-
  plot.high.cor <- plot.insitu <- matrix(NA,30000,100)

for(j in 1:100){
  #bioclim 5
  #order: "Low,Incu","low,Cor,PAIS", "Low,Cor,SPI","High,Incu","High,Cor","In situ")
  plot.low.incu[,j] = 1/(1+exp(-(out[,659]+ 
                                   out[,599]*seq.plot[j])))
  
  plot.low.cor.pais[,j] = 1/(1+exp(-(out[,660]+ 
                                   out[,600]*seq.plot[j])))

  plot.low.cor.spi[,j] = 1/(1+exp(-(out[,661]+ 
                                   out[,601]*seq.plot[j])))
  
  plot.high.incu[,j] = 1/(1+exp(-(out[,662]+ 
                                   out[,602]*seq.plot[j])))
  
  plot.high.cor[,j] = 1/(1+exp(-(out[,663]+ 
                                   out[,603]*seq.plot[j])))
  
  plot.insitu[,j] = 1/(1+exp(-(out[,664]+ 
                                   out[,604]*seq.plot[j])))
  }

plot.low.incu2.low = apply(plot.low.incu, 2, function(x) quantile(x, probs = c(0.025)))
plot.low.incu2 = apply(plot.low.incu, 2, function(x) quantile(x, probs = c(0.5)))
plot.low.incu2.high = apply(plot.low.incu, 2, function(x) quantile(x, probs = c(0.975)))

plot.low.cor.pais2.low = apply(plot.low.cor.pais, 2, function(x) quantile(x, probs = c(0.025)))
plot.low.cor.pais2 = apply(plot.low.cor.pais, 2, function(x) quantile(x, probs = c(0.5)))
plot.low.cor.pais2.high = apply(plot.low.cor.pais, 2, function(x) quantile(x, probs = c(0.975)))

plot.low.cor.spi2.low = apply(plot.low.cor.spi,2,function(x) quantile(x,probs = 0.025))
plot.low.cor.spi2 = apply(plot.low.cor.spi,2,function(x) quantile(x,probs = 0.5))
plot.low.cor.spi2.high = apply(plot.low.cor.spi,2,function(x) quantile(x,probs = 0.975))

plot.high.incu2.low = apply(plot.high.incu,2,function(x) quantile(x,probs = 0.025))
plot.high.incu2 = apply(plot.high.incu,2,function(x) quantile(x,probs = 0.5))
plot.high.incu2.high = apply(plot.high.incu,2,function(x) quantile(x,probs = 0.975))

plot.high.cor2.low = apply(plot.high.cor,2,function(x) quantile(x,probs = 0.025))
plot.high.cor2 = apply(plot.high.cor,2,function(x) quantile(x,probs = 0.5))
plot.high.cor2.high = apply(plot.high.cor,2,function(x) quantile(x,probs = 0.975))

plot.insitu2.low = apply(plot.insitu,2,function(x) quantile(x,probs = 0.025))
plot.insitu2 = apply(plot.insitu,2,function(x) quantile(x,probs = 0.5))
plot.insitu2.high = apply(plot.insitu,2,function(x) quantile(x,probs = 0.975))

cov.mean = mean(c(bio5$mean,bio5.fore.245$mean[1]))
cov.sd = sd(c(bio5$mean,bio5.fore.245$mean[1]))
plot.tmp = (seq.plot*cov.sd)+cov.mean 
plot.x = seq(min(plot.tmp),max(plot.tmp),length.out = 100)

col.plot = viridis(6)

par(
  mfrow=c(2,3), # panels will plot in 1 row with 3 columns
  oma=c(3,3.25,3,3), # outer margins (bottom, left, top, right)
  mai=c(0.5,0.1,0.1,0.1), # inner margins (between panels)
  cex=1, 	# make sure the text and points are a normal size even if your plot window is an unusual size (R likes to resize spontaneously)
  mgp=c(3,0.8,0))

plot(plot.low.incu2,type="l",x=plot.x,cex=1.2,
     ylim=c(0,1),xaxt="n",yaxt="n",col=col.plot[1],
     lwd=2)
lines(plot.low.incu2.low,lty=2,x=plot.x,col=col.plot[1],lwd=2)
lines(plot.low.incu2.high,lty=2,x=plot.x,col=col.plot[1],lwd=2)
mtext("A) Low Risk, Incubator", side=3, line=.5, adj=0)
axis(side=2, labels=TRUE, las=1)
axis(side=1, labels=FALSE, las=1)

plot(plot.low.cor.pais2,type="l",x=plot.x,col=col.plot[2],lwd=2,ylim=c(0,1),
     xaxt="n",yaxt="n")
lines(plot.low.cor.pais2.low,lty=2,x=plot.x,col=col.plot[2],lwd=2)
lines(plot.low.cor.pais2.high,lty=2,x=plot.x,col=col.plot[2],lwd=2)
mtext("B) Low Risk, Corral, PAIS", side=3, line=.5, adj=0)
axis(side=2, labels=FALSE, las=1)
axis(side=1, labels=FALSE, las=1)

plot(plot.high.incu2,type="l",x=plot.x,
     col=col.plot[4],lwd=2,ylim=c(0,1),
     xaxt="n",yaxt="n")
lines(plot.high.incu2.low,lty=2,x=plot.x,col=col.plot[4],lwd=2)
lines(plot.high.incu2.high,lty=2,x=plot.x,col=col.plot[4],lwd=2)
mtext("D) High Risk, Incubator, PAIS and North", side=3, line=.5, adj=0)
axis(side=1, at=NULL, labels=FALSE)
axis(side=2, labels=FALSE, las=1)

plot(plot.low.cor.spi2,type="l",x=plot.x,col=col.plot[3],lwd=2,ylim=c(0,1),
     xaxt="n",yaxt="n")
lines(plot.low.cor.spi2.high,lty=2,x=plot.x,col=col.plot[3],lwd=2)
lines(plot.low.cor.spi2.low,lty=2,x=plot.x,col=col.plot[3],lwd=2)
mtext("C) Low Risk, Corral, SPI", side=3, line=.5, adj=0)
axis(side=1, labels=TRUE, las=1)
axis(side=2, labels=TRUE, las=1)

plot(plot.high.cor2,type="l",x=plot.x,
     col=col.plot[5],lwd=2,ylim=c(0,1),
     xaxt="n",yaxt="n")
lines(plot.high.cor2.low,lty=2,x=plot.x,col=col.plot[5],lwd=2)
lines(plot.high.cor2.high,lty=2,x=plot.x,col=col.plot[5],lwd=2)
mtext("E) High Risk, Corral, SPI", side=3, line=.5, adj=0)
axis(side=1, at=NULL, labels=TRUE)
axis(side=2, labels=FALSE, las=1)

plot(plot.insitu2,type="l",x=plot.x,
     col=col.plot[6],lwd=2,ylim=c(0,1),
     xaxt="n",yaxt="n")
lines(plot.insitu2.low,lty=2,x=plot.x,col=col.plot[6],lwd=2)
lines(plot.insitu2.high,lty=2,x=plot.x,col=col.plot[6],lwd=2)
mtext("F) In situ", side=3, line=.5, adj=0)
axis(side=1, at=NULL, labels=TRUE)
axis(side=2, labels=FALSE, las=1)

mtext("Max Temp (C) of Warmest Month", 
      side=1, line=1, at=0.5,
      outer=TRUE)
mtext("Egg Success", side=2, line=1.5, 
      at=0.5,outer=TRUE,las=0)


#nest management effects
#659-664 - not correct names below!
good.corr.pais.plot = 1/(1+exp(-(out[,324])))
good.corr.spi.plot = 1/(1+exp(-(out[,325])))
good.incu.plot = 1/(1+exp(-(out[,323])))
ngood.corr.plot = 1/(1+exp(-(out[,327])))
ngood.incu.plot = 1/(1+exp(-(out[,326])))
insitu.plot = 1/(1+exp(-(out[,328])))

par(mfrow = c(1,1))

boxplot(good.corr.pais.plot,good.corr.spi.plot,
        good.incu.plot,ngood.corr.plot,
        ngood.incu.plot,insitu.plot,
        names = c("Low,Cor,PAIS", "Low,Cor,SPI",
                  "Low,Incu",
                  "High,Incu","High,Cor","In situ"),
        ylab = "Survival Probability")

points(x=5,y=0.63,col="red",pch=16)

#Plots of interactions

#betas for each management*may tmax
betas = cbind(out[,311:316])
vioplot(betas, names = c("LR,PAIS","LR,SPI",
                               "LR,Incu","HR,Incu",
                               "HR,Cor","In situ"))
abline(h = 0)

#good spi and good pais w/ sea level rise
seq.plot = seq(min(slr[25:31]),
                     max(slr[25:31]),
                     length.out=100)

plot.good.corr.spi <- plot.good.corr.pais <- 
  plot.good.incu <- plot.ngood.corr <- 
  plot.ngood.incu <- matrix(NA,nrow=150000,ncol=100)

for(j in 1:100){
  #low precip
  plot.good.corr.spi[,j]=1/(1+(exp(-(out[,351] + 
                  (out[,394]*seq.plot[j]) +
                  (out[,345]) + 
                  (out[,346]*seq.plot[j])))))
  
  plot.good.corr.pais[,j]=1/(1+(exp(-(out[,351] + 
                  (out[,394]*seq.plot[j]) +
                  (out[,343]) + 
                  (out[,344]*seq.plot[j])))))
  
  plot.good.incu[,j] = 1/(1+(exp(-(out[,351] + 
                  (out[,394]*seq.plot[j]) +
                  (out[,347]) + 
                  (out[,348]*seq.plot[j])))))
  
  plot.ngood.corr[,j] = 1/(1+(exp(-(out[,351] + 
                  (out[,394]*seq.plot[j]) +
                  (out[,354]) + 
                  (out[,355]*seq.plot[j])))))
  
  plot.ngood.incu[,j] = 1/(1+(exp(-(out[,351] + 
                  (out[,394]*seq.plot[j])))))
  }

plot.good.spi2.low = apply(plot.good.corr.spi, 2, function(x) quantile(x, probs = c(0.025)))
plot.good.spi2 = apply(plot.good.corr.spi, 2, function(x) quantile(x, probs = c(0.5)))
plot.good.spi2.high = apply(plot.good.corr.spi, 2, function(x) quantile(x, probs = c(0.975)))

plot.pais2.low = apply(plot.good.corr.pais, 2, function(x) quantile(x, probs = c(0.025)))
plot.pais2 = apply(plot.good.corr.pais, 2, function(x) quantile(x, probs = c(0.5)))
plot.pais2.high = apply(plot.good.corr.pais, 2, function(x) quantile(x, probs = c(0.975)))

plot.good.incu2.low = apply(plot.good.incu,2,function(x) quantile(x,probs = 0.025))
plot.good.incu2 = apply(plot.good.incu,2,function(x) quantile(x,probs = 0.5))
plot.good.incu2.high = apply(plot.good.incu,2,function(x) quantile(x,probs = 0.975))

plot.ngood.corr2.low = apply(plot.ngood.corr,2,function(x) quantile(x,probs = 0.025))
plot.ngood.corr2 = apply(plot.ngood.corr,2,function(x) quantile(x,probs = 0.5))
plot.ngood.corr2.high = apply(plot.ngood.corr,2,function(x) quantile(x,probs = 0.975))

plot.ngood.incu2.low = apply(plot.ngood.incu,2,function(x) quantile(x,probs = 0.025))
plot.ngood.incu2 = apply(plot.ngood.incu,2,function(x) quantile(x,probs = 0.5))
plot.ngood.incu2.high = apply(plot.ngood.incu,2,function(x) quantile(x,probs = 0.975))

cov.mean = mean(sl.summary$msl.mean)
cov.sd = sd(sl.summary$msl.mean)
plot.tmp = (seq.plot*cov.sd)+cov.mean 
plot.x = seq(min(plot.tmp),max(plot.tmp),length.out = 100)

col.plot = viridis(5)

plot(plot.good.spi2,type="l",x=plot.x,cex=1.2,
     ylim=c(0,1),xaxt="n",yaxt="n",col=col.plot[1],
     xlab="SLR",ylab="Nest success",
     lwd=2)
lines(plot.good.spi2.low,lty=2,x=plot.x,col=col.plot[1],lwd=2)
lines(plot.good.spi2.high,lty=2,x=plot.x,col=col.plot[1],lwd=2)

lines(plot.pais2,type="l",x=plot.x,col=col.plot[2],lwd=2)
lines(plot.pais2.low,lty=2,x=plot.x,col=col.plot[2],lwd=2)
lines(plot.pais2.high,lty=2,x=plot.x,col=col.plot[2],lwd=2)

lines(plot.good.incu2,type="l",x=plot.x,col=col.plot[3],lwd=2)
lines(plot.good.incu2.high,lty=2,x=plot.x,col=col.plot[3],lwd=2)
lines(plot.good.incu2.low,lty=2,x=plot.x,col=col.plot[3],lwd=2)

lines(plot.ngood.corr2,type="l",x=plot.x,col=col.plot[4],lwd=2)
lines(plot.ngood.corr2.low,lty=2,x=plot.x,col=col.plot[4],lwd=2)
lines(plot.ngood.corr2.high,lty=2,x=plot.x,col=col.plot[4],lwd=2)

lines(plot.ngood.incu2,type="l",x=plot.x,col=col.plot[5],lwd=2)
lines(plot.ngood.incu2.low,lty=2,x=plot.x,col=col.plot[5],lwd=2)
lines(plot.ngood.incu2.high,lty=2,x=plot.x,col=col.plot[5],lwd=2)

axis(side=1, at=NULL, labels=TRUE)
axis(side=2, labels=TRUE, las=1)	# will add tick marks in every case, but labels only on the first panel
legend("bottomright",title="Management Type",
       legend=c("SPI Good, Cor","PAIS Good, Cor",
                "Good, Incu","Not Good, Cor",
                "Not Good, Incu"), 
       pch=16,col=col.plot,bty="n",
       x.intersp=0.7, y.intersp=1.0, cex=1, pt.cex=1)
