library(viridis)
library(coda)
library(vioplot)
library(gplots)

load("~/KRST/Analysis/results_currentnm_5852_2100.RData")
load("~/KRST/Analysis/results_currentnm_2455_2100.RData")
load("~/KRST/Analysis/results_lowis_2455_2100.RData")
load("~/KRST/Analysis/results_noincu_2455_2100.RData")
load("~/KRST/Analysis/results_currentnm_2455_is1_2100.RData")
load("~/KRST/Analysis/results_currentnm_2455_is0_2100.RData")
load("~/KRST/Analysis/results_5852_is25_2100.RData")
load("~/KRST/Analysis/results_currentnm_2455_is5_2100.RData")
load("~/KRST/Analysis/results_currentnm_2455_is25_2100.RData")

out.2455 = as.mcmc(rbind(out.current.2455[[1]], out.current.2455[[2]], out.current.2455[[3]]))
out.5852 = as.mcmc(rbind(out.current.5852[[1]], out.current.5852[[2]], out.current.5852[[3]]))

out.is0 = as.mcmc(rbind(out.current.2455.is0[[1]], out.current.2455.is0[[2]], out.current.2455.is0[[3]]))
out.is1 =as.mcmc(rbind(out.current.2455.is1[[1]], out.current.2455.is1[[2]], out.current.2455.is1[[3]]))
out.is255 = as.mcmc(rbind(out.current.5852.is25[[1]], out.current.5852.is25[[2]], out.current.5852.is25[[3]]))
out.is252 = as.mcmc(rbind(out.current.2455.is25[[1]], out.current.2455.is25[[2]], out.current.2455.is25[[3]]))
out.is5 = as.mcmc(rbind(out.current.2455.is5[[1]], out.current.2455.is5[[2]], out.current.2455.is5[[3]]))

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

name.imm.low = which(varnames(out.5852)=="N.i[1]")
name.imm.high = which(varnames(out.5852)=="N.i[107]")

imm = apply(out.5852[,(name.imm.low):(name.imm.high)],2,quantile,.5)
imm.l = apply(out.5852[,(name.imm.low):(name.imm.high)],2,quantile,0.025)
imm.h = apply(out.5852[,(name.imm.low):(name.imm.high)],2,quantile,0.975)
Nnb = apply(out.5852[,which(varnames(out.5852)=="N.nb[1]"):(which(varnames(out.5852)=="N.nb[107]"))],2,mean)
Nnb.low = apply(out.5852[,which(varnames(out.5852)=="N.nb[1]"):which(varnames(out.5852)=="N.nb[107]")],2,quantile,0.025)
Nnb.high = apply(out.5852[,which(varnames(out.5852)=="N.nb[1]"):which(varnames(out.5852)=="N.nb[107]")],2,quantile,0.975)
#Np1 = apply(out[,63:93],2,mean)
#Np2 = apply(out[,94:124],2,mean)
#Np3 = apply(out[,125:155],2,mean)
#Np4 = apply(out[,156:186],2,mean)
#Np5 = apply(out[,187:217],2,mean)
#Np6 = apply(out[,218:248],2,mean)
#Np7 = apply(out[,249:279],2,mean)

#Np = rowSums(cbind(Np1,Np2,Np3,Np4,Np5,Np6,Np7))
Nad = rowSums(cbind(Nb.5852,Nnb,imm))

Nad.low <- Nad.high <- rep(NA,(Time+Time.fore))
for(i in 1:(Time+Time.fore)){
Nad.low[i]=quantile(rowSums(cbind(out.5852[,i],out.5852[,i+(2*(Time+Time.fore))],out.5852[,(i+name.imm.low-1)])),0.025)
Nad.high[i] =quantile(rowSums(cbind(out.5852[,i],out.5852[,i+(2*(Time+Time.fore))],out.5852[,(i+name.imm.low-1)])),0.975)}

par(mfrow=c(2,1))

## Plot of total TX adults (breed & nb) to current"

plot(Nad[1:Time],x=seq(1994,(1993+Time)), type="l", 
     ylab = "Abundance",
     xlab= "",col=col.plot[1],ylim=c(0,max(Nad.high[1:28])))
polygon(c(seq(1994,2021),rev(seq(1994,2021))),
        c(Nad.low[1:28],rev(Nad.high[1:28])),col=col.plot[1])
mtext("A) Total Abundance", side=3, line=.5, adj=0)
#axis(side=1, at=seq(1994,2021), labels=FALSE)

plot(Nb.5852[1:28], x = seq(1994,2021), type="l", ylab="Abundance",
     xlab="Year",ylim=c(0,max(Nnb.high[1:28])))
polygon(c(seq(1994,2021), rev(seq(1994,2021))),
        c(Nb.low.5852[1:28],rev(Nb.high.5852[1:28])),col=col.plot[2])


lines(Nnb[1:28], x= seq(1994,2021),type="l")
polygon(c(seq(1994,2021),rev(seq(1994,2021))),
        c(Nnb.low[1:28],rev(Nnb.high[1:28])),col=col.plot[3])

lines(imm[1:28], x = seq(1994,2021), type="l", ylab="Abundance",
      xlab="Year",ylim=c(0,max(imm.h)))
polygon(c(seq(1994,2021), rev(seq(1994,2021))),
        c(imm.l[1:28],rev(imm.h[1:28])),col=col.plot[5])

legend("topleft", 
       legend=c("Breeders","Non-Breeders","Immigrants"), pch=c(16,16),
       bty="n",col=c(col.plot[2],col.plot[3],col.plot[5]),
       x.intersp=0.7, y.intersp=1.0, cex=1.1, pt.cex=1.5)

mtext("B) Abundance by category", side=3, line=.5, adj=0)
dev.print(tiff,"tx_ad.tiff",res=600,width=9,units="in")


#forecast adults

plot(Nad,x=seq(1994,plot.year), type="l", ylab = "Abundance", 
     xlab= "Year",col=col.plot[1],ylim=c(0,max(Nad.high)))
polygon(c(seq(1994,plot.year),rev(seq(1994,plot.year))),
     c(Nad.low,rev(Nad.high)),col=col.plot[1])
#mtext("A) Total Abundance", side=3, line=.5, adj=0)
#axis(side=1, at=seq(1991,2021), labels=FALSE)

plot(Nb.5852, x = seq(1994,plot.year), type="l", ylab="Abundance",
     xlab="Year",ylim=c(0,max(Nb.high.5852)))
polygon(c(seq(1994,plot.year), rev(seq(1994,plot.year))),
        c(Nb.low.5852,rev(Nb.high.5852)),col=col.plot[2])


lines(Nnb, x= seq(1994,plot.year),type="l")
polygon(c(seq(1994,plot.year),rev(seq(1994,plot.year))),
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
  mfrow=c(2,1), # panels will plot in 1 row with 3 columns
  oma=c(3,3.25,3,3), # outer margins (bottom, left, top, right)
  mai=c(0.75,0.1,0.1,0.1), # inner margins (between panels)
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
lines(imm[1:28], x = seq(1994,2021), type="l", ylab="Abundance",
     xlab="Year",ylim=c(0,max(imm.h)))
polygon(c(seq(1994,2021), rev(seq(1994,2021))),
        c(imm.l[1:28],rev(imm.h[1:28])),col=col.plot[5])
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
Nb.is255 = apply(out.is255[,1:(Time+Time.fore)],2,quantile,.5)
Nb.low.is255 = apply(out.is255[,1:(Time+Time.fore)],2,quantile,0.025)
Nb.high.is255 = apply(out.is255[,1:(Time+Time.fore)],2,quantile,0.975)
Nb.is252 = apply(out.is252[,1:(Time+Time.fore)],2,quantile,.5)
Nb.low.is252 = apply(out.is252[,1:(Time+Time.fore)],2,quantile,0.025)
Nb.high.is252 = apply(out.is252[,1:(Time+Time.fore)],2,quantile,0.975)
Nb.is5 = apply(out.is5[,1:(Time+Time.fore)],2,quantile,.5)
Nb.low.is5 = apply(out.is5[,1:(Time+Time.fore)],2,quantile,0.025)
Nb.high.is5 = apply(out.is5[,1:(Time+Time.fore)],2,quantile,0.975)


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
lam.b.is255 = matrix(NA,30000, (Time+Time.fore-1))
lam.b.is252 = matrix(NA,30000, (Time+Time.fore-1))
lam.b.is5 = matrix(NA,30000, (Time+Time.fore-1))
lam.b.lis = matrix(NA,30000, (Time+Time.fore-1))
lam.b.noin = matrix(NA,30000, (Time+Time.fore-1))

for(t in 1:(Time+Time.fore-1)){ 
  lam.b.2455[,t] = out.2455[,t+1]/out.2455[,t]
  lam.b.5852[,t] = out.5852[,t+1]/out.5852[,t]
  lam.b.is0[,t] = out.is0[,t+1]/out.is0[,t]
  lam.b.is1[,t] = out.is1[,t+1]/out.is1[,t]
  lam.b.is255[,t] = out.is255[,t+1]/out.is255[,t]
  lam.b.is252[,t] = out.is252[,t+1]/out.is252[,t]
  lam.b.is5[,t] = out.is5[,t+1]/out.is5[,t]
  lam.b.lis[,t] = out.lis[,t+1]/out.lis[,t]
  lam.b.noin[,t] = out.noincu[,t+1]/out.noincu[,t]
  }
#remove Inf from lam
lam.b.2455[!is.finite(lam.b.2455)] = 0
lam.b.5852[!is.finite(lam.b.5852)] = 0
lam.b.is0[!is.finite(lam.b.is0)] = 0
lam.b.is1[!is.finite(lam.b.is1)]= 0
lam.b.is255[!is.finite(lam.b.is255)]=0
lam.b.is252[!is.finite(lam.b.is252)]=0
lam.b.is5[!is.finite(lam.b.is5)]=0
lam.b.lis[!is.finite(lam.b.lis)] = 0
lam.b.noin[!is.finite(lam.b.noin)] = 0

#calculate average year population is extripated
extrip.year = matrix(NA,30000,5)
for(i in 1:30000){
  extrip.year[i,1] = max(which(out.is0[i,1:107]>0))
  extrip.year[i,2] = max(which(out.is252[i,1:107]>0))
  extrip.year[i,3] = max(which(out.is255[i,1:107]>0))
  extrip.year[i,4] = max(which(out.is5[i,1:107]>0))
  extrip.year[i,5] = max(which(out.is1[i,1:107]>0))}

idx = extrip.year<107

extirp.out = matrix(NA,3,5)
for(i in 1:5){
extirp.out[1,i]=quantile(extrip.year[idx[,i],i],.025)+(2021-28)
extirp.out[2,i]=quantile(extrip.year[idx[,i],i],.5)+(2021-28)
extirp.out[3,i]=quantile(extrip.year[idx[,i],i],.975)+(2021-28)
}

#geo mean of lambda through future
geo_mean <- function(data){
  gm = rep(NA,30000)
  for(i in 1:length(data[,1])){
  gm[i] <- exp(mean(log(data[i,])))
  }
  return(gm)
  }

gm.lam.2455 = geo_mean(lam.b.2455[,29:106])
gm.lam.5852 = geo_mean(lam.b.5852[,29:106])
gm.lam.is1 = geo_mean(lam.b.is1[,29:106])
gm.lam.is0 = geo_mean(lam.b.is0[,29:106])
gm.lam.is255 = geo_mean(lam.b.is255[,29:106])
gm.lam.is252=geo_mean(lam.b.is252[,29:106])
gm.lam.is5=geo_mean(lam.b.is5[,29:106])
gm.lam.lis = geo_mean(lam.b.lis[,29:106])
gm.lam.noin = geo_mean(lam.b.noin[,29:106])

#add zeros to infinite values
gm.lam.is255[!is.finite(gm.lam.is255)] = 0
gm.lam.is252[!is.finite(gm.lam.is252)] = 0
gm.lam.is5[!is.finite(gm.lam.is5)] = 0
gm.lam.is0[!is.finite(gm.lam.is0)] =0
gm.lam.is1[!is.finite(gm.lam.is1)] = 0

#boxplot of not all in situ management
boxplot(cbind(gm.lam.2455,gm.lam.5852,gm.lam.lis,
              gm.lam.noin),
        ylab = expression(lambda),
        names = c("Current,SSP2","Current,SSP5",
                  "Limited In Situ","All Corrals"),
        outline=FALSE)
abline(h=1)

dev.print(tiff,"geo_lam.tiff",res=600,width=9,units="in")

zero.gm.2455 = length(which(gm.lam.2455==0))/30000
zero.gm.5852 = length(which(gm.lam.5852==0))/30000
zero.gm.is1 = length(which(gm.lam.is1==0))/30000
zero.gm.is0 = length(which(gm.lam.is0==0))/30000
zero.gm.is5 = length(which(gm.lam.is5==0))/30000
zero.gm.is252 = length(which(gm.lam.is252==0))/30000
zero.gm.is255 = length(which(gm.lam.is255==0))/30000
zero.gm.lis = length(which(gm.lam.lis==0))/30000
zero.gm.noin = length(which(gm.lam.noin==0))/30000

#boxplot of all in situ management
boxplot(cbind(gm.lam.is0, 
              gm.lam.is252,gm.lam.is255,
              gm.lam.is5,gm.lam.is1),
        ylab = expression(lambda),
        names = c("Effect = 0", 
                  "Effect = -0.25, SSP2",
                  "Effect = -0.25, SSP5",
                  "Effect =  -0.5",
                  "Effect = -1"),outline=FALSE)
mtext(expression(lambda), side=2, line=2.5,las=1)

text(x=1,y=0.5,labels=paste0("P(N = 0) = ",round(zero.gm.is0,2)))
text(x=2,y=0.5,labels=paste0("P(N = 0) = ",round(zero.gm.is252,2)))
text(x=3,y=0.5,labels=paste0("P(N = 0) = ",round(zero.gm.is255,2)))
text(x=4,y=0.5,labels=paste0("P(N = 0) = ",round(zero.gm.is5,2)))
text(x=5,y=0.5,labels=paste0("P(N = 0) = ",round(zero.gm.is1,2)))

abline(h=1)

dev.print(tiff,"geo_lam_insitu.tiff",res=600,width=11,units="in")

lam.low.fore.2455 = apply(lam.b.2455[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.5852 = apply(lam.b.5852[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.is0 = apply(lam.b.is0[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.is1 = apply(lam.b.is1[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.is255 = apply(lam.b.is255[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.is252 = apply(lam.b.is252[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.is5 = apply(lam.b.is5[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.lis = apply(lam.b.lis[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)
lam.low.fore.noin = apply(lam.b.noin[,(Time):(Time+Time.fore-1)],2,quantile,0.025,na.rm=TRUE)

lam.high.fore.2455 = apply(lam.b.2455[,(Time):(Time+Time.fore-1)],2,quantile,0.975)
lam.high.fore.5852 = apply(lam.b.5852[,(Time):(Time+Time.fore-1)],2,quantile,0.975)

lam.b.is0[!is.finite(lam.b.is0)] = 0
lam.b.is1[!is.finite(lam.b.is1)] = 0
lam.b.is255[!is.finite(lam.b.is255)] = 0
lam.b.is252[!is.finite(lam.b.is252)] = 0
lam.b.is5[!is.finite(lam.b.is5)] = 0
lam.high.fore.is0 = apply(lam.b.is0[,(Time):(Time+Time.fore-1)],2,quantile,0.975,na.rm=TRUE)
lam.high.fore.is1 = apply(lam.b.is1[,(Time):(Time+Time.fore-1)],2,quantile,0.975,na.rm=TRUE)
lam.high.fore.is255 = apply(lam.b.is255[,(Time):(Time+Time.fore-1)],2,quantile,0.975,na.rm=TRUE)
lam.high.fore.is252 = apply(lam.b.is252[,(Time):(Time+Time.fore-1)],2,quantile,0.975,na.rm=TRUE)
lam.high.fore.is5 = apply(lam.b.is5[,(Time):(Time+Time.fore-1)],2,quantile,0.975,na.rm=TRUE)

lam.high.fore.lis = apply(lam.b.lis[,(Time):(Time+Time.fore-1)],2,quantile,0.975,na.rm=TRUE)
lam.high.fore.noin = apply(lam.b.noin[,(Time):(Time+Time.fore-1)],2,quantile,0.975,na.rm=TRUE)

lam.med.all.2455 = apply(lam.b.2455[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.5852 = apply(lam.b.5852[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.is0 = apply(lam.b.is0[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.is1 = apply(lam.b.is1[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.is255 = apply(lam.b.is255[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.is252 = apply(lam.b.is252[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.is5 = apply(lam.b.is5[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.lis = apply(lam.b.lis[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)
lam.med.all.noin = apply(lam.b.noin[,1:(Time+Time.fore-1)],2,quantile,0.5,na.rm=TRUE)

par(mfrow=c(1,1))
#vioplot(lam[,-1],ylim=c(0,3),names=seq(1993,(1991+Time+Time.fore)))
boxplot(cbind(lam.b.2455[,1:(Time-1)]),
        ylab = expression(lambda),
        names = seq(1995,(1993+Time)),outline=FALSE)

#Not all in situ scenarios
par(
  mfrow=c(2,2), # panels will plot in 1 row with 3 columns
  oma=c(3,3.25,3,3), # outer margins (bottom, left, top, right)
  mai=c(0.5,0.1,0.1,0.1), # inner margins (between panels)
  cex=1, 	# make sure the text and points are a normal size even if your plot window is an unusual size (R likes to resize spontaneously)
  mgp=c(3,0.8,0))

lam.not.insitu = array(data = cbind(lam.low.fore.2455,
                              lam.low.fore.5852,
                              lam.low.fore.lis,
                              lam.low.fore.noin,
                              lam.med.all.2455[Time:(Time+Time.fore-1)],
                              lam.med.all.5852[Time:(Time+Time.fore-1)],
                              lam.med.all.lis[Time:(Time+Time.fore-1)],
                              lam.med.all.noin[Time:(Time+Time.fore-1)],
                              lam.high.fore.2455,
                              lam.high.fore.5852,
                              lam.high.fore.lis,
                              lam.high.fore.noin),
                       dim=c(79,4,3))

scenarios = c("Current, SSP2","Current, SSP5",
              "Limited In Situ","All Corrals")

for(i in 1:4){
  
  boxplot(#cbind(lam.b.2455[,(Time-2):(Time-1)],
          matrix(NA,30000,(Time.fore+1)),
          ylab = expression(lambda),outline=FALSE,
          xaxt='n',yaxt="n",
          ylim = c(0.5, max(lam.not.insitu[,,3])))

  polygon(c(seq(3,(Time.fore+2)), 
          rev(seq(3,(Time.fore+2)))),
        c(lam.not.insitu[,i,1],rev(lam.not.insitu[,i,3])),
          col=col.plot[3])

  lines(y=lam.not.insitu[,i,2],
      x=seq(3,(Time.fore+2)),lwd=2)

  abline(h=1)
  
  if(i==1 | i==2){
  axis(side=1, at=seq(1,81,by=10),labels=FALSE,las=1)}
  if(i==3 | i==4){
  axis(side=1, at=seq(1,81,by=10),
       labels=seq(2020,2100,by=10),las=1)}
  
  axis(side=2, labels=(i==1|i==3), las=1)	# will add tick marks in every case, but labels only on the first panel
  if(i==1|i==3) mtext(expression(lambda), side=2, line=2.5,las=1)
  if(i==3|i==4) mtext("Year", side=1, line=2)
  
  # add panel title:
  mtext(paste0("Scenario ",i,": ", scenarios[i]), side=3,  adj=0)
  }

dev.print(tiff,"lam_2100_notinsitu.tiff",res=600,width=9,units="in")

## in situ only

par(
  mfrow=c(3,2),
  oma=c(3,3.25,3,3), # outer margins (bottom, left, top, right)
  mai=c(0.5,0.1,0.1,0.1), # inner margins (between panels)
  cex=1, 	# make sure the text and points are a normal size even if your plot window is an unusual size (R likes to resize spontaneously)
  mgp=c(3,0.8,0))

lam.all.insitu = array(cbind(lam.low.fore.is0,
                     lam.low.fore.is252,
                     lam.low.fore.is255,
                     lam.low.fore.is5,
                     lam.low.fore.is1,
                     lam.med.all.is0[Time:(Time+Time.fore-1)],
                     lam.med.all.is252[Time:(Time+Time.fore-1)],
                     lam.med.all.is255[Time:(Time+Time.fore-1)],
                     lam.med.all.is5[Time:(Time+Time.fore-1)],
                     lam.med.all.is1[Time:(Time+Time.fore-1)],
                     lam.high.fore.is0,
                     lam.high.fore.is252,
                     lam.high.fore.is255,
                     lam.high.fore.is5,
                     lam.high.fore.is1),
                     dim=c(79,5,3))

scenarios = c("All In Situ 0","All In Situ -0.25,SSP2",
              "All In Situ -0.25,SSP5","All In Situ -0.5",
              "All In Situ -1")

for(i in 1:5){
  
  boxplot(#cbind(lam.b.2455[,(Time-2):(Time-1)],
          matrix(NA,30000,(Time.fore+1)),
          ylab = expression(lambda),outline=FALSE,
          xaxt='n',yaxt="n",
          ylim = c(0, 3))

  polygon(c(seq(3,(Time.fore+2)), 
          rev(seq(3,(Time.fore+2)))),
        c(lam.all.insitu[,i,1],rev(lam.all.insitu[,i,3])),
        col=col.plot[3])
  lines(y=lam.all.insitu[,i,2],
        x=seq(3,(Time.fore+2)),lwd=2)
  abline(h=1)

if(i==1|i==2|i==3){
  axis(side=1, at=seq(1,81,by=10),labels=FALSE,las=1)}
if(i==4 | i==5){
  axis(side=1, at=seq(1,81,by=10),
       labels=seq(2020,2100,by=10),las=1)}

axis(side=2, labels=(i==1|i==3|i==5), las=1)	# will add tick marks in every case, but labels only on the first panel
if(i==1|i==3|i==5) mtext(expression(lambda), side=2, line=2.5,las=1)
if(i==4|i==5) mtext("Year", side=1, line=2)

# add panel title:
mtext(paste0("Scenario ",i+4,": ", scenarios[i]), side=3,  adj=0)
}

dev.print(tiff,"lam_2100_insitu.tiff",res=600,width=9,units="in")

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


lines(y=lam.med.all.is0[Time:(Time+Time.fore-1)],
      x=seq(Time,(Time+Time.fore-1)),lty=1,lwd=2)
lines(y=lam.med.all.is1[Time:(Time+Time.fore-1)],
      x=seq(Time,(Time+Time.fore-1)),lty=2,lwd=2)
lines(y=lam.med.all.lis[Time:(Time+Time.fore-1)],
      x=seq(Time,(Time+Time.fore-1)),lty=2,lwd=2)
lines(y=lam.med.all.noin[Time:(Time+Time.fore-1)],
      x=seq(Time,(Time+Time.fore-1)),lty=2,lwd=2)


abline(h=1)

dev.print(tiff,"lam_2021.tiff",res=600,width=9,units="in")


lines(x=c(20,20),y=c(0,4))
lines(x=c(81,81),y=c(0,4))

#probability that lambda is > 1 at certain time periods
#20 years in future
lam.name = array(c(lam.b.2455,lam.b.5852,
                   lam.b.lis,lam.b.noin,lam.b.is0,
                   lam.b.is255,lam.b.is252,
                   lam.b.is5,lam.b.is1),
                 dim=c(30000,106,9))
lam.2042 = rep(NA,9)
lam.2100 = rep(NA,9)

for(i in 1:9){
prob.less.0 = ecdf(lam.name[,48,i])
lam.2042[i]=1-prob.less.0(.99)

#2100
prob.less.0=ecdf(lam.name[,106,i])
lam.2100[i]=1-prob.less.0(.99)
}

text(x=53,y=4,labels=expression(P(lambda > 1) == 0.53))
text(x=101,y=4,lables=expression(P(lambda > 1) == 0.60))

legend("topright", 
       legend=c("SSP2-4.5, 0.5 m", "SSP5-8.5, 2 m"), pch=c(16,16),
       bty="n",col=c(col.plot[3],col.plot[5]),
       x.intersp=0.7, y.intersp=1.0, cex=1.1, pt.cex=1.5)

#Detection probability
det.prob = out.5852[,1415]
plot(density(det.prob))

boxplot(det.prob,outline=FALSE,
        names="",
        ylab = "Detection Probability",
        col="white",ylim = c(.5,1))

#transition probabilities
boxplot(out.5852[,1522],out.5852[,1523],
        names=c("Breeder to Non","Non to Breeder"),
        outline=FALSE)

#### survival rates through time ####

phi.plot = as.matrix(out.2455[,1416:1442])
boxplot(phi.plot,names=seq(1995,(1993+Time)),
        ylab = "Adult Survival",outline=FALSE)
#abline(h=(1/(1+exp(-(mean(out.5852[,1407]))))))

dev.print(tiff,"ad_surv.tiff",res=600,width=9,units="in")

##### covs plots ####
#climate effects
#clim.cov.plot = out[,311]

bio5 = subset(hist.bioc.means,bioc == "bio5")
bio5.fore = subset(bioc.forecast,bioc=="bio5")
ssp2 = subset(bio5.fore, ssp == "245")
ssp5 = subset(bio5.fore,ssp == "585")

#plot of climate covariate past and future
plot(c(bio5$mean[4:31],rep(NA,79)),x=seq(1994,2100),type="l",
     xlab = "Year",ylab="Max Temperature (C)",ylim=c(31,42))
lines(rep(ssp2$mean,each=20),x=seq(2021,2100),col=col.plot[3])
lines(rep(ssp5$mean,each=20),x=seq(2021,2100),col=col.plot[2])

legend("topleft",legend=c("SSP2","SSP5"),pch=c(16,16),
       col=c(col.plot[3],col.plot[2]))

#yrs 1-24
seq.plot = seq(min(c(bio5.245.5$bio5.ann,bio5.245.5$bio5.fore.cov[1:79])),
               max(c(bio5.245.5$bio5.ann,bio5.245.5$bio5.fore.cov[1:79])),
               length.out=100)

#demonstration of effect of -1 on prob of nest success
demo.plot = array(NA,dim=c(30000,100,3))
for(j in 1:100){
  demo.plot[,j,1] = 1/(1+exp(-(out.is252[,1414]+out.is252[,1082]*seq.plot[j])))
  demo.plot[,j,2] = 1/(1+exp(-(out.is5[,1414]+out.is5[,1082]*seq.plot[j])))
  demo.plot[,j,3] = 1/(1+exp(-(out.is1[,1414]+out.is1[,1082]*seq.plot[j])))
}
clim = seq(min(subset(hist.bioc.means,bioc == "bio5")$mean),41,length.out=100)
plot(y = apply(demo.plot[,,1],2,quantile,.5), x = clim,
     type= "l", ylab = "Hatching Probability", 
     xlab = expression(paste("Max Temp ",degree,"C")),
     ylim = c(0,1),col=col.plot[1])
lines(y = apply(demo.plot[,,2],2,quantile,.5), x = clim,
      col=col.plot[2])
lines(y = apply(demo.plot[,,3],2,quantile,.5), x = clim,
      col=col.plot[3])
lines(y=rep(.63,100),x=clim, col=col.plot[4])

legend("bottomleft", 
       legend=c("No Effect","-0.25 Effect","-0.5 Effect",
                "-1 Effect"), pch=c(16,16,16,16),
       bty="n",col=c(col.plot[4],col.plot[1:3]),
       x.intersp=0.7, y.intersp=1.0, cex=1.1, pt.cex=1.5)

dev.print(tiff,"insitu_effect.tiff",res=600,width=9,units="in")

plot.low.incu <- plot.low.cor.pais <-
  plot.low.cor.spi <- plot.high.incu <-
  plot.high.cor <- plot.insitu <- matrix(NA,30000,100)

for(j in 1:100){
  #bioclim 5
  #order: "Low,Incu","low,Cor,PAIS", "Low,Cor,SPI","High,Incu","High,Cor","In situ")
  plot.low.incu[,j] = 1/(1+exp(-(out.2455[,which(varnames(out.2455)=="nest.cov[1]")]+ 
                                   out.2455[,which(varnames(out.2455)=="beta.temp[1]")]*seq.plot[j])))
  
  plot.low.cor.pais[,j] = 1/(1+exp(-(out.2455[,which(varnames(out.2455)=="nest.cov[2]")]+ 
                                   out.2455[,which(varnames(out.2455)=="beta.temp[2]")]*seq.plot[j])))

  plot.low.cor.spi[,j] = 1/(1+exp(-(out.2455[,which(varnames(out.2455)=="nest.cov[3]")]+ 
                                   out.2455[,which(varnames(out.2455)=="beta.temp[3]")]*seq.plot[j])))
  
  plot.high.incu[,j] = 1/(1+exp(-(out.2455[,which(varnames(out.2455)=="nest.cov[4]")]+ 
                                   out.2455[,which(varnames(out.2455)=="beta.temp[4]")]*seq.plot[j])))
  
  plot.high.cor[,j] = 1/(1+exp(-(out.2455[,which(varnames(out.2455)=="nest.cov[5]")]+ 
                                   out.2455[,which(varnames(out.2455)=="beta.temp[5]")]*seq.plot[j])))
  
  plot.insitu[,j] = 1/(1+exp(-(out.2455[,which(varnames(out.2455)=="nest.cov[6]")]+ 
                                   out.2455[,which(varnames(out.2455)=="beta.temp[6]")]*seq.plot[j])))
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

bio5.fore.245 = subset(bio5.fore,ssp == "245")

cov.mean = mean(c(bio5$mean,bio5.fore.245$mean))
cov.sd = sd(c(bio5$mean,bio5.fore.245$mean))
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
     ylim=c(0,1),xaxt="n",yaxt="n",#col=col.plot[1],
     lwd=2)
lines(plot.low.incu2.low,lty=2,x=plot.x,lwd=2)
lines(plot.low.incu2.high,lty=2,x=plot.x,lwd=2)
mtext("A) Low Risk, Incubator", side=3, line=.5, adj=0)
axis(side=2, labels=TRUE, las=1)
axis(side=1, labels=FALSE, las=1)

plot(plot.low.cor.pais2,type="l",x=plot.x,lwd=2,ylim=c(0,1),
     xaxt="n",yaxt="n")
lines(plot.low.cor.pais2.low,lty=2,x=plot.x,lwd=2)
lines(plot.low.cor.pais2.high,lty=2,x=plot.x,lwd=2)
mtext("B) Low Risk, Corral, PAIS", side=3, line=.5, adj=0)
axis(side=2, labels=FALSE, las=1)
axis(side=1, labels=FALSE, las=1)

plot(plot.high.incu2,type="l",x=plot.x,
     lwd=2,ylim=c(0,1),
     xaxt="n",yaxt="n")
lines(plot.high.incu2.low,lty=2,x=plot.x,lwd=2)
lines(plot.high.incu2.high,lty=2,x=plot.x,lwd=2)
mtext("D) High Risk, Incubator, PAIS and North", side=3, line=.5, adj=0)
axis(side=1, at=NULL, labels=FALSE)
axis(side=2, labels=FALSE, las=1)

plot(plot.low.cor.spi2,type="l",x=plot.x,lwd=2,ylim=c(0,1),
     xaxt="n",yaxt="n")
lines(plot.low.cor.spi2.high,lty=2,x=plot.x,lwd=2)
lines(plot.low.cor.spi2.low,lty=2,x=plot.x,lwd=2)
mtext("C) Low Risk, Corral, SPI", side=3, line=.5, adj=0)
axis(side=1, labels=TRUE, las=1)
axis(side=2, labels=TRUE, las=1)

plot(plot.high.cor2,type="l",x=plot.x,
     lwd=2,ylim=c(0,1),
     xaxt="n",yaxt="n")
lines(plot.high.cor2.low,lty=2,x=plot.x,lwd=2)
lines(plot.high.cor2.high,lty=2,x=plot.x,lwd=2)
mtext("E) High Risk, Corral, SPI", side=3, line=.5, adj=0)
axis(side=1, at=NULL, labels=TRUE)
axis(side=2, labels=FALSE, las=1)

plot(plot.insitu2,type="l",x=plot.x,
     lwd=2,ylim=c(0,1),
     xaxt="n",yaxt="n")
lines(plot.insitu2.low,lty=2,x=plot.x,lwd=2)
lines(plot.insitu2.high,lty=2,x=plot.x,lwd=2)
mtext("F) In situ", side=3, line=.5, adj=0)
axis(side=1, at=NULL, labels=TRUE)
axis(side=2, labels=FALSE, las=1)

mtext("Max Temp (C) of Warmest Month", 
      side=1, line=1, at=0.5,
      outer=TRUE)
mtext("Hatching Probability", side=2, line=1.5, 
      at=0.5,outer=TRUE,las=0)

dev.print(tiff,"climate_nests.tiff",res=600,width=9,units="in")

#nest management effects
#order: "Low,Incu","low,Cor,PAIS", "Low,Cor,SPI",
# "High,Incu","High,Cor","In situ")
good.corr.pais.plot = 1/(1+exp(-(out.2455[,1410])))
good.corr.spi.plot = 1/(1+exp(-(out.2455[,1411])))
good.incu.plot = 1/(1+exp(-(out.2455[,1409])))
ngood.corr.plot = 1/(1+exp(-(out.2455[,1413])))
ngood.incu.plot = 1/(1+exp(-(out.2455[,1412])))
insitu.plot = 1/(1+exp(-(out.2455[,1414])))

quantile(1/(1+exp(-(out.2455[,which(varnames(out.2455)=="nest.cov[3]")]))),c(0.025,0.5,0.975))

par(mfrow = c(1,1))

boxplot(good.corr.pais.plot,good.corr.spi.plot,
        good.incu.plot,ngood.corr.plot,
        ngood.incu.plot,insitu.plot,
        names = c("Low,Cor,PAIS", "Low,Cor,SPI",
                  "Low,Incu",
                  "High,Cor","High,Incu","In situ"),
        ylab = "Hatching Probability",outline = FALSE)

#points(x=5,y=0.63,col="red",pch=16)
dev.print(tiff,"nest_mange.tiff",res=600,width=9,units="in")

