library(viridis)
library(coda)
library(vioplot)
library(gplots)

load("~/KRST/Analysis/KRST_out.RData")

out = as.mcmc(rbind(out.up[[1]], out.up[[2]], out.up[[3]]))

col.plot=viridis(5,alpha=.6)

Time = 31
Time.fore = 5
plot.year = 1990+Time+Time.fore

## Figure for abundance
Nb = apply(out[,1:(Time+Time.fore)],2,mean)
Nb.low = apply(out[,1:(Time+Time.fore)],2,quantile,0.025)
Nb.high = apply(out[,1:(Time+Time.fore)],2,quantile,0.975)
imm = apply(out[,33:64],2,mean)
imm.l = apply(out[,33:64],2,quantile,0.025)
imm.h = apply(out[,33:64],2,quantile,0.975)
Nnb = apply(out[,65:96],2,mean)
Nnb.low = apply(out[,65:96],2,quantile,0.025)
Nnb.high = apply(out[,65:96],2,quantile,0.975)
#Np1 = apply(out[,63:93],2,mean)
#Np2 = apply(out[,94:124],2,mean)
#Np3 = apply(out[,125:155],2,mean)
#Np4 = apply(out[,156:186],2,mean)
#Np5 = apply(out[,187:217],2,mean)
#Np6 = apply(out[,218:248],2,mean)
#Np7 = apply(out[,249:279],2,mean)

#Np = rowSums(cbind(Np1,Np2,Np3,Np4,Np5,Np6,Np7))
Nad = rowSums(cbind(Nb,Nnb))

Nad.low <- Nad.high <- rep(NA,32)
for(i in 1:32){
Nad.low[i]=quantile(rowSums(cbind(out[,i],out[,i+64])),0.025)
Nad.high[i] =quantile(rowSums(cbind(out[,i],out[,i+64])),0.975)}

par(mfrow=c(2,1))

plot(Nad,x=seq(1991,plot.year), type="l", ylab = "Abundance", 
     xlab= "Year",col=col.plot[1],ylim=c(0,max(Nad.high)))
polygon(c(seq(1991,plot.year),rev(seq(1991,plot.year))),
     c(Nad.low,rev(Nad.high)),col=col.plot[1])
#mtext("A) Total Abundance", side=3, line=.5, adj=0)
#axis(side=1, at=seq(1991,2021), labels=FALSE)

plot(Nb, x = seq(1991,plot.year), type="l", ylab="Abundance",
     xlab="Year",ylim=c(0,max(Nad.high)))
polygon(c(seq(1991,plot.year), rev(seq(1991,plot.year))),
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

N.t = rowSums(cbind(Nb,imm))

Nt.low <- Nt.high <- rep(NA,(Time+Time.fore))
for(i in 1:(Time+Time.fore)){
  Nt.low[i]=quantile(rowSums(cbind(out[,i],out[,i+(Time+Time.fore)])),0.025)
  Nt.high[i] =quantile(rowSums(cbind(out[,i],out[,i+(Time+Time.fore)])),0.975)}

par(
  mfrow=c(3,1), # panels will plot in 1 row with 3 columns
  oma=c(3,3.25,3,3), # outer margins (bottom, left, top, right)
  mai=c(0.5,0.1,0.1,0.1), # inner margins (between panels)
  cex=1, 	# make sure the text and points are a normal size even if your plot window is an unusual size (R likes to resize spontaneously)
  mgp=c(3,0.8,0))

plot(N.t,x=seq(1991,plot.year), type="l", 
     xlab= "Year",col="black",
     ylab = "Abundance",
     ylim=c(0,max(Nt.high)))
polygon(c(seq(1991,plot.year),rev(seq(1991,plot.year))),
        c(Nt.low,rev(Nt.high)),col=col.plot[4])


#axis(side=1, at=seq(1991,2021), labels=FALSE)


plot(Nb, x= seq(1991,plot.year),type="l",
     ylim=c(0,max(imm.h)),
     xlab = "Year", ylab = "Abundance")
polygon(c(seq(1991,plot.year),rev(seq(1991,plot.year))),
        c(Nb.low,rev(Nb.high)),col=col.plot[3])
lines(imm, x = seq(1991,plot.year), type="l", ylab="Abundance",
     xlab="Year",ylim=c(0,max(imm.h)))
polygon(c(seq(1991,plot.year), rev(seq(1991,plot.year))),
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

#### LAMBDA ####
lam = matrix(NA,150000,30)
n.tot = matrix(NA,150000,31)

for(t in 1:31){
  n.tot[,t] = rowSums(cbind(out[,t], out[,62+t]))
}
for(t in 1:30){ 
  lam[,t] = n.tot[,t+1]/n.tot[,t]
}
#remove Inf from lam
lam[!is.finite(lam)] = NA

par(mfrow=c(1,1))
vioplot(lam[,-1],ylim=c(0,3),names=seq(1993,2021))
abline(h=1)

det.prob = out[,385]
plot(density(det.prob))

boxplot(det.prob,outline=FALSE,
        names="",
        ylab = "Detection Probability",
        col="white",ylim = c(.5,1))

#transition probabilities
boxplot(out[,416],out[,417],
        names=c("Breeder to Non","Non to Breeder"),
        outline=FALSE)

#survival rates through time
phi.plot = as.matrix(out[,342:372])
boxplot(phi.plot,names=seq(1991,2021),
        ylab = "Adult Survival",outline=FALSE)


##### covs plots ####
#climate effects
clim.cov.plot = out[,311]

#yrs 1-24
seq.plot = seq(min(bio5.ann[1:24]),
               max(bio5.ann[1:24]),
               length.out=100)

all.plot = matrix(NA,30000,100)
for(j in 1:100){
all.plot[,j] = 1/(1+exp(-(out[,330]+out[,337]*seq.plot[j])))
}

plot(all.plot[1,],x=seq.plot,xlab="Bio5",type="l",
     ylim=c(0,1))

for(j in 2:100){
  lines(all.plot[j,],x=seq.plot)
}


#nest management effects
#360-380
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
