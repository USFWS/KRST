#make SLR and climate data
make.fore.dat <- function(hist.bioc.means,
                      clim.var, ssp.in,slrp.summary,
                      sl.summary,slr.scenario){

bio5 = subset(hist.bioc.means,bioc == clim.var)
bio5.out = rep(NA,length(krst_data$J.ind))
bio5.fore.245 = subset(bioc.forecast,bioc == clim.var & ssp == ssp.in)
bio5.scale = as.numeric(scale(c(bio5$mean,bio5.fore.245$mean)))
bio5.tmp = data.frame(mean = c(rep(bio5.scale[32],length.out=19),
                               rep(bio5.scale[33],length.out=20),
                               rep(bio5.scale[34],length.out=20),
                               rep(bio5.scale[35],length.out=20)),
                      year = seq(2022,2100))

bio5$mean.scale = bio5.scale[1:31]

slr.cov = rep(NA,length(krst_data$J.ind))
slr.fore.med = subset(slrp.summary, Scenario == slr.scenario)
slr.scale = as.numeric(scale(c(sl.summary$msl.mean[2:8],slr.fore.med$msl.mean)))
sl.summary$mean.scale = c(NA,slr.scale[1:7])
slr.tmp = data.frame(mean = c(rep(slr.scale[8],length.out=19),
                              rep(slr.scale[9],length.out=20),
                              rep(slr.scale[10],length.out=20),
                              rep(slr.scale[11],length.out=20)),
                     year = seq(2022,2100))


#may.temp = subset(hist.temp.means,month == 5)
#may.tmax = rep(NA,length(nest.pro$Year)-1)
for(i in 1:1597){
  #may.tmax[i] = may.temp$tmax[which(may.temp$year==(year+1990)[i])]
  bio5.out[i] = bio5$mean.scale[which(bio5$year==(year+1993)[i])]
  slr.cov[i] = sl.summary$mean.scale[which(sl.summary==(year+1993)[i])]
}
#may.tmax.ann = as.numeric(scale(may.temp$tmax))
bio5.ann = bio5$mean.scale
#slr.cov = as.numeric(scale(slr.cov))

bio5.fore.cov = bio5.tmp$mean[1:(Time+Time.fore)]
slr.cov.fore = slr.tmp$mean[1:(Time+Time.fore)]

return(list(bio5.ann = bio5.ann, slr.cov = slr.cov,
       bio5.fore.cov = bio5.fore.cov, 
       slr.cov.fore = slr.cov.fore,
       bio5.out = bio5.out))
}
