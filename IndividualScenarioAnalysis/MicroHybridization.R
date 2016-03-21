source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
setwd("C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/mh")

WrapMtFunc(Simuls =  c("mh1","mh2"),FilePrefix ="mh",Param = c(0.9,0.3))

read.table(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/mh/mhresults.txt",header = T)

mhDistriMt<-read.table(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/mh/mhDistriMtMax.txt",header = T)

plot(density(mhDistriMt$DistriMtMax[mhDistriMt$Param==0.3]-mhDistriMt$DistriNuMtMax[mhDistriMt$Param==0.3]),col="red")
lines(density(mhDistriMt$DistriMtMax[mhDistriMt$Param==0.9]-mhDistriMt$DistriNuMtMax[mhDistriMt$Param==0.9]),col="blue")

hixtx<-hist(x = mhDistriMt$DistriMtMax[mhDistriMt$Param==0.3]-mhDistriMt$DistriNuMtMax[mhDistriMt$Param==0.3],breaks = 50,plot=F)
hixtx$counts<-log(hixtx$counts+1)
plot(hixtx,xlim=c(-0.4,1),col=rgb(red = 1,0,0,0.5))
hixty<-hist(x = mhDistriMt$DistriMtMax[mhDistriMt$Param==0.9]-mhDistriMt$DistriNuMtMax[mhDistriMt$Param==0.9],breaks = 50,plot = F)
hixty$counts<-log(hixty$counts+1)
plot(hixty,add=T,col=rgb(red = 0,0,1,0.5))

hist(x = c(mhDistriMt$DistriMtMax[mhDistriMt$Param==0.3]-mhDistriMt$DistriNuMtMax[mhDistriMt$Param==0.3],mhDistriMt$DistriMtMax[mhDistriMt$Param==0.9]-mhDistriMt$DistriNuMtMax[mhDistriMt$Param==0.9]),
     col=c(rgb(red = 1,0,0,0.5),rgb(0,0,blue=1,0.5)),breaks = 50,xlim=c(-0.4,1))

out<-mtfunc("mh1")
results<-out$results
resultsNoFix<-out$resultsNoFix
write.table(results,file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/mh/mhresults.txt", sep="\t",quote=F,append=F,col.names=T)
hist(out$DistribAuto$Introg10)

hist(out$DistribAuto$Introg)

m0000<-lm(out$DistribAuto$Introg~1)
summary(m0000)
plot(m0000)

min(out$DistribAuto$Introg10)
min(out$DistribAuto$Introg)
min(out$DistribAuto$Fix)
min(out$DistribAuto$MeanExo)


out<-mtfunc("mh2")
results<-out$results
resultsNoFix<-out$resultsNoFix
write.table(results,file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/mh/mhresults.txt", sep="\t",quote=F,append=T,col.names=F)
hist(out$DistribAuto$Introg10)

par(mfrow=c(1,1))
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/mh/")
simul<-read.table(file="mhresults.txt",header=T)
simul

Simul="mh1"
