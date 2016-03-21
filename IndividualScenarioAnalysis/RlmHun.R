setwd("C:/Users/Timothée/Dropbox/PodarcisCBGP/Renew simulations/lmHun/")
source(file="C:/Users/Timothée/Dropbox/PodarcisCBGP/New simulations/FunctionsSimul.R")
Param<-matrix(data=c(1,0.0015/0.0005,0.0005/0.0015,0.0019/0.0001,0.0001/0.0019,0.00199/0.00001,0.00001/0.00199),nrow=7,ncol=1)

Simul<-"lmHun0"
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax")
DistriMtMaxDF$Param<-Param[0+1]
write.table(results,file="lmHunresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="lmHunresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="lmHunDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)


for (nbsimul in 1:6)
{
  Simul<-paste("lmHun",nbsimul,sep="")
  outR<-mtfunc(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
  names(DistriMtMaxDF)<-c("Simul","DistriMtMax")
  DistriMtMaxDF$Param<-Param[nbsimul+1]
  write.table(results,file="lmHunresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="lmHunresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(DistriMtMaxDF,file="lmHunDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}


simul<-read.table(file="lmHunresults.txt",header=T)
simulnofix<-read.table(file="lmHunresultsnofix.txt",header=T)


######## Only the invaded part of the grid

Simul<-"lmHun0"
outR<-mtfuncInvasion(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax")
DistriMtMaxDF$Param<-Param[0+1]
write.table(results,file="lmHunPrresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="lmHunPrresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="lmHunPrDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)


for (nbsimul in 1:6)
{
  Simul<-paste("lmHun",nbsimul,sep="")
  outR<-mtfuncInvasion(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
  names(DistriMtMaxDF)<-c("Simul","DistriMtMax")
  DistriMtMaxDF$Param<-Param[nbsimul+1]
  write.table(results,file="lmHunPrresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="lmHunPrresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(DistriMtMaxDF,file="lmHunPrDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}


simulPr<-read.table(file="lmHunPrresults.txt",header=T)
simulnofixPr<-read.table(file="lmHunresultsnofix.txt",header=T)



Param<-matrix(data=c(1,1.6,0.626,19.11,0.052,181.25,0.006),nrow=7,ncol=1)
simul$Param<-log(Param)
simul<-simul[order(simul$Param),]
plot(simul$Param,simul$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n En fonction du logarithme du rapport de dispersion entre sexe (Femelle/Male)",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(simul$Param,simul$IntMt,col="blue",pch=2,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Simulations avec Introgression Mitochondriale"),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)

simulnofix$Param<-log(Param)
simulnofix<-simulnofix[order(simulnofix$Param),]

plot(simulnofix$Param,simulnofix$IntAut,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),pch=3,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion",main="Force de l'introgression Autosomale lors de captures mitochondriales,\n En fonction de la diff?rence de dispersion entre sexe (Femelle - Male)")
points(jitter(simulnofix$Param),simulnofix$IntAut10,col="green",pch=4,lwd=2,type="b")
points(simulnofix$Param,simulnofix$FixAut,col="orange",pch=5,lwd=2,type="b")
points(simulnofix$Param,simulnofix$MeanExoAut,lwd=2,pch=6,type="b")
legend(x="topleft",legend=c("Proportion moyenne de loci autosomaux introgressés","Proportion moyenne de loci autosomaux dont plus de 10% des copies sont introgressées","Proportion de loci autosomaux capturés","Proportion moyenne de copies introgressés"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)

plot(y=simulnofix$MaxMtI,x=simulnofix$Param)
