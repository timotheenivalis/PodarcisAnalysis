source(file="C:/Users/Timothée/Dropbox/PodarcisCBGP/New simulations/FunctionsSimul.R")
setwd("C:/Users/Timothée/Dropbox/PodarcisCBGP/Renew simulations/lmBC/")

Param<-matrix(data=c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01,0.001),nrow=11,ncol=1)


Simul<-"lmBC0"
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
DistriMtMaxDF$Param<-Param[0+1]
write.table(results,file="lmBCresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="lmBCresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="lmBCDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)


for (nbsimul in 1:10)
{
  Simul<-paste("lmBC",nbsimul,sep="")
  outR<-mtfunc(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
  names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
  DistriMtMaxDF$Param<-Param[nbsimul+1]
  write.table(results,file="lmBCresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="lmBCresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)  
  write.table(DistriMtMaxDF,file="lmBCDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}

MtMaxima<-read.table("lmBCDistriMtMax.txt",header=T)
plot(MtMaxima$DistriMtMax,MtMaxima$DistriNuMtMax,pch=16,
     col=rgb(red = as.integer(unique(MtMaxima$Simul)),0,max(as.integer(unique(MtMaxima$Simul)))-as.integer(unique(MtMaxima$Simul)),max(as.integer(unique(MtMaxima$Simul))),maxColorValue =max(as.integer(unique(MtMaxima$Simul)))))
abline(a = 0,b = 1,col="red")
m0<-lm(DistriNuMtMax~DistriMtMax*I(log(Param)),data = MtMaxima)
summary(m0)
MtMaxima$Disc<-MtMaxima$DistriMtMax/MtMaxima$DistriNuMtMax
plot(x=log(MtMaxima$Param),y=MtMaxima$Disc)
boxplot(log(MtMaxima$Disc+1) ~ log(MtMaxima$Param))
quantile(MtMaxima$Disc,probs = c(0.75,0.9,0.925,0.95),na.rm=T)

simul<-read.table(file="lmBCresults.txt",header=T)
simulnofix<-read.table(file="lmBCresultsnofix.txt",header=T)

Param<-matrix(data=c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01,0.001),nrow=11,ncol=1)
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


plot(simulnofix$Param,simulnofix$FstAut,col="dark green",type="p",ylim=c(0,1),pch=1,lwd=2,xlab=c(expression(paste("D",Delta,sigma^2))),ylab="Fst entre les deux habitats",main="Fst en fonction de la différence de dispersion entre sexe (Femelle - Male)")
points(simulnofix$Param,simulnofix$FstZ,col="green",pch=2,lwd=2)
points(jitter(simulnofix$Param),simulnofix$FstW,col="blue",pch=3,lwd=2)
points(jitter(simulnofix$Param),simulnofix$FstMt,col="red",pch=4,lwd=2)
legend(x="topleft",legend=c("Autosomes","Z","W","Mt"),col=c("dark green","green","blue","red"),pt.lwd=2,pch=c(1,2,3,4))
