###################################################################
source(file="C:/Users/Timothée/Dropbox/PodarcisCBGP/New simulations/FunctionsSimul.R")
setwd("C:/Users/Timothée/Dropbox/PodarcisCBGP/Renew simulations/lmBH/")

Param<-matrix(data=c(0.2/0.4,0.05/0.55,0.4/0.2,0.55/0.05,0.001/0.599,0.599/0.001),nrow=6,ncol=1)


Simul<-"lmBH0"
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
DistriMtMaxDF$Param<-Param[0+1]
write.table(results,file="lmBHresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="lmBHresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="lmBHDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)


for (nbsimul in 1:5)
{
  Simul<-paste("lmBH",nbsimul,sep="")
  outR<-mtfunc(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
  names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
  DistriMtMaxDF$Param<-Param[nbsimul+1]
  write.table(results,file="lmBHresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="lmBHresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(DistriMtMaxDF,file="lmBHDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}

MtMaxima<-read.table("lmBHDistriMtMax.txt",header=T)
MtMaxima$Disc<-MtMaxima$DistriMtMax/MtMaxima$DistriNuMtMax
plot(x=log(MtMaxima$Param),y=MtMaxima$Disc)
boxplot(log(MtMaxima$Disc+1) ~ log(MtMaxima$Param))

plot(MtMaxima$DistriMtMax,MtMaxima$DistriNuMtMax,pch=16,
     col=rgb(red = as.integer(unique(MtMaxima$Simul)),0,max(as.integer(unique(MtMaxima$Simul)))-as.integer(unique(MtMaxima$Simul)),max(as.integer(unique(MtMaxima$Simul))),maxColorValue =max(as.integer(unique(MtMaxima$Simul)))))

abline(a = 0,b = 1,col="red")
m0<-lm(DistriNuMtMax~DistriMtMax*I(log(Param)),data = MtMaxima)
summary(m0)



simul<-read.table(file="lmBHresults.txt",header=T)
simulnofix<-read.table(file="lmBHresultsnofix.txt",header=T)
Param<-c(0.5,0.09090909,2,11,0.001788909,559)

simul$Param<-log(Param)
simul<-simul[order(simul$Param),]
plot(simul$Param,simul$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab="log du rapport des fitness hybrides female/male",ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n En fonction du logarithme du rapport de fitness hybride entre sexe (Femelle/Male)",cex.lab=1.8,cex.axis=1.5,cex.main=1.3)
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
