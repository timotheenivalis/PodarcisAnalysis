source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
setwd("C:/Users/Timothée/Documents/Studies/PodarcisCBGP/Renew simulations/LowI/")


#ds2<-c(4.89,4.34,3.8,3.26,2.72,2.17,1.63,1.09,0.54,0.27,0.05)
migrate<-c(0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01,0.005,0.001)

Simul<-"LI0"
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
results$Param<-migrate[0+1]
resultsnofix$Param<-migrate[0+1]
DistriMtMaxDF$Param<-migrate[0+1]

write.table(results,file="LowIresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="LowIresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="LowIDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)

for (nbsimul in 1:10)
{
  Simul<-paste("LI",nbsimul,sep="")
  outR<-mtfunc(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
  names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
  results$Param<-migrate[nbsimul+1]
  resultsnofix$Param<-migrate[nbsimul+1]
  DistriMtMaxDF$Param<-migrate[nbsimul+1]
  
  write.table(results,file="LowIresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="LowIresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(DistriMtMaxDF,file="LowIDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}
simul<-read.table(file="LowIresults.txt",header=T)
simulNF<-read.table(file="LowIresultsnofix.txt",header=T)

simul$ds2<-ds2
simulNF$ds2<-ds2


plot(simul$ds2,simul$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab="ds2",ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n ",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(simul$ds2,simul$IntMt,col="blue",pch=2,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Simulations avec Introgression Mitochondriale"),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)

plot(simulNF$ds2,simulNF$NoFixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab="ds2",ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n ",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(simulNF$ds2,simulNF$IntMt,col="blue",pch=2,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Simulations avec Introgression Mitochondriale"),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)


plot(simul$ds2,simul$FixMt,col="dark blue",type="b",yaxp=c(0,1,10),ylim=c(0,1.2),pch=1,lwd=2,xlab="ds2",ylab="",main="Proportion de simulations avec Capture Mitochondriale, et parmi elles, \n proportion d'Introgression Nucléaire",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(jitter(simul$ds2),simul$IntAut10,col="light green",pch=4,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Proportion de loci nucléaires introgressés"),col=c("dark blue","light green"),pch=c(1,4),ncol=1,pt.lwd=2,cex=1.2)

plot(simul$ds2,simul$IntAut,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),pch=3,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion",main="Force de l'introgression Autosomale lors de captures mitochondriales,\n En fonction de la diff?rence de dispersion entre sexe (Femelle - Male)")
points(jitter(simul$ds2),simul$IntAut10,col="green",pch=4,lwd=2,type="b")
points(simul$ds2,simul$FixAut,col="orange",pch=5,lwd=2,type="b")
points(simul$ds2,simul$MeanExoAut,lwd=2,pch=6,type="b")
points(simul$ds2,simul$MtExo,lwd=2,pch=7,type="b",col="red")
legend(x="topleft",legend=c("Proportion moyenne de loci autosomaux introgress?s","Proportion moyenne de loci autosomaux dont plus de 10% des copies sont introgress?es","Proportion de loci autosomaux captur?s","Proportion moyenne de copies introgresses"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)

plot(simulNF$ds2,simulNF$IntAut,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),pch=3,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion",main="Force de l'introgression Autosomale lors de captures mitochondriales,\n En fonction de la diff?rence de dispersion entre sexe (Femelle - Male)")
points(jitter(simulNF$ds2),simulNF$IntAut10,col="green",pch=4,lwd=2,type="b")
points(simulNF$ds2,simulNF$FixAut,col="orange",pch=5,lwd=2,type="b")
points(simulNF$ds2,simulNF$MeanExoAut,lwd=2,pch=6,type="b")
points(simulNF$ds2,simulNF$MtExo,lwd=2,pch=7,type="b",col="red")
legend(x="topleft",legend=c("Proportion moyenne de loci autosomaux introgress?s","Proportion moyenne de loci autosomaux dont plus de 10% des copies sont introgress?es","Proportion de loci autosomaux captur?s","Proportion moyenne de copies introgresses"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)

plot(simul$ds2,simul$FstAut,col="dark green",type="p",ylim=c(0,1),pch=1,lwd=2,xlab="sélection mt",ylab="Fst entre les deux habitats",main="Fst en fonction de la sélection mt")
points(simul$ds2,simul$FstZ,col="green",pch=2,lwd=2)
points(jitter(simul$ds2),simul$FstW,col="blue",pch=3,lwd=2)
points(jitter(simul$ds2),simul$FstMt,col="red",pch=4,lwd=2)
legend(x="topleft",legend=c("Autosomes","Z","W","Mt"),col=c("dark green","green","blue","red"),pt.lwd=2,pch=c(1,2,3,4))

plot(simulNF$ds2,simulNF$FstAut,col="dark green",type="p",ylim=c(0,1),pch=1,lwd=2,xlab="sélection mt",ylab="Fst entre les deux habitats",main="Fst en fonction de la sélection mt")
points(simulNF$ds2,simulNF$FstZ,col="green",pch=2,lwd=2)
points(jitter(simulNF$ds2),simulNF$FstW,col="blue",pch=3,lwd=2)
points(jitter(simulNF$ds2),simulNF$FstMt,col="red",pch=4,lwd=2)
legend(x="topleft",legend=c("Autosomes","Z","W","Mt"),col=c("dark green","green","blue","red"),pt.lwd=2,pch=c(1,2,3,4))

