source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
#setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBD/")
setwd("C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/InvasionTotale/")

Param<-matrix(data=c(1,0.15/0.05,0.19/0.01,0.199/0.001,0.05/0.15,0.01/0.19,0.001/0.199),nrow=7,ncol=1)

wanted<-c(1:4,13:15)


nbsimul<-wanted[1]
Simul<-paste("Hun",nbsimul,sep="")
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
results$Param<-Param[1]
resultsnofix$Param<-Param[1]
DistriMtMaxDF$Param<-Param[1]


write.table(results,file="Hunresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="Hunresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="HunDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)

for (i in 2:length(wanted))
{
  nbsimul<-wanted[i]
  Simul<-paste("Hun",nbsimul,sep="")
  outR<-mtfunc(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  results$Param<-Param[i]
  resultsnofix$Param<-Param[i]
  write.table(results,file="Hunresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="Hunresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(cbind(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax,rep(Param[i],nrow(outR$DistriMtMax))),file="HunDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}


##############################################################################################
#############################################################################################
par(mfrow=c(1,1))
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/InvasionTotale/")
simul<-read.table(file="Hunresults.txt",header=T)
Param<-matrix(data=c(1:21),nrow=21,ncol=1)
simul$Param<-Param
plot(simul$Param,simul$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab="sélection mt",ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n En fonction de la sélecion mt",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(simul$Param,simul$IntMt,col="blue",pch=2,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Simulations avec Introgression Mitochondriale"),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)

plot(simul$Param,simul$FixMt,col="dark blue",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="",main="Proportion de simulations avec Capture ,\n En fonction du logarithme du rapport de dispersion entre sexe (Femelle/Male)",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(jitter(simul$Param),simul$IntAut10,col="light green",pch=4,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Proportion de loci nucl?aires introgress?s"),col=c("dark blue","light green"),pch=c(1,4),ncol=1,pt.lwd=2,cex=1.6)


plot(simul$Param,simul$IntAut,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),pch=3,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion",main="Force de l'introgression Autosomale lors de captures mitochondriales,\n En fonction de la diff?rence de dispersion entre sexe (Femelle - Male)")
points(jitter(simul$Param),simul$IntAut10,col="green",pch=4,lwd=2,type="b")
points(simul$Param,simul$FixAut,col="orange",pch=5,lwd=2,type="b")
points(simul$Param,simul$MeanExoAut,lwd=2,pch=6,type="b")
legend(x="topleft",legend=c("Proportion moyenne de loci autosomaux introgress?s","Proportion moyenne de loci autosomaux dont plus de 10% des copies sont introgress?es","Proportion de loci autosomaux captur?s","Proportion moyenne de copies introgresses"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)