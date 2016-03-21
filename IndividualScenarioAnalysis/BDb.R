source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
#setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBD/")
setwd("C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/BDb/")

Param<-matrix(data=c(1,0.15/0.05,0.05/0.15,0.19/0.01,0.01/0.19,0.199/0.001,0.001/0.199),nrow=7,ncol=1)

nbsimul<-0
Simul<-paste("BDb",nbsimul,sep="")
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
results$Param<-Param[nbsimul+1]
resultsnofix$Param<-Param[nbsimul+1]
DistriMtMaxDF$Param<-Param[nbsimul+1]


write.table(results,file="BDbresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="BDbresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="BDbDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)

for (nbsimul in 1:6)
{
  Simul<-paste("BDb",nbsimul,sep="")
  outR<-mtfunc(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  results$Param<-Param[nbsimul+1]
  resultsnofix$Param<-Param[nbsimul+1]
  write.table(results,file="BDbresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="BDbresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(cbind(rep(Simul,nrow(outR$DistriMtMax)),
                    outR$DistriMtMax,rep(Param[nbsimul+1],nrow(outR$DistriMtMax))),
              file="BDbDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}


##############################################################################################
#############################################################################################
par(mfrow=c(1,1))
simul<-read.table(file="BDbresults.txt",header=T)
simulnofix<-read.table(file="BDbresultsnofix.txt",header=T)
DistriMtMax<-read.table(file="BDbDistriMtMax.txt",header=T)
simul$Param<-log(Param)
simul<-simul[order(simul$Param),]

DistriMtMax$Param<-log(DistriMtMax$Param)
DistriMtMax<-DistriMtMax[order(DistriMtMax$Param),]
boxplot(DistriMtMax~Param, at=simul$Param, data=DistriMtMax,add=T,axes=F,border=gray(level=0.3)
        ,boxwex=0.3,boxlwd=1.5,whisklwd=1.5,staplelwd=1.5,outlwd=1.5,medlwd=1.5)

plot(simul$Param,simul$FixMt,col="black",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=21,bg="black",lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n En fonction du logarithme du rapport de dispersion entre sexe (Femelle/Male)",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(simul$Param,simul$IntMt,col="blue",pch=2,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Simulations avec Introgression Mitochondriale"),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)


simul<-read.table(file=paste(diri,"/BDbresults.txt",sep=""),header=T)
Param<-matrix(data=c(0.006,0.052,0.626,1,1.6,19.11,181.25),nrow=7,ncol=1)
simul$Param<-log(Param)
plot(simul$Param,simul$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n En fonction du logarithme du rapport de dispersion entre sexe (Femelle/Male)",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(simul$Param,simul$IntMt,col="blue",pch=2,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Simulations avec Introgression Mitochondriale"),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)

plot(simul$Param,simul$FixMt,col="dark blue",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="",main="Proportion de simulations avec Capture ,\n En fonction du logarithme du rapport de dispersion entre sexe (Femelle/Male)",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(jitter(simul$Param),simul$IntAut10,col="light green",pch=4,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Proportion de loci nucléaires introgressés"),col=c("dark blue","light green"),pch=c(1,4),ncol=1,pt.lwd=2,cex=1.6)


plot(simul$Param,simul$IntAut,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),pch=3,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion",main="Force de l'introgression Autosomale lors de captures mitochondriales,\n En fonction de la diff?rence de dispersion entre sexe (Femelle - Male)")
points(jitter(simul$Param),simul$IntAut10,col="green",pch=4,lwd=2,type="b")
points(simul$Param,simul$FixAut,col="orange",pch=5,lwd=2,type="b")
points(simul$Param,simul$MeanExoAut,lwd=2,pch=6,type="b")
legend(x="topleft",legend=c("Proportion moyenne de loci autosomaux introgressés","Proportion moyenne de loci autosomaux dont plus de 10% des copies sont introgressées","Proportion de loci autosomaux captur?s","Proportion moyenne de copies introgressés"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)

plot(simul$Param,simul$FstAut,col="dark green",type="p",ylim=c(0,1),pch=1,lwd=2,xlab=c(expression(paste("D",Delta,sigma^2))),ylab="Fst entre les deux habitats",main="Fst en fonction de la différence de dispersion entre sexe (Femelle - Male)")
points(simul$Param,simul$FstZ,col="green",pch=2,lwd=2)
points(jitter(simul$Param),simul$FstW,col="blue",pch=3,lwd=2)
points(jitter(simul$Param),simul$FstMt,col="red",pch=4,lwd=2)
legend(x="topleft",legend=c("Autosomes","Z","W","Mt"),col=c("dark green","green","blue","red"),pt.lwd=2,pch=c(1,2,3,4))
