
###################################################################
source(file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/FunctionsSimul.R")
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun")
Simul<-"gHun1"
outR<-mtfunc(Simul)
results<-outR$results
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHunresults.txt", sep="\t",quote=F,append=F,col.names=T)


for (nbsimul in 2:36)
{
  #if (nbsimul!=2 & nbsimul!=6 & nbsimul!=18)
    {
      Simul<-paste("gHun",nbsimul,sep="")
      outR<-mtfunc(Simul)
      results<-outR$results
      write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHunresults.txt", sep="\t",quote=F,append=T,col.names=F)
    }
}

##############################################################################################
#############################################################################################
par(mfrow=c(1,1))
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/")
simul<-read.table(file="gHunresults.txt",header=T)
Param<-matrix(data=c(1:dim(simul)[1]),nrow=dim(simul)[1],ncol=1)
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

plot(simul$Param,simul$FstAut,col="dark green",type="p",ylim=c(0,1),pch=1,lwd=2,xlab=c(expression(paste("D",Delta,sigma^2))),ylab="Fst entre les deux habitats",main="Fst en fonction de la diff?rence de dispersion entre sexe (Femelle - Male)")
points(simul$Param,simul$FstZ,col="green",pch=2,lwd=2)
points(jitter(simul$Param),simul$FstW,col="blue",pch=3,lwd=2)
points(jitter(simul$Param),simul$FstMt,col="red",pch=4,lwd=2)
legend(x="topleft",legend=c("Autosomes","Z","W","Mt"),col=c("dark green","green","blue","red"),pt.lwd=2,pch=c(1,2,3,4))


simul2<-read.table(file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/InvasionTotale/Hunresults.txt",header=T)
interesting<-simul2[c(1,13,14,15),]
interesting$ds2<-c(1,0.626,0.052,0.006)
interesting$G<-6000
interesting$r<-0.5
interesting$multiS<-F

simul<-read.table(file="gHunresults.txt",header=T)
simul$ds2<-rep(c(1,0.626,0.052,0.006),times=9)
simul$G<-c(rep(12000,times=4),rep(18000,times=4),rep(6000,times=4),rep(12000,times=4),rep(18000,times=4),rep(6000,times=16))
simul$r<-c(rep(0.5,times=8),rep(0.1,times=12),rep(0.01,times=4),rep(0.5,times=4),rep(0.1,times=4),rep(0.01,times=4))
simul$multiS<-c(rep(F,times=24),rep(T,times=12))

simult<-rbind(simul,interesting)

plot(simult$IntAut,x=simult$ds2)
library(ggplot2)
qplot(x=jitter(ds2),y=IntAut10,data=simult[ simult$multiS==F, ],facets=~G,color=as.factor(r),size=3,ylim=c(0,1),xlab="S² female/male")
qplot(x=ds2,y=IntAut,data=simult[ simult$multiS==F, ],facets=~G,color=as.factor(r),size=3,ylim=c(0,1),xlab="S² female/male")
qplot(x=ds2,y=FixAut,data=simult[ simult$multiS==F, ],facets=~G,color=as.factor(r),size=3,ylim=c(0,1),xlab="S² female/male")
qplot(x=ds2,y=MeanExoAut,data=simult[ simult$multiS==F, ],facets=~G,color=as.factor(r),size=3,ylim=c(0,1),xlab="S² female/male")

qplot(x=ds2,y=IntAut10,data=simult[ simult$multiS==T, ],facets=~G,color=as.factor(r),size=3,ylim=c(0,1),xlab="S² female/male")
qplot(x=jitter(ds2),y=IntAut,data=simult[ simult$multiS==T, ],facets=~G,color=as.factor(r),size=3,ylim=c(0,1),xlab="S² female/male")
qplot(x=ds2,y=FixAut,data=simult[ simult$multiS==T, ],facets=~G,color=as.factor(r),size=3,ylim=c(0,1),xlab="S² female/male")
qplot(x=ds2,y=MeanExoAut,data=simult[ simult$multiS==T, ],facets=~G,color=as.factor(r),size=3,ylim=c(0,1),xlab="S² female/male")

log(c(1,0.626,0.052,0.006))
