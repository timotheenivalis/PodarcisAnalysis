#only BDWR0 is interesting here
source(file="C:/Users/Timoth�e/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
setwd("C:/Users/Timoth�e/Documents/Studies/PodarcisCBGP/New simulations/BDWR/")

Simul<-paste("BDWR0")
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
results$Param<-1
resultsnofix$Param<-1
DistriMtMaxDF$Param<-1


write.table(results,file="BDWRresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="BDWRresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="BDWRDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)
results




###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWR5")
Simul<-"BDWR5"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWRresults.txt", sep="\t",quote=F,append=F,col.names=T)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWR4")
Simul<-"BDWR4"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWRresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWR2")
Simul<-"BDWR2"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWRresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWR0")
Simul<-"BDWR0"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWRresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWR1")
Simul<-"BDWR1"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWRresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWR3")
Simul<-"BDWR3"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWRresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWR6")
Simul<-"BDWR6"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR/BDWRresults.txt", sep="\t",quote=F,append=T,col.names=F)


##############################################################################################
#############################################################################################
par(mfrow=c(1,1))
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BDWR")
simul<-read.table(file="BDWRresults.txt",header=T)
Param<-matrix(data=c(1:7),nrow=7,ncol=1)
simul$Param<-Param
plot(simul$Param,simul$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab="s�lection mt",ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n En fonction de la s�lecion mt",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
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
