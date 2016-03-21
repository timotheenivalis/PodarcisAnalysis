source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/InvasionFunctionsSimul.R")
#setwd("C:/Users/Thimothee Admin/Dropbox/PodarcisCBGP/New simulations/InvasionTotale/")
setwd("C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/InvasionTotale/")

Param<-matrix(data=c(1,0.15/0.05,0.19/0.01,0.199/0.001,0.05/0.15,0.01/0.19,0.001/0.199),nrow=7,ncol=1)
Param<-matrix(data=c(1,0.19/0.01,0.199/0.001,0.01/0.19,0.001/0.199),nrow=5,ncol=1)
wanted<-c(1,3:4,14:15)# miss 2 and 13


nbsimul<-wanted[1]
Simul<-paste("Hun",nbsimul,sep="")
outR<-mtfuncInvasion(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
results$Param<-Param[1]
resultsnofix$Param<-Param[1]
DistriMtMaxDF$Param<-Param[1]

write.table(results,file="HunresultsProf.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="HunresultsnofixProf.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="HunDistriMtMaxProf.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)

for (i in 2:length(wanted))
{
  nbsimul<-wanted[i]
  Simul<-paste("Hun",nbsimul,sep="")
  outR<-mtfuncInvasion(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  results$Param<-Param[i]
  resultsnofix$Param<-Param[i]
  write.table(results,file="HunresultsProf.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="HunresultsnofixProf.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(cbind(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax,
                    rep(Param[i],nrow(outR$DistriMtMax))),
              file="HunDistriMtMaxProf.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}

##Now fetch the missing profiles in the old files Attila/Hun2 and Attila/
Simul<-"Attila/Hun2/"
outR<-mtfuncInvasion(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
results$Param<-0.05/0.15
resultsnofix$Param<-0.05/0.15
write.table(results,file="HunresultsProf.txt", sep="\t",quote=F,append=T,col.names=F)
write.table(resultsnofix,file="HunresultsnofixProf.txt", sep="\t",quote=F,append=T,col.names=F)
write.table(cbind(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax,
                  rep(0.05/0.15,nrow(outR$DistriMtMax))),
            file="HunDistriMtMaxProf.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)


simul<-read.table(file="HunresultsProf.txt",header=T)
simulnofix<-read.table(file="HunresultsnofixProf.txt",header=T)
DMT<-read.table(file="HunDistriMtMaxProf.txt",header=T)

simul2<-read.table(file="Hunresults.txt",header=T)
simulnofix2<-read.table(file="Hunresultsnofix.txt",header=T)
DMT2<-read.table(file="HunDistriMtMax.txt",header=T)
missingDMT<-DMT2[which(DMT2$Param==3),]

boxplot(DistriMtMax~Param,rbind(DMT,missingDMT))


simul<-rbind(simul,simul2[which(simul2$Param==3),])
simul<-simul[order(simul$Param),]
simulnofix<-rbind(simulnofix,simulnofix2[which(simul2$Param==3),])
simulnofix<-simulnofix[order(simulnofix$Param),]

simul["Hun2",]<-c(100,0.21,0.52,0.35,0.984,0.954,0.09,0.623,0.278,NA,NA,NA,NA,0.90,0.95,1,1,1,0.85,0.8775,0.95,1,1,0,0.05,0.05,0.2,0.2,0.5001,0.516,0.629,0.729,0.751,3)
simulnofix["Hun2",]<-c(100,0.79,0.31,0.18,0.984,0.954,0.101,0.609,0.278,NA,NA,NA,NA,0.45,0.75,0.95,1,1,0.45,0.65,0.9,1,1,0,0.05,0.2,0.5,0.85,0.33,0.44,0.6109,0.745,0.865,0.975,3)

write.table(x = simul,file = "HunresultsProf.txt",quote = F,col.names = T,append=F)
write.table(x = simulnofix,file = "HunresultsnofixProf.txt",quote = F,col.names = T,append=F)
write.table(x=rbind(DMT,missingDMT),file = "HunDistriMtMaxProf.txt",quote=F,col.names = T)


#########################################
# On all the grid
#########################################


nbsimul<-wanted[1]
Simul<-paste("Hun",nbsimul,sep="")
outR<-mtfuncInvasionAllGrid(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax[,1],outR$DistriMtMax[,2])
names(DistriMtMaxDF)<-c("Simul","DistriMtMax")
results$Param<-Param[1]
resultsnofix$Param<-Param[1]
DistriMtMaxDF$Param<-Param[1]

write.table(results,file="HunresultsProfGrid.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="HunresultsnofixProfGrid.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="HunDistriMtMaxProfGrid.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)

for (i in 2:length(wanted))
{
  nbsimul<-wanted[i]
  Simul<-paste("Hun",nbsimul,sep="")
  outR<-mtfuncInvasionAllGrid(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  results$Param<-Param[i]
  resultsnofix$Param<-Param[i]
  write.table(results,file="HunresultsProfGrid.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="HunresultsnofixProfGrid.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(cbind(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax[,1],outR$DistriMtMax[,2],rep(Param[i],length(outR$DistriMtMax))),file="HunDistriMtMaxProfGrid.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}

##Now fetch the missing profiles in the old files Attila/Hun2 and Attila/
Simul<-"Attila/Hun2/"
outR<-mtfuncInvasionAllGrid(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
results$Param<-0.05/0.15
resultsnofix$Param<-0.05/0.15
write.table(results,file="HunresultsProfGrid.txt", sep="\t",quote=F,append=T,col.names=F)
write.table(resultsnofix,file="HunresultsnofixProfGrid.txt", sep="\t",quote=F,append=T,col.names=F)
write.table(cbind(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax[,1],outR$DistriMtMax[,2],rep(0.05/0.15,length(outR$DistriMtMax))),file="HunDistriMtMaxProfGrid.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)


simul<-read.table(file="HunresultsProfGrid.txt",header=T)
simulnofix<-read.table(file="HunresultsnofixProfGrid.txt",header=T)
DMT<-read.table(file="HunDistriMtMaxProfGrid.txt",header=T)

boxplot(DistriMtMax~Param,DMT)
simul2<-read.table(file="Hunresults.txt",header=T)
simulnofix2<-read.table(file="Hunresultsnofix.txt",header=T)
DMT2<-read.table(file="HunDistriMtMax.txt",header=T)
missingDMT<-DMT2[which(DMT2$Param==3),]

names(DMT)<- names(missingDMT)
boxplot(DistriMtMax~Param,rbind(DMT,missingDMT))

simul<-rbind(simul,simul2[which(simul2$Param==3),])
simul<-simul[order(simul$Param),]
simulnofix<-rbind(simulnofix,simulnofix2[which(simul2$Param==3),])
simulnofix<-simulnofix[order(simulnofix$Param),]

simul["Hun2",]<-c(100,0.16,0.52,0.39,0.984,0.954,0.09,0.619,0.278,NA,NA,NA,NA,0.90,0.95,1,1,1,0.85,0.8775,0.95,1,1,0,0.05,0.05,0.2,0.2,0.5001,0.516,0.629,0.729,0.751,3)
simulnofix["Hun2",]<-c(100,0.79,0.31,0.18,0.984,0.954,0.101,0.609,0.278,NA,NA,NA,NA,0.45,0.75,0.95,1,1,0.45,0.65,0.9,1,1,0,0.05,0.2,0.5,0.85,0.33,0.44,0.6109,0.745,0.865,0.975,3)

write.table(x = simul,file = "HunresultsProfGrid.txt",quote = F,col.names = T,append=F)
write.table(x = simulnofix,file = "HunresultsnofixProfGrid.txt",quote = F,col.names = T,append=F)
write.table(x=rbind(DMT,missingDMT),file = "HunDistriMtMaxProfGrid.txt",quote=F,col.names = T)


##############################################################################################
#############################################################################################
# par(mfrow=c(1,1))
# setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/InvasionTotale/")
# simul<-read.table(file="Hunresults.txt",header=T)
# Param<-matrix(data=c(1:21),nrow=21,ncol=1)
# simul$Param<-Param
# simul<-simul[order(simul$Param),]
# plot(log(simul$Param),simul$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab="sélection mt",ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n En fonction de la sélecion mt",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
# points(log(simul$Param),simul$IntMt,col="blue",pch=2,lwd=2,type="b")
# legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Simulations avec Introgression Mitochondriale"),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)
# 
# plot(log(simul$Param),simul$FixMt,col="dark blue",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="",main="Proportion de simulations avec Capture ,\n En fonction du logarithme du rapport de dispersion entre sexe (Femelle/Male)",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
# points(log(jitter(simul$Param)),simul$IntAut10,col="light green",pch=4,lwd=2,type="b")
# legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Proportion de loci nucl?aires introgress?s"),col=c("dark blue","light green"),pch=c(1,4),ncol=1,pt.lwd=2,cex=1.6)
# 
# 
# plot(simul$Param,simul$IntAut,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),pch=3,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion",main="Force de l'introgression Autosomale lors de captures mitochondriales,\n En fonction de la diff?rence de dispersion entre sexe (Femelle - Male)")
# points(jitter(simul$Param),simul$IntAut10,col="green",pch=4,lwd=2,type="b")
# points(simul$Param,simul$FixAut,col="orange",pch=5,lwd=2,type="b")
# points(simul$Param,simul$MeanExoAut,lwd=2,pch=6,type="b")
# legend(x="topleft",legend=c("Proportion moyenne de loci autosomaux introgress?s","Proportion moyenne de loci autosomaux dont plus de 10% des copies sont introgress?es","Proportion de loci autosomaux captur?s","Proportion moyenne de copies introgresses"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)
