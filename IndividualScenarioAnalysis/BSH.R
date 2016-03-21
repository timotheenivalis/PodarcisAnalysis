source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
#setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBD/")
setwd("C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/BSH/")

Param<-matrix(data=round(c(0.2/0.4,0.05/0.55,0.4/0.2,0.55/0.05,0.001/0.559,0.559/0.001),3),nrow=6,ncol=1)

nbsimul<-0
Simul<-paste("BSH",nbsimul,sep="")
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
results$Param<-Param[nbsimul+1]
resultsnofix$Param<-Param[nbsimul+1]
DistriMtMaxDF$Param<-Param[nbsimul+1]


write.table(results,file="BSHresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="BSHresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="BSHDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)

for (nbsimul in 1:5)
{
  Simul<-paste("BSH",nbsimul,sep="")
  outR<-mtfunc(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  results$Param<-Param[nbsimul+1]
  resultsnofix$Param<-Param[nbsimul+1]
  write.table(results,file="BSHresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="BSHresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(cbind(rep(Simul,nrow(outR$DistriMtMax)),
                    outR$DistriMtMax,rep(Param[nbsimul+1],nrow(outR$DistriMtMax))),
              file="BSHDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}

##############################################################################################
#############################################################################################
par(mfrow=c(1,1))
simul<-read.table(file="BSHresults.txt",header=T)
simulnofix<-read.table(file="BSHresultsnofix.txt",header=T)
DistriMtMax<-read.table(file="BSHDistriMtMax.txt",header=T)
simul$Param<-log(Param)
simul<-simul[order(simul$Param),]

DistriMtMax$Param<-log(DistriMtMax$Param)
DistriMtMax<-DistriMtMax[order(DistriMtMax$Param),]
boxplot(DistriMtMax~Param, at=simul$Param, data=DistriMtMax,add=T,axes=F,border=gray(level=0.3)
        ,boxwex=0.3,boxlwd=1.5,whisklwd=1.5,staplelwd=1.5,outlwd=1.5,medlwd=1.5)

plot(simul$Param,simul$FixMt,col="black",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=21,bg="black",lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n En fonction du logarithme du rapport de dispersion entre sexe (Femelle/Male)",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(simul$Param,simul$IntMt,col="blue",pch=2,lwd=2,type="b")

