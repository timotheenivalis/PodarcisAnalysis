source(file="D:/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
#setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBD/")
setwd("D:/Documents/Studies/PodarcisCBGP/New simulations/HomogamyMulti/")

Simuls<-c("HM0","HM1","HM6","HM7","HM8","HBM2","HBM4","HBM5")

Param<-matrix(data=1-c(0.5 ,0.99,0.1 ,0.3 ,0,0.3 ,0.5 ,0.7),nrow=8,ncol=1)


Simul<-Simuls[1]
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
results$Param<-Param[1]
resultsnofix$Param<-Param[1]
DistriMtMaxDF$Param<-Param[1]

write.table(results,file="HomogamyMultiresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="HomogamyMultiresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="HomogamyMultiDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)

for (i in 2:8)
{
  Simul<-Simuls[i]
  outR<-mtfunc(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  results$Param<-Param[i]
  resultsnofix$Param<-Param[i]
  write.table(results,file="HomogamyMultiresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="HomogamyMultiresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(cbind(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax,rep(Param[i],nrow(outR$DistriMtMax))),file="HomogamyMultiDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}

simul<-read.table(file="HomogamyMultiresults.txt",header=T)
simulnofix<-read.table(file="HomogamyMultiresultsnofix.txt",header=T)
DistriMtMax<-read.table(file="HomogamyMultiDistriMtMax.txt",header=T)
simul<-simul[order(simul$Param),]

DistriMtMax<-DistriMtMax[order(DistriMtMax$Param),]
boxplot(DistriMtMax~Param, at=1:9, data=DistriMtMax,add=F,axes=F,border=gray(level=0.3)
        ,boxwex=0.03,boxlwd=1.5,whisklwd=1.5,staplelwd=1.5,outlwd=1.5,medlwd=1.5)

plot(simul$Param,simul$FixMt,col="black",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=21,bg="black",lwd=2,xlab='',ylab="Proportion de simulations",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)

DistriMtMax$Disc <- DistriMtMax$DistriMtMax - DistriMtMax$DistriNuMtMax
plot(DistriMtMax$Disc, x=DistriMtMax$Param,ylim=c(-0.4,1))
abline(h=0.8)
