source(file="D:/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
#setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBD/")
setwd("D:/Documents/Studies/PodarcisCBGP/New simulations/Homogamy/")

Simuls<-c("H0","H1","H6","H7","HB6")

Param<-matrix(data=c(0.5,0.99,0.1,0.3,0.7,0),nrow=6,ncol=1)

Simul<-Simuls[1]
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
results$Param<-Param[1]
resultsnofix$Param<-Param[1]
DistriMtMaxDF$Param<-Param[1]

write.table(results,file="Homogamyresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="Homogamyresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="HomogamyDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)

for (i in 2:5)
{
  Simul<-Simuls[i]
  outR<-mtfunc(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  results$Param<-Param[i]
  resultsnofix$Param<-Param[i]
  write.table(results,file="Homogamyresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="Homogamyresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(cbind(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax,rep(Param[i],nrow(outR$DistriMtMax))),file="HomogamyDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}

simul<-read.table(file="Homogamyresults.txt",header=T)
simulnofix<-read.table(file="Homogamyresultsnofix.txt",header=T)
DistriMtMax<-read.table(file="HomogamyDistriMtMax.txt",header=T)
simul<-simul[order(simul$Param),]

DistriMtMax<-DistriMtMax[order(DistriMtMax$Param),]
boxplot(DistriMtMax~Param, at=simul$Param, data=DistriMtMax,add=T,axes=F,border=gray(level=0.3)
        ,boxwex=0.03,boxlwd=1.5,whisklwd=1.5,staplelwd=1.5,outlwd=1.5,medlwd=1.5)

plot(simul$Param,simul$FixMt,col="black",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=21,bg="black",lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n En fonction du logarithme du rapport de dispersion entre sexe (Femelle/Male)",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)


###direction introg
dirInt<-vector(length=5)
dirCapt<-vector(length=5)
for (i in 1:5)
{
  Simul<-Simuls[i]
  dirwork<-paste(getwd(),Simul,sep="/")
  results<-as.data.frame(matrix(nrow=1,ncol=13,dimnames=list(Simul,c("NbRun","FixMt","IntMt","MtExo","IntAut","IntAut10","FixAut","MeanExoAut","SDExoAut","FstAut","FstZ","FstW","FstMt"))))
  data<-read.table(paste(dirwork,"IntrogStats.txt",sep="/"),header=T) 
  dirCapt[i]<-mean(data$Sp[which(data$Mt==1)])
  dirInt[i]<-mean(data$Sp*data$Mt)
  }
plot(Param,c(1-dirInt,0.5))
plot(Param,c(1-dirCapt,0.5))
