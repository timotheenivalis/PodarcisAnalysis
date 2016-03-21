###################################################################
source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/InvasionFunctionsSimul.R")
setwd("C:/Users/Timothée/Documents/Studies/PodarcisCBGP/Renew simulations/lmHunS/")
Param<-matrix(data=c(1,0.0015/0.0005,0.0005/0.0015,0.0019/0.0001,0.0001/0.0019,0.00199/0.00001,0.00001/0.00199),nrow=7,ncol=1)

Simul<-"lmHunS0"
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
DistriMtMaxDF$Param<-Param[0+1]
write.table(results,file="lmHunSresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="lmHunSresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="lmHunSDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)


for (nbsimul in 1:6)
{
  Simul<-paste("lmHunS",nbsimul,sep="")
  outR<-mtfunc(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
  names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
  DistriMtMaxDF$Param<-Param[nbsimul+1]
  write.table(results,file="lmHunSresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="lmHunSresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(DistriMtMaxDF,file="lmHunSDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}


Simul<-"lmHunS0"
outR<-mtfuncInvasion(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
results$Param<-Param[0+1]
resultsnofix$Param<-Param[0+1]
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
DistriMtMaxDF$Param<-Param[0+1]
write.table(results,file="lmHunSresultsProf.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="lmHunSresultsnofixProf.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="lmHunSDistriMtMaxProf.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)


for (i in 1:6)
{
  nbsimul<-i
  Simul<-paste("lmHunS",nbsimul,sep="")
  outR<-mtfuncInvasion(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  results$Param<-Param[i+1]
  resultsnofix$Param<-Param[i+1]
  DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
  names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
  DistriMtMaxDF$Param<-Param[i+1]
  write.table(results,file="lmHunSresultsProf.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="lmHunSresultsnofixProf.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(DistriMtMaxDF,file="lmHunSDistriMtMaxProf.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
  #write.table(cbind(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax,rep(Param[i+1],length(outR$DistriMtMax))),file="lmHunSDistriMtMaxProf.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}

simul<-read.table("lmHunSresultsProf.txt")
simul<-simul[order(simul$Param),]

prDMT<-read.table("lmHunSDistriMtMaxProf.txt",header=T)
plot(prDMT$DistriNuMtMax)

##################################"
### whole grid
##################################"
Simul<-"lmHunS0"
outR<-mtfuncInvasionAllGrid(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
results$Param<-Param[0+1]
resultsnofix$Param<-Param[0+1]
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
DistriMtMaxDF$Param<-Param[0+1]
write.table(results,file="lmHunSresultsProfGrid.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="lmHunSresultsnofixProfGrid.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="lmHunSDistriMtMaxProfGrid.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)


for (i in 1:6)
{
  nbsimul<-i
  Simul<-paste("lmHunS",nbsimul,sep="")
  outR<-mtfuncInvasionAllGrid(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  results$Param<-Param[i+1]
  resultsnofix$Param<-Param[i+1]
  DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
  names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
  DistriMtMaxDF$Param<-Param[i+1]
  write.table(results,file="lmHunSresultsProfGrid.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="lmHunSresultsnofixProfGrid.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(DistriMtMaxDF,file="lmHunSDistriMtMaxProfGrid.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}
