###################################################################
source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
setwd("C:/Users/Timothée/Documents/Studies/PodarcisCBGP/Renew simulations/lmBHCD/")

Simul<-"lmBHCD0"
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
write.table(results,file="lmBHCDresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="lmBHCDresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="lmBHCDDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)


for (nbsimul in 1:17)
{
  Simul<-paste("lmBHCD",nbsimul,sep="")
  outR<-mtfunc(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
  names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
  write.table(results,file="lmBHCDresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="lmBHCDresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(DistriMtMaxDF,file="lmBHCDDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}
