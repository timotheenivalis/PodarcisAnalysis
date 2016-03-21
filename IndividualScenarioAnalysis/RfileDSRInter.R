source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
setwd("C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/DSRInter/")

Simul<-paste("DSRInter1")
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,nrow(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
results$Param<-1
resultsnofix$Param<-1
DistriMtMaxDF$Param<-1


write.table(results,file="DSRIresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="DSRIresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="DSRIDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)
results
plot(DistriMtMaxDF$DistriMtMax-DistriMtMaxDF$DistriNuMtMax)
