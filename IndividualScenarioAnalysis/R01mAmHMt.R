source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
setwd("C:/Users/Timothée//Documents/Studies/PodarcisCBGP/New simulations/R01mAmHMt/")

Simuls<-"R01mAmHMt0"
for (nbsimul in 1:17)
{
  Simuls<-c(Simuls,paste("R01mAmHMt",nbsimul,sep=""))
}

WrapMtFunc(Simuls = Simuls,FilePrefix = "R01mAmHMt",Param = rep(x = 1,times=length(Simuls)))


# 
# outR<-mtfunc(Simul)
# results<-outR$results
# resultsnofix<-outR$resultsNoFix
# DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
# names(DistriMtMaxDF)<-c("Simul","DistriMtMax")
# results$Param<-1
# resultsnofix$Param<-1
# DistriMtMaxDF$Param<-1
# 
# 
# write.table(results,file="R01mAmHMtresults.txt", sep="\t",quote=F,append=F,col.names=T)
# write.table(resultsnofix,file="R01mAmHMtresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
# write.table(DistriMtMaxDF,file="R01mAmHMtDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)
# 
# 
# for (nbsimul in 1:17)
# {
#   Simul<-paste("R01mAmHMt",nbsimul,sep="")
#   outR<-mtfunc(Simul)
#   results<-outR$results
#   resultsnofix<-outR$resultsNoFix
#   DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
#   names(DistriMtMaxDF)<-c("Simul","DistriMtMax")
#   results$Param<-1
#   resultsnofix$Param<-1
#   DistriMtMaxDF$Param<-1
#   
#   write.table(results,file="R01mAmHMtresults.txt", sep="\t",quote=F,append=T,col.names=F)
#   write.table(resultsnofix,file="R01mAmHMtresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
#   write.table(DistriMtMaxDF,file="R01mAmHMtDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
# }
