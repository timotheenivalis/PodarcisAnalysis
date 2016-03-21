###################################################################
source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
setwd("C:/Users/Timothée/Documents/Studies/PodarcisCBGP/Renew simulations/lmMtA/")
Param<-matrix(data=c(1,0.9975,0.995,0.9925,0.99,0.975,0.95,0.925,0.9,0.8,0.7),nrow=11,ncol=1)


Simul<-"MtA0"
outR<-mtfunc(Simul)
results<-outR$results
resultsnofix<-outR$resultsNoFix
DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
DistriMtMaxDF$Param<-Param[0+1]
write.table(results,file="lmMtAresults.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(resultsnofix,file="lmMtAresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
write.table(DistriMtMaxDF,file="lmMtADistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)

plot(outR$DistribAutoNoFix$MeanExo[seq(from = 1,to = length(outR$DistribAutoNoFix$MeanExo),by = 2)],x=outR$DistribMtnoFix[seq(from = 1,to = length(outR$DistribAutoNoFix$MeanExo),by = 2)])

for (nbsimul in 1:10)
{
  Simul<-paste("MtA",nbsimul,sep="")
  outR<-mtfunc(Simul)
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
  names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
  DistriMtMaxDF$Param<-Param[nbsimul+1]
  write.table(results,file="lmMtAresults.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(resultsnofix,file="lmMtAresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
  write.table(DistriMtMaxDF,file="lmMtADistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
}

ResultsSim<-read.table("lmMtAresults.txt",header=T)
MtMaxima<-read.table("lmMtADistriMtMax.txt",header=T)
plot(MtMaxima$DistriMtMax,MtMaxima$DistriNuMtMax,pch=16,col=rgb(red = exp(MtMaxima$Param),0,0,1,maxColorValue = exp(1)))

##Now reading Arlequin
ArlSum<-data.frame(Param,ResultsSim$FixMt,ResultsSim$MtExo*2)
names(ArlSum)<-c("Param","FixMt","CopiesMtPop2")
ArlSum$PowSimulPop1<-NA
ArlSum$PowProbSD0Pop1<-NA
ArlSum$PowSimulPop2<-NA
ArlSum$PowProbSD0Pop2<-NA
ArlSum$ThetaSPop1<-NA
ArlSum$ThetaSPop2<-NA
ArlSum$ThetaPiPop1<-NA
ArlSum$ThetaPiPop2<-NA
ArlSum$DPop1<-NA
ArlSum$DPop2<-NA

for (i in 0:10)
{
  Sim<-paste("MtA",i,sep="")
  tajima<-read.table(file=paste(Sim,"/tajima.sum",sep=""),skip=6,header=F)
  tajima<-tajima[,-1]
  names(tajima)<-c("pop_name","Ncopies","NpolymSites","theta_S","theta_pi","TajimaD","PvalueBetaApprox","Pvaluesimul","ProbSimulatedD0","NumSimulations")  
  tajima<-tajima[which(substr(x=tajima$pop_name,start=1,stop=2)=="Gr"),]
  tajima$pop_name<-as.character(tajima$pop_name)
  ArlSum$PowSimulPop1[i+1]<-length(which(tajima$Pvaluesimul[which(tajima$pop_name=="Group1")]<0.05))/length(tajima$Pvaluesimul[which(tajima$pop_name=="Group1")])
  ArlSum$PowProbSD0Pop1[i+1]<-length(which(tajima$ProbSimulatedD0[which(tajima$pop_name=="Group1")]<0.05))/length(tajima$ProbSimulatedD0[which(tajima$pop_name=="Group1")])
  ArlSum$PowSimulPop2[i+1]<-length(which(tajima$Pvaluesimul[which(tajima$pop_name=="Group2")]<0.05))/length(tajima$Pvaluesimul[which(tajima$pop_name=="Group2")])
  ArlSum$PowProbSD0Pop2[i+1]<-length(which(tajima$ProbSimulatedD0[which(tajima$pop_name=="Group2")]<0.05))/length(tajima$ProbSimulatedD0[which(tajima$pop_name=="Group2")])
  ArlSum$ThetaSPop1[i+1]<-mean(tajima$theta_S[which(tajima$pop_name=="Group1")])
  ArlSum$ThetaSPop2[i+1]<-mean(tajima$theta_S[which(tajima$pop_name=="Group2")])
  ArlSum$ThetaPiPop1[i+1]<-mean(tajima$theta_pi[which(tajima$pop_name=="Group1")])
  ArlSum$ThetaPiPop2[i+1]<-mean(tajima$theta_pi[which(tajima$pop_name=="Group2")])
  ArlSum$DPop1[i+1]<-mean(tajima$TajimaD[which(tajima$pop_name=="Group1")])
  ArlSum$DPop2[i+1]<-mean(tajima$TajimaD[which(tajima$pop_name=="Group2")])
}

FuSum<-data.frame(Param,ResultsSim$FixMt,ResultsSim$MtExo*2)
names(FuSum)<-c("Param","FixMt","CopiesMtPop2")
FuSum$ProbSimFSPop1<-NA
FuSum$ProbSinfo0Pop1<-NA
FuSum$ProbSimFSPop2<-NA
FuSum$ProbSinfo0Pop2<-NA
FuSum$thetaPop1<-NA
FuSum$thetaPop2<-NA
FuSum$FSPop1<-NA
FuSum$FSPop2<-NA
for (i in 0:10)
{
  Sim<-paste("MtA",i,sep="")
  fu<-read.table(file=paste(Sim,"/fu_fs.sum",sep=""),skip=6,header=F)
  fu<-fu[,-1]
  names(fu)<-c("pop_name","Ncopies","NHaplo","theta","FS","ProbsimFS","PFSinf0","Nsimul")
  fu<-fu[which(substr(x=fu$pop_name,start=1,stop=2)=="Gr"),]                 
  fu$FS[which(fu$FS=="340282346638528859800000000000000000000.000000")]<-NA
  fu$FS[which(fu$FS>1000)]<-NA
  fu$FS<-as.numeric(as.character(fu$FS))
  FuSum$ProbSimFSPop1[i+1]<-length(which(fu$ProbsimFS[which(fu$pop_name=="Group1")]<0.05))/length(fu$ProbsimFS[which(fu$pop_name=="Group1")])
  FuSum$ProbSinfo0Pop1[i+1]<-length(which(fu$PFSinf0[which(fu$pop_name=="Group1")]<0.05))/length(fu$PFSinf0[which(fu$pop_name=="Group1")])
  FuSum$ProbSimFSPop2[i+1]<-length(which(fu$ProbsimFS[which(fu$pop_name=="Group2")]<0.05))/length(fu$ProbsimFS[which(fu$pop_name=="Group2")])
  FuSum$ProbSinfo0Pop2[i+1]<-length(which(fu$PFSinf0[which(fu$pop_name=="Group2")]<0.05))/length(fu$PFSinf0[which(fu$pop_name=="Group2")])
  FuSum$thetaPop1[i+1]<-mean(fu$theta[which(fu$pop_name=="Group1")],na.rm=T)
  FuSum$thetaPop2[i+1]<-mean(fu$theta[which(fu$pop_name=="Group2")],na.rm=T)
  FuSum$FSPop1[i+1]<-mean(fu$FS[which(fu$pop_name=="Group1")],na.rm=T)
  FuSum$FSPop2[i+1]<-mean(fu$FS[which(fu$pop_name=="Group2")],na.rm=T)
}

ArlSum<-cbind(ArlSum,FuSum[,-c(1,2)])

write.table(x=ArlSum,file="MtSelectionTests",quote=F,row.names=F)
