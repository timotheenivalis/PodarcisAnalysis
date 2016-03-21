mtfunc<-function(Simul)
{
  dirwork<-paste(getwd(),Simul,sep="/")
  results<-as.data.frame(matrix(nrow=1,ncol=13,dimnames=list(Simul,c("NbRun","FixMt","IntMt","MtExo","IntAut","IntAut10","FixAut","MeanExoAut","SDExoAut","FstAut","FstZ","FstW","FstMt"))))
  data<-read.table(paste(dirwork,"IntrogStats.txt",sep="/"),header=T) 
  fixedMt<-data[data$Mt==1,]

  distriMtMax<-tapply(X=data$Mt,INDEX=data$Run,FUN=function(x){max(x,na.rm=T)})
  
  MtCop<-tapply(X=data$Mt,INDEX=data$Run,FUN=function(x){which(x==max(x,na.rm=T))[1]})
  MtCopMax<-data$Mt[as.numeric(as.character(names(MtCop)))*2-2+MtCop]
  NuCop<-apply(X=data[,3:22],MARGIN=1,FUN=function(x){mean(x=x,na.rm=T)})[as.numeric(as.character(names(MtCop)))*2-2+MtCop]
  distriMtMax<-cbind(MtCopMax,NuCop)
  
  results$NbRun<-max(data$Run)
  results$FixMt<-length(fixedMt[,1])/max(data$Run)
  results$IntMt<-length(unique(data[data$Mt>0,1]))/(max(data$Run))
  results$MtExo<-mean(data$Mt)
  results$IntAut<-mean((apply(X=fixedMt[,3:22],MARGIN=1,FUN=function(x){length(x[x!=0])}))/20)
  results$IntAut10<-mean((apply(X=fixedMt[,3:22],MARGIN=1,FUN=function(x){length(x[x>=0.1])}))/20)
  results$FixAut<-mean((apply(X=fixedMt[,3:22],MARGIN=1,FUN=function(x){length(x[x==1])}))/20)
  results$MeanExoAut<-mean(apply(X=fixedMt[,3:22],MARGIN=1,FUN=function(x){mean(x=x,na.rm=T)}))
  results$SDExoAut<-mean(apply(X=fixedMt[,3:22],MARGIN=1,FUN=sd))
  FstHe<-read.table(paste(dirwork,"FstHeFile.txt",sep="/"),header=T)
  FixFstHe<-FstHe[data$Run[data$Mt==1],]
  results$FstAut<-mean(colMeans(FixFstHe[,3:21],na.rm=TRUE))
  results$FstZ<-mean(mean(FixFstHe[,22],na.rm=TRUE))
  results$FstW<-mean(mean(FixFstHe[,23],na.rm=TRUE))
  results$FstMt<-mean(mean(FixFstHe[,24],na.rm=TRUE))
  
  
  dIntAut<-as.vector((apply(X=fixedMt[,3:22],MARGIN=1,FUN=function(x){length(x[x!=0])}))/20)
  dIntAut10<-as.vector((apply(X=fixedMt[,3:22],MARGIN=1,FUN=function(x){length(x[x>=0.1])}))/20)
  dIntFixAut<-as.vector((apply(X=fixedMt[,3:22],MARGIN=1,FUN=function(x){length(x[x==1])}))/20)
  dIntExoAut<-as.vector(apply(X=fixedMt[,3:22],MARGIN=1,FUN=mean))
  
  results[,c("Introgq0","Introgq025","Introgq50","Introgq975","Introgq100")]<-quantile(dIntAut,probs=c(0,0.025,0.5,0.975,1))
  results[,c("Introg10q0","Introg10q025","Introg10q50","Introg10q975","Introg10q100")]<-quantile(dIntAut10,probs=c(0,0.025,0.5,0.975,1))
  results[,c("Fixq0","Fixq025","Fixq50","Fixq975","Fixq100")]<-quantile(dIntFixAut,probs=c(0,0.025,0.5,0.975,1))
  results[,c("MeanExoq0","MeanExoq025","MeanExoq50","MeanExoq975","MeanExoq100")]<-quantile(dIntExoAut,probs=c(0,0.025,0.5,0.975,1))
  
  
  resultsNoFix<-as.data.frame(matrix(nrow=1,ncol=13,dimnames=list(Simul,c("NbRun","NoFixMt","IntMt","MtExo","IntAut","IntAut10","FixAut","MeanExoAut","SDExoAut","FstAut","FstZ","FstW","FstMt"))))
  
  NofixedMt<-data[data$Mt<1,]
  resultsNoFix$NbRun<-max(data$Run)
  resultsNoFix$NoFixMt<-1-length(fixedMt[,1])/max(data$Run)
  resultsNoFix$MtExo<-mean(NofixedMt$Mt)
  resultsNoFix$IntMt<-length(unique(data[data$Mt>0 & data$Mt<1,1]))/(max(data$Run))
  resultsNoFix$IntAut<-mean((apply(X=NofixedMt[,3:22],MARGIN=1,FUN=function(x){length(x[x!=0])}))/20)
  resultsNoFix$IntAut10<-mean((apply(X=NofixedMt[,3:22],MARGIN=1,FUN=function(x){length(x[x>=0.1])}))/20)
  resultsNoFix$FixAut<-mean((apply(X=NofixedMt[,3:22],MARGIN=1,FUN=function(x){length(x[x==1])}))/20)
  resultsNoFix$MeanExoAut<-mean(apply(X=NofixedMt[,3:22],MARGIN=1,FUN=function(x){mean(x=x,na.rm=T)}),na.rm=T)
  resultsNoFix$SDExoAut<-mean(apply(X=NofixedMt[,3:22],MARGIN=1,FUN=sd))
  FstHe<-read.table(paste(dirwork,"FstHeFile.txt",sep="/"),header=T)
  NoFixFstHe<-FstHe[data$Run[data$Mt<1],]
  resultsNoFix$FstAut<-mean(colMeans(NoFixFstHe[,3:21],na.rm=TRUE))
  resultsNoFix$FstZ<-mean(mean(NoFixFstHe[,22],na.rm=TRUE))
  resultsNoFix$FstW<-mean(mean(NoFixFstHe[,23],na.rm=TRUE))
  resultsNoFix$FstMt<-mean(mean(NoFixFstHe[,24],na.rm=TRUE))
  dMtnoFix<-as.vector(NofixedMt$Mt)

  dIntAutNF<-as.vector((apply(X=NofixedMt[,3:22],MARGIN=1,FUN=function(x){length(x[x!=0])}))/20)
  dIntAut10NF<-as.vector((apply(X=NofixedMt[,3:22],MARGIN=1,FUN=function(x){length(x[x>=0.1])}))/20)
  dIntFixAutNF<-as.vector((apply(X=NofixedMt[,3:22],MARGIN=1,FUN=function(x){length(x[x==1])}))/20)
  dIntExoAutNF<-as.vector(apply(X=NofixedMt[,3:22],MARGIN=1,FUN=mean))
  
  resultsNoFix[,c("Introgq0","Introgq025","Introgq50","Introgq975","Introgq100")]<-quantile(dIntAutNF,probs=c(0,0.025,0.5,0.975,1),na.rm=T)
  resultsNoFix[,c("Introg10q0","Introg10q025","Introg10q50","Introg10q975","Introg10q100")]<-quantile(dIntAut10NF,probs=c(0,0.025,0.5,0.975,1),na.rm=T)
  resultsNoFix[,c("Fixq0","Fixq025","Fixq50","Fixq975","Fixq100")]<-quantile(dIntFixAutNF,probs=c(0,0.025,0.5,0.975,1),na.rm=T)
  resultsNoFix[,c("MeanExoq0","MeanExoq025","MeanExoq50","MeanExoq975","MeanExoq100")]<-quantile(dIntExoAutNF,probs=c(0,0.025,0.5,0.975,1),na.rm=T)
  
  resultsNoFix$MaxMtI<-max(dMtnoFix)
  
  
  toreturn<-list(results,list(dIntAut,dIntAut10,dIntFixAut,dIntExoAut),resultsNoFix,list(dIntAutNF,dIntAut10NF,dIntFixAutNF,dIntExoAutNF),dMtnoFix,distriMtMax)
  
  names(toreturn)<-list("results","DistribAuto","resultsNoFix","DistribAutoNoFix","DistribMtnoFix","DistriMtMax")
  names(toreturn$DistribAuto)<-list("Introg","Introg10","Fix","MeanExo")
  names(toreturn$DistribAutoNoFix)<-list("Introg","Introg10","Fix","MeanExo")
  return(toreturn)
}

WrapMtFunc<-function(Simuls, FilePrefix,Param)
{
  outR<-mtfunc(Simuls[1])
  results<-outR$results
  resultsnofix<-outR$resultsNoFix
  DistriMtMaxDF<-data.frame(rep(Simuls[1],nrow(outR$DistriMtMax)),outR$DistriMtMax)
  names(DistriMtMaxDF)<-c("Simul","DistriMtMax","DistriNuMtMax")
  results$Param<-Param[1]
  resultsnofix$Param<-Param[1]
  DistriMtMaxDF$Param<-Param[1]
  
  write.table(results,file=paste(FilePrefix,"results.txt",sep=""), sep="\t",quote=F,append=F,col.names=T)
  write.table(resultsnofix,file=paste(FilePrefix,"resultsnofix.txt",sep=""), sep="\t",quote=F,append=F,col.names=T)
  write.table(DistriMtMaxDF,file=paste(FilePrefix,"DistriMtMax.txt",sep=""), sep="\t",quote=F,append=F,col.names=T,row.names=F)
  
  if(length(Simuls>1))
    {
      for (i  in 2:length(Simuls))
        {
          outR<-mtfunc(Simuls[i])
          results<-outR$results
          resultsnofix<-outR$resultsNoFix
          results$Param<-Param[i]
          resultsnofix$Param<-Param[i]
          write.table(results,file=paste(FilePrefix,"results.txt",sep=""), sep="\t",quote=F,append=T,col.names=F)
          write.table(resultsnofix,file=paste(FilePrefix,"resultsnofix.txt",sep=""), sep="\t",quote=F,append=T,col.names=F)
          write.table(cbind(rep(Simuls[i],nrow(outR$DistriMtMax)),
                            outR$DistriMtMax,rep(Param[i],nrow(outR$DistriMtMax))),
                      file=paste(FilePrefix,"DistriMtMax.txt",sep=""), sep="\t",quote=F,append=T,col.names=F,row.names=F)
        }
  }
}