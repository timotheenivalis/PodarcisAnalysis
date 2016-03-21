mtfuncInvasion<-function(Simul)
{
  dirwork<-paste(getwd(),Simul,sep="/")
  results<-as.data.frame(matrix(nrow=1,ncol=13,dimnames=list(Simul,c("NbRun","FixMt","IntMt","MtExo","IntAut","IntAut10","FixAut","MeanExoAut","SDExoAut","FstAut","FstZ","FstW","FstMt"))))
  data<-read.table(paste(dirwork,"IntrogStats.txt",sep="/"),header=T) 
  profile<-read.table(paste(dirwork,"IntrogProfile.txt",sep="/"),header=T) 
  
  mtdistri<-1-tapply(X = profile[which(profile$x>17),"LocusMt"],INDEX = profile[which(profile$x>17),"Run"],FUN = mean)
  captures<-which(mtdistri==1)
  notcaptures<-which(mtdistri!=1)
  
  MtCopMax<-as.vector(mtdistri)
  
  x<-cbind(rowMeans(profile[which(profile$x>17),3:22]),INDEX=profile[which(profile$x>17),"Run"])
  NuCop<-1-tapply(X =x[,1],INDEX=x[,2], FUN=function(x){mean(x=x,na.rm=T)})
  distriMtMax<-cbind(MtCopMax,NuCop)
  
  fixedMt<-profile[which(profile$Run %in% captures),]
      
  autD<-vector(mode = "list")
  if(length(captures)>0)
    {
      for (i in 1:length(captures))
        {
          autD[[i]]<-1-apply(profile[which(profile$Run==captures[i] & profile$x>17),3:22],MARGIN = 2,FUN = mean)
        }
    }else{autD[[1]]<-NA}
  
  results$NbRun<-max(profile$Run)
  results$FixMt<-length(captures)/max(profile$Run)
  results$IntMt<-length(which(mtdistri>0))/(max(profile$Run))
  results$MtExo<-mean(mtdistri)
  results$IntAut<-mean(unlist(lapply(X = autD,FUN = function(x){length(which(x!=0))/20})))
  results$IntAut10<-mean(unlist(lapply(X = autD,FUN = function(x){length(which(x>=0.1))/20})))
  results$FixAut<-mean(unlist(lapply(X = autD,FUN = function(x){length(which(x==1))/20})))
  results$MeanExoAut<-mean(unlist(lapply(autD,FUN = function(x){mean(x,na.rm=T)})))
  results$SDExoAut<-mean(unlist(lapply(autD,FUN = function(x){sd(x,na.rm=T)})))
  results$FstAut<-NA
  results$FstZ<-NA
  results$FstW<-NA
  results$FstMt<-NA
  
  dIntAut<-unlist(lapply(X = autD,FUN = function(x){length(which(x!=0))/20}))
  dIntAut10<-unlist(lapply(X = autD,FUN = function(x){length(which(x>=0.1))/20}))
  dIntFixAut<-unlist(lapply(X = autD,FUN = function(x){length(which(x==1))/20}))
  dIntExoAut<-unlist(lapply(autD,FUN = function(x){mean(x,na.rm=T)}))
  
  results[,c("Introgq0","Introgq025","Introgq50","Introgq975","Introgq100")]<-quantile(dIntAut,probs=c(0,0.025,0.5,0.975,1),na.rm=T)
  results[,c("Introg10q0","Introg10q025","Introg10q50","Introg10q975","Introg10q100")]<-quantile(dIntAut10,probs=c(0,0.025,0.5,0.975,1),na.rm=T)
  results[,c("Fixq0","Fixq025","Fixq50","Fixq975","Fixq100")]<-quantile(dIntFixAut,probs=c(0,0.025,0.5,0.975,1),na.rm=T)
  results[,c("MeanExoq0","MeanExoq025","MeanExoq50","MeanExoq975","MeanExoq100")]<-quantile(dIntExoAut,probs=c(0,0.025,0.5,0.975,1),na.rm=T)
  
  
  
  autDnofix<-vector(mode = "list")
  if(length(notcaptures)>0)
    {
      for (i in 1:length(notcaptures))
        {
          autDnofix[[i]]<-1-apply(profile[which(profile$Run==notcaptures[i] & profile$x>17),3:22],MARGIN = 2,FUN = mean)
        }
    }else{autDnofix[[1]]<-NA}
  resultsNoFix<-as.data.frame(matrix(nrow=1,ncol=13,dimnames=list(Simul,c("NbRun","NoFixMt","IntMt","MtExo","IntAut","IntAut10","FixAut","MeanExoAut","SDExoAut","FstAut","FstZ","FstW","FstMt"))))
  NofixedMt<-profile[which(! profile$Run %in% captures),]  
  resultsNoFix$NbRun<-max(profile$Run)
  resultsNoFix$NoFixMt<-1-length(captures)/max(profile$Run)
  resultsNoFix$MtExo<-mean(mtdistri[which(mtdistri<1)])
  resultsNoFix$IntMt<-length(which(mtdistri>0 & mtdistri<1))/max(profile$Run)
  resultsNoFix$IntAut<-mean(unlist(lapply(X = autDnofix,FUN = function(x){length(which(x!=0))/20})))
  resultsNoFix$IntAut10<-mean(unlist(lapply(X = autDnofix,FUN = function(x){length(which(x>=0.1))/20})))
  resultsNoFix$FixAut<-mean(unlist(lapply(X = autDnofix,FUN = function(x){length(which(x==1))/20})))
  resultsNoFix$MeanExoAut<-mean(unlist(lapply(autDnofix,FUN = function(x){mean(x,na.rm=T)})))
  resultsNoFix$SDExoAut<-mean(unlist(lapply(autDnofix,FUN = function(x){sd(x,na.rm=T)})))
  resultsNoFix$FstAut<-NA
  resultsNoFix$FstZ<-NA
  resultsNoFix$FstW<-NA
  resultsNoFix$FstMt<-NA
  dMtnoFix<-as.vector(mtdistri[which(mtdistri<1)])
  
  dIntAutNF<-unlist(lapply(X=autDnofix,FUN=function(x){length(x[x!=0])}))/20
  dIntAut10NF<-unlist(lapply(X=autDnofix,FUN=function(x){length(x[x>=0.1])}))/20
  dIntFixAutNF<-unlist(lapply(X=autDnofix,FUN=function(x){length(x[x==0])}))/20
  dIntExoAutNF<-unlist(lapply(X=autDnofix,FUN=function(x){mean(x)}))
  
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


mtfuncInvasionAllGrid<-function(Simul)
{
  dirwork<-paste(getwd(),Simul,sep="/")
  results<-as.data.frame(matrix(nrow=1,ncol=13,dimnames=list(Simul,c("NbRun","FixMt","IntMt","MtExo","IntAut","IntAut10","FixAut","MeanExoAut","SDExoAut","FstAut","FstZ","FstW","FstMt"))))
  data<-read.table(paste(dirwork,"IntrogStats.txt",sep="/"),header=T) 
  profile<-read.table(paste(dirwork,"IntrogProfile.txt",sep="/"),header=T) 
  
  mtdistri<-1-tapply(X = profile[,"LocusMt"],INDEX = profile[,"Run"],FUN = mean)
  captures<-which(mtdistri==1)
  notcaptures<-which(mtdistri!=1)
  
  distriMtMax<-as.vector(mtdistri)
  
  fixedMt<-profile[which(profile$Run %in% captures),]
  
  distriMtMax<-tapply(X=data$Mt,INDEX=data$Run,FUN=function(x){max(x,na.rm=T)})
  
  MtCopMax<-as.vector(mtdistri)
  
  x<-cbind(rowMeans(profile[,3:22]),INDEX=profile[,"Run"])
  NuCop<-1-tapply(X =x[,1],INDEX=x[,2], FUN=function(x){mean(x=x,na.rm=T)})
  distriMtMax<-cbind(MtCopMax,NuCop)

  
  autD<-vector(mode = "list")
  if(length(captures)>0)
  {
    for (i in 1:length(captures))
    {
      autD[[i]]<-1-apply(profile[which(profile$Run==captures[i]),3:22],MARGIN = 2,FUN = mean)
    }
  }else{autD[[1]]<-NA}
  
  results$NbRun<-max(profile$Run)
  results$FixMt<-length(captures)/max(profile$Run)
  results$IntMt<-length(which(mtdistri>0))/(max(profile$Run))
  results$MtExo<-mean(mtdistri)
  results$IntAut<-mean(unlist(lapply(X = autD,FUN = function(x){length(which(x!=0))/20})))
  results$IntAut10<-mean(unlist(lapply(X = autD,FUN = function(x){length(which(x>=0.1))/20})))
  results$FixAut<-mean(unlist(lapply(X = autD,FUN = function(x){length(which(x==1))/20})))
  results$MeanExoAut<-mean(unlist(lapply(autD,FUN = function(x){mean(x,na.rm=T)})))
  results$SDExoAut<-mean(unlist(lapply(autD,FUN = function(x){sd(x,na.rm=T)})))
  results$FstAut<-NA
  results$FstZ<-NA
  results$FstW<-NA
  results$FstMt<-NA
  
  
  dIntAut<-unlist(lapply(X = autD,FUN = function(x){length(which(x!=0))/20}))
  dIntAut10<-unlist(lapply(X = autD,FUN = function(x){length(which(x>=0.1))/20}))
  dIntFixAut<-unlist(lapply(X = autD,FUN = function(x){length(which(x==1))/20}))
  dIntExoAut<-unlist(lapply(autD,FUN = function(x){mean(x,na.rm=T)}))
  
  results[,c("Introgq0","Introgq025","Introgq50","Introgq975","Introgq100")]<-quantile(dIntAut,probs=c(0,0.025,0.5,0.975,1),na.rm=T)
  results[,c("Introg10q0","Introg10q025","Introg10q50","Introg10q975","Introg10q100")]<-quantile(dIntAut10,probs=c(0,0.025,0.5,0.975,1),na.rm=T)
  results[,c("Fixq0","Fixq025","Fixq50","Fixq975","Fixq100")]<-quantile(dIntFixAut,probs=c(0,0.025,0.5,0.975,1),na.rm=T)
  results[,c("MeanExoq0","MeanExoq025","MeanExoq50","MeanExoq975","MeanExoq100")]<-quantile(dIntExoAut,probs=c(0,0.025,0.5,0.975,1),na.rm=T)
  
  
  
  autDnofix<-vector(mode = "list")
  if(length(notcaptures)>0)
  {
    for (i in 1:length(notcaptures))
    {
      autDnofix[[i]]<-1-apply(profile[which(profile$Run==notcaptures[i]),3:22],MARGIN = 2,FUN = mean)
    }
  }else{autDnofix[[1]]<-NA}
  resultsNoFix<-as.data.frame(matrix(nrow=1,ncol=13,dimnames=list(Simul,c("NbRun","NoFixMt","IntMt","MtExo","IntAut","IntAut10","FixAut","MeanExoAut","SDExoAut","FstAut","FstZ","FstW","FstMt"))))
  NofixedMt<-profile[which(! profile$Run %in% captures),]  
  resultsNoFix$NbRun<-max(profile$Run)
  resultsNoFix$NoFixMt<-1-length(captures)/max(profile$Run)
  resultsNoFix$MtExo<-mean(mtdistri[which(mtdistri<1)])
  resultsNoFix$IntMt<-length(which(mtdistri>0 & mtdistri<1))/max(profile$Run)
  resultsNoFix$IntAut<-mean(unlist(lapply(X = autDnofix,FUN = function(x){length(which(x!=0))/20})))
  resultsNoFix$IntAut10<-mean(unlist(lapply(X = autDnofix,FUN = function(x){length(which(x>=0.1))/20})))
  resultsNoFix$FixAut<-mean(unlist(lapply(X = autDnofix,FUN = function(x){length(which(x==1))/20})))
  resultsNoFix$MeanExoAut<-mean(unlist(lapply(autDnofix,FUN = function(x){mean(x,na.rm=T)})))
  resultsNoFix$SDExoAut<-mean(unlist(lapply(autDnofix,FUN = function(x){sd(x,na.rm=T)})))
  resultsNoFix$FstAut<-NA
  resultsNoFix$FstZ<-NA
  resultsNoFix$FstW<-NA
  resultsNoFix$FstMt<-NA
  dMtnoFix<-as.vector(mtdistri[which(mtdistri<1)])
  
  dIntAutNF<-unlist(lapply(X=autDnofix,FUN=function(x){length(x[x!=0])}))/20
  dIntAut10NF<-unlist(lapply(X=autDnofix,FUN=function(x){length(x[x>=0.1])}))/20
  dIntFixAutNF<-unlist(lapply(X=autDnofix,FUN=function(x){length(x[x==0])}))/20
  dIntExoAutNF<-unlist(lapply(X=autDnofix,FUN=function(x){mean(x)}))
  
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