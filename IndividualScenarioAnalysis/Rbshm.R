mtfunc<-function(Sce)
{
  results<-matrix(nrow=1,ncol=12,dimnames=list(Simul,c("NbRun","FixMt","IntMt","IntAut","IntAut10","FixAut","MeanExoAut","SDExoAut","FstAut","FstZ","FstW","FstMt")))
  data<-read.table(paste(Sce,"/IntrogStats.txt",sep=""),header=T) 
  mt<-data[,c(1,2,26)]
  condfix<-vector() ### recupere les numero de lignes ou il y a de la fixation
  condfixrun<-vector() ###recupere les numero de run ou il y a fixation
  compteurF=0
  compteurI=0
  for (run in 1:max(mt$Run))
  {
    mtrun<-mt[which(mt$Run==run),]
    if ((mtrun[1,3])==1||(mtrun[2,3]==1))
    {
      compteurF=compteurF+1   
      condfixrun<-c(condfixrun,run)
      if (mtrun[1,3]==1)
      {
        condfix<-c(condfix,2*run-1)
      }
      if(mtrun[2,3]==1)
      {
        condfix<-c(condfix,2*run)
      }
    }
    if((mtrun[1,3]!=0)||(mtrun[2,3]!=0))
    {
      compteurI=compteurI+1
    }
  }
  Fix=compteurF/max(mt$Run)
  Introg=compteurI/max(mt$Run)
  results[1,1]<-max(mt$Run)
  results[1,2]<-Fix  ### % de simul avec capture mt
  results[1,3]<-Introg ### % de simul avec introg mt
  Fixdata<-data[condfix,] ### sur les simul avec capture, prop de loci nucleaires introg, dans chaque espece (pas globalement)
  
  comptIntro<-vector(l=22) # prop de simul avec intro de chaque locus
  comptIntro10p<-vector(l=22) #prop de simul avec plus de 10% des echantillons introgresses pour chaque locus
  comptFix<-vector(l=22) #prop de simul avec fixation de chaque locus
  vectorFix<-vector(l=length(levels(as.factor(Fixdata$Run))))
  vectorIntro10p<-vector(l=length(levels(as.factor(Fixdata$Run))))
  c<-0
  b<-0
  
  if (length(levels(as.factor(Fixdata$Run)))!=0)
  {
    for (r in 1:length(levels(as.factor(Fixdata$Run))))
    {
      run<-as.integer(levels(as.factor(Fixdata$Run))[r])
      autrun<-Fixdata[which(Fixdata$Run==run),]
      for (aut in 3:22)
      {
        if(autrun[1,aut]==1)
        {
          comptFix[aut-2]=1+comptFix[aut-2]
          b=b+1/(20)
        }
        if(autrun[1,aut]!=0)
        {
          comptIntro[aut-2]=1+comptIntro[aut-2]
        }
        if(autrun[1,aut]>0.1)
        {
          comptIntro10p[aut-2]=1+comptIntro10p[aut-2]
          c=c+1/(20)
        }
      }
      vectorIntro10p[r]<-c
      c<-0
      vectorFix[r]<-b
      b<-0
    }
    results[1,4]<-mean(comptIntro[1:19]/length(levels(as.factor(Fixdata$Run))))
    results[1,5]<-mean(comptIntro10p[1:19]/length(levels(as.factor(Fixdata$Run))))
    results[1,6]<-mean(comptFix[1:19]/length(levels(as.factor(Fixdata$Run)))) #doit valloir 1 pour le W (comme il suit la mt)
    results[1,7]<-mean(colMeans(Fixdata[,3:22]))
    results[1,8]<-mean(sapply((Fixdata[,3:22]),sd))
  }
  FstHe<-read.table(paste(Sce,"/FstHeFile.txt",sep=""),header=T)
  
  FixFstHe<-FstHe[condfixrun,]
  #colMeans(FixFstHe[,3:24],na.rm=TRUE)
  #mean(colMeans(FixFstHe[,3:21],na.rm=TRUE))
  results[1,9]<-mean(colMeans(FixFstHe[,3:21],na.rm=TRUE))
  results[1,10]<-mean(mean(FixFstHe[,22],na.rm=TRUE))
  results[1,11]<-mean(mean(FixFstHe[,23],na.rm=TRUE))
  results[1,12]<-mean(mean(FixFstHe[,24],na.rm=TRUE))
  results<-round(results,3)
  return (results)
}

setwd("C:/Users/Timothée/Dropbox/PodarcisCBGP/BSHM")
diri<-getwd()

###################################################################
Sce<-"BSHM0"
results<-mtfunc(Sce)
write.table(results,file=paste(diri,"/BSHMresults.txt",sep=""), sep="\t",quote=F,append=F,col.names=T)

for (i in 1:5)
{
  Sce<-paste("BSHM",i, sep="")
  results<-mtfunc(Sce)
  write.table(results,file=paste(diri,"/BSHMresults.txt",sep=""), sep="\t",quote=F,append=T,col.names=F)
}

simul<-read.table(file="BSHMresults.txt",header=T)
