MeanAutLociIntro<-list() #Loci nucl?aires introgresses
MeanAutCopyIntro<-list() #Proportions de copies d'origine introgressee sur les loci nucl?aires
MeanAutLociI10<-list() #proportion de loci nucl?aires avec + de 10% d'introgression
MeanAutLociFix<-list() #proportion de loci nucl?aires avec capture

mtfunc<-function()
{
  results<-matrix(nrow=1,ncol=12,dimnames=list(Simul,c("NbRun","FixMt","IntMt","IntAut","IntAut10","FixAut","MeanExoAut","SDExoAut","FstAut","FstZ","FstW","FstMt")))
  data<-read.table("IntrogStats.txt",header=T) 
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
  FstHe<-read.table("FstHeFile.txt",header=T)
  
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

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/TestHybridRate/HR0/")
Simul<-"HR0"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/TestHybridRate/HRresults.txt", sep="\t",quote=F,append=F,col.names=T)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations//TestHybridRate/HR1/")
Simul<-"HR1"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/TestHybridRate/HRresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations//TestHybridRate/HR2/")
Simul<-"HR2"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/TestHybridRate/HRresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations//TestHybridRate/HR3/")
Simul<-"HR3"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/TestHybridRate/HRresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations//TestHybridRate/HR4/")
Simul<-"HR4"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/TestHybridRate/HRresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations//TestHybridRate/HR5/")
Simul<-"HR5"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/TestHybridRate/HRresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations//TestHybridRate/HR6/")
Simul<-"HR6"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/TestHybridRate/HRresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations//TestHybridRate/HR7/")
Simul<-"HR7"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/TestHybridRate/HRresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations//TestHybridRate/HR8/")
Simul<-"HR8"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/TestHybridRate/HRresults.txt", sep="\t",quote=F,append=T,col.names=F)

##############################################################################################
#############################################################################################
par(mfrow=c(1,1))
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/TestHybridRate/")
simul<-read.table(file="HRresults.txt",header=T)
Param<-matrix(data=c(0.3,0.25,0.2,0.15,0.1,0.05,0.01,0.005,0.001),nrow=9,ncol=1)
simul$Param<-Param
simul<-simul[sort.list(simul$Param,decreasing=F),]
simultitle<-"Hybrid fitness"

plot(simul$Param,simul$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab=simultitle,ylab="Proportion de simulations",main=paste("Proportion de simulations avec Capture ou Introgression Mitochondriale,\n",simultitle),cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(simul$Param,simul$IntMt,col="blue",pch=2,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Simulations avec Introgression Mitochondriale"),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)

plot(simul$Param,simul$FixMt,col="dark blue",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab=simultitle,ylab="",main=paste("Proportion de simulations avec Capture Mitochondriale, et parmi elles, proportion d'Introgression Nucl?aire,\n",simultitle),cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(jitter(simul$Param),simul$IntAut10,col="light green",pch=4,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Proportion de loci nucl?aires introgress?s"),col=c("dark blue","light green"),pch=c(1,4),ncol=1,pt.lwd=2,cex=1.2)


plot(simul$Param,simul$IntAut,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),pch=3,lwd=2,xlab=simultitle,ylab="Proportion",main=paste("Force de l'introgression Autosomale lors de captures mitochondriales,\n ",simultitle))
points(jitter(simul$Param),simul$IntAut10,col="green",pch=4,lwd=2,type="b")
points(simul$Param,simul$FixAut,col="orange",pch=5,lwd=2,type="b")
points(simul$Param,simul$MeanExoAut,lwd=2,pch=6,type="b")
legend(x="topleft",legend=c("Proportion moyenne de loci autosomaux introgress?s","Proportion moyenne de loci autosomaux dont plus de 10% des copies sont introgress?es","Proportion de loci autosomaux captur?s","Proportion moyenne de copies introgresses"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)

plot(simul$Param,simul$FstAut,col="dark green",type="p",ylim=c(0,1),pch=1,lwd=2,xlab=simultitle,ylab="Fst entre les deux habitats",main=paste("Fst en fonction de ",simultitle))
points(simul$Param,simul$FstZ,col="green",pch=2,lwd=2)
points(jitter(simul$Param),simul$FstW,col="blue",pch=3,lwd=2)
points(jitter(simul$Param),simul$FstMt,col="red",pch=4,lwd=2)
legend(x="topleft",legend=c("Autosomes","Z","W","Mt"),col=c("dark green","green","blue","red"),pt.lwd=2,pch=c(1,2,3,4))
