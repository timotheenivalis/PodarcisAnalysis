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
  for (run1 in 1:(length(mt$Run)/2))
  {
    run<-mt$Run[run1]
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
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt0")
Simul<-"HMt0"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=F,col.names=T)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt1")
Simul<-"HMt1"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt2")
Simul<-"HMt2"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt3")
Simul<-"HMt3"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt4")
Simul<-"HMt4"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt5")
Simul<-"HMt5"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt6")
Simul<-"HMt6"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt7")
Simul<-"HMt7"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt8")
Simul<-"HMt8"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt9")
Simul<-"HMt9"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt10")
Simul<-"HMt10"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt11")
Simul<-"HMt11"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt12")
Simul<-"HMt12"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt13")
Simul<-"HMt13"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt14")
Simul<-"HMt14"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt15")
Simul<-"HMt15"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt16")
Simul<-"HMt16"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

###################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt17")
Simul<-"HMt17"
results<-mtfunc()
write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMtresults.txt", sep="\t",quote=F,append=T,col.names=F)

##############################################################################################
#############################################################################################
par(mfrow=c(1,1))
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt")
simul<-read.table(file="HMtresults.txt",header=T)
Param<-matrix(data=0:17,nrow=18,ncol=1)
simul$Param<-Param
simul<-simul[sort.list(simul$Param,decreasing=F),]

plot(simul$Param,simul$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab="sélection mt",ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n En fonction de la sélecion mt",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(simul$Param,simul$IntMt,col="blue",pch=2,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Simulations avec Introgression Mitochondriale"),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)

plot(simul$Param,simul$FixMt,col="dark blue",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab="sélection mt",ylab="",main="Proportion de simulations avec Capture Mitochondriale, et parmi elles, proportion d'Introgression Nucléaire,\n En fonction de la sélecion mt",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
points(jitter(simul$Param),simul$IntAut10,col="light green",pch=4,lwd=2,type="b")
legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Proportion de loci nucléaires introgressés"),col=c("dark blue","light green"),pch=c(1,4),ncol=1,pt.lwd=2,cex=1.2)


plot(simul$Param,simul$IntAut,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),pch=3,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion",main="Force de l'introgression Autosomale lors de captures mitochondriales,\n En fonction de la diff?rence de dispersion entre sexe (Femelle - Male)")
points(jitter(simul$Param),simul$IntAut10,col="green",pch=4,lwd=2,type="b")
points(simul$Param,simul$FixAut,col="orange",pch=5,lwd=2,type="b")
points(simul$Param,simul$MeanExoAut,lwd=2,pch=6,type="b")
legend(x="topleft",legend=c("Proportion moyenne de loci autosomaux introgress?s","Proportion moyenne de loci autosomaux dont plus de 10% des copies sont introgress?es","Proportion de loci autosomaux captur?s","Proportion moyenne de copies introgresses"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)

plot(simul$Param,simul$FstAut,col="dark green",type="p",ylim=c(0,1),pch=1,lwd=2,xlab="sélection mt",ylab="Fst entre les deux habitats",main="Fst en fonction de la sélection mt")
points(simul$Param,simul$FstZ,col="green",pch=2,lwd=2)
points(jitter(simul$Param),simul$FstW,col="blue",pch=3,lwd=2)
points(jitter(simul$Param),simul$FstMt,col="red",pch=4,lwd=2)
legend(x="topleft",legend=c("Autosomes","Z","W","Mt"),col=c("dark green","green","blue","red"),pt.lwd=2,pch=c(1,2,3,4))




