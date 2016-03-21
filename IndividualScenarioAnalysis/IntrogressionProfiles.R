setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/TestAnalyses/Disp/BD3")
Introg<-read.table(file="IntrogProfile.txt", header=TRUE)
par(mfrow=c(3,3))
for (i in 2:(dim(Introg)[2]))
  {
    plot(Introg[,2],Introg[,i],xlab=names(Introg)[i],ylab="% de genes de type 1")
    abline(v=15, col="red")
  }

#############Introgression stat now

MeanAutLociIntro<-list() #Loci nucl?aires introgresses
MeanAutCopyIntro<-list() #Proportions de copies d'origine introgressee sur les loci nucl?aires
MeanAutLociI10<-list() #proportion de loci nucl?aires avec + de 10% d'introgression
MeanAutLociFix<-list() #proportion de loci nucl?aires avec capture
# Merge the data-frame
data<-read.table("IntrogStats.txt",header=T)
Simul<-"BD5"
results<-matrix(nrow=1,ncol=12,dimnames=list(Simul,c("NbRun","FixMt","IntMt","IntAut","IntAut10","FixAut","MeanExoAut","SDExoAut","FstAut","FstZ","FstW","FstMt")))

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

mean(data$Nucl_Intro) ### sur toutes les simul, prop de loci nucleaires introg, dans chaque espece (pas globalement)
var(data$Nucl_Intro)

Fixdata<-data[condfix,] ### sur les simul avec capture, prop de loci nucleaires introg, dans chaque espece (pas globalement)
mean(Fixdata$Nucl_Intro)
var(Fixdata$Nucl_Intro)

mean(colMeans(data[which(data$Sp==0),3:22]))
mean(colMeans(data[which(data$Sp==1),3:22]))
mean(data[which(data$Sp==0),23])
mean(data[which(data$Sp==1),23])

mean(data[which(data$Sp==0),24])
mean(data[which(data$Sp==1),24])

mean(data[which(data$Sp==0),25])
mean(data[which(data$Sp==1),25])

mean(data[which(data$Sp==0),27])
mean(data[which(data$Sp==1),27])

MeanAutLociIntro$BD5<-Fixdata$Nucl_Intro

comptIntro<-vector(l=22) # prop de simul avec intro de chaque locus
comptIntro10p<-vector(l=22) #prop de simul avec plus de 10% des echantillons introgresses pour chaque locus
comptFix<-vector(l=22) #prop de simul avec fixation de chaque locus
vectorFix<-vector(l=length(levels(as.factor(Fixdata$Run))))
vectorIntro10p<-vector(l=length(levels(as.factor(Fixdata$Run))))
c<-0
b<-0

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
MeanAutLociI10$BD5<-vectorIntro10p
MeanAutLociFix$BD5<-vectorFix
results[1,4]<-mean(comptIntro[1:19]/length(levels(as.factor(Fixdata$Run))))
results[1,5]<-mean(comptIntro10p[1:19]/length(levels(as.factor(Fixdata$Run))))
results[1,6]<-mean(comptFix[1:19]/length(levels(as.factor(Fixdata$Run)))) #doit valloir 1 pour le W (comme il suit la mt)
results[1,7]<-mean(mean(Fixdata[,3:22]))
results[1,8]<-mean(sd(Fixdata[,3:22]))

FstHe<-read.table("FstHeFile.txt",header=T)

FixFstHe<-FstHe[condfixrun,]
mean(FixFstHe[,3:22],na.rm=TRUE)
mean(mean(FixFstHe[,3:21],na.rm=TRUE))

results[1,9]<-mean(mean(FixFstHe[,3:21],na.rm=TRUE))
results[1,10]<-mean(mean(FixFstHe[,22],na.rm=TRUE))
results[1,11]<-mean(mean(FixFstHe[,23],na.rm=TRUE))
results[1,12]<-mean(mean(FixFstHe[,24],na.rm=TRUE))
results<-round(results,3)

write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/TestAnalyses/Disp/BDresults.txt", sep="\t",quote=F,append=F,col.names=T)
MeanAutCopyIntro$BD5<-as.double(rowMeans(Fixdata[,3:22],na.rm=T))

###################################################################

setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/TestAnalyses/Disp/BD3")
Introg<-read.table(file="IntrogProfile.txt", header=TRUE)
par(mfrow=c(2,4))
for (i in 2:(dim(Introg)[2]))
{
  plot(Introg[,2],Introg[,i],xlab=names(Introg)[i],ylab="% de genes de type 1")
  abline(v=15, col="red")
}

#############Introgression stat now

MeanAutLociIntro<-list() #Loci nucl?aires introgresses
MeanAutCopyIntro<-list() #Proportions de copies d'origine introgressee sur les loci nucl?aires
MeanAutLociI10<-list() #proportion de loci nucl?aires avec + de 10% d'introgression
MeanAutLociFix<-list() #proportion de loci nucl?aires avec capture
# Merge the data-frame
data<-read.table("IntrogStats.txt",header=T)
Simul<-"BD3"

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

results[2,1]<-max(mt$Run)
results[2,2]<-Fix  ### % de simul avec capture mt
results[2,3]<-Introg ### % de simul avec introg mt

mean(data$Nucl_Intro) ### sur toutes les simul, prop de loci nucleaires introg, dans chaque espece (pas globalement)
var(data$Nucl_Intro)

Fixdata<-data[condfix,] ### sur les simul avec capture, prop de loci nucleaires introg, dans chaque espece (pas globalement)
mean(Fixdata$Nucl_Intro)
var(Fixdata$Nucl_Intro)

MeanAutLociIntro$BD3<-Fixdata$Nucl_Intro

comptIntro<-vector(l=22) # prop de simul avec intro de chaque locus
comptIntro10p<-vector(l=22) #prop de simul avec plus de 10% des echantillons introgresses pour chaque locus
comptFix<-vector(l=22) #prop de simul avec fixation de chaque locus
vectorFix<-vector(l=length(levels(as.factor(Fixdata$Run))))
vectorIntro10p<-vector(l=length(levels(as.factor(Fixdata$Run))))
c<-0
b<-0

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
MeanAutLociI10$BD3<-vectorIntro10p
MeanAutLociFix$BD3<-vectorFix
results[2,4]<-mean(comptIntro[1:19]/length(levels(as.factor(Fixdata$Run))))
results[2,5]<-mean(comptIntro10p[1:19]/length(levels(as.factor(Fixdata$Run))))
results[2,6]<-mean(comptFix[1:19]/length(levels(as.factor(Fixdata$Run)))) #doit valloir 1 pour le W (comme il suit la mt)
results[2,7]<-mean(mean(Fixdata[,3:22]))
results[2,8]<-mean(sd(Fixdata[,3:22]))

FstHe<-read.table("FstHeFile.txt",header=T)

FixFstHe<-FstHe[condfixrun,]
mean(FixFstHe[,3:22],na.rm=TRUE)
mean(mean(FixFstHe[,3:21],na.rm=TRUE))

results[2,9]<-mean(mean(FixFstHe[,3:21],na.rm=TRUE))
results[2,10]<-mean(mean(FixFstHe[,22],na.rm=TRUE))
results[2,11]<-mean(mean(FixFstHe[,23],na.rm=TRUE))
results[2,12]<-mean(mean(FixFstHe[,24],na.rm=TRUE))
results<-round(results,3)

write.table(results,file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/TestAnalyses/Disp/BDresults.txt", sep="\t",quote=F,append=F,col.names=T)
MeanAutCopyIntro$BD3<-as.double(rowMeans(Fixdata[,3:22],na.rm=T))