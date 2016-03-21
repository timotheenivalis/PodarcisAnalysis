MergeFolders<-function(prim,sec)
{
  filesPrim<-list.files(path=prim)
  MaxGPrim1<-max(as.numeric(substr(filesPrim[substr(filesPrim,start=1,stop=11)=="GenepopFile"],start=13,stop=14)))
  MaxGPrim2<-max(as.numeric(substr(filesPrim[substr(filesPrim,start=1,stop=11)=="GenepopFile"],start=13,stop=15)),na.rm=T)
  MaxGPrim<-max(MaxGPrim1,MaxGPrim2)
  
  filesSec<-list.files(path=sec)
  GenepopSec<-filesSec[substr(filesSec,start=1,stop=11)=="GenepopFile"]
  
  for (i in 1:length(GenepopSec))
    {  
    stopping<-14
    if(substr(GenepopSec[i],16,16)=="."){stopping<-15}
    newIndex<-as.numeric(substr(GenepopSec[i],13,stopping))+MaxGPrim
    newName<-paste("BenepopFile_",newIndex,".txt",sep="")
    newDir<-paste(sec,newName,sep="")
    oldDir<-paste(sec,GenepopSec[i],sep="")
    
    file.rename(oldDir,newDir)
    renewDir<-paste(prim,paste("GenepopFile_",newIndex,".txt",sep=""),sep="")
    file.copy(from=newDir,to=renewDir)
    
  }
  fstSecDir<-paste(sec,"FstHeFile.txt",sep="")
  fstSec<-read.table(fstSecDir,header=T)
  fstSec$Run<-fstSec$Run+MaxGPrim
  fstPrimDir<-paste(prim,"FstHeFile.txt",sep="")
  write.table(x=fstSec,file=fstPrimDir,append=TRUE,quote=F,row.names=F,col.names=F,sep="\t")
  
  StatsSecDir<-paste(sec,"IntrogStats.txt",sep="")
  StatsSec<-read.table(StatsSecDir,header=T)
  StatsSec$Run<-StatsSec$Run+MaxGPrim
  StatsPrimDir<-paste(prim,"IntrogStats.txt",sep="")
  write.table(x=StatsSec,file=StatsPrimDir,append=TRUE,quote=F,row.names=F,col.names=F,sep="\t")
  
  ProfSecDir<-paste(sec,"IntrogProfile.txt",sep="")
  ProfSec<-read.table(ProfSecDir,header=T)
  ProfSec$Run<-ProfSec$Run+MaxGPrim
  ProfPrimDir<-paste(prim,"IntrogProfile.txt",sep="")
  write.table(x=ProfSec,file=ProfPrimDir,append=TRUE,quote=F,row.names=F,col.names=F,sep="\t")

}

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt3/"
#sec<-("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/multiHMt3bis/")

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/HMt14/"
#sec<-("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/multiHMt14bis/")

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MAMHMt/MAMHMt4/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MAMHMt/MAMHMt4bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MAMHMt/MAMHMt15/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MAMHMt/MAMHMt15bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MAMHMt/MAMHMt16/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MAMHMt/MAMHMt16bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt0/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt0bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt0/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt0tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt1/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt1bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt2/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt2bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt2/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt2tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt3/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt3bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt3/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt3tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt4/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt4bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt4/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt4tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt6/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt6bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt6/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt6tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt7/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt7bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt7/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt7tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt8/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt8bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt9/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt9bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt9/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt9tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt11/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt11bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt11/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt11tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt12/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt12bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt12/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt12tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt13/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt13bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt13/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt13tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt15/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt15bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt15/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/wRmAmHMt15tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun12/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun12bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun12/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun12tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun9/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun9bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun9/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun9tris/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun17/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun17bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun19/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun19bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun20/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/gHun/gHun20bis/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun0/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun0b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun1/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun1b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun2/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun2b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun3/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun3b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun4/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun4b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun5/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun5b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun6/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmHun/lmHun6b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC0/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC0b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC1/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC1b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC2/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC2b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC3/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC3b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC4/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC4b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC5/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC5b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC6/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC6b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC7/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC7b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC8/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDC/lmBDC8b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC0/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC0b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC1/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC1b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC2/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC2b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC3/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC3b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC4/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC4b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC5/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC5b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC6/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC6b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC7/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC7b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC8/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBDCC/lmBDCC8b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC0/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC0b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC1/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC1b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC2/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC2b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC3/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC3b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC4/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC4b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC5/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC5b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC6/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC6b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC7/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC7b/"

#prim<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC8/"
#sec<-"C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBHC/lmBHC8b/"

MergeFolders(prim,sec)
