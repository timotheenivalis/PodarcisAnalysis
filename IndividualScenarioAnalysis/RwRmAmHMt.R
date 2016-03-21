source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
setwd("C:/Users/Timothée//Documents/Studies/PodarcisCBGP/New simulations/wRmAmHMt/")

Simuls<-"wRmAmHMt0"
for (nbsimul in 1:17)
{
  Simuls<-c(Simuls,paste("wRmAmHMt",nbsimul,sep=""))
}

WrapMtFunc(Simuls = Simuls,FilePrefix = "wRmAmHMt",Param = 1:length(Simuls))


# source(file="C:/Users/Timothée/Dropbox/PodarcisCBGP/New simulations/FunctionsSimul.R")
# setwd("C:/Users/Timothée/Dropbox/PodarcisCBGP/New simulations/wRmAmHMt/")
# 
# Simul<-"wRmAmHMt0"
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
# write.table(results,file="wRmAmHMtresults.txt", sep="\t",quote=F,append=F,col.names=T)
# write.table(resultsnofix,file="wRmAmHMtresultsnofix.txt", sep="\t",quote=F,append=F,col.names=T)
# write.table(DistriMtMaxDF,file="wRmAmHMtDistriMtMax.txt", sep="\t",quote=F,append=F,col.names=T,row.names=F)
# 
# 
# for (nbsimul in 1:17)
# {
#   Simul<-paste("wRmAmHMt",nbsimul,sep="")
#   outR<-mtfunc(Simul)
#   results<-outR$results
#   resultsnofix<-outR$resultsNoFix
#   DistriMtMaxDF<-data.frame(rep(Simul,length(outR$DistriMtMax)),outR$DistriMtMax)
#   names(DistriMtMaxDF)<-c("Simul","DistriMtMax")
#   results$Param<-1
#   resultsnofix$Param<-1
#   DistriMtMaxDF$Param<-1
#   
#   write.table(results,file="wRmAmHMtresults.txt", sep="\t",quote=F,append=T,col.names=F)
#   write.table(resultsnofix,file="wRmAmHMtresultsnofix.txt", sep="\t",quote=F,append=T,col.names=F)
#   write.table(DistriMtMaxDF,file="wRmAmHMtDistriMtMax.txt", sep="\t",quote=F,append=T,col.names=F,row.names=F)
# }
# 
# ##############################################################################################
# #############################################################################################
# par(mfrow=c(1,1))
# setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/")
# simul<-read.table(file="wRmAmHMtresults.txt",header=T)
# Param<-matrix(data=c(0:17),nrow=18,ncol=1)
# simul$Param<-Param
# plot(simul$Param,simul$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab="sélection mt",ylab="Proportion de simulations",main="Proportion de simulations avec Capture ou Introgression Mitochondriale,\n En fonction de la sélecion mt",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
# points(simul$Param,simul$IntMt,col="blue",pch=2,lwd=2,type="b")
# legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Simulations avec Introgression Mitochondriale"),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)
# 
# plot(simul$Param,simul$FixMt,col="dark blue",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),pch=1,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="",main="Proportion de simulations avec Capture ,\n En fonction du logarithme du rapport de dispersion entre sexe (Femelle/Male)",cex.lab=1.8,cex.axis=1.5,cex.main=1.6)
# points(jitter(simul$Param),simul$IntAut10,col="light green",pch=4,lwd=2,type="b")
# legend(x="topleft",legend=c("Simulations avec Capture Mitochondriale","Proportion de loci nucl?aires introgress?s"),col=c("dark blue","light green"),pch=c(1,4),ncol=1,pt.lwd=2,cex=1.6)
# 
# 
# plot(simul$Param,simul$IntAut,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),pch=3,lwd=2,xlab=c(expression(paste("log du rapport des ",sigma^2))),ylab="Proportion",main="Force de l'introgression Autosomale lors de captures mitochondriales,\n En fonction de la diff?rence de dispersion entre sexe (Femelle - Male)")
# points(jitter(simul$Param),simul$IntAut10,col="green",pch=4,lwd=2,type="b")
# points(simul$Param,simul$FixAut,col="orange",pch=5,lwd=2,type="b")
# points(simul$Param,simul$MeanExoAut,lwd=2,pch=6,type="b")
# legend(x="topleft",legend=c("Proportion moyenne de loci autosomaux introgress?s","Proportion moyenne de loci autosomaux dont plus de 10% des copies sont introgress?es","Proportion de loci autosomaux captur?s","Proportion moyenne de copies introgresses"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)
# 
# plot(simul$Param,simul$FstAut,col="dark green",type="p",ylim=c(0,1),pch=1,lwd=2,xlab=c(expression(paste("D",Delta,sigma^2))),ylab="Fst entre les deux habitats",main="Fst en fonction de la diff?rence de dispersion entre sexe (Femelle - Male)")
# points(simul$Param,simul$FstZ,col="green",pch=2,lwd=2)
# points(jitter(simul$Param),simul$FstW,col="blue",pch=3,lwd=2)
# points(jitter(simul$Param),simul$FstMt,col="red",pch=4,lwd=2)
# legend(x="topleft",legend=c("Autosomes","Z","W","Mt"),col=c("dark green","green","blue","red"),pt.lwd=2,pch=c(1,2,3,4))
# 
# 
# simul$Nselection<-matrix(data=c(rep(x=0.9,6),rep(0.5,6),rep(0.1,6)),nrow=18,ncol=1)
# simul$Homogamy<-matrix(data=rep(c(0.1,0.2,0.3,0.4,0.45,0.49),3),nrow=18,ncol=1)
# 
# library("reshape")#for melt
# forgrid<-melt.data.frame(simul,id.vars=c(14,15),measure.vars="IntAut10")
# forgrid<-forgrid[order(forgrid$Nselection),]
# 
# forim<-matrix(forgrid[,4],nrow=6,ncol=3)
# bonsens<-seq(min(forim),max(forim),by=1/15)
# #bonsens<-forim[order(forim)]
# image(x=sort(unique(simul[,15])),y=sort(unique(simul[,14])),z=forim,xlab="homogamy",ylab="nuclear selection",col=rgb(red=bonsens,green=0,blue=1-bonsens), main="Homogamy, selecion N et selection mt \n recombi faible")
# col=rgb(red=bonsens,green=0,blue=1-bonsens)
# legend(title="introgression nucléaire >10%",legend=c(max(forim),min(forim)),col=c(col[length(col)],col[1]),x="topleft",lwd=10)
# image(x=sort(unique(simul[,15])),y=sort(unique(simul[,14])),z=forim,xlab="homogamy",ylab="nuclear selection",col=heat.colors(30), main="Homogamy, selecion N et selection mt \n Recombi faible")
# col=heat.colors(30)
# legend(title="introgression nucléaire >10%",legend=c(max(forim),min(forim)),col=c(col[30],col[1]),x="topleft",lwd=10)