par(mfrow=c(1,1))
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MAMHMt/")
simulM<-read.table(file="MAMHMtresults.txt",header=T)
Param<-matrix(data=c(0:17),nrow=18,ncol=1)
simulM$Param<-Param

simulM$Nselection<-matrix(data=c(rep(x=0.9,6),rep(0.5,6),rep(0.1,6)),nrow=18,ncol=1)
simulM$Homogamy<-matrix(data=rep(c(0.1,0.2,0.3,0.4,0.45,0.49),3),nrow=18,ncol=1)

library("reshape")#for melt
forgridM<-melt.data.frame(simulM,id.vars=c(14,15),measure.vars="IntAut10")
forgridM<-forgridM[order(forgridM$Nselection),]

forimM<-matrix(forgridM[,4],nrow=6,ncol=3)
bonsensM<-seq(min(forimM),max(forimM),by=1/15)

ColorRamp <- rgb( seq(0,1,length=256),  # Red
                  seq(0,0,length=256),  # Green
                  seq(1,0,length=256))  # Blue
ColorLevels <- seq(0, 1, length=length(ColorRamp))
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
par(mar=c(5.5, 5, 4, 1.5))
image(x=sort(unique(simulM[,15])),y=sort(unique(simulM[,14])),z=forimM,cex.axis=1.2,cex.lab=1.2,xlab="homogamy",ylab="nuclear selection",col=rgb(red=bonsensM,green=0,blue=1-bonsensM))
for (i in 1:dim(forimM)[1])
  for (j in 1:dim(forimM)[2])
  {
    text(x=sort(unique(simulM[,15]))[i],y=sort(unique(simulM[,14]))[j], labels=round(forimM[i,j],digits=3),col="white")
  }
par(mar = c(5.5,0.5,4,5))
image(1, ColorLevels,axes=FALSE,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n",cex.lab=1.2)
axis(RIGHT <-4, at=seq(from=0,to=1,by=0.2), labels=seq(from=0,to=1,by=0.2), las= HORIZONTAL<-1, cex.axis=1)
mtext(text="IAL10",side=4,line=3,cex=1.2)
layout(1)
################################################################################################################

setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/sRmAmHMt/")
simuls<-read.table(file="sRmAmHMtresults.txt",header=T)
Param<-matrix(data=c(0:17),nrow=18,ncol=1)
simuls$Param<-Param
simuls$Nselection<-matrix(data=c(rep(x=0.9,6),rep(0.5,6),rep(0.1,6)),nrow=18,ncol=1)
simuls$Homogamy<-matrix(data=rep(c(0.1,0.2,0.3,0.4,0.45,0.49),3),nrow=18,ncol=1)

library("reshape")#for melt
forgrids<-melt.data.frame(simuls,id.vars=c(14,15),measure.vars="IntAut10")
forgrids<-forgrids[order(forgrids$Nselection),]

forims<-matrix(forgrids[,4],nrow=6,ncol=3)
bonsenss<-seq(min(forims),max(forims),by=1/15)

ColorRamp <- rgb( seq(0,1,length=256),  # Red
                  seq(0,0,length=256),  # Green
                  seq(1,0,length=256))  # Blue
ColorLevels <- seq(0, 1, length=length(ColorRamp))
ColorLevelsMMs <- seq(0, max(forims), length=length(ColorRamps))
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
par(mar=c(5.5, 5, 4, 1.5))
image(x=sort(unique(simuls[,15])),y=sort(unique(simuls[,14])),z=forims,cex.axis=1.2,cex.lab=1.2,xlab="homogamy",ylab="nuclear selection",col=rgb(red=bonsenss,green=0,blue=1-bonsenss))
for (i in 1:dim(forims)[1])
  for (j in 1:dim(forims)[2])
  {
    text(x=sort(unique(simuls[,15]))[i],y=sort(unique(simuls[,14]))[j], labels=round(forims[i,j],digits=3),col="white")
  }
par(mar = c(5.5,0.5,4,5))
image(1, ColorLevels,axes=FALSE,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n",cex.lab=1.2)
axis(RIGHT <-4, at=seq(from=0,to=1,by=0.2), labels=seq(from=0,to=1,by=0.2), las= HORIZONTAL<-1, cex.axis=1)
mtext(text="IAL10",side=4,line=3,cex=1.2)
layout(1)

################################################################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/R01mAmHMt/")
simulR<-read.table(file="R01mAmHMtresults.txt",header=T)
Param<-matrix(data=c(0:17),nrow=18,ncol=1)
simulR$Param<-Param
simulR$Nselection<-matrix(data=c(rep(x=0.9,6),rep(0.5,6),rep(0.1,6)),nrow=18,ncol=1)
simulR$Homogamy<-matrix(data=rep(c(0.1,0.2,0.3,0.4,0.45,0.49),3),nrow=18,ncol=1)

library("reshape")#for melt
forgridR<-melt.data.frame(simulR,id.vars=c(14,15),measure.vars="IntAut10")
forgridR<-forgridR[order(forgridR$Nselection),]

forimR<-matrix(forgridR[,4],nrow=6,ncol=3)
bonsensR<-seq(min(forimR),max(forimR),by=1/15)

ColorRamp <- rgb( seq(0,1,length=256),  # Red
                  seq(0,0,length=256),  # Green
                  seq(1,0,length=256))  # Blue
ColorLevels <- seq(0, 1, length=length(ColorRamp))
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
par(mar=c(5.5, 5, 4, 1.5))
image(x=sort(unique(simulR[,15])),y=sort(unique(simulR[,14])),z=forimR,cex.axis=1.2,cex.lab=1.2,xlab="homogamy",ylab="nuclear selection",col=rgb(red=bonsensR,green=0,blue=1-bonsensR))
for (i in 1:dim(forimR)[1])
  for (j in 1:dim(forimR)[2])
  {
    text(x=sort(unique(simulR[,15]))[i],y=sort(unique(simulR[,14]))[j], labels=round(forimR[i,j],digits=3),col="white")
  }
par(mar = c(5.5,0.5,4,5))
image(1, ColorLevels,axes=FALSE,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n",cex.lab=1.2)
axis(RIGHT <-4, at=seq(from=0,to=1,by=0.2), labels=seq(from=0,to=1,by=0.2), las= HORIZONTAL<-1, cex.axis=1)
mtext(text="IAL10",side=4,line=3,cex=1.2)
layout(1)


################################################################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/wRmAmHMt/")
simulw<-read.table(file="wRmAmHMtresults.txt",header=T)
Param<-matrix(data=c(0:17),nrow=18,ncol=1)
simulw$Param<-Param
simulw$Nselection<-matrix(data=c(rep(x=0.9,6),rep(0.5,6),rep(0.1,6)),nrow=18,ncol=1)
simulw$Homogamy<-matrix(data=rep(c(0.1,0.2,0.3,0.4,0.45,0.49),3),nrow=18,ncol=1)

library("reshape")#for melt
forgridw<-melt.data.frame(simulw,id.vars=c(14,15),measure.vars="IntAut10")
forgridw<-forgridw[order(forgridw$Nselection),]

forimw<-matrix(forgridw[,4],nrow=6,ncol=3)
bonsensw<-seq(min(forimw),max(forimw),by=1/15)

ColorRamp <- rgb( seq(0,1,length=256),  # Red
                  seq(0,0,length=256),  # Green
                  seq(1,0,length=256))  # Blue
ColorLevels <- seq(0, 1, length=length(ColorRamp))
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
par(mar=c(5.5, 5, 4, 1.5))
image(x=sort(unique(simulw[,15])),y=sort(unique(simulw[,14])),z=forimw,cex.axis=1.2,cex.lab=1.2,xlab="homogamy",ylab="nuclear selection",col=rgb(red=bonsensw,green=0,blue=1-bonsensw))
for (i in 1:dim(forimw)[1])
  for (j in 1:dim(forimw)[2])
  {
    text(x=sort(unique(simulw[,15]))[i],y=sort(unique(simulw[,14]))[j], labels=round(forimw[i,j],digits=3),col="white")
  }
par(mar = c(5.5,0.5,4,5))
image(1, ColorLevels,axes=FALSE,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n",cex.lab=1.2)
axis(RIGHT <-4, at=seq(from=0,to=1,by=0.2), labels=seq(from=0,to=1,by=0.2), las= HORIZONTAL<-1, cex.axis=1)
mtext(text="IAL10",side=4,line=3,cex=1.2)
layout(1)

################################################################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/UniHomogMt/")
simuluni<-read.table(file="HMtresults.txt",header=T)
Param<-matrix(data=c(0:17),nrow=18,ncol=1)
simuluni$Param<-Param
simuluni$Nselection<-matrix(data=c(rep(x=1,6),rep(0.99,6),rep(0.9,6)),nrow=18,ncol=1)
simuluni$Homogamy<-matrix(data=rep(c(0.1,0.2,0.3,0.4,0.45,0.49),3),nrow=18,ncol=1)

library("reshape")#for melt
forgriduni<-melt.data.frame(simuluni,id.vars=c(14,15),measure.vars="IntAut10")
forgriduni<-forgriduni[order(forgriduni$Nselection),]

forimuni<-matrix(forgriduni[,4],nrow=6,ncol=3)
bonsensuni<-seq(min(forimuni),max(forimuni),by=1/15)

ColorRamp <- rgb( seq(0,1,length=256),  # Red
                   seq(0,0,length=256),  # Green
                   seq(1,0,length=256))  # Blue
ColorLevels <- seq(0, 1, length=length(ColorRamp))
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
par(mar=c(5.5, 5, 4, 1.5))
image(x=sort(unique(simuluni[,15])),y=sort(unique(simuluni[,14])),z=forimuni,cex.axis=1.2,cex.lab=1.2,xlab="homogamy",ylab="mt selection",col=rgb(red=bonsensuni,green=0,blue=1-bonsensuni))
for (i in 1:dim(forimuni)[1])
  for (j in 1:dim(forimuni)[2])
  {
    text(x=sort(unique(simuluni[,15]))[i],y=sort(unique(simuluni[,14]))[j], labels=round(forimuni[i,j],digits=3),col="white")
  }
par(mar = c(5.5,0.5,4,5))
image(1, ColorLevels,axes=FALSE,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n",cex.lab=1.2)
axis(RIGHT <-4, at=seq(from=0,to=1,by=0.2), labels=seq(from=0,to=1,by=0.2), las= HORIZONTAL<-1, cex.axis=1)
mtext(text="IAL10",side=4,line=3,cex=1.2)
layout(1)

################################################################################################################
setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/MultiHomogMt/")
simulmul<-read.table(file="HMtresults.txt",header=T)
Param<-matrix(data=c(0:17),nrow=18,ncol=1)
simulmul$Param<-Param
simulmul$Nselection<-matrix(data=c(rep(x=1,6),rep(0.99,6),rep(0.9,6)),nrow=18,ncol=1)
simulmul$Homogamy<-matrix(data=rep(c(0.1,0.2,0.3,0.4,0.45,0.49),3),nrow=18,ncol=1)

library("reshape")#for melt
forgridmul<-melt.data.frame(simulmul,id.vars=c(14,15),measure.vars="IntAut10")
forgridmul<-forgridmul[order(forgridmul$Nselection),]

forimmul<-matrix(forgridmul[,4],nrow=6,ncol=3)
bonsensmul<-seq(min(forimmul),max(forimmul),by=1/15)
ColorRamp <- rgb( seq(0,1,length=256),  # Red
                  seq(0,0,length=256),  # Green
                  seq(1,0,length=256))  # Blue
ColorLevels <- seq(0, 1, length=length(ColorRamp))
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
par(mar=c(5.5, 5, 4, 1.5))
image(x=sort(unique(simulmul[,15])),y=sort(unique(simulmul[,14])),z=forimmul,cex.axis=1.2,cex.lab=1.2,xlab="homogamy",ylab="mt selection",col=rgb(red=bonsensmul,green=0,blue=1-bonsensmul))
for (i in 1:dim(forimmul)[1])
  for (j in 1:dim(forimmul)[2])
  {
    text(x=sort(unique(simulmul[,15]))[i],y=sort(unique(simulmul[,14]))[j], labels=round(forimmul[i,j],digits=3),col="white")
  }
par(mar = c(5.5,0.5,4,5))
image(1, ColorLevels,axes=FALSE,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n",cex.lab=1.2)
axis(RIGHT <-4, at=seq(from=0,to=1,by=0.2), labels=seq(from=0,to=1,by=0.2), las= HORIZONTAL<-1, cex.axis=1)
mtext(text="IAL10",side=4,line=3,cex=1.2)
layout(1)
