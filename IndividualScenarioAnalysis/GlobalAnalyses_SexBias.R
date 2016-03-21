#######Sex-Biased hybrid fitness###########
setwd(dir="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations")
simulBSH<-read.table(file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/BSH/BSHresults.txt",header=T)
simulDisp<-read.table(file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/Disp/BDresults.txt")
simulHomog<-read.table(file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/Homogamy/Homogamyresults.txt")

simulDispNofix<-read.table(file="C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/New simulations/Disp/BDresultsnofix.txt")


simulBSHc<-rbind(simulBSH[1:3,],simulDisp["BD0",],simulBSH[4:6,])

Param<-matrix(data=c(0.001,0.05,0.2,0.3,0.4,0.55,0.599),nrow=7,ncol=1)
Param<-Param/(0.6-Param)
logParam<-log(Param)
simulBSHc$Bsurv<-Param
simulBSHc$logBsurv<-logParam

plot(simulBSHc$logBsurv,simulBSHc$IntAut,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),xaxt="n",pch=3,lwd=2,xlab="female hybrid fitness over male hybrid fitness",ylab="Mean proportion",main="")
points(jitter(simulBSHc$logBsurv),simulBSHc$IntAut10,col="green",pch=4,lwd=2,type="b")
points(simulBSHc$logBsurv,simulBSHc$FixAut,col="orange",pch=5,lwd=2,type="b")
points(simulBSHc$logBsurv,simulBSHc$MeanExoAut,lwd=2,pch=6,type="b")
axis(BOTTOM <-1, at=simulBSHc$logBsurv, labels=round(simulBSHc$Bsurv,digits=3), las= HORIZONTAL<-1, cex.axis=1)

#legend(x="topleft",legend=c("ALI_0","ALI_0.1","ALI_1","AI"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)

c(expression(paste("log du rapport des ",sigma^2)))

plot(x=simulBSHc$logBsurv,y=simulBSHc$Introgq50,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),xaxt="n",pch=3,lwd=2,xlab="female hybrid fitness over male hybrid fitness",ylab="Median proportion and CI",main="")
points(jitter(simulBSHc$logBsurv),simulBSHc$Introg10q50,col="green",pch=4,lwd=2,type="b")
points(simulBSHc$logBsurv,simulBSHc$Fixq50,col="orange",pch=5,lwd=2,type="b")
points(simulBSHc$logBsurv,simulBSHc$MeanExoq50,lwd=2,pch=6,type="b")
axis(BOTTOM <-1, at=simulBSHc$logBsurv, labels=round(simulBSHc$Bsurv,digits=3), las= HORIZONTAL<-1, cex.axis=1)

for(i in 1:dim(simulBSHc)[1])
{
  points(x=c(simulBSHc$logBsurv[i],simulBSHc$logBsurv[i]),y=c(simulBSHc$Introgq025[i],simulBSHc$Introgq975[i]),type="l",col="dark green")
  points(x=c(simulBSHc$logBsurv[i],simulBSHc$logBsurv[i]),y=c(simulBSHc$Introg10q025[i],simulBSHc$Introg10q975[i]),type="l",col="green")
  points(x=c(simulBSHc$logBsurv[i],simulBSHc$logBsurv[i]),y=c(simulBSHc$Fixq025[i],simulBSHc$Fixq975[i]),type="l",col="orange")
  points(x=c(simulBSHc$logBsurv[i],simulBSHc$logBsurv[i]),y=c(simulBSHc$MeanExoq025[i],simulBSHc$MeanExoq975[i]),type="l",col="black")
}
#legend(x="topleft",legend=c("ALI_0","ALI_0.1","ALI_1","AI"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)

plot(simulBSHc$logBsurv,simulBSHc$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),xaxt="n",pch=1,lwd=2,xlab="female hybrid fitness over male hybrid fitness",ylab="Simulation proportion",cex.lab=1,cex.axis=1)
points(simulBSHc$logBsurv,simulBSHc$IntMt,col="blue",pch=2,lwd=2,type="b")
#legend(x="topleft",legend=c("Mt_1    ","Mt_0    "),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)
axis(BOTTOM <-1, at=simulBSHc$logBsurv, labels=round(simulBSHc$Bsurv,digits=3), las= HORIZONTAL<-1, cex.axis=1)


#######Sex-Biased hybrid dispersal###########
Param<-matrix(data=c(0.006,0.052,0.626,1,1.6,19.11,181.25),nrow=7,ncol=1)
logParam<-log(Param)
simulDisp$Bdisp<-Param
simulDisp$logBdisp<-logParam
simulDispNofix$Bdisp<-Param
simulDispNofix$logBdisp<-logParam

plot(x=simulDisp$logBdisp,y=simulDisp$Introgq50,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),xaxt="n",pch=3,lwd=2,xlab=c(expression(paste(sigma^2,"ratio female/male"))),ylab="Median proportion",main="")
points(simulDisp$logBdisp,simulDisp$Introg10q50,col="green",pch=4,lwd=2,type="b")
points(simulDisp$logBdisp,simulDisp$Fixq50,col="orange",pch=5,lwd=2,type="b")
points(simulDisp$logBdisp,simulDisp$MeanExoq50,lwd=2,pch=6,type="b")
axis(BOTTOM <-1, at=simulDisp$logBdisp, labels=round(simulDisp$Bdisp,digits=3), las= HORIZONTAL<-1, cex.axis=1)

for(i in 1:dim(simulDisp)[1])
{
  points(x=c(simulDisp$logBdisp[i],simulDisp$logBdisp[i]),y=c(simulDisp$Introgq025[i],simulDisp$Introgq975[i]),type="l",col="dark green")
  points(x=c(simulDisp$logBdisp[i],simulDisp$logBdisp[i]),y=c(simulDisp$Introg10q025[i],simulDisp$Introg10q975[i]),type="l",col="green")
  points(x=c(simulDisp$logBdisp[i],simulDisp$logBdisp[i]),y=c(simulDisp$Fixq025[i],simulDisp$Fixq975[i]),type="l",col="orange")
  points(x=c(simulDisp$logBdisp[i],simulDisp$logBdisp[i]),y=c(simulDisp$MeanExoq025[i],simulDisp$MeanExoq975[i]),type="l",col="black")
}
legend(x="topleft",legend=c("ALI_0","ALI_0.1","ALI_1","AI"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)

plot(simulDisp$logBdisp,simulDisp$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),xaxt="n",pch=1,lwd=2,xlab=c(expression(paste(sigma^2,"ratio female/male"))),cex.lab=1.2,ylab="Simulation proportion")
points(simulDisp$logBdisp,simulDisp$IntMt,col="blue",pch=2,lwd=2,type="b")
legend(x="topleft",legend=c("Mt_1    ","Mt_0    "),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)
axis(BOTTOM <-1, at=simulDisp$logBdisp, labels=round(simulDisp$Bdisp,digits=3), las= HORIZONTAL<-1, cex.axis=1)

plot(x=simulDisp$logBdisp,y=simulDisp$IntAut10,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),xaxt="n",pch=3,lwd=2,xlab=c(expression(paste(sigma^2,"ratio female/male"))),ylab="Median proportion and CI",main="")
points(x=simulDispNofix$logBdisp,y=simulDispNofix$IntAut10,type="b")

plot(x=simulDisp$logBdisp,y=simulDisp$MeanExoAut,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),xaxt="n",pch=3,lwd=2,xlab=c(expression(paste(sigma^2,"ratio female/male"))),ylab="Median proportion and CI",main="")
points(x=simulDispNofix$logBdisp,y=simulDispNofix$MeanExoAut,type="b")

#######Sex-Biased inter-taxa crosses###########
simulHomogc1<-rbind(simulDisp["BD0",-c(25,26)],simulHomog[1:5,])
Param<-matrix(data=c(0,0.05,0.1,0.3,0.5,0.99),nrow=6,ncol=1)
logParam<-log(Param+1)
simulHomogc1$Bmat<-Param
plot(x=simulHomogc1$Bmat,y=simulHomogc1$Introgq50,col="dark green",type="b",yaxp=c(0,1,10),ylim=c(0,1.4),xaxt="n",pch=3,lwd=2,xlab="Male of taxa 1 / female of taxa 2 mating refusal rate",ylab="Median proportion and CI",main="")
points(simulHomogc1$Bmat,simulHomogc1$Introg10q50,col="green",pch=4,lwd=2,type="b")
points(simulHomogc1$Bmat,simulHomogc1$Fixq50,col="orange",pch=5,lwd=2,type="b")
points(simulHomogc1$Bmat,simulHomogc1$MeanExoq50,lwd=2,pch=6,type="b")
axis(BOTTOM <-1, at=simulHomogc1$Bmat, labels=round(simulHomogc1$Bmat,digits=3), las= HORIZONTAL<-1, cex.axis=1)

for(i in 1:dim(simulHomogc1)[1])
{
  points(x=c(simulHomogc1$Bmat[i],simulHomogc1$Bmat[i]),y=c(simulHomogc1$Introgq025[i],simulHomogc1$Introgq975[i]),type="l",col="dark green")
  points(x=c(simulHomogc1$Bmat[i],simulHomogc1$Bmat[i]),y=c(simulHomogc1$Introg10q025[i],simulHomogc1$Introg10q975[i]),type="l",col="green")
  points(x=c(simulHomogc1$Bmat[i],simulHomogc1$Bmat[i]),y=c(simulHomogc1$Fixq025[i],simulHomogc1$Fixq975[i]),type="l",col="orange")
  points(x=c(simulHomogc1$Bmat[i],simulHomogc1$Bmat[i]),y=c(simulHomogc1$MeanExoq025[i],simulHomogc1$MeanExoq975[i]),type="l",col="black")
}
legend(x="topleft",legend=c("ALI_0","ALI_0.1   ","ALI_1","AI"),col=c("dark green","green","orange","black"),pch=c(3,4,5,6),ncol=1,pt.lwd=2,cex=1)

plot(simulHomogc1$Bmat,simulHomogc1$FixMt,col="red",type="b",yaxp=c(0,1,10),ylim=c(0,1.3),xaxt="n",pch=1,lwd=2,xlab="Male of taxa 1 / female of taxa 2 mating refusal rate",cex.lab=1.2,ylab="Simulation proportion")
points(simulHomogc1$Bmat,simulHomogc1$IntMt,col="blue",pch=2,lwd=2,type="b")
legend(x="topleft",legend=c("Mt_1    ","Mt_0    "),col=c("red","blue"),pch=c(1,2),ncol=1,pt.lwd=2,cex=1.2)
axis(BOTTOM <-1, at=simulHomogc1$Bmat, labels=round(simulHomogc1$Bmat,digits=3), las= HORIZONTAL<-1, cex.axis=1)
