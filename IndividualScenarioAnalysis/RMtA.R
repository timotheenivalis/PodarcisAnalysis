source(file="C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/FunctionsSimul.R")
#setwd("C:/Users/Thimothee Admin/Documents/external research/PodarcisCBGP/Renew simulations/lmBD/")
setwd("C:/Users/Timothée/Documents/Studies/PodarcisCBGP/New simulations/MtA/")

Param<-c(1,0.9975,0.995,0.9925,0.99,0.975,0.95)
Simuls<-c("MtA0","MtA1","MtA2","MtA3","MtA4","MtA5","MtA6")


WrapMtFunc(Simuls = Simuls,FilePrefix = "MtA",Param = Param)
