EV1<-model$biplot[,3]^2
EV2<-model$biplot[,4]^2
EV3<-model$biplot[,5]^2
EV<-EV1+EV2+EV3
rEV<-rank(EV)
rEV<-data.frame(EV1,EV2,EV3,EV,rEV)
rEV
