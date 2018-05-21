Za1 = model$biplot[,3]*model$analysis[1,1]
Za1<-Za1[-c(53:56)]
Za2 = model$biplot[,4]*model$analysis[2,1]
Za2<-Za2[-c(53:56)]
Za3 = model$biplot[,5]*model$analysis[3,1]
Za3<-Za3[-c(53:56)]
Za = Za1+Za2+Za3
rZ<-rank(Z)
rZi<-data.frame(Za1,Za2,Za3,Za,rZa)
rZi
