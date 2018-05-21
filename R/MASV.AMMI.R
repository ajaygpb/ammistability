MASV.AMMI<-function(model)
{
  A<-model$biplot
  A<-A[A[,1]=="GEN",-c(1,2)]
  pc1<-model$analysis[1,4]/model$analysis[2,4]
  pc2<-model$analysis[2,4]/model$analysis[3,4]
  MASV<-apply(A,1,function(x) sqrt((pc1*(x[1])^2+(x[2])^2)+(pc2*(x[2])^2+(x[3])^2)))
  rk<-rank(MASV)
  B<-model$means
  W<-tapply.stat(B[,3],B[,2],function(x) mean(x,rm.na=TRUE))
  Rx<-rank(-W[,2])
  YSI<-rk+Rx
  ranking<-data.frame(MASV,YSI,rMASV=rk,rYSI=Rx,means=W[,2])
  invisible(ranking)
}
