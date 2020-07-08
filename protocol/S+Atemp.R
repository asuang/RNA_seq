# find the max of row(in treatment group) and means of con1 and con2

Sbeta$max = apply(Sbeta[,3:7], 1, max)
Sbeta<-transform(Sbeta,means=(con1+con2)/2)
Sbeta$min = apply(Sbeta[,3:7], 1, min)

#filter
Sbeta1<-Sbeta[Sbeta$max!=0 | Sbeta$means!=0,]
Sbeta2<-Sbeta[Sbeta$min!=0 | Sbeta$means!=0,]
Sbeta1<-transform(Sbeta1,min=NULL)
Sbeta2<-transform(Sbeta2,max=NULL)
Sbeta1<-transform(Sbeta1,ratio=log2((max+1)/(means+1)))
Sbeta2<-transform(Sbeta2,ratio=log2((min+1)/(means+1)))
SbetaU<-Sbeta1[Sbeta1$ratio>1,]
SbetaD<-Sbeta2[-1>Sbeta2$ratio,]
write.csv(SbetaU,"1 SbetaU.csv")
write.csv(SbetaD,"1 SbetaD.csv")
for (i  in seq(3,7)) {
  SbetaU[,i]=log2((SbetaU[,i]+1)/(SbetaU[,9]+1))
}
for (i  in seq(3,7)) {
  SbetaD[,i]=log2((SbetaD[,i]+1)/(SbetaD[,8]+1))
}

SbetaUL<-SbetaU[,3:7]
write.csv(SbetaUL,"2 SbetaU.csv")
SbetaUL<-as.matrix(SbetaUL)
SbetaDL<-SbetaD[,3:7]
write.csv(SbetaDL,"2 SbetaD.csv")
SbetaDL<-as.matrix(SbetaDL)

pheatmap(AalphaUL,cluster_cols = F,cluster_rows = F,color=colorRampPalette(c("navy","white","red"))(10))
pheatmap(AbetaDL,cluster_cols = F,cluster_rows = T,color=colorRampPalette(c("navy","white","red"))(10))
