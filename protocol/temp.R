# find the max of row(in treatment group) and means of con1 and con2

Gbeta$max = apply(Gbeta[,2:6], 1, max)
Gbeta<-transform(Gbeta,means=(con1+con2)/2)
Gbeta$min = apply(Gbeta[,2:6], 1, min)

#filter
Gbeta1<-Gbeta[Gbeta$max!=0 | Gbeta$con!=0,]
Gbeta2<-Gbeta[Gbeta$min!=0 | Gbeta$con!=0,]
Gbeta1<-transform(Gbeta1,min=NULL)
Gbeta2<-transform(Gbeta2,max=NULL)
Gbeta1<-transform(Gbeta1,ratio=log2((max+1)/(con+1)))
Gbeta2<-transform(Gbeta2,ratio=log2((min+1)/(con+1)))
GbetaU<-Gbeta1[Gbeta1$ratio>1,]
GbetaD<-Gbeta2[-2>Gbeta2$ratio,]
write.csv(GbetaU,"1 GBU.csv")
write.csv(GbetaD,"1 GBD.csv")
for (i  in seq(2,6)) {
  GbetaU[,i]=log2((GbetaU[,i]+1)/(GbetaU[,1]+1))
}
for (i  in seq(2,6)) {
  GbetaD[,i]=log2((GbetaD[,i]+1)/(GbetaD[,1]+1))
}

GbetaUL<-GbetaU[,2:6]
write.csv(GbetaUL,"2 GBU.csv")
SbetaUL<-as.matrix(SbetaUL)
GbetaDL<-GbetaD[,2:6]
write.csv(GbetaDL,"2 GBD.csv")
SbetaDL<-as.matrix(SbetaDL)

pheatmap(AalphaUL,cluster_cols = F,cluster_rows = F,color=colorRampPalette(c("navy","white","red"))(10))
pheatmap(GbetaDL,cluster_cols = F,cluster_rows = T,color=colorRampPalette(c("navy","white","red"))(10))
