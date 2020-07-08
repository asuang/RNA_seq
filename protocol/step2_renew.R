setwd("/media/song/_dde_data/20200626_RNAseq_result/G_sample")

install.packages('dplyr')
install.packages("devtools")
install.packages("ggplot2")
BiocManager::install("genefilter")
BiocManager::install("ballgown")
install.packages("pheatmap")

library(dplyr)
library(devtools)
library(genefilter)
library(ballgown)
devtools::install_github('alyssafrazee/RSkittleBrewer')
library(RSkittleBrewer)
library(ggplot2)
library(pheatmap)

phenotype<-read.csv("phenotype.csv",header = T)
bg<-ballgown(dataDir = 'ballgown2',samplePattern = 'G',pData=phenotype)

#merge into a table(gene-level)
G<-gexpr(bg)
Galpha<-G[,c(1,5,6,7,8,9)]
Galpha<-as.data.frame(Galpha)
colnames(Galpha)<-c('con','A1','A2','A3','A4','A5')
Galpha<-G[,c(10,11,12,2,3,4)]
Galpha<-as.data.frame(Galpha)
colnames(Galpha)<-c('con','B1','B2','B3','B4','B5')

# find the max of row(in treatment group) and means of con1 and con2
Galpha$max = apply(Galpha[,2:6], 1, max)
Galpha<-transform(Galpha,means=(con1+con2)/2)
Galpha$min = apply(Galpha[,2:6], 1, min)

#filter
Galpha1<-Galpha[Galpha$max!=0 | Galpha$con!=0,]
Galpha2<-Galpha[Galpha$min!=0 | Galpha$con!=0,]
Galpha1<-transform(Galpha1,min=NULL)
Galpha2<-transform(Galpha2,max=NULL)
Galpha1<-transform(Galpha1,ratio=log2((max+1)/(con+1)))
Galpha2<-transform(Galpha2,ratio=log2((min+1)/(con+1)))
GalphaU<-Galpha1[Galpha1$ratio>1,]
GalphaD<-Galpha2[-1>Galpha2$ratio,]
write.csv(GalphaU,"1 SAU.csv")
write.csv(GalphaD,"1 SAD.csv")
for (i  in seq(2,6)) {
  GalphaU[,i]=log2((GalphaU[,i]+1)/(GalphaU[,1]+1))
}
for (i  in seq(2,6)) {
  GalphaD[,i]=log2((GalphaD[,i]+1)/(GalphaD[,1]+1))
}

GalphaUL<-GalphaU[,2:6]
write.csv(GalphaUL,"2 SAU.csv")
GalphaUL<-as.matrix(GalphaUL)
GalphaDL<-GalphaD[,2:6]
write.csv(GalphaDL,"2 SAD.csv")
AalphaDL<-as.matrix(AalphaDL)

pheatmap(GalphaUL,cluster_cols = F,cluster_rows = T,show_rownames = F,treeheight_row=0,
         cellwidth=20,cellheight = 3,main="B  Up",main.fontsize=20,color=colorRampPalette(c("navy","white","red"))(6))
pheatmap(GalphaDL,cluster_cols = F,cluster_rows = T,show_rownames = F,treeheight_row=0,
         cellwidth=20,cellheight = 4,main="B Down",main.fontsize=20,color=colorRampPalette(c("navy","white","red"))(6))

#filt non-coding RNA
list<-c("RNA","MIR","SNO","LINC","SNHG","LOC","RNVU","RNU","AS")
for (i in list) {
  Galpha1<-Galpha1[-grep(i, row.names(Galpha1)),]
}
for (i in list) {
  Galpha2<-Galpha2[-grep(i, row.names(Galpha2)),]
}
write.csv(SalphaUL,"3 SalphaU.csv")
write.csv(SalphaDL,"3 SalphaD.csv")
