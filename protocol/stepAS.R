setwd("/media/song/_dde_data/20200626_RNAseq_result/A_sample")

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
bg<-ballgown(dataDir = 'ballgown2',samplePattern = 'A',pData=phenotype)

#merge into a table(gene-level)
G<-gexpr(bg)
Galpha<-G[,c(1,4,5,6,7,8,9)]
Galpha<-as.data.frame(Galpha)
colnames(Galpha)<-c('con1','con2','alpha15','alpha30','alpha60','alpha90','alpha120')
#note!!!here using "alpha" to replace beta!!
Galpha<-G[,c(1,4,10,11,12,2,3)]
Galpha<-as.data.frame(Galpha)
colnames(Galpha)<-c('con1','con2','beta15','beta30','beta60','beta90','beta120')

# find the max of row(in treatment group) and means of con1 and con2
Galpha$max = apply(Galpha[,3:7], 1, max)
Galpha<-transform(Galpha,means=(con1+con2)/2)
Galpha$min = apply(Galpha[,3:7], 1, min)

#filter
Galpha1<-Galpha[Galpha$max!=0 | Galpha$means!=0,]
Galpha2<-Galpha[Galpha$min!=0 | Galpha$means!=0,]
Galpha1<-transform(Galpha1,min=NULL)
Galpha2<-transform(Galpha2,max=NULL)
Galpha1<-transform(Galpha1,ratio=log2((max+1)/(means+1)))
Galpha2<-transform(Galpha2,ratio=log2((min+1)/(means+1)))
GalphaU<-Galpha1[Galpha1$ratio>3,]
GalphaD<-Galpha2[-2>Galpha2$ratio,]
write.csv(GalphaU,"1 AalphaU.csv")
write.csv(GalphaD,"1 AalphaD.csv")
for (i  in seq(3,7)) {
  GalphaU[,i]=log2((GalphaU[,i]+1)/(GalphaU[,9]+1))
}
for (i  in seq(3,7)) {
  GalphaD[,i]=log2((GalphaD[,i]+1)/(GalphaD[,8]+1))
}

GalphaUL<-GalphaU[,3:7]
write.csv(GalphaUL,"2 AalphaU.csv")
GalphaUL<-as.matrix(GalphaUL)
GalphaDL<-GalphaD[,3:7]
write.csv(GalphaDL,"2 AalphaD.csv")
GalphaDL<-as.matrix(GalphaDL)

pheatmap(GalphaUL,cluster_cols = F,cluster_rows = T,show_rownames = F,treeheight_row=0,
         cellwidth=20,cellheight = 4,main="A  Up",main.fontsize=20,color=colorRampPalette(c("navy","white","red"))(6))
pheatmap(GalphaDL,cluster_cols = F,cluster_rows = T,show_rownames = F,treeheight_row=0,
         cellwidth=20,cellheight = 3,main="A Down",main.fontsize=20,color=colorRampPalette(c("navy","white","red"))(6))

#filt non-coding RNA
list<-c("RNA","MIR","SNO","LINC","SNHG","LOC","RNVU","RNU","AS")
for (i in list) {
  Galpha1<-Galpha1[-grep(i, row.names(Galpha1)),]
}
for (i in list) {
  Galpha2<-Galpha2[-grep(i, row.names(Galpha2)),]
}

write.csv(GalphaUL,"3 AalphaU.csv")
write.csv(GalphaDL,"3 AalphaD.csv")
