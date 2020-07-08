setwd("/media/song/_dde_data/20200626_RNAseq_result/S_sample")
#
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
#set path(next command the option dataDir set filename or data_dir,but when setting data_dir, it cannot work)
data_dir<-system.file('ballgown',package = 'ballgown')

bg<-ballgown(dataDir = 'ballgown2',samplePattern = 'S',pData=phenotype)
#if warnings:Rows of pData did not seem to be in the same order as the columns of the expression data. Attempting to rearrange pData...
#using list.files("ballgown") to know the order of ballgown files,and writing the phenotype.txt/csv
all(phenotype$Ids==list.files("ballgown"))

#merge into a table(transcription-level,contain transcription variance)
S<-texpr(bg,meas = 'all')
A_filt<-A[,c(10,12,14,16,18,20,22,24,26,28,30,32,34)]
#merge into a table(gene-level)
A<-gexpr(bg)
Aalpha<-A[,c(1,4,5,6,7,8,9)]
Aalpha<-as.data.frame(Aalpha)
colnames(Aalpha)<-c('con1','con2','alpha15','alpha30','alpha60','alpha90','alpha120')
Abeta<-A[,c(1,4,2,3,10,11,12)]

#filter(ask for whether to do this step) 
bg_filt<-subset(bg,"rowVars(texpr(bg))>0.1",genomesubset=T)
#statistibal test differential expression
result_transcript <- stattest(bg,feature = "transcript",getFC = F,meas = "FPKM",covariate = "treat")
result_transcripts<-data.frame(geneNames=ballgown::geneNames(bg), geneIDs = ballgown::geneIDs(bg), result_transcript)
sig_transcripts=subset(result_transcripts,result_transcript$pval<0.05)

# find the max of row

Aalpha$max = apply(Aalpha[,3:7], 1, max)
Aalpha<-transform(Aalpha,means=(con1+con2)/2)
Aalpha$min = apply(Aalpha[,1:7],1,min)

Aalpha<-transform(Aalpha,min=NULL)
#filter
Aalpha1<-Aalpha[Aalpha$max!=0&Aalpha$means!=0,]
Aalpha1<-transform(Aalpha1,ratio=log2((max+1)/(means+1)))
Aalphaex<-Aalpha1[-grep("MSTRG", row.names(Aalpha1)),]
AalphaU<-Aalpha1[Aalpha1$ratio>2,]
AalphaD<-Aalpha1[Aalpha1$ratio<-1,]
for (i  in seq(3,7)) {
  AalphaU[,i]=log2((AalphaU[,i]+1)/(AalphaU[,9]+1))
}
AalphaU<-Aalpha1[Aalpha1$max>100,]
Aalphaex<-Aalphaex[Aalphaex$ratio!=0,]
Aalphaex<-Aalphaex[order(-Aalphaex[,10]),]

Aalpha4<-Aalpha1[order(Aalpha1[,10]),]
AalphaU<-Aalpha2[c(1:200),]
AalphaD<-Aalpha4[c(1:200),]
AalphaM<-rbind(AalphaU,AalphaU)
#Aalpha1<-Aalpha[Aalpha$max>10&Aalpha$min>5,]
AalphaM<-AalphaM[,c(3:7,9)]
AalphaUL<-AalphaU[,3:7]
AalphaUL<-as.matrix(AalphaUL)
AalphaM<-log2(AalphaM+1)
pheatmap(AalphaUL,cluster_cols = T,cluster_rows = T,color=colorRampPalette(c("navy","white","red"))(10))



