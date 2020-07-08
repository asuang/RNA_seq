setwd("/media/song/_dde_data/20200626_RNAseq_result/G_sample/G")
AD<-read.csv("2 GAD.csv",header = T,row.names=1 )
#install packages
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationHub")

library(AnnotationHub)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
#gene symbol convert to Entrez ID
hs<-org.Hs.eg.db
gene<-row.names(AU)
list<-AnnotationDbi::select(hs, keys = gene,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
geneid<-list[,2]
select(hs, keys = gene,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
ego<-enrichGO(gene = geneid,OrgDb = "org.Hs.eg.db",ont = "ALL",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T)
dotplot(ego,showCategory=20,font.size=10,title="GO analysis   shsy5y+beta(up-regulated genes)") #泡泡图
barplot(ego, showCategory=20,title="EnrichmentGO")  #柱状图
summary(ego)
write.csv(as.data.frame(summary(ego)),"GBD.csv")
