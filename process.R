setwd("~/Data/RNA_Seq/chendan/de_ribo/ballgown2_total_new") 
setwd("~/Data/RNA_Seq/chendan/de_ribo/ballgown2") 
library(ballgown)
phenotype<-read.csv("sta/phenotype.csv",header = T)

b<-ballgown(dataDir = "sta",samples=c("sta/N2_3sta","sta/N2_6ctr","sta/N2_6sta","sta/N2_12ctr","sta/N2_12sta")
            ,pData=phenotype)
g<-gexpr(b)
rm(b)
#BiocManager::install("org.Ce.eg.db")
library(org.Ce.eg.db)
col<-colnames(g)
new_cols<-c(col[c(2,4)],col[c(1,3,5)])
g<-g[,new_cols]
#screen genes  the FPKM more than 1 in wt(no repeat in samples)
g<-g[g[,1]>1 & g[,2]>1,]
#only choose the 6h control
g<-transform(g,FPKM.N2_12ctr=NULL)
#choose two control samples
g<-transform(g,con_avg=(FPKM.N2_6ctr+FPKM.N2_12ctr)/2)

# if g is the transcitions,you need do follow steps 
g<-texpr(b,meas="all")
rm(b)
g<-g[,-c(11,13,15,17,18,19)]
library(org.Ce.eg.db)
col<-colnames(g)
new_cols<-c(col[c(1:10)],col[c(1,3,5)])
g<-g[,new_cols]

## only choose 6h crl ,3h sta ,6h sta 
g<-data.frame(g)
## 3 time points ,6h control,3h starvation,6h starvation 
g<-g[,c(1,3,4)]
g<-transform(g,fold_change=FPKM.3_N2_6_h_Sta/FPKM.2_N2_6_h_Ctr)
up<-g[g$fold_change>1,]
down<-g[g$fold_change<1,]
name<-read.csv("/home/song/Data/Genomes/celegans/c.elegans_gene_name.csv",header = T)
up$Wb_id<-rownames(up)
down$Wb_id<-rownames(down)
up<-merge(up,name,by.x="Wb_id",by.y="GeneID",all.x=T)
down<-merge(down,name,by.x="Wb_id",by.y="GeneID",all.x=T)
col<-colnames(up)
up<-up[,c(col[1],col[6:9],col[2:5])]
col<-colnames(down)
down<-down[,c(col[1],col[6:9],col[2:5])]
up<-up[order(up$fold_change,decreasing = T),]
down<-down[order(down$fold_change),]
write.csv(up,"celegans_sta_up_636.csv")
write.csv(down,"celegans_sta_down_636.csv")
##function : screen the differential expression genes
differential<-function(x){
  #the columns number must match the treatment group
  x$min<-apply(x[,2:4],1,min)
  x$max<-apply(x[,2:4],1,max)
  x<-transform(x,down=log2((FPKM.N2_6ctr+1)/(max+1)))
  x<-transform(x,up=log2((min+1)/(FPKM.N2_6ctr+1)))
  #x<-x[x$max>1,] #threshold value use max value
  return(x)
}
G<-differential(g)
G<-transform(G,tem=log2((FPKM.N2_3sta+1)/(FPKM.N2_6ctr+1)))
#for the down-regulated genes ,wt_RPKM must be higher than 10
down<-G[G$FPKM.N2_6ctr>10 &G$down >0.5,]
#down<-G[G$down>1,]
down1<-G[G$FPKM.N2_6ctr>10 &-0.5>G$tem & G$down<=0.5, ] #the gene expressiong is up-regulated in 3h and down-regulated in 6h and 12h
#up<-G[G$up>1, ]  #normally screen 
#for the up-regulated genes , max in treatments mut be higher than 10
up<-G[G$max>10 & G$up>1,]
down<-rbind(down,down1)
write.csv(G,"raw_FPKM.csv",row.names = T)

library(clusterProfiler)
library(AnnotationHub)
library(AnnotationDbi)
library(ggplot2)
#function annotation
go<-function(x){
  db<-org.Ce.eg.db
  gene<-rownames(x)
  gene_list<-AnnotationDbi::select(db,keytype = "WORMBASE",keys = gene,columns = c("WORMBASE","ENTREZID"))
  ego<-enrichGO(gene_list[,2],db,keyType = "ENTREZID",ont = "ALL",pvalueCutoff = 0.05,readable = T)
  return(ego)
}
#go for down-regulated genes
fa<-go(down)
jpeg(filename = paste0("celegans_","starvation_down_genes(GO analysis).jpeg"),width = 1200,height = 1200,units = "px", bg="white",quality = 100)
d1<-dotplot(fa,showCategory=20,font.size=20)
d1<-d1+scale_color_continuous(low="dodgerblue",high="brown1")+scale_size(range=c(2,15))+
  labs(title = paste("GO analysis for down-regulated genes after starvation","\n","n =",nrow(down)))+
  theme(plot.title = element_text(face = "bold",hjust=0.5,size = 30),
        axis.title.x = element_text(face="italic", color="grey39", size=30),
        legend.title=element_text(size=28),legend.text=element_text(size=20))
print(d1)
dev.off()
tem<-data.frame(fa)
write.csv(tem,paste0("celegans_","starvation_","down(GO analysis).csv"))
#go for up-regulated genes
fa<-go(up)
jpeg(filename = paste0("celegans_","starvation_up_genes(GO analysis).jpeg"),width = 1200,height = 1200,units = "px", bg="white",quality = 100)
d1<-dotplot(fa,showCategory=20,font.size=20)
d1<-d1+scale_color_continuous(low="dodgerblue",high="brown1")+scale_size(range=c(2,15))+
  labs(title = paste("GO analysis for up-regulated genes after starvation","\n","n =",nrow(up)))+
  theme(plot.title = element_text(face = "bold",hjust=0.5,size = 30),
        axis.title.x = element_text(face="italic", color="grey39", size=30),
        legend.title=element_text(size=28),legend.text=element_text(size=20))
print(d1)
dev.off()
tem<-data.frame(fa)
write.csv(tem,paste0("celegans_","starvation_","up(GO analysis).csv"))
library(graph)
plotGOgraph(ego)

kegg<-function(x){
  gene<-rownames(x)
  db<-org.Ce.eg.db
  gene_list<-AnnotationDbi::select(db,keytype = "WORMBASE",keys = gene,columns = c("WORMBASE","ENTREZID"))
  ekg<-enrichKEGG(gene_list[,2],organism = "cel",keyType="ncbi-geneid",pvalueCutoff = 0.05)
  return(ekg)
}

#kegg for down-regulated genes
fa<-kegg(down)
jpeg(filename = paste0("celegans_","starvation_down_genes(KEGG analysis).jpeg"),width = 1200,height = 1200,units = "px", bg="white",quality = 100)
d1<-dotplot(fa,showCategory=20,font.size=20)
d1<-d1+scale_color_continuous(low="dodgerblue",high="brown1")+scale_size(range=c(2,15))+
  labs(title = paste("KEGG analysis for down-regulated genes after starvation","\n","n =",nrow(down)))+
  theme(plot.title = element_text(face = "bold",hjust=0.5,size = 30),
        axis.title.x = element_text(face="italic", color="grey39", size=30),
        legend.title=element_text(size=28),legend.text=element_text(size=20))
print(d1)
dev.off()
tem<-data.frame(fa)
write.csv(tem,paste0("celegans_","starvation_","down(KEGG analysis).csv"))

#kegg for up-regulated genes
fa<-kegg(up)
jpeg(filename = paste0("celegans_","starvation_up_genes(KEGG analysis).jpeg"),width = 1200,height = 1200,units = "px", bg="white",quality = 100)
d1<-dotplot(fa,showCategory=20,font.size=20)
d1<-d1+scale_color_continuous(low="dodgerblue",high="brown1")+scale_size(range=c(2,15))+
  labs(title = paste("KEGG analysis for up-regulated genes after starvation","\n","n =",nrow(up)))+
  theme(plot.title = element_text(face = "bold",hjust=0.5,size = 30),
        axis.title.x = element_text(face="italic", color="grey39", size=30),
        legend.title=element_text(size=28),legend.text=element_text(size=20))
print(d1)
dev.off()
tem<-data.frame(fa)
write.csv(tem,paste0("celegans_","starvation_","up(KEGG analysis).csv"))
browseKEGG(fa,"cel00350") 

#lysosome
list<-AnnotationDbi::select(org.Ce.eg.db,keytype = "WORMBASE",keys = rownames(G),columns = c("WORMBASE","SYMBOL"))
L<-read.table("~/Data/RNA_Seq/chendan/lysosome_gene.txt",header=F)
L<-merge(list,L,by.x="SYMBOL",by.y="V1")
l<-G[unique(L$WORMBASE),]
rownames(l)<-unique(L$SYMBOL)
#autophage
A<-read.table("~/Data/RNA_Seq/chendan/autophage_gene.txt",header=F)
A<-merge(list,A,by.x="SYMBOL",by.y="V1")
a<-G[unique(A$WORMBASE),]
rownames(a)<-unique(A$SYMBOL)

#par(mfrow = c(1,2))
#jpeg(filename = "class2.jpeg",width = 1000,height = 1200,units = "px", bg="white",quality = 100)
#gpar(mfrow=c(1,2))
down_order<-plot.sort.hm(down[,c(1:4)]) #the function in script ribo_code/plot_ss2m.R
up_order<-plot.sort.hm(up[,c(1:4)])
library(pheatmap)
jpeg(filename = paste0("celegans_","starvation_down_genes.jpeg"),width = 1200,height = 1200,units = "px", bg="white",quality = 100)
h1<-pheatmap(log2(as.matrix(down_order[,c(1:4)])+1),cluster_cols = F,cluster_rows = F,cellwidth=40,cellheight=2,show_rownames = F,
             angle_col = 45,fontsize = 20,fontsize_col = 20,
             border_color=NA,cutree_rows=5,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),
         main= paste("down-regulated genes after starvation","\n","n =",nrow(down_order)))
dev.off()

jpeg(filename = paste0("celegans_","starvation_up_genes.jpeg"),width = 1200,height = 1200,units = "px", bg="white",quality = 100)
h2<-pheatmap(log2(as.matrix(up_order[,c(1:4)])+1),cluster_cols = F,cluster_rows = F,cellwidth=40,cellheight=2,show_rownames = F,
             angle_col = 45,fontsize = 20,fontsize_col = 20,
             border_color=NA,cutree_rows=5,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),
             main=paste("up-regulated genes after starvation","\n","n =",nrow(up_order)))
dev.off()
#after clustering row,if you don't , skip two steps
up<-up[rownames(up[h2$tree_row[["order"]],]),]
down<-down[rownames(down[h1$tree_row[["order"]],]),]
#add the gene name 
name<-read.csv("/home/song/Data/Genomes/celegans/c.elegans_gene_name.csv",header = T)
up_order$gene<-rownames(up_order)
down_order$gene<-rownames(down_order)
gene_up<-rownames(up_order)
gene_down<-rownames(down_order)
up_order<-merge(name,up_order,by.x="GeneID",by.y="gene",all.y=T)
down_order<-merge(name,down_order,by.x="GeneID",by.y="gene",all.y=T)
rownames(up_order)<-up_order$GeneID
rownames(down_order)<-down_order$GeneID
up_order<-up_order[gene_up,]
down_order<-down_order[gene_down,]
write.csv(up_order,"sta_up.csv",row.names = F)
write.csv(down_order,"sta_down.csv",row.names = F)

#choose the non-coding rna 
ncrna<-function(x){
  n<-grep("MSTRG.*",rownames(x))
  new<-x[n,]
  return(new)
}
nc_up<-ncrna(up)
nc_down<-ncrna(down)
nc_down_order<-plot.sort.hm(nc_down[,c(1:4)])
nc_up_order<-plot.sort.hm(nc_up[,c(1:4)])
library(pheatmap)
jpeg(filename = paste0("celegans_","starvation_down_nc_genes(log2).jpeg"),width = 1200,height = 1200,units = "px", bg="white",quality = 100)
h1<-pheatmap(log2(as.matrix(nc_down_order[,c(1:4)])+1),cluster_cols = F,cluster_rows = T,cellwidth=40,cellheight=20,show_rownames = T,
             angle_col = 45,fontsize = 20,fontsize_col = 20,
             border_color=NA,cutree_rows=5,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),
             main= paste("down-regulated non-coding genes after starvation(log2)","\n","n =",nrow(nc_down_order)))
dev.off()

jpeg(filename = paste0("celegans_","starvation_up_nc_genes(log2).jpeg"),width = 1200,height = 1200,units = "px", bg="white",quality = 100)
h2<-pheatmap(log2(as.matrix(nc_up_order[,c(1:4)])+1),cluster_cols = F,cluster_rows = T,cellwidth=30,cellheight=18,show_rownames = T,
             angle_col = 45,fontsize = 20,fontsize_col = 20,
             border_color=NA,cutree_rows=5,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),
             main=paste("up-regulated non-coding genes after starvation(log2)","\n","n =",nrow(nc_up_order)))
dev.off()
nc_up_order<-nc_up_order[rownames(nc_up_order[h2$tree_row[["order"]],]),]
nc_down_order<-nc_down_order[rownames(nc_down_order[h1$tree_row[["order"]],]),]
write.csv(nc_up_order,"nc_sta_up.csv",row.names = T)
write.csv(nc_down_order,"nc_sta_down.csv",row.names = T)

readdata<-function(x){
  phenotype<-data.frame(ids=c("N2_12ctr",x),treat=c("con","mut"))
  b<-ballgown(dataDir = "mut",samples=c("mut/N2_12ctr",paste0("mut/",x)),pData=phenotype)
  g<-gexpr(b)
  rm(b)
  g<-g[g[,1]>1 & g[,2]>1,]
  g<-transform(g,foldchange_log2=log2(g[,2]/g[,1]))
  down<-g[-1>g$foldchange_log2,]
  up<-g[g$foldchange_log2>1,]
  write.csv(down,paste0("celegans","mutant_",x,"_down(log2FD<-1).csv"))
  write.csv(up,paste0("celegans","mutant_",x,"_up(log2FD>1).csv"))
  write.csv(g,paste0("celegans","mutant_",x,"_rawFPKM.csv"))
  lys<-read.table("~/Data/RNA_Seq/chendan/lysosome_gene.txt",header=F)
  lys<-merge(list,lys,by.x="SYMBOL",by.y="V1")
  Lys<-g[unique(lys$WORMBASE),]
  rownames(Lys)<-unique(lys$SYMBOL)
  jpeg(filename = paste(x,"_lysosome.jpeg"),width = 800,height = 900,units = "px", bg="white",quality = 100)
  h1<-pheatmap(log2(as.matrix(Lys[,c(1,2)])),cluster_cols = F,cluster_rows = T,cellwidth=20,cellheight=10,show_rownames = T,
               border_color=NA,cutree_rows=5,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),main="lysosomes genes")
  print(h1)
  dev.off()
  auto<-read.table("~/Data/RNA_Seq/chendan/autophage_gene.txt",header=F)
  auto<-merge(list,auto,by.x="SYMBOL",by.y="V1")
  Auto<-g[unique(auto$WORMBASE),]
  rownames(Auto)<-unique(auto$SYMBOL)
  jpeg(filename = paste(x,"_autophage.jpeg"),width = 800,height = 1200,units = "px", bg="white",quality = 100)
  h2<-pheatmap(log2(as.matrix(Auto[,c(1,2)])),cluster_cols = F,cluster_rows = T,cellwidth=20,cellheight=10,show_rownames = T,
               border_color=NA,cutree_rows=5,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),main="autophage genes")
  print(h2)
  dev.off()
  
}
samples<-c("6_qx666_Day_1","7_epg-5_Day_1","8_dpy-7_Day_1")
for (i in samples) {
  readdata(i)
}

#process mutant data
DA<-function(x){
  #note the sample name(line 264--line 266,line 265--line267)
  #phenotype<-data.frame(ids=c("4_N2_12_h_Ctr",x),treat=c("con","mut"))
  #b<-ballgown(dataDir = "mut",samples=c("mut/4_N2_12_h_Ctr",paste0("mut/",x)),pData=phenotype)
  phenotype<-data.frame(ids=c("N2_12ctr",x),treat=c("con","mut"))
  b<-ballgown(dataDir = "mut",samples=c("mut/N2_12ctr",paste0("mut/",x)),pData=phenotype)
  g<-gexpr(b)
  rm(b)
  g<-as.data.frame(g)
  #g$Max<-apply(g,1,max)
  #g<-g[g$Max>20,]
  #g<-transform(g,Max=NULL)
  #g<-g[g[,1]>1 & g[,2]>1,]
  g<-transform(g,foldchange_log2=log2(g[,2]/g[,1]))
  # down<-g[-1>g$foldchange_log2,]
  # up<-g[g$foldchange_log2>1,]
  #For down-regulated genes , RPKM wt must be higher than 10  
  down<-g[-0.5>g$foldchange_log2 & g$FPKM.N2_12ctr>10,]
  up<-g[g$foldchange_log2>1 & g$FPKM.dpy.7 >10,]
  #GO analysis
  pgo<-go(up)
  jpeg(filename = paste("celegans","mutant_",x,"_up_genes(GO analysis).jpeg"),width = 1800,height = 1200,units = "px", bg="white",quality = 100)
  p1<-dotplot(pgo,showCategory=20,font.size=20)
  p1<-p1+scale_color_continuous(low="dodgerblue",high="brown1")+scale_size(range=c(2,15))+
    labs(title = paste("GO analysis for up-regulated genes in",x,"\n","n =",nrow(up)))+
    theme(plot.title = element_text(face = "bold",hjust=0.5,size = 30),
          axis.title.x = element_text(face="italic", color="grey39", size=30),
          legend.title=element_text(size=28),legend.text=element_text(size=20))
  print(p1)
  dev.off()
  tem<-data.frame(pgo)
  write.csv(tem,paste0("celegans","mutant_",x,"_up(GO analysis).csv"))
  pgo<-go(down)
  jpeg(filename = paste("celegans","mutant_",x,"_down_genes(GO analysis).jpeg"),width = 2100,height = 1200,units = "px", bg="white",quality = 100)
  p2<-dotplot(pgo,showCategory=20,font.size=20)
  p2<-p2+scale_color_continuous(low="dodgerblue",high="brown1")+scale_size(range=c(2,15))+
    labs(title = paste("GO analysis for down-regulated genes in",x,"\n","n =",nrow(down)))+
    theme(plot.title = element_text(face = "bold",hjust=0.5,size = 30),
          axis.title.x = element_text(face="italic", color="grey39", size=30),
          legend.title=element_text(size=28),legend.text=element_text(size=20))
  print(p2)
  dev.off()
  tem<-data.frame(pgo)
  write.csv(tem,paste("celegans","mutant_",x,"_down(GO analysis).csv"))
  #KEGG analysis
  fa<-kegg(up)
  jpeg(filename = paste("celegans","mutant_",x,"_up_genes(kegg analysis).jpeg"),width = 1500,height = 1200,units = "px", bg="white",quality = 100)
  p3<-dotplot(fa,showCategory=20,font.size=20)
  p3<-p3+scale_color_continuous(low="dodgerblue",high="brown1")+scale_size(range=c(2,15))+
    labs(title = paste("KEGG analysis for up-regulated genes in",x,"\n","n =",nrow(up)))+
    theme(plot.title = element_text(face = "bold",hjust=0.5,size = 30),
          axis.title.x = element_text(face="italic", color="grey39", size=30),
          legend.title=element_text(size=28),legend.text=element_text(size=20))
  print(p3)
  dev.off()
  tem<-data.frame(fa)
  write.csv(tem,paste("celegans","mutant_",x,"_up(KEGG analysis).csv"))
  fa<-kegg(down)
  jpeg(filename = paste("celegans","mutant_",x,"_down_genes(kegg analysis).jpeg"),width = 1500,height = 1200,units = "px", bg="white",quality = 100)
  p4<-dotplot(fa,showCategory=20,font.size=20)
  p4<-p4+scale_color_continuous(low="dodgerblue",high="brown1")+scale_size(range=c(2,15))+
    labs(title = paste("KEGG analysis for down-regulated genes in",x,"\n","n =",nrow(down)))+
    theme(plot.title = element_text(face = "bold",hjust=0.5,size = 30),
          axis.title.x = element_text(face="italic", color="grey39", size=30),
          legend.title=element_text(size=28),legend.text=element_text(size=20))
  print(p4)
  dev.off()
  tem<-data.frame(fa)
  write.csv(tem,paste("celegans","mutant_",x,"_down(KEGG analysis).csv"))
  
  #plot heatmap(RPKM)
  jpeg(filename = paste("celegans","mutant_",x,"_up_genes(log2).jpeg"),width = 1000,height = 1300,units = "px", bg="white",quality = 100)
  h1<-pheatmap(log2(as.matrix(up[,c(1,2)])),cluster_cols = F,cluster_rows = T,cellwidth=30,cellheight=1,show_rownames = F,
               treeheight_row=200,angle_col = 45,fontsize = 20,fontsize_col = 20,
               border_color=NA,cutree_rows=4,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),main=paste("up-regulated genes(log2) in",x,"\n","n = ",nrow(up)))
  print(h1)
  dev.off()
  jpeg(filename = paste("celegans","mutant_",x,"_down_genes(log2).jpeg"),width =1000,height = 1300,units = "px", bg="white",quality = 100)
  h2<-pheatmap(log2(as.matrix(down[,c(1,2)])),cluster_cols = F,cluster_rows = T,cellwidth=30,cellheight=1,show_rownames = F,
               angle_col = 45,fontsize = 20,fontsize_col = 20,treeheight_row=200,
               border_color=NA,cutree_rows=4,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),main=paste("down-regulated genes(log2) in",x,"\n","n = ",nrow(down)))
  print(h2)
  dev.off()
  name<-read.csv("/home/song/Data/Genomes/celegans/c.elegans_gene_name.csv",header = T)
  g$gene<-rownames(g)
  up$gene<-rownames(up)
  down$gene<-rownames(down)
  g<-merge(name,g,by.x="GeneID",by.y="gene",all.y=T)
  up<-merge(name,up,by.x="GeneID",by.y="gene",all.y=T)
  down<-merge(name,down,by.x="GeneID",by.y="gene",all.y=T)
  rownames(up)<-up$GeneID
  rownames(down)<-down$GeneID
  up<-up[rownames(up[h1$tree_row[["order"]],]),]
  down<-down[rownames(down[h2$tree_row[["order"]],]),]
  write.csv(down,paste0("celegans ","mutant_",x,"_down(log2FD<-1).csv"),row.names = F)
  write.csv(up,paste0("celegans ","mutant_",x,"_up(log2FD>1).csv"),row.names = F)
  write.csv(g,paste0("celegans ","mutant_",x,"_rawFPKM.csv"))
}
samples<-c("6_qx666_Day_1","7_epg-5_Day_1","8_dpy-7_Day_1")
samples<-c("qx666","epg-5","dpy-7")
for (i in samples) {
  DA(i)
}
#non-coding genes 
nc_up<-ncrna(up)
nc_down<-ncrna(down)

library(pheatmap)
jpeg(filename = paste("celegans","mutant_",x,"_up_nc_genes(log2).jpeg"),width = 1000,height = 1300,units = "px", bg="white",quality = 100)
h1<-pheatmap(log2(as.matrix(nc_up[,c(1:2)])+1),cluster_cols = F,cluster_rows = T,cellwidth=40,cellheight=10,show_rownames = T,
             angle_col = 45,fontsize = 20,fontsize_col = 20,
             border_color=NA,cutree_rows=5,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),
             main=paste("up-regulated non-coding genes(log2) in",x,"\n","n = ",nrow(nc_up)))
dev.off()

jpeg(filename = paste("celegans","mutant_",x,"_down_nc_genes(log2).jpeg"),width =800,height = 1300,units = "px", bg="white",quality = 100)
h2<-pheatmap(log2(as.matrix(nc_down[,c(1:2)])+1),cluster_cols = F,cluster_rows = T,cellwidth=30,cellheight=20,show_rownames = T,
             angle_col = 45,fontsize = 20,fontsize_col = 20,
             border_color=NA,cutree_rows=3,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),
             main=paste("down-regulated non-coding genes(log2) in",x,"\n","n = ",nrow(nc_down)))
dev.off()
nc_up<-nc_up[rownames(nc_up[h1$tree_row[["order"]],]),]
nc_down<-nc_down[rownames(nc_down[h2$tree_row[["order"]],]),]
write.csv(nc_up,paste("celegans","mutant_",x,"_up_nc_genes(log2).csv"),row.names = T)
write.csv(nc_down,paste("celegans","mutant_",x,"_down_nc_genes(log2).csv"),row.names = T)

#merge the GO convertion name with RPKM table
x<-"celegans mutant_7_epg-5_Day_1_rawFPKM.csv"
tables<-read.csv(x,header = T)
tem<-read.table("tem.txt",sep="\t",header = T,quote = "\"\"")
tem<-tem[,-3]
colnames(tem)<-c("WB_ID","Symbol_ID","Gene_name")
all<-merge(tem,tables,by.x="WB_ID",by.y="X",all.y=T)
write.csv(all,x,row.names = F)
#g<-g[g$min>0,]
h1<-pheatmap(log(as.matrix(g[,c(1:5)])),cluster_cols = F,cluster_rows = T,cellwidth=40,cellheight=0.05,show_rownames = F,
             border_color=NA,cutree_rows=4,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),main="total genes")
y<-"celegans mutant_8_dpy-7_Day_1_up(log2FD>1).csv"
aa<-read.csv(y,header = T,row.names = 1)
aa<-aa[rownames(down),]
aa<-aa[rownames(up),]
write.csv(aa,y)
