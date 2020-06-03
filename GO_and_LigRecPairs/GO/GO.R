# GO enrichment for selected FAPs markers up and downregulated
# If you use clusterProfiler in published research, please cite:
#  Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.

#previously ==> BiocManager::install("org.Mm.eg.db") 
library("clusterProfiler")
library(AnnotationDbi)
library(cowplot)
library(ggplot2)
library(dplyr)
  library(org.Mm.eg.db) #mouse genbank annotations

prloc="~/reanalysesPubs/"
setwd(prloc)
# get markers
fullse <- read.table("Opres_DeMichFAPSonly.100markers.17.txt", sep="\t", header=T)
typeof(fullse$cluster) # integer
# exclude cluster "10" because its Tenocytes
senp <- fullse %>% subset(cluster!=10)

# log2fold-changes in geneList:
geneList <- senp$avg_logFC
# gene to entrezID:
#detach("package:dplyr", unload = TRUE) # detach dplyr then "select" or use ensembldb::select
mouseannot <- org.Mm.eg.db
entrezid = ensembldb::select(mouseannot,
                  keys = as.character(senp$gene),
                  columns = c("ENTREZID","SYMBOL"),
                  keytype= "SYMBOL")
de <- entrezid$ENTREZID
names(geneList) <- de
# count genes in each cluster:
c0 =length(senp$cluster[senp$cluster==0])
c1 = length(senp$cluster[senp$cluster==1])
c5 = length(senp$cluster[senp$cluster==5])
c12 =length(senp$cluster[senp$cluster==12])
c15 = length(senp$cluster[senp$cluster==15])
mydf <- data.frame(Entrez=de, 
                   group= c(rep("0:Smoc2+",c0),rep("1:Pi16",c1),rep("5:Smoc2+",c5),rep("12:Smoc++",c12),rep("15:Smoc2,FAPS??",c15))
                   )
xx.formula <- compareCluster(Entrez~group,data=mydf,
                             fun="enrichGO", ont="all", OrgDb="org.Mm.eg.db")

pdf("GO/BP_CC_MF_ensemble.pdf", width=20,height = 15)
dotplot(xx.formula, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~.,scale="free") +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 9)) +
  labs(x="cluster label (number of genes)")
dev.off()  #errors in graph :( so,

# version step by step:
xx.formula2 <- compareCluster(Entrez~group,data=mydf,
                             fun="enrichGO",ont="BP", OrgDb="org.Mm.eg.db")
xx.formula3 <- compareCluster(Entrez~group,data=mydf,
                              fun="enrichGO",ont="CC", OrgDb="org.Mm.eg.db")
xx.formula4 <- compareCluster(Entrez~group,data=mydf,
                              fun="enrichGO",ont="MF", OrgDb="org.Mm.eg.db")

pdf("GO/plots_GOenrichmentDay0_DMichOprescuXX.pdf", width=15 )
par(mfrow=c(3,1))
dotplot(xx.formula2)+
  theme(axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 11))+
  ggtitle("Biological Process")
dotplot(xx.formula3)+
  theme(axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 11))+
  ggtitle("Cellular component")
dotplot(xx.formula4)+
  theme(axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 11))+
  labs(x="cluster label (number of genes)")+
  ggtitle("Molecular Function")
dev.off()

li <- list()
fc <- list()
res <- list()
li[[1]] = mydf$Entrez[mydf$group=="0:Smoc2+"]
fc[[1]] = geneList[1:30]
li[[2]] = mydf$Entrez[mydf$group=="1:Pi16"]
fc[[2]] = geneList[31:60]
li[[3]] = mydf$Entrez[mydf$group=="5:Smoc2+"]
fc[[3]] = geneList[61:90]
li[[4]] = mydf$Entrez[mydf$group=="12:Smoc++"]
fc[[4]] = geneList[91:97]
li[[5]] = mydf$Entrez[mydf$group=="15:Smoc2,FAPS??"]
fc[[5]] = geneList[98:104]

for (i in 1:5){
  res[[i]] = enrichGO(li[[i]], ont="BP", OrgDb="org.Mm.eg.db", readable=TRUE)
}

ho <- enrichGO(li[[1]], ont="BP", OrgDb="org.Mm.eg.db", readable=TRUE)
he <- enrichGO(li[[1]], ont="BP", OrgDb="org.Mm.eg.db")
cnetplot(ho, foldChange=fc[[1]], circular=TRUE, colorEdge = TRUE)
cnetplot(he, foldChange=fc[[1]], circular=TRUE, colorEdge = TRUE)
#cnetplots:
pdf('GO/cnetsB.pdf',width=9)
cnetplot(res[[1]], foldChange = fc[[1]], circular=TRUE, colorEdge=TRUE)+
  ggtitle("0:Smoc2+")
cnetplot(res[[2]], foldChange = fc[[2]], circular=TRUE, colorEdge=TRUE)+
  ggtitle("1:Pi16")
cnetplot(res[[3]], foldChange = fc[[3]], circular=TRUE, colorEdge=TRUE)+
  ggtitle("5:Smoc2+")
cnetplot(res[[4]], foldChange = fc[[4]], circular=TRUE, colorEdge=TRUE)+
  ggtitle("12:Smoc2++")
cnetplot(res[[5]], foldChange = fc[[5]], circular=TRUE, colorEdge=TRUE)+
  ggtitle("15:Smoc2,FAPS??")
dev.off()


# le groupe 5 est étrange: analyse tout seul

entrezid5 <- mydf[mydf$group=="5:Smoc2+",]
ego5 <- enrichGO(entrezid5$Entrez, OrgDb = "org.Mm.eg.db", ont="all", readable=TRUE)
pdf("GO/only_group5.pdf", width=15)
dotplot(ego5, split="ONTOLOGY") + facet_grid(ONTOLOGY~.,scale="free") +
  theme(axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 11))+
  ggtitle("Only group 5 : Smoc2+")
dev.off()