# GO enrichment
# If you use clusterProfiler in published research, please cite:
#  Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.

# TOY EXEMPLE WITH ONLY 3 GENES

#previously ==> BiocManager::install("org.Mm.eg.db") 
library("clusterProfiler")
library(AnnotationDbi)
library(cowplot)
library(ggplot2)
#using genome annotation for mouse: org.Mm.eg.db
mouseannot <- org.Mm.eg.db
# log2fold-changes in geneList:
geneList <- c(2.36354, 2.04829, 1.80765, 2.2143, 2.11138, 2.1023)
genesymbols = c("Cxcl14", "Smoc2", "Hsd11b1", "Ugdh", "Fn1","Has1")
# convert symbols to entrezID:
entrezid = select(mouseannot,
                  keys=genesymbols,
                  columns = c("ENTREZID","SYMBOL"),
                  keytype= "SYMBOL")
de <- entrezid$ENTREZID
names(geneList) <- de
# initial testing with Biological Process "BP":
ego <- enrichGO(de, OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
# ego is an S4 object
hist(ego@result[["p.adjust"]])  #simple view of FDR values in results 
selected <- subset(ego, subset=ego@result[["p.adjust"]]<=0.03)
dim(selected) 
topcat = dim(selected)[1] #rows number 
barplot(ego,showCategory = topcat)

ego2 <- enrichGO(de, OrgDb = "org.Mm.eg.db", ont="all", readable=TRUE)

dotplot(ego2, split="ONTOLOGY") + facet_grid(ONTOLOGY~.,scale="free")
  cnetplot(ego2, foldChange=geneList, circular=TRUE, colorEdge = TRUE)
# === try to create  object including clusters id =====

# formula interface
mydf <- data.frame(Entrez=de, group= c(rep("0",3),rep("1",3)))
xx.formula <- compareCluster(Entrez~group,data=mydf,
                             fun="enrichGO", ont="all", OrgDb="org.Mm.eg.db")

dotplot(xx.formula, split="ONTOLOGY") + facet_grid(ONTOLOGY~.,scale="free")

## testing with bigger list repeating 3 first genes

de2 = c(de,de[1:3]) 
mydf <- data.frame(Entrez=de2, group= c(rep("0",3),rep("1",3), rep("fake2",3)))
xx.formula2 <- compareCluster(Entrez~group,data=mydf,
                             fun="enrichGO", ont="all", OrgDb="org.Mm.eg.db")

dotplot(xx.formula2, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~.,scale="free") +
  theme(axis.text.y = element_text(size = 9)) +
  labs(x="cluster label (number of genes)")
