#Ligand-Receptor analysis
#TUTO: https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
#devtools::install_github("saeyslab/nichenetr")

# ** this version is setting pi16 as ligand and not as target !!!
library(nichenetr)
library(dplyr)
library(tidyverse)

library(Seurat)

prloc="~/reanalysesPubs/"
setwd(prloc)

# ,,,,,,,,,,,,,,,,,,,,
# ,,,,,,, DATA ,,,,,,,

# Oprescu&DeMicheli integrated expression data
# load Seurat  object:
seuratObj = readRDS("Opres_DeMichseu_withNames.rds")
seuratObj@meta.data %>% head()
#                       orig.ident nCount_RNA nFeature_RNA percent.mt integrated_snn_res.0.4 seurat_clusters
# D0_A_AAACCTGAGCTATGCT  DeMicheli       3764         1370   5.340064                      9               9
# D0_A_AAACCTGAGGCCCGTT  DeMicheli       1961         1164   4.232534                      2               2

#in my object, cell types are in @active.ident:
head(seuratObj@active.ident)
seuratObj@active.ident %>% table()
# FAPS, Smoc2+          FAPs activated       Endothelial cells                 unknown   Endothelial_lymphatic            FAPs, Smoc2+ 
#   2408                    1746                    1342                     906                     809                     801 
# Myonuclei       B_lymphocytes.CD4           Smooth muscle                    MuSc               Tenocytes                 Early T 
# 747                     292                     271                     253                     239                     174 
# FAPs?, Smoc2++           Schwann cells Neutrophyls_macrophages                 ?FAPs?? 
#   173                     191                     118                      68
DimPlot(seuratObj, reduction="umap", label = T)

length(seuratObj@meta.data$seurat_clusters)
#[1] 10538
length(seuratObj@active.ident)
#[1] 10538

head(rownames(seuratObj@meta.data))
head(names(seuratObj@active.ident))

if (all( rownames(seuratObj@meta.data) == names(seuratObj@active.ident) )){
  print("copying active.ident into celltype metadata")
  seuratObj@meta.data$celltypes <- as.character(seuratObj@active.ident)
}else{ print("problem!!different order in metadata") }
#define "aggregate" in the metadata:

FAPNONpi16 <- c("FAPS, Smoc2+", "FAPs, Smoc2+", "FAPs?, Smoc2++")
FAPpi16 <- c("FAPs activated")
seuratObj@meta.data$aggregate <- "nonfaps"
# 2. define gene set of interest : FAPpi16 vs FAPNONpi16
seuratObj@meta.data <- seuratObj@meta.data %>% 
  mutate(
    aggregate = case_when(
      celltypes %in% FAPNONpi16 ~ "FAPNONpi16",  
      celltypes %in% FAPpi16 ~ "FAPpi16", 
      TRUE ~ "nonfaps"
    )
  )
# problem with my mutate step: I deleted rows !! OMG the cells barcodes
seuratObj@meta.data %>% head()
# orig.ident nCount_RNA nFeature_RNA percent.mt integrated_snn_res.0.4 seurat_clusters         celltypes aggregate
# 1  DeMicheli       3764         1370   5.340064                      9               9              MuSc   nonfaps
# 2  DeMicheli       1961         1164   4.232534                      2               2 Endothelial cells   nonfaps
# RE-integrate rownames :
rownames(seuratObj@meta.data) <- names(seuratObj@active.ident) 

# %%% DOWNLOAD MODEL %%%

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

# We are working in mouse, it is necessary to convert NicheNet's model
# from human to mouse symbols: one-to-one orthologs were gathered from NCBI HomoloGene and also from ENSEMBL 
# also explained at https://github.com/saeyslab/nichenetr/blob/master/vignettes/symbol_conversion.md

lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

# ,,,,,,,,,,,,,,,,,,,,,,,,
# ,,,,,,, ANALYSIS ,,,,,,,

## receiver = target, cell type expressing receptors : all except senders, see below
receivers = allCells <- c("FAPS, Smoc2+","Endothelial cells","unknown",
                          "Endothelial_lymphatic","FAPs, Smoc2+","Myonuclei","B_lymphocytes.CD4",
                          "Smooth muscle","MuSc","Tenocytes","Early T","FAPs?, Smoc2++",
                          "Schwann cells","Neutrophyls_macrophages","?FAPs??","Schwann cells")  
expressed_genes_receiver = get_expressed_genes(receivers[1], seuratObj, pct = 0.10)
for (k in 2:length(receivers)){
  tmpvec = get_expressed_genes(receivers[k], seuratObj, pct = 0.10)
  expressed_genes_receiver = c(expressed_genes_receiver, tmpvec )
}
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender  : sender cells express ligands
sender_celltypes = c("FAPs activated")  
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()


# receiver : define differentially expr in two conditions, 
# but here I have only one : day zero !
seurat_obj_receiver= subset(seuratObj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["aggregate"]])

# TODO: when two days available, uncomment and modify:
# condition_oi = "day0"
# condition_reference = "day2"
# 
# DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")
# 
# geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= 0.25) %>% pull(gene)
# geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

# TODO: comment this part when days available to compare (above)
savedmarkers <- read.table("Opres_DeMichFAPSonly.100markers.17.txt", sep="\t", header =T)

geneset_oi1 <- savedmarkers$gene[savedmarkers$cluster==1 & savedmarkers$avg_logFC > 0.7]
geneset_oi <- geneset_oi1 %>% .[. %in% rownames(ligand_target_matrix)] 

# 3. SET OF POTENTIAL LIGANDS:
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

# 4. NicheNet ligand activity analysis
ligand_activities = predict_ligand_activities(geneset = as.character(geneset_oi), background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities

best_upstream_ligands = ligand_activities %>% top_n(40, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
DotPlot(seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
# We see here that top-rankd ligands can predict reaonably? 

# 5. infer receptors:

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 300) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                 ligand_target_matrix = ligand_target_matrix, cutoff = 0.01) # changed from 0.33
#Quantile cutoff on the ligand-target scores of the input weighted ligand-target network

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.006,0.012))
p_ligand_target_network

#receptors of top ranked ligands
  
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()  
  
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
pdf("ligrec_pi16aslig_heatmap.pdf", width=15)
p_ligand_receptor_network  
dev.off()
  

# ---------END-------------


# the aim is : a) to assign ligands to one of the posible sender cell types,
# being the sender cell the one that most strongly expresses this ligand (expression higher than average),
# or b) to assign ligands to the group of 'generally' expressed ligands,
# it is, ligands having expression less or equal the average expression 
# *I am just rephrasing help from: https://github.com/saeyslab/nichenetr/issues/5



