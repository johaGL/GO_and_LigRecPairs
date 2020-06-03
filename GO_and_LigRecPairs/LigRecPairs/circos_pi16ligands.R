# ,,,,,,,,,,,,,,,,,,,,
# ,,,,,, CIRCOS ,,,,,,
# THIS IS VERSION SHOWS PI16 LIGANDS
# and predicted interactions with receptors(=target) expressed by the other cell types 

library(nichenetr)
library(dplyr)
library(tidyverse)
library(circlize)

library(Seurat)

prloc="~/reanalysesPubs/"
setwd(prloc)

# ** ATTENTION if objects already loaded, not to run this part **

# Oprescu&DeMicheli integrated expression data
# load Seurat  object:
seuratObj = readRDS("Opres_DeMichseu_withNames.rds")
seuratObj@meta.data %>% head()
#in my object, cell types are in @active.ident:
head(seuratObj@active.ident)
seuratObj@active.ident %>% table()

length(seuratObj@meta.data$seurat_clusters)
#[1] 10538
length(seuratObj@active.ident)
#[1] 10538
head(rownames(seuratObj@meta.data))
head(names(seuratObj@active.ident))
#check if these entire vectors are the same (the barcodes):
if (all( rownames(seuratObj@meta.data) == names(seuratObj@active.ident) )){
  print("copying active.ident into celltype metadata")
  seuratObj@meta.data$celltypes <- as.character(seuratObj@active.ident)
}else{ print("problem!!different order in metadata") }

rownames(seuratObj@meta.data) <- names(seuratObj@active.ident) # restore rownames

# Load the ligand-target model
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

# mouse transf:
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

  # not used ?:
  #weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

# ** end ATTENTION if objects already loaded, not to run this part **


## start here if everything ok
ALLMARKERS <- read.table("ALLMARKERSoprescuDmicheli.txt", sep="\t", header=T)
ALLMARKERS$celltypes = seuratObj@meta.data$celltypes[match(ALLMARKERS$cluster,seuratObj@meta.data$seurat_clusters)]


FAPSmoc2.main = "FAPS, Smoc2+" # cluster 0
FAPSmoc2.others = c("FAPs, Smoc2+",  "FAPs?, Smoc2++")
FAPpi16 = "FAPs activated" # cluster 1
tenocytes = "Tenocytes" # cluster 10
endothelial = c("Endothelial cells", "Endothelial_lymphatic")
neumacr = "Neutrophyls_macrophages"
smooth = "Smooth muscle"
schwann = "Schwann cells"
musc = "MuSc"


allCells <- c("FAPS, Smoc2+","FAPs activated","Endothelial cells","unknown",
               "Endothelial_lymphatic","FAPs, Smoc2+","Myonuclei","B_lymphocytes.CD4",
               "Smooth muscle","MuSc","Tenocytes","Early T","FAPs?, Smoc2++",
               "Schwann cells","Neutrophyls_macrophages","?FAPs??","Schwann cells")  # background :all except target pi16

# expressed_genes_FAPSmoc2.main <- ALLMARKERS %>% 
#   filter(celltypes %in% FAPSmoc2.main & avg_logFC > 0.3) %>% pull(gene)
# # is a factor vector, convert into char vector:
# expressed_genes_FAPSmoc2.main <- as.character(expressed_genes_FAPSmoc2.main)
# 
# expressed_genes_FAPSmoc2.others <- ALLMARKERS %>% 
#   filter(celltypes %in% FAPSmoc2.others & avg_logFC > 0.3) %>% pull(gene)  
# expressed_genes_FAPSmoc2.others <- as.character(expressed_genes_FAPSmoc2.others)

expressed_genes_FAPpi16 <- ALLMARKERS %>% 
  filter(celltypes %in% FAPpi16 & avg_logFC > 0.3) %>% pull(gene)
expressed_genes_FAPpi16 <- as.character(expressed_genes_FAPpi16)

expressed_genes_tenocytes <- ALLMARKERS %>%
  filter(celltypes %in% tenocytes & avg_logFC > 0.3) %>% pull(gene)
expressed_genes_tenocytes <- as.character(expressed_genes_tenocytes)
# 
# expressed_genes_endothelial <- ALLMARKERS %>% 
#   filter(celltypes %in% endothelial & avg_logFC > 0.3) %>% pull(gene)
# expressed_genes_endothelial <- as.character(expressed_genes_endothelial)
# 
# expressed_genes_neumacr <- ALLMARKERS %>% 
#   filter(celltypes %in% neumacr & avg_logFC > 0.3) %>% pull(gene)  
# expressed_genes_neumacr <- as.character(expressed_genes_neumacr)
# 
# expressed_genes_smooth <- ALLMARKERS %>% 
#   filter(celltypes %in% smooth & avg_logFC > 0.3) %>% pull(gene)  
# expressed_genes_smooth <- as.character(expressed_genes_smooth)
# 
# expressed_genes_schwann <- ALLMARKERS %>% 
#   filter(celltypes %in% schwann & avg_logFC > 0.3) %>% pull(gene)  
# expressed_genes_schwann <- as.character(expressed_genes_schwann)
# 
# expressed_genes_musc <- ALLMARKERS %>% 
#   filter(celltypes %in% musc & avg_logFC > 0.3) %>% pull(gene)  
# expressed_genes_musc <- as.character(expressed_genes_musc)

expressed_genes_allCells <- ALLMARKERS %>% 
  filter(celltypes %in% allCells & avg_logFC > 0.7) %>% pull(gene)
expressed_genes_allCells <- as.character(expressed_genes_allCells)


# Load geneS setS of interest i.e. receptors at each celltype; and set background cells Other
# my set of interest is pi16 highest expressed:
FAPpi16_geneset = ALLMARKERS %>% 
  filter(celltypes %in% FAPpi16 & avg_logFC > 0.7) %>% pull(gene) %>%
   .[. %in% rownames(ligand_target_matrix)]
FAPpi16_geneset <- as.character(FAPpi16_geneset)

FAPSmoc2.main_geneset <- ALLMARKERS %>% 
  filter(celltypes %in% FAPSmoc2.main & avg_logFC > 0.7) %>% pull(gene)%>%
  .[. %in% rownames(ligand_target_matrix)]
FAPSmoc2.main_geneset <- as.character(FAPSmoc2.main_geneset)

FAPSmoc2.others_geneset <- ALLMARKERS %>% 
  filter(celltypes %in% FAPSmoc2.others & avg_logFC > 0.7) %>% pull(gene)  %>%
  .[. %in% rownames(ligand_target_matrix)]
FAPSmoc2.others_geneset <- as.character(FAPSmoc2.others_geneset)

tenocytes_geneset <- ALLMARKERS %>% 
  filter(celltypes %in% tenocytes & avg_logFC > 0.7) %>% pull(gene)%>%
  .[. %in% rownames(ligand_target_matrix)]
tenocytes_geneset <- as.character(tenocytes_geneset)

endothelial_geneset <- ALLMARKERS %>% 
  filter(celltypes %in% endothelial & avg_logFC > 0.7) %>% pull(gene) %>%
  .[. %in% rownames(ligand_target_matrix)]
endothelial_geneset <- as.character(endothelial_geneset)

neumacr_geneset <- ALLMARKERS %>% 
  filter(celltypes %in% neumacr & avg_logFC > 0.7) %>% pull(gene) %>%
  .[. %in% rownames(ligand_target_matrix)]
neumacr_geneset <- as.character(neumacr_geneset)

smooth_geneset <- ALLMARKERS %>% 
  filter(celltypes %in% smooth & avg_logFC > 0.7) %>% pull(gene) %>%
  .[. %in% rownames(ligand_target_matrix)]
smooth_geneset <- as.character(smooth_geneset)

schwann_geneset <- ALLMARKERS %>% 
  filter(celltypes %in% schwann & avg_logFC > 0.7) %>% pull(gene) %>%
  .[. %in% rownames(ligand_target_matrix)] 
schwann_geneset <- as.character(schwann_geneset)

musc_geneset <- ALLMARKERS %>% 
  filter(celltypes %in% musc & avg_logFC > 0.7) %>% pull(gene) %>%
  .[. %in% rownames(ligand_target_matrix)]
musc_geneset <- as.character(musc_geneset)

uni_geneset <- union(FAPpi16_geneset,union(FAPSmoc2.main_geneset,
                                         union(FAPSmoc2.others_geneset,union(tenocytes_geneset,
                                         union(endothelial_geneset,union(neumacr_geneset,
                                          union(smooth_geneset, union(schwann_geneset, musc_geneset)))))))
                   )

background_expressed_genes = expressed_genes_allCells %>% .[. %in% rownames(ligand_target_matrix)]


# perform ligand activity analysis on set of interest

ligands = lr_network %>% pull(from) %>% unique()

# expressed_ligands_FAPSmoc2.main = intersect(ligands,expressed_genes_FAPSmoc2.main)
# expressed_ligands_FAPSmoc2.others = intersect(ligands,expressed_genes_FAPSmoc2.others)
expressed_ligands_FAPpi16 = intersect(ligands,expressed_genes_FAPpi16)
expressed_ligands_tenocytes = intersect(ligands,expressed_genes_tenocytes)
# expressed_ligands_endothelial = intersect(ligands,expressed_genes_endothelial)
# expressed_ligands_neumacr = intersect(ligands,expressed_genes_neumacr)
# expressed_ligands_smooth = intersect(ligands,expressed_genes_smooth)
# expressed_ligands_schwann = intersect(ligands,expressed_genes_schwann)
# expressed_ligands_musc = intersect(ligands,expressed_genes_musc)
# expressed_ligands = union(expressed_ligands_FAPSmoc2.main,
#                           union(expressed_ligands_FAPSmoc2.others,
#                           union(expressed_ligands_FAPpi16,
#                           union(expressed_ligands_tenocytes,
#                           union(expressed_ligands_endothelial,
#                           union(expressed_ligands_neumacr,
#                           union(expressed_ligands_smooth,
#                           union(expressed_ligands_schwann,
#                           expressed_ligands_musc)))))))
#                           )
expressed_ligands = union(expressed_ligands_FAPpi16, expressed_ligands_tenocytes)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,background_expressed_genes)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & 
                                            to %in% expressed_receptors) %>%
  pull(from) %>% unique()

# if you want ot make circos plots for multiple gene sets, 
#combine the different data frames and differentiate which target 
#belongs to which gene set via the target type

ligand_activities = predict_ligand_activities(geneset = uni_geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

best_upstream_ligands = ligand_activities %>% top_n(50, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

#check overlaps and assign those "common" specifically:
# this part is not developped yet, general idea:
# lfc.common_i.group.a > lfc.common_i.group.b + 1
# and if not satisfying this, dropout the common gene
# preparation, likely:
# pi16_lig_df <- ALLMARKERS %>% filter(celltypes %in% FAPpi16 & avg_logFC>0) %>%
#  +   filter(gene %in% best_upstream_ligands)

# for now, continue as ligands were specific:
# FAPSmoc2.main_specific_ligands <- expressed_ligands_FAPSmoc2.main %>% intersect(best_upstream_ligands)
# FAPSmoc2.others_specific_ligands <- expressed_ligands_FAPSmoc2.others %>% intersect(best_upstream_ligands)
FAPpi16_specific_ligands <- expressed_ligands_FAPpi16 %>% intersect(best_upstream_ligands)
tenocytes_specific_ligands <- expressed_ligands_tenocytes %>% intersect(best_upstream_ligands)
# endothelial_specific_ligands <- expressed_ligands_endothelial %>% intersect(best_upstream_ligands)
# neumacr_specific_ligands <- expressed_ligands_neumacr %>% intersect(best_upstream_ligands)
# smooth_specific_ligands <- expressed_ligands_smooth %>% intersect(best_upstream_ligands)
# schwann_specific_ligands <- expressed_ligands_schwann %>% intersect(best_upstream_ligands)
# musc_specific_ligands <- expressed_ligands_musc  %>% intersect(best_upstream_ligands)


general_ligands = expressed_ligands %>% intersect(best_upstream_ligands)
  
# ligand_type_indication_df = tibble(
#   ligand_type = c(rep("FAPSmoc2.main", times = length(FAPSmoc2.main_specific_ligands)),
#                   rep("FAPSmoc2.others", times = length(FAPSmoc2.others_specific_ligands)),
#                   rep("FAPpi16", times = length(FAPpi16_specific_ligands)),
#                   rep("Tenocytes", times = length(tenocytes_specific_ligands)),
#                   rep("Endothelial", times = length(endothelial_specific_ligands)),
#                   rep("NeutrophMacroph", times = length(neumacr_specific_ligands)),
#                   rep("SmoothMus", times = length(smooth_specific_ligands)),
#                   rep("SchwannCells",times = length(schwann_specific_ligands)),
#                   rep("MuSC", times = length(musc_specific_ligands))
#                   #,rep("general", times= length(general_ligands))
#                   ),
#   ligand = c(FAPSmoc2.main_specific_ligands, FAPSmoc2.others_specific_ligands, FAPpi16_specific_ligands,
#              tenocytes_specific_ligands, endothelial_specific_ligands, 
#              neumacr_specific_ligands,smooth_specific_ligands,
#              schwann_specific_ligands, musc_specific_ligands ) #, general_ligands
#   )#end tibble

ligand_type_indication_df = tibble(
  ligand_type = c(    rep("FAPpi16", times = length(FAPpi16_specific_ligands))
                  #rep("Tenocytes", times = length(tenocytes_specific_ligands))
                  ),
  ligand = c( FAPpi16_specific_ligands ) #,tenocytes_specific_ligands
  )#end tibble


#**** infer target genes of top-ranked ligands and visualize in a circos plot
ligtarLIST = list()
ligtarLIST[[1]] = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = FAPpi16_geneset, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows()
ligtarLIST[[2]] = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = FAPSmoc2.main_geneset, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows()
ligtarLIST[[3]] = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = FAPSmoc2.others_geneset, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows()
ligtarLIST[[4]] = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = tenocytes_geneset, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows()
ligtarLIST[[5]] = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = endothelial_geneset, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows()
ligtarLIST[[6]] = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = neumacr_geneset, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows()
ligtarLIST[[7]] = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = smooth_geneset, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows()
ligtarLIST[[8]] = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = schwann_geneset, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows()
ligtarLIST[[9]] = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = musc_geneset, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows()

targets = c("FAPpi16","FAPSmoc2main", "FAPSmoc2others","Tenocytes","Endothelial","NeutrophMacroph","SmoothMus","SchwannCells","MuSC")

# if you want ot make circos plots for multiple gene sets, 
#combine the different data frames and differentiate which target 
#belongs to which gene set via the target type
active_ligand_target_links_df = ligtarLIST[[1]] %>% mutate(target_type = "FAPpi16") %>% inner_join(ligand_type_indication_df) 
for (i in 2:length(targets)){
  tmpdf = ligtarLIST[[i]] %>% mutate(target_type = targets[i]) %>% inner_join(ligand_type_indication_df)
  print(is.data.frame(tmpdf))
  active_ligand_target_links_df = rbind(active_ligand_target_links_df, tmpdf)
}  # other loong ugly way to do this: see at the END

# show only links with weight higher than cutoff
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.66, na.rm = TRUE)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

grid_col_target =c("FAPSmoc2.main" = "darkorange",
                   "FAPSmoc2.others" = "gold",
                   "FAPpi16" = "orangered",
                   "Tenocytes" = "deeppink4",
                   "Endothelial" = "steelblue2",
                   "NeutrophMacroph" = "wheat2",
                  "SmoothMus" = "salmon", 
                  "SchwannCells" = "turquoise4" ,
                  "MuSC" =  "purple"  )  # other colors: "saddlebrown"  ,  
grid_col_ligand =c("FAPpi16" = "lightblue")

# HOWEVER: see below that order of appeareance was changed , neutrophyls are moved last
# *-*

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

target_order = circos_links$target %>% unique()
ligand_order = c(FAPpi16_specific_ligands ) %>% 
            c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)


#svg saving commands do not work, better do manually from plots
#svg("cellsPi16receive_EndothSmocTeno.svg",width=10, height = 10)
circos.par(gap.degree = 1)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.7)
}, bg.border = NA) 

#dev.off()

# legend manually added, alphabetic order:
pdf("LigRecPairs/circosExtendedLegend_ligp16.pdf",  width=15)
  ggplot(active_ligand_target_links_df_circos, aes(x=target_type, y=weight, fill=target_type))+ geom_violin()+
    scale_fill_manual(values= c("steelblue2", "orangered","darkorange","gold", 
                                "purple", "wheat2","turquoise4" , "salmon", "deeppink4"))+
    ggtitle("Ligand-receptor weights as calculated by NicheNet")
dev.off()

#circos without transparency
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) 



circos.clear()


# END

#looong ugly ways :
# predf.1 = ligtarLIST[[1]] %>% mutate(target_type = "FAPpi16") %>% inner_join(ligand_type_indication_df) 
# predf.2 = ligtarLIST[[2]] %>% mutate(target_type = "FAPSmoc2main") %>% inner_join(ligand_type_indication_df) 
# predf.3 = ligtarLIST[[3]] %>% mutate(target_type = "FAPSmoc2others") %>% inner_join(ligand_type_indication_df) 
# predf.4 = ligtarLIST[[... %>% mutate(target_type = "Tenocytes") %>% inner_join(ligand_type_indication_df) 
# predf.5 = ligtarLIST[[... %>% mutate(target_type = "Endothelial") %>% inner_join(ligand_type_indication_df) 
# predf.6 = ligtarLIST[[... %>% mutate(target_type = "NeutrophMacroph") %>% inner_join(ligand_type_indication_df) 
# predf.7 = ligtarLIST[[... %>% mutate(target_type = "SmoothMus") %>% inner_join(ligand_type_indication_df) 
# predf.8 = ligtarLIST[[... %>% mutate(target_type = "SchwannCells") %>% inner_join(ligand_type_indication_df) 
# predf.9 = ligtarLIST[[... %>% mutate(target_type = "MuSC") %>% inner_join(ligand_type_indication_df) 


