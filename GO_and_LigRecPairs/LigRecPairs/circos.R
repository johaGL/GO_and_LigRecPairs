# ,,,,,,,,,,,,,,,,,,,,
# ,,,,,, CIRCOS ,,,,,,

library(nichenetr)
library(dplyr)
library(tidyverse)
library(circlize)

library(Seurat)

prloc="~/reanalysesPubs/"
setwd(prloc)


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
#check if these entire vectors are the same (the barcodes)
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

FAPpi16 = "FAPs activated" # cluster 1
FAPSmoc2.main = "FAPS, Smoc2+" # cluster 0
#FAPSmoc2.others: "FAPs, Smoc2+" "FAPs?, Smoc2++"
# "Tenocytes"  # test also
otherCells <- c("Endothelial cells","unknown",
               "Endothelial_lymphatic","Myonuclei","B_lymphocytes.CD4",
               "Smooth muscle","MuSc","Early T",
               "Schwann cells","Neutrophyls_macrophages","?FAPs??","Schwann cells") 

expressed_genes_FAPpi16 <- ALLMARKERS %>% 
  filter(celltypes %in% FAPpi16 & avg_logFC > 0) %>% pull(gene)
# is a factor vector, convert into char vector:
expressed_genes_FAPpi16 <- as.character(expressed_genes_FAPpi16)

expressed_genes_FAPSmoc2.main <- ALLMARKERS %>% 
  filter(celltypes %in% FAPSmoc2.main & avg_logFC > 0) %>% pull(gene)
# is a factor vector, convert into char vector:
expressed_genes_FAPSmoc2.main <- as.character(expressed_genes_FAPSmoc2.main)

expressed_genes_otherCells <- ALLMARKERS %>% 
  filter(celltypes %in% otherCells & avg_logFC > 0) %>% pull(gene)
# is a factor vector, convert into char vector:
expressed_genes_otherCells <- as.character(expressed_genes_otherCells)


# Load gene set of interest and set background cells other than faps
# my set of interest is pi16 highest expressed:
i_geneset = ALLMARKERS %>% 
  filter(celltypes %in% FAPpi16 & avg_logFC > 1) %>% pull(gene) %>%
   .[. %in% rownames(ligand_target_matrix)]
# is a factor vector, convert into char vector:
i_geneset <- as.character(i_geneset)

background_expressed_genes = expressed_genes_otherCells %>% .[. %in% rownames(ligand_target_matrix)]


# perform ligand activity analysis on set of interest

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands_FAPpi16 = intersect(ligands,expressed_genes_FAPpi16)
expressed_ligands_FAPSmoc2.main = intersect(ligands,expressed_genes_FAPSmoc2.main)
expressed_ligands = union(expressed_genes_FAPpi16, expressed_genes_FAPSmoc2.main)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_otherCells)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & 
                                            to %in% expressed_receptors) %>%
  pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = i_geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

best_upstream_ligands = ligand_activities %>% top_n(30, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

#check overlaps:
a <- best_upstream_ligands %>% intersect(expressed_ligands_FAPpi16)
b <- best_upstream_ligands %>% intersect(expressed_ligands_FAPSmoc2.main)
intersect(a,b) # "Bmp2"    "Rarres2" "Cxcl1"  # good news , little overlapping
common <- intersect(a,b)

# assign those "common" specifically:
# this part is not developped yet, general idea:
# lfc.common_i.group.a > lfc.common_i.group.b + 1
# and if not satisfying this, dropout the common gene
# preparation, likely:
# pi16_lig_df <- ALLMARKERS %>% filter(celltypes %in% FAPpi16 & avg_logFC>0) %>%
#  +   filter(gene %in% best_upstream_ligands)

# for now, continue as ligands were specific:
FAPpi16_specific_ligands <- expressed_ligands_FAPpi16 %>% intersect(best_upstream_ligands)
FAPSmoc2.main_specific_ligands <- expressed_ligands_FAPSmoc2.main %>% intersect(best_upstream_ligands)
general_ligands = setdiff(FAPpi16_specific_ligands,FAPSmoc2.main_specific_ligands)
  
ligand_type_indication_df = tibble(
  ligand_type = c(rep("FAPpi16", times = length(FAPpi16_specific_ligands)),
                  rep("General", times = length(general_ligands)),
                  rep("FAPSmoc2.main", times = length(FAPSmoc2.main_specific_ligands))
                  ),
  ligand = c(FAPpi16_specific_ligands, general_ligands, FAPSmoc2.main_specific_ligands)
)

#**** infer target genes of top-ranked ligands and visualize in a circos plot

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = i_geneset, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows()

active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "i_fapPi16") %>% inner_join(ligand_type_indication_df) 
# if you want ot make circos plots for multiple gene sets, 
#combine the different data frames and differentiate which target 
#belongs to which gene set via the target type

# show only links with weight higher than cutoff
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.15, na.rm = TRUE)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

grid_col_ligand =c("General" = "lawngreen",
                   "FAPpi16" = "royalblue",
                   "FAPSmoc2.main" = "gold")
grid_col_target =c("i_fapPi16" = "tomato")

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
ligand_order = c(FAPpi16_specific_ligands,general_ligands,FAPSmoc2.main_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

# DEFINING GAPS FOR FIGURE: DID NOT WORK
# width_same_cell_same_ligand_type = 0.1
# width_different_cell = 1
# width_ligand_receptor = 1
# width_same_cell_same_receptor_type = 0.1
# 
# gaps = c(
#   # width_ligand_receptor,
#   rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "FAPpi16") %>% distinct(ligand) %>% nrow() -1)),
#   width_different_cell,
#   rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
#   width_different_cell,
#   rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type =="FAPSmoc2.main" ) %>% distinct(ligand) %>% nrow() -1)), 
#   width_ligand_receptor,
#   rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(target_type == "i_fapPi16") %>% distinct(target) %>% nrow() -1)),
#   width_ligand_receptor
# )

# circos.par(gap.degree = gaps) NOT WORKING :(!!!!!
circos.par(gap.degree = 1) # OK
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
circos.clear()


# slight improvement setting widths 0.1  1  1  0.1
 
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
               preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #  

