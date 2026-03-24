message("\n\n##########################################################################\n",
        "# Start B02a replot from saved seur.qs ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n")

library(qs)
library(tidyverse)
library(Seurat)
library(colorRamps)
library(viridis)
library(ComplexHeatmap)
library(patchwork)

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

script_ind = "LD_B02a_"
out_dir = paste0(main_dir,"LD_B_AST_analysis_output/")

gr_tab = read_csv("data_TREM2_michael/A_input/group_tab.csv")
seur = qread(file = paste0(out_dir, "LD_B02a_seur.qs"))

#get marker gene panels 
GOI = list() # creates an empty list
t1 = read_csv("data_TREM2_michael/A_input/cell_type_markers_241219_w_astr_subtype_markers.csv")
GOI$cell_type_markers = t1$gene[t1$level %in% c("cell_types", "neuronal_lineage")] # stores genes markingbroad cell types 
GOI$subtype_markers = t1$gene[t1$level %in% c("Astrocyte_subtypes")] # stores genes marking astrocyte cell types 

#replace APOE metadata column (else complications with plotting APOE feature vs metadata col)
names(seur@meta.data)[names(seur@meta.data) == "APOE"] = "APOE_genotype"


####################################
#Functions
####################################

#custom colour palette for variable values defined in vector v
pal = function(v){
  v2 = length(unique(v))
  if (v2 == 2){
    p2 = c("grey20", "dodgerblue")
  } else if (v2 ==3){
    p2 = c("dodgerblue", "grey20", "orange")
  } else if (v2 ==4){
    p2 = c("dodgerblue", "green4","grey20", "orange")
  } else if (v2<6){
    p2 = matlab.like(6)[1:v2]
  } else {
    p2 = matlab.like(v2)
  }
  return(p2)
}


#distinct scale (larger number of colours, colour vector gets shuffled)
pal_dist = function(v){
  v2 = length(unique(v))
  if (v2 < 6){p2 = matlab.like(6)[1:v2]} else {
    p2 = matlab.like(v2)
    set.seed(12)
    p2 = sample(p2)
  }
  return(p2)
}



#############################################
#Plot UMAP and marker expression (test different clustering resolutions)
#############################################

#seur = qread(file = paste0(out_dir, script_ind, "seur.qs"))

set.seed(42)

#define resolutions for clustering tests
test_res = c(0.1, 0.15, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5)


#subsample seurat object for plotting for large datasets (else plots become too large (pdf with hundreds of MB))

seur0 = seur

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)] # randomly subsample 100,000 cells for plotting if dataset is very large
} 


#create lists for plot objects
pl_umap = list()
pl_dotplots_cell_type_markers = list()
pl_dotplots_subtype_markers = list()
pl_umap_expr = list()


### generate umap plot by group and sample


#define grouping variables and colors for plotting
gr = unique(gr_tab$group)
samples = unique(gr_tab$sample)


pl_umap[["group"]] = DimPlot(seur, group.by = "group", shuffle = TRUE,
                             label = FALSE, reduction = "umap", 
                             pt.size = 0.01)+
  scale_color_manual(limits = gr, values = pal(gr))+
  labs(title = "UMAP by group")

pl_umap[["sample"]] = DimPlot(seur, group.by = "sample", shuffle = TRUE,
                              label = FALSE, reduction = "umap", 
                              pt.size = 0.01)+
  scale_color_manual(limits = samples, values = pal(samples))+
  labs(title = "UMAP by sample")




### create plots for clusters and marker dot plot for different resolutions (each test_res)
j = test_res[1]
for (j in test_res){
  
  set.seed(1234)
  
  seur0 <- FindClusters(
    object = seur0,
    graph.name = "SCT_snn",
    algorithm = 1,
    resolution = j,
    verbose = FALSE
  )
  
  #subset for plotting
  
  seur = seur0
  
  cells = colnames(seur@assays$SCT@scale.data)
  
  if (length(cells)>100000){
    set.seed(1234)
    seur = seur[, sample(cells, size =100000, replace=F)]
  } 
  
  cl = unique(seur$seurat_clusters)
  pl_umap[[paste0("resol_", j)]] = DimPlot(seur, group.by = "seurat_clusters", shuffle = TRUE,
                                           label = TRUE, reduction = "umap", pt.size = 0.01)+
    scale_color_manual(limits = cl, values = pal(cl))+
    NoLegend()+labs(title = paste0("clusters resol_", j))
  
  
  #create dotplot with gene expression of marker genes
  
  DefaultAssay(seur0) = "SCT"
  
  pl_dotplots_cell_type_markers[[paste0("resol_", j)]] = 
    DotPlot(seur0, features = intersect(GOI$cell_type_markers, rownames(seur0)), 
            scale.by = "size") +RotatedAxis()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste0("resol_", j))
  
  pl_dotplots_subtype_markers[[paste0("resol_", j)]] = 
    DotPlot(seur0, features = intersect(GOI$subtype_markers, rownames(seur0)), 
            scale.by = "size") +RotatedAxis()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste0("resol_", j))
  
}


seur = seur0

#save dataset with cluster assignment for different cluster resolutions
qsave(seur, file = paste0(out_dir, script_ind, "seur.qs"))




### save umap clustering and dotplots

pdf(file = paste0(out_dir, script_ind, "UMAP_clustering_test_integr_dataset.pdf"), width = 6, height = 5)
lapply(pl_umap, function(x){x})
dev.off()

pdf(file = paste0(out_dir, script_ind, "Dotplots_clustering_test_integr_dataset_cell_type_markers.pdf"), width = 12, height = 10)
lapply(pl_dotplots_cell_type_markers, function(x){x})
dev.off()

pdf(file = paste0(out_dir, script_ind, "Dotplots_clustering_test_integr_dataset_subtype_markers.pdf"), width = 25, height = 10)
lapply(pl_dotplots_subtype_markers, function(x){x})
dev.off()



############################################################
# UMAP plots split by group vs repl/batch_seq
############################################################

seur = seur0

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


t1 = FetchData(seur, vars = c("umap_1", "umap_2", "group", "sample"))


#plot group vs sample

p1 = ggplot(data = t1, aes(x = umap_1, y = umap_2))+geom_point(size = 0.05, alpha = 0.5)+
  ggplot2::facet_wrap(facets = factor(t1$sample, levels = samples), nrow = 2)+
  theme_bw()+
  labs(title = "UMAP by group vs sample")+RestoreLegend()

pdf(file = paste0(out_dir, script_ind, "UMAP_split_by_group_vs_sample.pdf"), width = 10, height = 2.5)
plot(p1)
dev.off()


######################################################
# plots for marker expression on UMAP
######################################################

#subsample seurat object for large datasets (else plots become too large (pdf with hundreds of MB))

seur = seur0

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>10000){
  set.seed(1234)
  seur = seur[, sample(cells, size =10000, replace=F)]
} 

DefaultAssay(seur) = "SCT"

pl_umap_expr = FeaturePlot(seur, features = c(GOI$cell_type_markers, GOI$subtype_markers), order = TRUE,
                           reduction = "umap", pt.size = 0.01, ncol = 9)

pdf(file = paste0(out_dir, script_ind, "marker_expr_UMAP_integr_dataset.pdf"), width = 30, height = 15)
plot(pl_umap_expr)
dev.off()


###################################################################################
# save table for cluster assignment/labelling (with highest resolution) 
###################################################################################
# for downstream analyses assign cell types in csv; remove not needed clusters for lower resolutions

seur = seur0

v1 = unique(seur$seurat_clusters)
v2 = v1[order(v1)]
t1 = tibble(cluster = v2, cluster_name = paste0("c",v2), cell_type = "all", cell_class = "all")

if (file.exists(paste0(out_dir, script_ind, "cluster_assignment.csv"))){
  write_csv(t1, paste0(out_dir, script_ind, "cluster_assignment_new.csv"))
} else {write_csv(t1, paste0(out_dir, script_ind, "cluster_assignment.csv"))}



#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed B02 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")

