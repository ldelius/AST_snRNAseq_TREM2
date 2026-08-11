message("\n\n##########################################################################\n",
        "# Start LD_B04a: Astrocyte subcluster characterisation ", Sys.time(),
        "\n##########################################################################\n\n")

library(qs)
library(tidyverse)
library(Signac) # Seurat extension for chromatin accessibility data (not called in this file)
library(Seurat)
library(EnsDb.Hsapiens.v86) # used for gene things, not called in this file
library(BSgenome.Hsapiens.UCSC.hg38) # used for gene things, not called in this file
library(colorRamps)
library(viridis)
library(lmerTest) # useful for stats, not called in this script
library(pheatmap)
library(sccomp) # stats engine for differential abundance testing
library(clusterProfiler) # GO enrichment
library(DOSE) # GO enrichment
library(org.Hs.eg.db) # GO enrichment
library(ggrepel)
library(presto)
# the not called packages might be dependencies from other packages, therefore i keep them for now

### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "LD_B04a_v02_" # this is version 02 now, cause I changed cluster s1 from SLC1A2 to GFAP.

#specify output directory
out_dir = paste0(main_dir,"LD_B_AST_analysis_output/")

### load group and file info for analysis dataset
gr_tab = read_csv("data_TREM2_michael/A_input/group_tab.csv") # metadata table with group and sample information for all samples in the dataset, used for plotting and stats

###load seurat dataset --> cleaned astrocytes
seur = qread(file = paste0(out_dir, "LD_B03a_seur.qs"))

#add plaque density (Amyloid plaques)
t1 = read_csv("data_TREM2_michael/A_input/TREM2_plaque_data_Sam.csv")
gr_tab$plaque_dens = t1$TotalDensity[match(gr_tab$BrainBankNetworkIDFormatted, t1$BrainBankNetworkIDFormatted)]
seur$plaque_dens = t1$TotalDensity[match(seur$BrainBankNetworkIDFormatted, t1$BrainBankNetworkIDFormatted)]

#add tau pathology (PHF1) from group table
seur$pctPHF1PositiveArea = gr_tab$pctPHF1PositiveArea[match(seur$sample, gr_tab$sample)]

#load subset dataset
clust_tab = read_csv(paste0(out_dir,"LD_B03a_cluster_assignment.csv")) # adds the CSV with manually assigned cluster names

#select subcluster resolution for further analyses
clust_res = 0.3 # adjust with the resolution I think fits best for further analyses

### get marker gene panels

GOI = list()
t1 = read_csv("data_TREM2_michael/A_input/cell_type_markers_241219_w_astr_subtype_markers.csv")
GOI$cell_type_markers = t1$gene[t1$level %in% c("cell_types", "neuronal_lineage")]
GOI$subtype_markers = t1$gene[t1$level %in% c("Astrocyte_subtypes")]

t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol

### get all subtype markers from Gazestani et al., 2023, Manucso et al., 2024 --> later used to copare our astrocyte subclusters

subtype_markers = list()

t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/Green24_S2_subpopulation_markers.csv"))

for (cl in unique(t1$state)){
  t2 = t1[t1$state == cl & t1$avg_log2FC>log2(1.2) & t1$p_val_adj<0.05,]
  subtype_markers$Green24[[cl]] = t2$gene
}

lengths(subtype_markers$Green24)


t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/Gazestani23_TableS2_MIC_markers.csv"))

for (cl in unique(t1$cluster_name)){
  t2 = t1[t1$cluster_name == cl & t1$Biopsy_avg_logFC>log2(1.2),]
  subtype_markers$Gazestani23[[cl]] = t2$gene_short_name
}

lengths(subtype_markers$Gazestani23)


###define group comparisons for each cluster (for sccomp)

comp_groups = list(Control_R47H_vs_CV = c("Control_R47H", "Control_CV"),
                   Control_R62H_vs_CV = c("Control_R62H", "Control_CV"),
                   AD_R47H_vs_CV = c("AD_R47H", "AD_CV"),
                   AD_R62H_vs_CV = c("AD_R62H", "AD_CV"),
                   AD_R47H_vs_R62H = c("AD_R47H", "AD_R62H"),
                   CV_AD_vs_Control = c("AD_CV", "Control_CV"),
                   R47H_AD_vs_Control = c("AD_R47H", "Control_R47H"),
                   R62H_AD_vs_Control = c("AD_R62H", "Control_R62H")
                   )


write_csv(gr_tab, file = paste0(out_dir,script_ind,"gr_tab_updated.csv"))



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


#manual credible interval plot (replaces plot_1D_intervals)
#  shows posterior effect estimate per cluster with 95% credible intervals,
#  faceted by model parameter, coloured by FDR significance
plot_sccomp_intervals = function(sccomp_res, title_label = ""){
  t_plot = sccomp_res[!is.na(sccomp_res$factor),]
  if (nrow(t_plot) == 0) return(NULL)
  t_plot$sig = ifelse(t_plot$c_FDR < 0.05, "FDR < 0.05", "n.s.")
  ggplot(t_plot, aes(x = c_effect, y = cluster_name, colour = sig)) +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = c_lower, xmax = c_upper), height = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("FDR < 0.05" = "red3", "n.s." = "grey50")) +
    facet_wrap(~parameter, scales = "free_x") +
    labs(x = "Effect (logit scale)", y = "", colour = "", title = title_label) +
    theme_light() +
    theme(strip.text = element_text(size = 7))
}


###########################################################
# save cluster//cell_type/cell_class labels to dataset
###########################################################
# maps the Seurat clusters (at resolution 0.3) to the manually assigned cluster names and cell types using table from LD_B03a
DefaultAssay(seur) = "SCT"

#add clusters with default resolution
seur$seurat_clusters = seur@meta.data[,paste0("SCT_snn_res.", clust_res)]

seur$cluster_name = clust_tab$cluster_name[match(seur$seurat_clusters, clust_tab$cluster)]
seur$cell_type = clust_tab$cell_type[match(seur$seurat_clusters, clust_tab$cluster)]
seur$cell_class = clust_tab$cell_class[match(seur$seurat_clusters, clust_tab$cluster)]

qsave(seur, file = paste0(out_dir, script_ind, "seur.qs"))



###########################################################
# plot with cluster/cell_type/cell_class names
###########################################################

#define grouping variables
gr = unique(gr_tab$group)
samples = unique(gr_tab$sample)
clusters = clust_tab$cluster
cluster_names = clust_tab$cluster_name
cell_types = unique(clust_tab$cell_type)
cell_classes = unique(clust_tab$cell_class)


#subsample seurat object for plotting for large datasets (else plots become too large (pdf with hundreds of MB))

seur0 = seur

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


#umap plots for cell cluster and cell class (labelled)
pl = list()

pl[["group"]] = DimPlot(seur, group.by = "group", shuffle = TRUE,
                        label = FALSE, reduction = "umap", 
                        pt.size = 0.01)+
  scale_color_manual(limits = gr, values = pal(gr))+
  labs(title = "UMAP by group")

pl[["cluster"]] = DimPlot(seur, group.by = "seurat_clusters", shuffle = TRUE, repel = TRUE, 
                          label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = clusters, values = pal(clusters))+
  NoLegend()+labs(title = "clusters")

pl[["cluster_dist"]] = DimPlot(seur, group.by = "seurat_clusters", shuffle = TRUE, repel = TRUE, 
                               label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = clusters, values = pal_dist(clusters))+
  NoLegend()+labs(title = "clusters")

pl[["cluster_name_umap"]] = DimPlot(seur, group.by = "cluster_name", shuffle = TRUE, repel = TRUE, 
                                    label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cluster_names, values = pal(cluster_names))+
  NoLegend()+labs(title = "cluster_names")

pl[["cluster_name_umap_dist"]] = DimPlot(seur, group.by = "cluster_name", shuffle = TRUE, repel = TRUE, 
                                         label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  NoLegend()+labs(title = "cluster_names")

pl[["cell_type_umap_dist"]] = DimPlot(seur, group.by = "cell_type", shuffle = TRUE, repel = TRUE, 
                                      label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cell_types, values = pal_dist(cell_types))+
  NoLegend()+labs(title = "cell_types")

pl[["cell_class_umap_dist"]] = DimPlot(seur, group.by = "cell_class", shuffle = TRUE, repel = TRUE, 
                                       label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cell_classes, values = pal_dist(cell_classes))+
  NoLegend()+labs(title = "cell_classes")

pdf(file = paste0(out_dir,script_ind, "UMAP_clusters_labelled.pdf"), width = 6, height = 5)
lapply(pl, function(x){x})
dev.off()


########################################################################################
### add plots split by group vs repl colored by cluster_name
########################################################################################
# this adds UMAPS per sample to check whether clusters are consistent across individual samples
seur = seur0

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


t1 = FetchData(seur, vars = c("umap_1", "umap_2", "cluster_name","group", "sample"))


#plot group vs sample

p1 = ggplot(data = t1, aes(x = umap_1, y = umap_2, color = cluster_name))+geom_point(size = 0.05, alpha = 0.5)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  ggplot2::facet_wrap(facets = factor(t1$sample, levels = samples), nrow = 2)+
  theme_bw()+
  labs(title = "UMAP by group vs sample")+RestoreLegend()

pdf(file = paste0(out_dir, script_ind, "UMAP_split_by_group_vs_sample.pdf"), width = 30, height = 5)
plot(p1)
dev.off()


################################################################################
#labelled dotplot for marker gene expression by cluster (split by cell type markers and subtype markers)
################################################################################
seur = seur0

DefaultAssay(seur) = "SCT"

p1 = DotPlot(seur, features = intersect(GOI$cell_type_markers, rownames(seur)), 
             group.by = "cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names)

p2 = DotPlot(seur, features = intersect(GOI$subtype_markers, rownames(seur)), 
             group.by = "cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names)

pdf(file = paste0(out_dir,script_ind, "Cell_markers_dotplot_clusters_labelled.pdf"), 
    width = 10, height = 4)
plot(p1)
dev.off()

pdf(file = paste0(out_dir,script_ind, "Cell_subset_markers_dotplot_clusters_labelled.pdf"),
    width = 15, height = 4)
plot(p2)
dev.off()

# size-adjusted version: wider figure, smaller gene labels
p2_adj = DotPlot(seur, features = intersect(GOI$subtype_markers, rownames(seur)),
                 group.by = "cluster_name", scale.by = "size") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
  scale_y_discrete(limits = cluster_names)

pdf(file = paste0(out_dir, script_ind, "Cell_subset_markers_dotplot_clusters_labelled_size_adj.pdf"),
    width = 25, height = 5)
plot(p2_adj)
dev.off()



################################################################################
#labelled dotplot for marker module activity by cluster (subtype marker set)
################################################################################
# do the cluster resamble published subtypes based on the two publications? 
seur = seur0

for(subtype_dataset in names(subtype_markers)){
  
  seur <- AddModuleScore(seur,
                          features = subtype_markers[[subtype_dataset]],
                          name = subtype_dataset)
  
  colnames(seur@meta.data)[grepl(subtype_dataset, colnames(seur@meta.data))] = 
    paste0(subtype_dataset, "_", names(subtype_markers[[subtype_dataset]]))
}


p1 = DotPlot(seur, features = colnames(seur@meta.data)[grepl("Green24", colnames(seur@meta.data))],
             group.by = "cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names)

p2 = DotPlot(seur, features = colnames(seur@meta.data)[grepl("Gazestani", colnames(seur@meta.data))], 
             group.by = "cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names)

pdf(file = paste0(out_dir,script_ind, "Cell_markers_module_score_dotplot_Green24.pdf"),
    width = 6, height = 8)
plot(p1)
dev.off()

pdf(file = paste0(out_dir,script_ind, "Cell_markers_module_score_dotplot_Gazestani23.pdf"),
    width = 6, height = 8)
plot(p2)
dev.off()

# Green24: size-adjusted wide overview
green24_cols = colnames(seur@meta.data)[grepl("Green24", colnames(seur@meta.data))]

p_green24_adj = DotPlot(seur, features = green24_cols,
                        group.by = "cluster_name", scale.by = "size") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
  scale_y_discrete(limits = cluster_names)

pdf(file = paste0(out_dir, script_ind, "Cell_markers_module_score_dotplot_Green24_size_adj.pdf"),
    width = 25, height = 6)
plot(p_green24_adj)
dev.off()

# Green24: astrocyte subtypes only
green24_ast_cols = green24_cols[grepl("Green24_Ast", green24_cols)]

p_green24_ast = DotPlot(seur, features = green24_ast_cols,
                        group.by = "cluster_name", scale.by = "size") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
  scale_y_discrete(limits = cluster_names)

pdf(file = paste0(out_dir, script_ind, "Cell_markers_module_score_dotplot_Green24_AST_only.pdf"),
    width = 8, height = 6)
plot(p_green24_ast)
dev.off()



#############################################################################
# De novo cluster marker identification
#############################################################################

###identify markers
# tests every gene in the SCT assay across all clusters
seur = seur0

t1 <- FindAllMarkers(seur,group.by = "cluster_name", only.pos = TRUE)
t1 = as_tibble(t1[t1$p_val_adj<0.001 & t1$avg_log2FC > 0.5, ])
t1$cluster = as.character(t1$cluster)
t2 = t1[order(t1$cluster, -t1$avg_log2FC),]
seur_markers = t2

write_csv(seur_markers, file = paste0(out_dir,script_ind, "Seurat_markers.csv"))
# tells us for each of the clusters, which genes are specifically upregulated in that cluster relative to the other astrocytes, and by how much

###plot top10 markers per cluster (max 200 cells/cluster)

top_markers = seur_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

write_csv(top_markers, file = paste0(out_dir, script_ind, "Seurat_markers_top10_per_cluster.csv"))

#subsample seur (done for the following DoHeatmap)
meta = seur@meta.data
v1 = NULL

for (cl in cluster_names){
  v2 = rownames(meta[meta$cluster_name == cl,])
  if (length(v2)>200){
    set.seed(1234)
    v2 = sample(v2, 200)
    }
  v1 = c(v1,v2)
}

seur_plot = seur[, v1]


#plot
# this step serves as a visual quality check. It shows the top markers for 200 cells per cluster (one column per cell).
# if there is a biological meaningfuls distinction, we should see a clear pattern of marker expression across clusters.
p1 = DoHeatmap(seur_plot, group.by = "cluster_name", slot = "data", 
               features = top_markers$gene) + NoLegend()

pdf(file = paste0(out_dir,script_ind, "Seurat_markers_top10_heatmap.pdf"), 
    width = 10, height = 10)
plot(p1)
dev.off()

# One dot per gene per cluster with the full dataset
p1 = DotPlot(seur, features = unique(top_markers$gene), 
             group.by = "cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names)

pdf(file = paste0(out_dir,script_ind, "Seurat_markers_top10_dotplot.pdf"), 
    width = length(cluster_names)*2.5+1, height = length(cluster_names)/3+1.5)
plot(p1)
dev.off()


### Gene Ontology (GO) analysis
# standardised database for gene functions, and biological processes
# uses all markers from the filtered FindAllMarkers output to check "are there any biological processes that appear in this gene list more often than expected by chance?
GO_list = list()
GO_results_tab = NULL

for (cl in cluster_names){
  
  message("\n          *** GO analysis DEGs vs cluster ", cl, " - ", Sys.time(),"\n")
  ego = NULL

  ego = enrichGO(gene         = seur_markers$gene[seur_markers$cluster == cl],
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP", # restricts to biological processes
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
  
  if (!is.null(ego)){

    t1 = ego@result[ego@result$p.adjust<=0.05,]

      if (nrow(t1)>0){

        GO_list[[cl]] = ego
        t2 = cbind(cluster_name = cl, t1)
        GO_results_tab = rbind(GO_results_tab, t2)
    }
  }
}

write_csv(GO_results_tab, file = paste0(out_dir,script_ind, "Seurat_markers_GO_terms.csv"))

save(GO_list, file = paste0(out_dir,script_ind, "Seurat_markers_GO_analysis.rda"))


### visualise GO analysis (network plot)
# each clusters with >5 significantly enriched GO terms gets a network plot. 

pl = list()

for (cl in names(GO_list)){
  
  ego = GO_list[[cl]]
  t1 = ego@result

  if (nrow(t1[t1$p.adjust<0.01,])>5){
      edo = enrichplot::pairwise_termsim(GO_list[[cl]])
      pl[[cl]] = emapplot(edo, showCategory = 100, repel = TRUE)+labs(title = cl) 
    }
}

pdf(file = paste0(out_dir,script_ind, "Seurat_markers_GO_terms_network_plot.pdf"), 
    width = 20, height = 20)
lapply(pl, function(x){x})
dev.off()



### visualise GO analysis (dotplot by cluster for top 10 terms by cluster)

t1 = GO_results_tab

t2 = NULL

for (cl in unique(t1$cluster_name)){
  t3 = t1[t1$cluster_name == cl,]
  if (nrow(t3)>10){t3 = t3[1:10,]}
  t2 = rbind(t2, t3)
}

t2 = t1[t1$Description %in% t2$Description,]

p1 = ggplot(t2, aes(x = cluster_name, y = Description, size = Count, colour = -log10(p.adjust)))+
  geom_point()+
  scale_color_gradient(low = "blue", high = "red")+
  scale_x_discrete(limits = cluster_names)+
  scale_y_discrete(limits = unique(t2$Description))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


pdf(file = paste0(out_dir,script_ind, "Seurat_markers_GO_terms_dotplot_top10_terms.pdf"), 
    width = 9, height = 15)
plot(p1)
dev.off()

### visualise GO analysis (dotplot by cluster for top 5 terms by cluster)

t1 = GO_results_tab

t2 = NULL

for (cl in unique(t1$cluster_name)){
  t3 = t1[t1$cluster_name == cl,]
  if (nrow(t3)>5){t3 = t3[1:5,]}
  t2 = rbind(t2, t3)
}

t2 = t1[t1$Description %in% t2$Description,]

p1 = ggplot(t2, aes(x = cluster_name, y = Description, size = Count, colour = -log10(p.adjust)))+
  geom_point()+
  scale_color_gradient(low = "blue", high = "red")+
  scale_x_discrete(limits = cluster_names)+
  scale_y_discrete(limits = unique(t2$Description))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


pdf(file = paste0(out_dir,script_ind, "Seurat_markers_GO_terms_dotplot_top5_terms.pdf"), 
    width = 9, height = 8)
plot(p1)
dev.off()



#############################################################################
# cluster abundance quantification
#############################################################################
# quantifies how many cells for each sample end up in each ubcluster to answer "Are certain astrocyte subtypes more or less abundant in certain experimental groups?"
seur = seur0

meta = seur@meta.data

# create table with one row for each cluster for each sample

stat_tab = tibble(cluster = unlist(lapply(cluster_names, rep, length.out = length(samples))), 
                  sample = rep(samples, length(cluster_names)))
stat_tab$group = gr_tab$group[match(stat_tab$sample, gr_tab$sample)] #add group assignment

#count cells per cluster per sample, add to stat_tab
t1 = meta %>% group_by(cluster_name, sample) %>% summarize(N_cells = n())
stat_tab$N_cells = t1$N_cells[match( paste0(stat_tab$cluster,stat_tab$sample), paste0(t1$cluster_name,t1$sample) )]
stat_tab$N_cells[is.na(stat_tab$N_cells)] = 0

#add total cells per sample and per cluster, fraction of cluster, fraction of sample
N_sample = stat_tab %>% group_by(sample) %>% summarize(N = sum(N_cells))
stat_tab$N_sample = N_sample$N[match(stat_tab$sample, N_sample$sample)]
stat_tab$fract_sample = stat_tab$N_cells / stat_tab$N_sample
N_cluster = stat_tab %>% group_by(cluster) %>% summarize(N = sum(N_cells))
stat_tab$N_cluster = N_cluster$N[match(stat_tab$cluster, N_cluster$cluster)]
stat_tab$fract_cluster = stat_tab$N_cells / stat_tab$N_cluster

write_csv(stat_tab, file = paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster.csv"))


### crossbar-dotplot quantification of cluster contribution fraction of sample

t1 = stat_tab
t2 = t1 %>% group_by(cluster, group) %>% 
  summarise(mean_fract = mean(fract_sample), sd_fract = sd(fract_sample))

t1$group = factor(t1$group, levels = gr)
t2$group = factor(t2$group, levels = gr) #required to fix order of groups and samples
#t2$sample = factor(t2$sample, levels = samples)


p1 = ggplot()+
  geom_col(data = t2, aes(x = cluster, y = mean_fract, color = group, group = group), 
           fill = "grey90", position = position_dodge(), width = 0.7, lwd = 0.3)+
  geom_errorbar(data = t2, aes(x = cluster,
                               ymin = mean_fract-sd_fract,
                               y = mean_fract,
                               ymax = mean_fract+sd_fract,
                               color = group,
                               group = group),
                position = position_dodge(width = 0.7), width = 0.3, lwd = 0.2)+
  geom_point(data = t1, aes(x = cluster, y = fract_sample, color = group, group = group), 
             position = position_dodge(width = 0.7), size = 0.5, stroke = 0.3)+
  geom_hline(yintercept = 0)+
  scale_x_discrete(limits = cluster_names)+
  scale_color_manual(limits = gr, values = pal(gr))+
  scale_fill_manual(limits = gr, values = pal(gr))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster.pdf"), 
    width = 6, height = 3)
plot(p1)
dev.off()



### crossbar-dotplot quantification of cluster contribution fraction of sample (alternative plot test)

t1 = stat_tab
t2 = t1 %>% group_by(cluster, group) %>% 
  summarise(mean_fract = mean(fract_sample), sd_fract = sd(fract_sample))

p1 = ggplot()+
  geom_col(data = t2, aes(x = cluster, y = mean_fract, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.3)+
  geom_errorbar(data = t2, aes(x = cluster,
                               ymin = mean_fract-sd_fract,
                               y = mean_fract,
                               ymax = mean_fract+sd_fract,
                               color = group),
                position = position_dodge(width = 0.5), width = 0.3, lwd = 0.2)+
  geom_point(data = t1, aes(x = cluster, y = fract_sample, color = group), 
             position = position_dodge(width = 0.5), size = 0.5, stroke = 0.3)+
  geom_hline(yintercept = 0)+
  scale_x_discrete(limits = cluster_names)+
  scale_color_manual(limits = gr, values = pal(gr))+
  scale_fill_manual(limits = gr, values = pal(gr))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster2.pdf"), 
    width = 8, height = 5)
plot(p1)
dev.off()

gc()


###########################################################
# sccomp differential cell cluster abundance analysis
###########################################################
# adds statistics on the abundance testing form just before
seur$group = factor(seur$group, levels = gr) #required to fix order of groups

sccomp_result =
  seur |>
  sccomp_glm(
    formula_composition = ~ group,
    .sample =  sample,
    .cell_group = cluster_name,
    bimodal_mean_variability_association = TRUE,
    cores = 8
  )

message("sccomp_result columns: ", paste(names(sccomp_result), collapse = ", "))

# handle sccomp version differences: column may be "factor" or "parameter"
if (!"factor" %in% names(sccomp_result) & "parameter" %in% names(sccomp_result)){
  sccomp_result$factor = sccomp_result$parameter
}

write_csv(sccomp_result, paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group.csv"))


### credible interval estimates plot (originally plot_1D_intervals)
# try original function, save to _orig file
tryCatch({
  pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates_orig.pdf"),
      width = 8, height = 5)
  sccomp_result |> plot_1D_intervals()
  dev.off()
}, error = function(e) {
  message("plot_1D_intervals (group) failed: ", conditionMessage(e))
  try(dev.off(), silent = TRUE)
})

# manual ggplot version
pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates.pdf"),
    width = 10, height = 6)
p = plot_sccomp_intervals(sccomp_result, title_label = "sccomp: ~ group (estimates)")
if (!is.null(p)) print(p)
dev.off()


### boxplot for significant clusters (originally sccomp_boxplot)
t1 = sccomp_result
t2 = t1[!is.na(t1$factor),]
t3 = t2[t2$c_FDR<0.05,]

if (nrow(t3)>0){
  # try original function, save to _orig file
  tryCatch({
    pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot_orig.pdf"),
        width = 8, height = 8)
    sccomp_result |> sccomp_boxplot(factor = "group")
    dev.off()
  }, error = function(e) {
    message("sccomp_boxplot (group) failed: ", conditionMessage(e))
    try(dev.off(), silent = TRUE)
  })

  # manual ggplot version
  tryCatch({
    sig_clusters = unique(t3$cluster_name)
    plot_data = stat_tab[stat_tab$cluster %in% sig_clusters,]

    pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot.pdf"),
        width = 8, height = 8)
    p = ggplot(plot_data, aes(x = group, y = fract_sample, fill = group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 1) +
      facet_wrap(~cluster, scales = "free_y") +
      scale_fill_manual(values = pal(gr)) +
      labs(y = "Fraction of sample", title = "Significant clusters - abundance by group") +
      theme_light()
    print(p)
    dev.off()
  }, error = function(e) { message("Boxplot by group failed: ", conditionMessage(e)); try(dev.off(), silent=TRUE) })
}


### abundance barplot with sccomp significance asterisks
tryCatch({
  t1 = stat_tab
  t2 = t1 %>% group_by(cluster, group) %>%
    summarise(mean_fract = mean(fract_sample), sd_fract = sd(fract_sample), .groups = "drop")

  t1$group = factor(t1$group, levels = gr)
  t2$group = factor(t2$group, levels = gr)

  # extract significance from sccomp result
  # parameter values are like "groupControl_R47H"; strip prefix to get group name
  sig_tab = sccomp_result %>%
    filter(factor == "group") %>%
    mutate(
      group = sub("^group", "", parameter),
      sig_label = case_when(
        c_FDR < 0.05 ~ "**",
        c_FDR < 0.1  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    filter(sig_label != "") %>%
    select(cluster_name, group, sig_label, c_FDR)

  if (nrow(sig_tab) > 0){
    sig_tab = sig_tab %>%
      left_join(t2, by = c("cluster_name" = "cluster", "group" = "group")) %>%
      mutate(y_pos = mean_fract + sd_fract + max(t2$mean_fract, na.rm = TRUE) * 0.02)

    sig_tab$group = factor(sig_tab$group, levels = gr)

    n_groups = length(gr)
    dodge_width = 0.7
    sig_tab = sig_tab %>%
      mutate(
        x_num = as.numeric(factor(cluster_name, levels = cluster_names)),
        group_idx = as.numeric(group),
        x_dodge = x_num + dodge_width * (group_idx - (n_groups + 1) / 2) / n_groups
      )
  }

  p_sig = ggplot() +
    geom_col(data = t2, aes(x = cluster, y = mean_fract, color = group, group = group),
             fill = "grey90", position = position_dodge(width = 0.7), width = 0.7, lwd = 0.3) +
    geom_errorbar(data = t2, aes(x = cluster,
                                 ymin = mean_fract - sd_fract,
                                 y = mean_fract,
                                 ymax = mean_fract + sd_fract,
                                 color = group,
                                 group = group),
                  position = position_dodge(width = 0.7), width = 0.4, lwd = 0.3) +
    geom_point(data = t1, aes(x = cluster, y = fract_sample, color = group, group = group),
               position = position_dodge(width = 0.7), size = 0.3, stroke = 0.2) +
    geom_hline(yintercept = 0) +
    scale_x_discrete(limits = cluster_names) +
    scale_color_manual(limits = gr, values = pal(gr)) +
    scale_fill_manual(limits = gr, values = pal(gr)) +
    labs(caption = "** p<0.05, * p<0.1 (sccomp FDR vs Control_CV)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.caption = element_text(size = 7))

  if (nrow(sig_tab) > 0){
    p_sig = p_sig +
      geom_text(data = sig_tab, aes(x = x_dodge, y = y_pos, label = sig_label),
                size = 4, vjust = 0)
  }

  pdf(file = paste0(out_dir, script_ind, "cell_abundance_by_sample_cluster_sig.pdf"),
      width = 10, height = 4)
  print(p_sig)
  dev.off()
}, error = function(e) {
  message("Abundance significance plot failed: ", conditionMessage(e))
  try(dev.off(), silent = TRUE)
})



###########################################################
# plot cluster abundance vs plaque_dens by TREM2Variant
###########################################################
# this section asks whther the abundance of certain clusters correlates with plaque density, and whether this correlation is different between TREM2 variant carriers and non-carriers
# for each clsuter we see: as plaque pathology increases, does the cluster ge more or less abundant? and does this differ between the variants?
t1 = stat_tab

t1$plaque_dens = gr_tab$plaque_dens[match(t1$sample, gr_tab$sample)]
t1$TREM2Variant = gr_tab$TREM2Variant[match(t1$sample, gr_tab$sample)]

stat_tab = t1


###plot with linear trend

pl = list()

for (cl in cluster_names){
  
  t2 = t1[t1$cluster == cl,]
  
  pl[[paste0(cl)]] = ggplot(t2, aes(x = plaque_dens, y = fract_sample, color = TREM2Variant))+
    geom_hline(yintercept = 0)+
    geom_smooth(method = "lm")+
    geom_point()+
    scale_color_manual(limits = unique(t2$TREM2Variant), values = pal(unique(t2$TREM2Variant)))+
    labs(title = paste0(cl , " abundance"))+
    theme_light()
  
}

pdf(file = paste0(out_dir,script_ind,"cell_abundance_vs_plaque_dens_by_TREM2Variant_by_cluster.pdf"), 
    width = 4, height = 3)
lapply(pl, function(x){x})
dev.off()



###########################################################
# plot cluster abundance vs PHF1 by TREM2Variant
###########################################################
# same as just above but for tau pathology
t1 = stat_tab

t1$pctPHF1PositiveArea = gr_tab$pctPHF1PositiveArea[match(t1$sample, gr_tab$sample)]
t1$TREM2Variant = gr_tab$TREM2Variant[match(t1$sample, gr_tab$sample)]

stat_tab = t1


###plot with linear trend

pl = list()

for (cl in cluster_names){
  
  t2 = t1[t1$cluster == cl,]
  
  pl[[paste0(cl)]] = ggplot(t2, aes(x = pctPHF1PositiveArea, y = fract_sample, color = TREM2Variant))+
    geom_hline(yintercept = 0)+
    geom_smooth(method = "lm")+
    geom_point()+
    scale_color_manual(limits = unique(t2$TREM2Variant), values = pal(unique(t2$TREM2Variant)))+
    labs(title = paste0(cl , " abundance"))+
    theme_light()
  
}

pdf(file = paste0(out_dir,script_ind,"cell_abundance_vs_pctPHF1PositiveArea_by_TREM2Variant_by_cluster.pdf"),
    width = 4, height = 3)
lapply(pl, function(x){x})
dev.off()


###########################################################
# combined facet regression plots (all clusters in one figure)
###########################################################

# plaque density combined facet plot
p_plaque = ggplot(stat_tab, aes(x = plaque_dens, y = fract_sample, color = TREM2Variant)) +
  geom_hline(yintercept = 0) +
  geom_smooth(method = "lm") +
  geom_point(size = 0.8) +
  scale_color_manual(values = pal(unique(stat_tab$TREM2Variant))) +
  facet_wrap(~ factor(cluster, levels = cluster_names), scales = "free_y") +
  labs(x = "Plaque density", y = "Fraction of sample") +
  theme_light() +
  theme(strip.text = element_text(size = 7))

pdf(file = paste0(out_dir, script_ind, "cell_abundance_vs_plaque_dens_by_TREM2Variant_facet.pdf"),
    width = 14, height = 10)
plot(p_plaque)
dev.off()

# tau (PHF1) combined facet plot
t1_tau = stat_tab[!is.na(stat_tab$pctPHF1PositiveArea), ]

p_tau = ggplot(t1_tau, aes(x = pctPHF1PositiveArea, y = fract_sample, color = TREM2Variant)) +
  geom_hline(yintercept = 0) +
  geom_smooth(method = "lm") +
  geom_point(size = 0.8) +
  scale_color_manual(values = pal(unique(t1_tau$TREM2Variant))) +
  facet_wrap(~ factor(cluster, levels = cluster_names), scales = "free_y") +
  labs(x = "% PHF1 positive area", y = "Fraction of sample") +
  theme_light() +
  theme(strip.text = element_text(size = 7))

pdf(file = paste0(out_dir, script_ind, "cell_abundance_vs_pctPHF1_by_TREM2Variant_facet.pdf"),
    width = 14, height = 10)
plot(p_tau)
dev.off()


###########################################################
# sccomp differential cell cluster abundance analysis by TREM2Variant vs plaque_dens
###########################################################
# adds statistic to the questions jsut answered before
seur = seur0

seur$TREM2Variant = factor(seur$TREM2Variant, levels = unique(gr_tab$TREM2Variant))
seur = seur[,!is.na(seur$plaque_dens)]

sccomp_result =
  seur |>
  sccomp_glm(
    formula_composition = ~ TREM2Variant*plaque_dens,
    .sample =  sample,
    .cell_group = cluster_name,
    bimodal_mean_variability_association = TRUE,
    cores = 8
  )

message("sccomp_result columns: ", paste(names(sccomp_result), collapse = ", "))

if (!"factor" %in% names(sccomp_result) & "parameter" %in% names(sccomp_result)){
  sccomp_result$factor = sccomp_result$parameter
}

write_csv(sccomp_result, paste0(out_dir,script_ind,"sccomp_cell_abundance_by_TREM2Variant_vs_plaque_dens.csv"))


### credible interval estimates plot (originally plot_1D_intervals)
# try original function, save to _orig file
tryCatch({
  pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates_TREM2Variant_vs_plaque_dens_orig.pdf"),
      width = 8, height = 5)
  sccomp_result |> plot_1D_intervals()
  dev.off()
}, error = function(e) {
  message("plot_1D_intervals (TREM2*plaque) failed: ", conditionMessage(e))
  try(dev.off(), silent = TRUE)
})

# manual ggplot version
pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates_TREM2Variant_vs_plaque_dens.pdf"),
    width = 10, height = 6)
p = plot_sccomp_intervals(sccomp_result, title_label = "sccomp: ~ TREM2Variant * plaque_dens (estimates)")
if (!is.null(p)) print(p)
dev.off()


### boxplot for significant clusters (originally sccomp_boxplot)
t1 = sccomp_result
t2 = t1[!is.na(t1$factor),]
t3 = t2[t2$c_FDR<0.05,]

if (nrow(t3)>0){
  # try original function, save to _orig file
  tryCatch({
    pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot_TREM2Variant_vs_plaque_dens_orig.pdf"),
        width = 8, height = 8)
    sccomp_result |> sccomp_boxplot(factor = "TREM2Variant")
    dev.off()
  }, error = function(e) {
    message("sccomp_boxplot (TREM2*plaque) failed: ", conditionMessage(e))
    try(dev.off(), silent = TRUE)
  })

  # manual ggplot version
  tryCatch({
    sig_clusters = unique(t3$cluster_name)
    plot_data = stat_tab[stat_tab$cluster %in% sig_clusters,]

    pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot_TREM2Variant_vs_plaque_dens.pdf"),
        width = 8, height = 8)
    p = ggplot(plot_data, aes(x = TREM2Variant, y = fract_sample, fill = TREM2Variant)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 1) +
      facet_wrap(~cluster, scales = "free_y") +
      scale_fill_manual(values = pal(unique(plot_data$TREM2Variant))) +
      labs(y = "Fraction of sample", title = "Significant clusters - abundance by TREM2Variant") +
      theme_light()
    print(p)
    dev.off()
  }, error = function(e) { message("Boxplot by TREM2Variant failed: ", conditionMessage(e)); try(dev.off(), silent=TRUE) })
}


###########################################################
# sccomp differential cell cluster abundance analysis by TREM2Variant vs plaque_dens corrected for APOE and CD33
###########################################################
# same as sccomp just before but adds correction for APOE and CD33 genotype
seur = seur0

seur$TREM2Variant = factor(seur$TREM2Variant, levels = unique(gr_tab$TREM2Variant))
seur$APOEgroup = factor(seur$APOEgroup, levels = unique(gr_tab$APOEgroup))
seur$CD33Group = factor(seur$CD33Group, levels = unique(gr_tab$CD33Group))

seur = seur[,!is.na(seur$plaque_dens)&!is.na(seur$APOEgroup)&!is.na(seur$CD33Group)]

sccomp_result =
  seur |>
  sccomp_glm(
    formula_composition = ~ APOEgroup + CD33Group + TREM2Variant*plaque_dens,
    .sample =  sample,
    .cell_group = cluster_name,
    bimodal_mean_variability_association = TRUE,
    cores = 8
  )

message("sccomp_result columns: ", paste(names(sccomp_result), collapse = ", "))

if (!"factor" %in% names(sccomp_result) & "parameter" %in% names(sccomp_result)){
  sccomp_result$factor = sccomp_result$parameter
}

write_csv(sccomp_result, paste0(out_dir,script_ind,"sccomp_cell_abundance_by_TREM2Variant_vs_plaque_dens_corr_APOE_CD33.csv"))


### credible interval estimates plot (originally plot_1D_intervals)
# try original function, save to _orig file
tryCatch({
  pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates_TREM2Variant_vs_plaque_dens_corr_APOE_CD33_orig.pdf"),
      width = 8, height = 5)
  sccomp_result |> plot_1D_intervals()
  dev.off()
}, error = function(e) {
  message("plot_1D_intervals (TREM2*plaque corr.) failed: ", conditionMessage(e))
  try(dev.off(), silent = TRUE)
})

# manual ggplot version
pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates_TREM2Variant_vs_plaque_dens_corr_APOE_CD33.pdf"),
    width = 10, height = 6)
p = plot_sccomp_intervals(sccomp_result, title_label = "sccomp: ~ APOEgroup + CD33Group + TREM2Variant * plaque_dens (estimates)")
if (!is.null(p)) print(p)
dev.off()


### boxplot for significant clusters (originally sccomp_boxplot)
t1 = sccomp_result
t2 = t1[!is.na(t1$factor),]
t3 = t2[t2$c_FDR<0.05,]

if (nrow(t3)>0){
  # try original function, save to _orig file
  tryCatch({
    pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot_TREM2Variant_vs_plaque_dens_corr_APOE_CD33_orig.pdf"),
        width = 8, height = 8)
    sccomp_result |> sccomp_boxplot(factor = "TREM2Variant")
    dev.off()
  }, error = function(e) {
    message("sccomp_boxplot (TREM2*plaque corr.) failed: ", conditionMessage(e))
    try(dev.off(), silent = TRUE)
  })

  # manual ggplot version
  tryCatch({
    sig_clusters = unique(t3$cluster_name)
    plot_data = stat_tab[stat_tab$cluster %in% sig_clusters,]

    pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot_TREM2Variant_vs_plaque_dens_corr_APOE_CD33.pdf"),
        width = 8, height = 8)
    p = ggplot(plot_data, aes(x = TREM2Variant, y = fract_sample, fill = TREM2Variant)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 1) +
      facet_wrap(~cluster, scales = "free_y") +
      scale_fill_manual(values = pal(unique(plot_data$TREM2Variant))) +
      labs(y = "Fraction of sample", title = "Significant clusters - abundance by TREM2Variant (corr. APOE, CD33)") +
      theme_light()
    print(p)
    dev.off()
  }, error = function(e) { message("Boxplot by TREM2Variant (corr.) failed: ", conditionMessage(e)); try(dev.off(), silent=TRUE) })
}


###########################################################
# sccomp differential cell cluster abundance analysis by TREM2Variant vs tau (PHF1)
###########################################################
# mirrors the plaque_dens sccomp analyses but uses tau pathology (pctPHF1PositiveArea) instead
seur = seur0

seur$TREM2Variant = factor(seur$TREM2Variant, levels = unique(gr_tab$TREM2Variant))
seur = seur[,!is.na(seur$pctPHF1PositiveArea)]

sccomp_result =
  seur |>
  sccomp_glm(
    formula_composition = ~ TREM2Variant*pctPHF1PositiveArea,
    .sample =  sample,
    .cell_group = cluster_name,
    bimodal_mean_variability_association = TRUE,
    cores = 8
  )

message("sccomp_result columns: ", paste(names(sccomp_result), collapse = ", "))

if (!"factor" %in% names(sccomp_result) & "parameter" %in% names(sccomp_result)){
  sccomp_result$factor = sccomp_result$parameter
}

write_csv(sccomp_result, paste0(out_dir,script_ind,"sccomp_cell_abundance_by_TREM2Variant_vs_PHF1.csv"))


### credible interval estimates plot
tryCatch({
  pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates_TREM2Variant_vs_PHF1_orig.pdf"),
      width = 8, height = 5)
  sccomp_result |> plot_1D_intervals()
  dev.off()
}, error = function(e) {
  message("plot_1D_intervals (TREM2*PHF1) failed: ", conditionMessage(e))
  try(dev.off(), silent = TRUE)
})

pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates_TREM2Variant_vs_PHF1.pdf"),
    width = 10, height = 6)
p = plot_sccomp_intervals(sccomp_result, title_label = "sccomp: ~ TREM2Variant * pctPHF1PositiveArea (estimates)")
if (!is.null(p)) print(p)
dev.off()


### boxplot for significant clusters
t1 = sccomp_result
t2 = t1[!is.na(t1$factor),]
t3 = t2[t2$c_FDR<0.05,]

if (nrow(t3)>0){
  tryCatch({
    pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot_TREM2Variant_vs_PHF1_orig.pdf"),
        width = 8, height = 8)
    sccomp_result |> sccomp_boxplot(factor = "TREM2Variant")
    dev.off()
  }, error = function(e) {
    message("sccomp_boxplot (TREM2*PHF1) failed: ", conditionMessage(e))
    try(dev.off(), silent = TRUE)
  })

  tryCatch({
    sig_clusters = unique(t3$cluster_name)
    plot_data = stat_tab[stat_tab$cluster %in% sig_clusters,]

    pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot_TREM2Variant_vs_PHF1.pdf"),
        width = 8, height = 8)
    p = ggplot(plot_data, aes(x = TREM2Variant, y = fract_sample, fill = TREM2Variant)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 1) +
      facet_wrap(~cluster, scales = "free_y") +
      scale_fill_manual(values = pal(unique(plot_data$TREM2Variant))) +
      labs(y = "Fraction of sample", title = "Significant clusters - abundance by TREM2Variant (PHF1)") +
      theme_light()
    print(p)
    dev.off()
  }, error = function(e) { message("Boxplot by TREM2Variant (PHF1) failed: ", conditionMessage(e)); try(dev.off(), silent=TRUE) })
}


###########################################################
# sccomp differential cell cluster abundance analysis by TREM2Variant vs tau (PHF1) corrected for APOE and CD33
###########################################################
seur = seur0

seur$TREM2Variant = factor(seur$TREM2Variant, levels = unique(gr_tab$TREM2Variant))
seur$APOEgroup = factor(seur$APOEgroup, levels = unique(gr_tab$APOEgroup))
seur$CD33Group = factor(seur$CD33Group, levels = unique(gr_tab$CD33Group))

seur = seur[,!is.na(seur$pctPHF1PositiveArea)&!is.na(seur$APOEgroup)&!is.na(seur$CD33Group)]

sccomp_result =
  seur |>
  sccomp_glm(
    formula_composition = ~ APOEgroup + CD33Group + TREM2Variant*pctPHF1PositiveArea,
    .sample =  sample,
    .cell_group = cluster_name,
    bimodal_mean_variability_association = TRUE,
    cores = 8
  )

message("sccomp_result columns: ", paste(names(sccomp_result), collapse = ", "))

if (!"factor" %in% names(sccomp_result) & "parameter" %in% names(sccomp_result)){
  sccomp_result$factor = sccomp_result$parameter
}

write_csv(sccomp_result, paste0(out_dir,script_ind,"sccomp_cell_abundance_by_TREM2Variant_vs_PHF1_corr_APOE_CD33.csv"))


### credible interval estimates plot
tryCatch({
  pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates_TREM2Variant_vs_PHF1_corr_APOE_CD33_orig.pdf"),
      width = 8, height = 5)
  sccomp_result |> plot_1D_intervals()
  dev.off()
}, error = function(e) {
  message("plot_1D_intervals (TREM2*PHF1 corr.) failed: ", conditionMessage(e))
  try(dev.off(), silent = TRUE)
})

pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates_TREM2Variant_vs_PHF1_corr_APOE_CD33.pdf"),
    width = 10, height = 6)
p = plot_sccomp_intervals(sccomp_result, title_label = "sccomp: ~ APOEgroup + CD33Group + TREM2Variant * pctPHF1PositiveArea (estimates)")
if (!is.null(p)) print(p)
dev.off()


### boxplot for significant clusters
t1 = sccomp_result
t2 = t1[!is.na(t1$factor),]
t3 = t2[t2$c_FDR<0.05,]

if (nrow(t3)>0){
  tryCatch({
    pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot_TREM2Variant_vs_PHF1_corr_APOE_CD33_orig.pdf"),
        width = 8, height = 8)
    sccomp_result |> sccomp_boxplot(factor = "TREM2Variant")
    dev.off()
  }, error = function(e) {
    message("sccomp_boxplot (TREM2*PHF1 corr.) failed: ", conditionMessage(e))
    try(dev.off(), silent = TRUE)
  })

  tryCatch({
    sig_clusters = unique(t3$cluster_name)
    plot_data = stat_tab[stat_tab$cluster %in% sig_clusters,]

    pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot_TREM2Variant_vs_PHF1_corr_APOE_CD33.pdf"),
        width = 8, height = 8)
    p = ggplot(plot_data, aes(x = TREM2Variant, y = fract_sample, fill = TREM2Variant)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 1) +
      facet_wrap(~cluster, scales = "free_y") +
      scale_fill_manual(values = pal(unique(plot_data$TREM2Variant))) +
      labs(y = "Fraction of sample", title = "Significant clusters - abundance by TREM2Variant (PHF1, corr. APOE, CD33)") +
      theme_light()
    print(p)
    dev.off()
  }, error = function(e) { message("Boxplot by TREM2Variant (PHF1 corr.) failed: ", conditionMessage(e)); try(dev.off(), silent=TRUE) })
}


#get info on version of R, used packages etc
sessionInfo()


message("\n\n##########################################################################\n",
        "# Completed C03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")
