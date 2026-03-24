message("\n\n##########################################################################\n",
        "# Start C03: subcluster characterisation ", Sys.time(),
        "\n##########################################################################\n",
        "\n   ",
        "\n##########################################################################\n\n")

# Open packages necessary for analysis.
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
# the not called packages might be dependencies from other packages, therefore i keep them for now

### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "LD_B04a_"

#specify output directory
out_dir = paste0(main_dir,"LD_B_AST_analysis_output/")

### load group and file info for analysis dataset
gr_tab = read_csv("data_TREM2_michael/A_input/group_tab.csv") # metadata table with group and sample information for all samples in the dataset, used for plotting and stats

###load seurat dataset --> cleaned astrocytes
seur = qread(file = paste0(main_dir, "data_TREM2_michael/B_load_from_scflow_subcluster/LD_B03a_seur.qs"))

#add plaque density (Amyloid plaques)
t1 = read_csv("data_TREM2_michael/A_input/TREM2_plaque_data_Sam.csv")
gr_tab$plaque_dens = t1$TotalDensity[match(gr_tab$BrainBankNetworkIDFormatted, t1$BrainBankNetworkIDFormatted)]
seur$plaque_dens = t1$TotalDensity[match(seur$BrainBankNetworkIDFormatted, t1$BrainBankNetworkIDFormatted)]

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

### get all subtype markers from Gazestani et al., 2023, Manucso et al., 2024 --> 

subtype_markers = list()

t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/Green24_S2_subpopulation_markers.csv"))

for (cl in unique(t1$cluster)){
  t2 = t1[t1$cluster == cl & t1$avg_log2FC>log2(1.2) & t1$p_val_adj<0.05,]
  subtype_markers$Mancuso24[[cl]] = t2$gene
}

lengths(subtype_markers$Mancuso24)


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


###########################################################
# save cluster//cell_type/cell_class labels to dataset
###########################################################

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



################################################################################
#labelled dotplot for marker module activity by cluster (subtype marker set)
################################################################################

seur = seur0

for(subtype_dataset in names(subtype_markers)){
  
  seur <- AddModuleScore(seur,
                          features = subtype_markers[[subtype_dataset]],
                          name = subtype_dataset)
  
  colnames(seur@meta.data)[grepl(subtype_dataset, colnames(seur@meta.data))] = 
    paste0(subtype_dataset, "_", names(subtype_markers[[subtype_dataset]]))
}


p1 = DotPlot(seur, features = colnames(seur@meta.data)[grepl("Mancuso", colnames(seur@meta.data))], 
             group.by = "cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names)

p2 = DotPlot(seur, features = colnames(seur@meta.data)[grepl("Gazestani", colnames(seur@meta.data))], 
             group.by = "cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names)

pdf(file = paste0(out_dir,script_ind, "Cell_markers_module_score_dotplot_Mancuso24.pdf"), 
    width = 6, height = 8)
plot(p1)
dev.off()

pdf(file = paste0(out_dir,script_ind, "Cell_markers_module_score_dotplot_Gazestani23.pdf"), 
    width = 6, height = 8)
plot(p2)
dev.off()



#############################################################################
# De novo cluster marker identification
#############################################################################

###identify markers

seur = seur0

t1 <- FindAllMarkers(seur,group.by = "cluster_name", only.pos = TRUE)
t1 = as_tibble(t1[t1$p_val_adj<0.001 & t1$avg_log2FC > 0.5, ])
t1$cluster = as.character(t1$cluster)
t2 = t1[order(t1$cluster, -t1$avg_log2FC),]
seur_markers = t2

write_csv(seur_markers, file = paste0(out_dir,script_ind, "Seurat_markers.csv"))


###plot top10 markers per cluster (max 200 cells/cluster)

top_markers = seur_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

#subsample seur
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
p1 = DoHeatmap(seur_plot, group.by = "cluster_name", slot = "data", 
               features = top_markers$gene) + NoLegend()

pdf(file = paste0(out_dir,script_ind, "Seurat_markers_top10_heatmap.pdf"), 
    width = 10, height = 10)
plot(p1)
dev.off()

p1 = DotPlot(seur, features = unique(top_markers$gene), 
             group.by = "cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names)

pdf(file = paste0(out_dir,script_ind, "Seurat_markers_top10_dotplot.pdf"), 
    width = length(cluster_names)*2.5+1, height = length(cluster_names)/3+1.5)
plot(p1)
dev.off()


### GO analysis

GO_list = list()
GO_results_tab = NULL

for (cl in cluster_names){
  
  message("\n          *** GO analysis DEGs vs cluster ", cl, " - ", Sys.time(),"\n")
  ego = NULL

  ego = enrichGO(gene         = seur_markers$gene[seur_markers$cluster == cl],
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
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

seur$group = factor(seur$group, levels = gr) #required to fix order of groups

sccomp_result = 
  seur |>
  sccomp_estimate( 
    formula_composition = ~ group, 
    .sample =  sample, 
    .cell_group = cluster_name, 
    bimodal_mean_variability_association = TRUE,
    cores = 8 
  ) |> 
  #sccomp_remove_outliers(cores = 8) |> # Optional
  sccomp_test()

write_csv(sccomp_result, paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group.csv"))


pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates.pdf"), 
    width = 8, height = 5)
sccomp_result |> 
  plot_1D_intervals()
dev.off()

### plot boxplot if any comparison significant difference (else creates error)
t1 = sccomp_result
t2 = t1[!is.na(t1$factor),]
t3 = t2[t2$c_FDR<0.05,]

if (nrow(t3)>0){
  pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot.pdf"), 
      width = 8, height = 8)
  sccomp_result |> 
    sccomp_boxplot(factor = "group")
  dev.off()
}



###########################################################
# plot cluster abundance vs plaque_dens by TREM2Variant
###########################################################

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
# sccomp differential cell cluster abundance analysis by TREM2Variant vs plaque_dens
###########################################################

seur = seur0

seur$TREM2Variant = factor(seur$TREM2Variant, levels = unique(gr_tab$TREM2Variant))
seur = seur[,!is.na(seur$plaque_dens)]

sccomp_result = 
  seur |>
  sccomp_estimate( 
    formula_composition = ~ TREM2Variant*plaque_dens, 
    .sample =  sample, 
    .cell_group = cluster_name, 
    bimodal_mean_variability_association = TRUE,
    cores = 8 
  ) |> 
  #sccomp_remove_outliers(cores = 8) |> # Optional
  sccomp_test()

write_csv(sccomp_result, paste0(out_dir,script_ind,"sccomp_cell_abundance_by_TREM2Variant_vs_plaque_dens.csv"))


pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates_TREM2Variant_vs_plaque_dens.pdf"), 
    width = 8, height = 5)
sccomp_result |> 
  plot_1D_intervals()
dev.off()

### plot boxplot if any comparison significant difference (else creates error)
t1 = sccomp_result
t2 = t1[!is.na(t1$factor),]
t3 = t2[t2$c_FDR<0.05,]

if (nrow(t3)>0){
  pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot_TREM2Variant_vs_plaque_dens.pdf"), 
      width = 8, height = 8)
  sccomp_result |> 
    sccomp_boxplot(factor = "TREM2Variant")
  dev.off()
}


###########################################################
# sccomp differential cell cluster abundance analysis by TREM2Variant vs plaque_dens corrected for APOE and CD33
###########################################################

seur = seur0

seur$TREM2Variant = factor(seur$TREM2Variant, levels = unique(gr_tab$TREM2Variant))
seur$APOEgroup = factor(seur$APOEgroup, levels = unique(gr_tab$APOEgroup))
seur$CD33Group = factor(seur$CD33Group, levels = unique(gr_tab$CD33Group))

seur = seur[,!is.na(seur$plaque_dens)&!is.na(seur$APOEgroup)&!is.na(seur$CD33Group)]

sccomp_result = 
  seur |>
  sccomp_estimate( 
    formula_composition = ~ APOEgroup + CD33Group + TREM2Variant*plaque_dens, 
    .sample =  sample, 
    .cell_group = cluster_name, 
    bimodal_mean_variability_association = TRUE,
    cores = 8 
  ) |> 
  #sccomp_remove_outliers(cores = 8) |> # Optional
  sccomp_test()

write_csv(sccomp_result, paste0(out_dir,script_ind,"sccomp_cell_abundance_by_TREM2Variant_vs_plaque_dens_corr_APOE_CD33.csv"))


pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates_TREM2Variant_vs_plaque_dens_corr_APOE_CD33.pdf"), 
    width = 8, height = 5)
sccomp_result |> 
  plot_1D_intervals()
dev.off()

### plot boxplot if any comparison significant difference (else creates error)
t1 = sccomp_result
t2 = t1[!is.na(t1$factor),]
t3 = t2[t2$c_FDR<0.05,]

if (nrow(t3)>0){
  pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot_TREM2Variant_vs_plaque_dens_corr_APOE_CD33.pdf"), 
      width = 8, height = 8)
  sccomp_result |> 
    sccomp_boxplot(factor = "TREM2Variant")
  dev.off()
}


#get info on version of R, used packages etc
sessionInfo()


message("\n\n##########################################################################\n",
        "# Completed C03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


