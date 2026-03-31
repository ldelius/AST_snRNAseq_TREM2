message("\n\n##########################################################################\n",
        "# Start LD_F02a1 Pseudobulk DESeq2 analysis: ", Sys.time(),
        "\n##########################################################################\n",
        "\n   LRT test by cluster with correction for covariates",
        "\n##########################################################################\n\n")

# what this script does:
# - runs DESeq2 LRT test by cluster with correction for covariates (cohort, BrainRegion, APOEgroup, CD33Group)
# - extracts DESeq2 results and DEGs
# - calculates vst norm expression matrix (based on all clusters combined), corrects for covariates and calculates gene z-scores
# - plots PCA of samples by cluster based on uncorrected and corrected vst matrix
# - plots heatmap of top 3000 variable genes and DEGs by cluster_sample based on vst matrix corrected for covariates  

# Open packages necessary for analysis.
library(qs)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Seurat)
library(DESeq2)
library(colorRamps)
library(viridis)
library(pheatmap)
library(ggrepel)
library(WGCNA)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify output directory
out_dir = paste0(main_dir,"LD_F_DESeq_pseudobulk_WGCNA/")

#specify script/output index as prefix for file names
script_ind = "LD_F02a1_"


### load dataset; remove Control samples

bulk_data = qread(file = paste0(out_dir, "LD_F01_bulk_data.qs"))

t1 = bulk_data$gr_tab

#remove Control samples
t1 = t1[t1$NeuropathologicalDiagnosis != "Control" &
          !is.na(t1$CD33Group) & !is.na(t1$APOEgroup) & !is.na(t1$BrainRegion),]

gr_tab = t1

bulk_data$gr_tab = gr_tab


### adapt metadata for modelling

t1 = bulk_data$meta
t1 = t1[t1$sample %in% gr_tab$sample,]

t1$cluster_name = factor(t1$cluster_name, levels = unique(t1$cluster_name))
t1$TREM2Variant = factor(t1$TREM2Variant, levels = unique(t1$TREM2Variant))
t1$APOEgroup = factor(t1$APOEgroup, levels = c("APOE4-neg", "APOE4-pos"))
t1$CD33Group = factor(t1$CD33Group, levels = c("CV","CD33var"))
t1$BrainRegion = factor(t1$BrainRegion, levels = c("MTG","SSC"))
t1$cohort = factor(t1$cohort, levels = unique(t1$cohort))

bulk_data$meta = t1

m1 = bulk_data$counts
m2 = m1[,bulk_data$meta$cluster_sample]
bulk_data$counts = m2


### define formula for extracting differential genes to use for WGCNA 
form_full = "~cohort + BrainRegion + APOEgroup  + CD33Group + cluster_name  + TREM2Variant"
form_red = "~cohort + BrainRegion + APOEgroup  + CD33Group  + cluster_name"

### define covariates to correct vst matrix for
covars_corr = c("cohort", "BrainRegion", "APOEgroup", "CD33Group")


###get marker gene panels
GOI = list()
t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol




###########################################################
# functions
###########################################################

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


###function: plot bulk gene expression heatmap with annotations

bulkdata_heatmap = function(pl_mat, pl_meta, x_col, 
                            p_mat = FALSE,
                            pl_genes = NULL, 
                            meta_annot_cols = NULL,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            show_rownames = TRUE, show_colnames = TRUE,
                            color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                            lims = NULL,  cellwidth = 15, cellheight = 10, 
                            fontsize = 10, number_color = "grey70",
                            title = "Z-score vst-norm gene expression"){
  
  if (is.null(pl_genes)){pl_genes = rownames(pl_mat)}
  
  pl_mat = pl_mat[match(pl_genes, rownames(pl_mat), nomatch = 0),]
  pl_meta = pl_meta[match(pl_meta[[x_col]], colnames(pl_mat), nomatch = 0),]
  
  #define annotation bars
  
  if (!is.null(meta_annot_cols)){
    
    annot_col = data.frame(row.names = pl_meta[[x_col]]) #if not defined cbind converts factor values to factor levels
    
    for (col1 in meta_annot_cols){
      v1 = pl_meta[match(colnames(pl_mat), pl_meta[[x_col]]),][[col1]]
      if (is.numeric(v1)){v1 = factor(v1, levels = unique(pl_meta[[col1]][order(pl_meta[[col1]])]))} else {
        v1 = factor(v1, levels = unique(pl_meta[[col1]]))
      }
      annot_col = as.data.frame(cbind(annot_col, v1))
    }
    colnames(annot_col) = meta_annot_cols
    rownames(annot_col) = colnames(pl_mat)
    
    annot_colors = lapply(meta_annot_cols, function(x){
      v1 = pal(unique(annot_col[[x]]) )
      names(v1) = levels(annot_col[[x]])
      return(v1)
    })
    names(annot_colors) = meta_annot_cols
    
  }else{
    annot_col = NULL
    annot_colors = NULL
  }
  
  # define plot limits
  
  if (is.null(lims)){lims = c(-0.7*max(abs(na.omit(pl_mat))), 0.7*max(abs(na.omit(pl_mat))))}
  
  #create plot
  
  p1 = pheatmap::pheatmap(pl_mat, 
                          display_numbers = p_mat,
                          cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                          show_rownames = show_rownames, show_colnames = show_colnames,
                          color = color,
                          breaks = seq(lims[1], lims[2], length.out = length(color)+1),
                          annotation_col = annot_col, annotation_colors = annot_colors,
                          border_color = NA, cellwidth = cellwidth, cellheight = cellheight, 
                          fontsize = fontsize, number_color = number_color, fontsize_number = fontsize, 
                          main = title
  )
  
  return(p1)
  
}




###########################################################
# Preparing dataset for DESeq2 analysis
###########################################################

message("\n\n          *** Preparing dataset for DESeq2 analysis... ", Sys.time(),"\n\n")

t1 = bulk_data$meta

#remove clusters present in <4 samples
t2 = t1 %>% dplyr::group_by(cluster_name) %>% dplyr::summarise(N_samples = n())

comp_clusters = unique(t2$cluster_name[t2$N_samples>3])

comp_meta = t1[t1$cluster_name %in% comp_clusters,]

bulk_data$meta = comp_meta

#remove genes from countmatrix with <0.1 counts/cell in all pseudobulks  

m1 = bulk_data$counts

keep_genes = rownames(m1)[apply(m1, 1, max)>0.1]

comp_counts = m1[keep_genes, comp_meta$cluster_sample]

bulk_data$counts = comp_counts


#######################################
# Run DESeq2, calculate uncorrected and corrected vst matrix and z-scores
#######################################

message("\n\n          *** Running DESeq2 analysis groups and clusters combined - ", Sys.time(),"\n\n")

meta_cl = comp_meta
counts_cl = comp_counts


#run deseq2

dds = DESeqDataSetFromMatrix(counts_cl, colData = meta_cl, 
                             design = as.formula(form_full))
dds = DESeq(dds, test = "LRT", reduced = as.formula(form_red)) 

bulk_data$deseq_dataset_groups_clusters_combined = dds
bulk_data$deseq_formula_groups_clusters_combined = list(form_full = form_full,
                                                        form_red = form_red)

### extract vst norm expression matrix (based on all clusters combined), correct for covariates, calculate gene z-scores 

dds = bulk_data$deseq_dataset_groups_clusters_combined

vst_mat = assay(vst(dds))

bulk_data$vst_mat_uncorr = vst_mat

#batch-correct vst matrix
vst_mat_corr = vst_mat

for (cov1 in covars_corr){
  if (length(unique(comp_meta[[cov1]]))>1){
    vst_mat_corr = limma::removeBatchEffect(vst_mat_corr, batch = comp_meta[[cov1]], group = comp_meta$group)
  } else {vst_mat_corr = vst_mat_corr}
}

bulk_data$vst_mat = vst_mat_corr

#calculate Z-score per gene by pseudobulk (cluster_sample) for all clusters combined (from uncorrected and corrected vst matrix)
cluster_sample_mat = t(apply(vst_mat, 1, scale))
colnames(cluster_sample_mat) = colnames(vst_mat)
bulk_data$gene_Z_scores_uncorr[["clusters_combined"]] = cluster_sample_mat

cluster_sample_mat = t(apply(vst_mat_corr, 1, scale))
colnames(cluster_sample_mat) = colnames(vst_mat_corr)
bulk_data$gene_Z_scores[["clusters_combined"]] = cluster_sample_mat


### extract DESeq results and DEGs

t0 = as.data.frame(results(dds))
t1 = cbind(gene = rownames(t0), t0)
bulk_data$deseq_results[["clusters_combined"]] = t1
t2 = t1[!is.na(t1$pval) & t1$pval<0.05,]
bulk_data$DEGs[["clusters_combined"]] = t2$gene


#save bulk_dataset with DESeq results

qsave(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.rda")) 



#######################################
# plot PCA plot samples (by cluster, based on uncorrected vst matrix)
#######################################

meta = bulk_data$meta

samples = unique(bulk_data$gr_tab$sample)
gr = unique(bulk_data$gr_tab$group)

vst_mat = bulk_data$vst_mat_uncorr

#identify top 3000 variable genes
v1 = rowVars(vst_mat)
var_genes = names(v1[order(-v1)][1:3000])


### plot PCA by cluster

pl = list()

cl = "AST_SLC1A2_s4"  # changed from HOM_s16 for MIC

for (cl in comp_clusters){
  
  meta_cl = meta[meta$cluster_name == cl,]
  
  m1 = vst_mat[var_genes,meta_cl$cluster_sample]
  z_mat_cl = scale(m1)
  
  pc_analysis = prcomp(z_mat_cl)
  
  meta_pca = cbind(meta_cl, pc_analysis$rotation)
  meta_pca$sample_label = paste0(meta_pca$sample)
  
  pl[[paste0(cl, "_PC1_PC2")]] = ggplot(meta_pca)+geom_point(aes(x = PC1, y = PC2, color = group))+
    geom_text_repel(aes(x = PC1, y = PC2, label = sample_label, color = group))+
    scale_color_manual(limits = gr, values = pal(gr))+
    theme_minimal()+labs(title = paste0(cl, "- PC1 vs PC2 "))
  
}


### plot PCA all clusters

meta_cl = meta

m1 = vst_mat[var_genes,]
z_mat_cl = scale(m1)

pc_analysis = prcomp(z_mat_cl)

meta_pca = cbind(meta_cl, pc_analysis$rotation)
meta_pca$sample_label = paste0(meta_pca$cluster_sample)

pl[["all_clusters"]] = ggplot(meta_pca)+geom_point(aes(x = PC1, y = PC2, color = group))+
  geom_text_repel(aes(x = PC1, y = PC2, label = sample_label, color = group))+
  scale_color_manual(limits = gr, values = pal(gr))+
  theme_minimal()+labs(title = paste0(" all_clusters - PCA "))


pdf(file = paste0(out_dir,script_ind, "PCA_plot_samples_by_cluster_uncorr.pdf"), 
    width = 7, height = 6)
{
  lapply(pl, function(x){x})
}
dev.off()



#######################################
# plot PCA plot samples (by cluster, based on corrected vst matrix)
#######################################

meta = bulk_data$meta

samples = unique(bulk_data$gr_tab$sample)
gr = unique(bulk_data$gr_tab$group)

vst_mat = bulk_data$vst_mat

#identify top 3000 variable genes
v1 = rowVars(vst_mat)
var_genes = names(v1[order(-v1)][1:3000])


### plot PCA by cluster

pl = list()

for (cl in comp_clusters){
  
  meta_cl = meta[meta$cluster_name == cl,]
  
  m1 = vst_mat[var_genes,meta_cl$cluster_sample]
  z_mat_cl = scale(m1)
  
  pc_analysis = prcomp(z_mat_cl)
  
  meta_pca = cbind(meta_cl, pc_analysis$rotation)
  meta_pca$sample_label = paste0(meta_pca$sample)
  
  pl[[paste0(cl, "_PC1_PC2")]] = ggplot(meta_pca)+geom_point(aes(x = PC1, y = PC2, color = group))+
    geom_text_repel(aes(x = PC1, y = PC2, label = sample_label, color = group))+
    scale_color_manual(limits = gr, values = pal(gr))+
    theme_minimal()+labs(title = paste0(cl, "- PC1 vs PC2 "))
  
}


### plot PCA all clusters

meta_cl = meta

m1 = vst_mat[var_genes,]
z_mat_cl = scale(m1)

pc_analysis = prcomp(z_mat_cl)

meta_pca = cbind(meta_cl, pc_analysis$rotation)
meta_pca$sample_label = paste0(meta_pca$cluster_sample)

pl[["all_clusters"]] = ggplot(meta_pca)+geom_point(aes(x = PC1, y = PC2, color = group))+
  geom_text_repel(aes(x = PC1, y = PC2, label = sample_label, color = group))+
  scale_color_manual(limits = gr, values = pal(gr))+
  theme_minimal()+labs(title = paste0(" all_clusters - PCA "))


pdf(file = paste0(out_dir,script_ind, "PCA_plot_samples_by_cluster_corr.pdf"), 
    width = 7, height = 6)
{
  lapply(pl, function(x){x})
}
dev.off()


#################################################
# summarise all DEGs (all clusters combined)
#################################################

v1 = unique(unlist(bulk_data$DEGs))

t1 = tibble(N_DEGs_all = length(v1),
            N_DEGs_all_TF = length(intersect(v1, GOI$TF)),)

write_csv(t1, file = paste0(out_dir,script_ind, "DEGs_all_clusters_combined_N_genes.csv"))



#################################################
# plot heatmap top3000 variable genes by cluster_sample
#################################################

vst_mat = bulk_data$vst_mat

#identify top 3000 variable genes
v1 = rowVars(vst_mat)
pl_genes = names(v1[order(-v1)][1:3000])

pl_mat_X = bulk_data$gene_Z_scores$clusters_combined[pl_genes,]
lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

names(meta)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_var_genes_top3000.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = meta,
                        pl_genes = pl_genes,
                        x_col = "cluster_sample", 
                        meta_annot_cols = c("TREM2Variant","NeuropathologicalDiagnosis", 
                                            "cluster_name"),
                        show_rownames = FALSE, show_colnames = FALSE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 2, cellheight = 0.2, 
                        fontsize = 5, title = paste0("Top3000 variable genes - Z-score"))
  
}
dev.off()



#################################################
# plot heatmap DEGs combined by cluster_sample
#################################################

vst_mat = bulk_data$vst_mat

pl_genes = unique(unlist(bulk_data$DEGs))

pl_mat_X = bulk_data$gene_Z_scores$clusters_combined[pl_genes,]
lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

names(meta)

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_comb.pdf"), 
    width = 15, height = 30)
{
  p1 = bulkdata_heatmap(pl_mat = pl_mat_X, 
                        pl_meta = meta,
                        pl_genes = pl_genes,
                        x_col = "cluster_sample", 
                        meta_annot_cols = c("TREM2Variant","NeuropathologicalDiagnosis", 
                                            "cluster_name"),
                        show_rownames = FALSE, show_colnames = FALSE,
                        cluster_rows = TRUE, cluster_cols = FALSE,
                        color = viridis(250),
                        lims = lims_X,  cellwidth = 2, cellheight = 0.2, 
                        fontsize = 5, title = paste0("Combined DEGs - Z-score"))
  
}
dev.off()





#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed C01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


