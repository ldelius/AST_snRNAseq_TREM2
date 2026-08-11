message("\n\n##########################################################################\n",
        "# Start LD_E02c: Five-covariate pseudobulk DESeq2 analysis ", Sys.time(),
        "\n##########################################################################\n\n")

# Adds cohort and Sex to both the E02a2 VST correction and DESeq2 design to
# test robustness to expanded adjustment.
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
library(matrixStats)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "LD_E02c_v01_"

#specify output directory (input pseudobulk object lives here too)
out_dir = paste0(main_dir,"LD_E_DESeq_pseudobulk/")

#dedicated sub-folder for all LD_E02c outputs
sub_out_dir = paste0(out_dir, "LD_E02c/")
dir.create(sub_out_dir, showWarnings = FALSE, recursive = TRUE)



### load dataset; define variables required for DESeq2 models
# keep clusters with >= 2 samples per level and at least 2 levels for each categorical variable

bulk_data = qread(file = paste0(out_dir, "LD_E01_v02_bulk_data.qs"))

t1 = bulk_data$meta
names(t1)

#define categorical model variables (used for NA removal AND the per-cluster
# level-based filter below).
model_vars_cat = c("cohort", "APOEgroup", "CD33Group", "BrainRegion", "Sex",
                   "NeuropathologicalDiagnosis", "TREM2Variant")

t1$NeuropathologicalDiagnosis = factor(t1$NeuropathologicalDiagnosis, levels = c("Control", "AD"))

t1$APOEgroup   = factor(t1$APOEgroup, levels = c("APOE4-neg", "APOE4-pos"))
t1$CD33Group   = factor(t1$CD33Group, levels = c("CV", "CD33var"))
t1$BrainRegion = factor(t1$BrainRegion, levels = c("SSC", "MTG"))
t1$cohort      = factor(t1$cohort, levels = c("BiogenInitial", "BiogenExtra"))

# NEW covariate: Sex (categorical factor)
t1$Sex = factor(t1$Sex)

# define covariates to correct the vst matrix for
#   - categorical covariates go through removeBatchEffect(batch = ...)
covars_corr_cat  = c("cohort", "APOEgroup", "CD33Group", "BrainRegion", "Sex")

# define reference cluster for between cluster "response" comparisons
ref_cluster = "AST_SLC1A2_s0" # changed from HOM_s0 for MIC

#remove samples with missing covariate values
for (var1 in model_vars_cat){
  t1 = t1[!is.na(t1[[var1]]),]
}

#filter for clusters containing at least 2 samples for 2 different levels for all
# analysed CATEGORICAL variables -->  DESeq2 can't estimate effects without variation

for (var1 in model_vars_cat){
  t1$model_var = t1[[var1]]
  t2 = t1 %>% group_by(cluster_name, model_var) %>% summarise(N_samples = n())
  t3 = t2[t2$N_samples>=2,] %>% group_by(cluster_name) %>% summarise(N_levels = n())
  t1 = t1[t1$cluster_name %in% t3$cluster_name[t3$N_levels>1],]
}

bulk_data$meta = t1


### update bulkdata object, extract group tab
m1 = bulk_data$counts
m2 = m1[,bulk_data$meta$cluster_sample]
bulk_data$counts = m2

gr_tab = bulk_data$gr_tab[bulk_data$gr_tab$sample %in% bulk_data$meta$sample,]
bulk_data$gr_tab = gr_tab



###get marker gene panels
# reads a list of genes known to encode transcriton factors, to be used for plotting later on.
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

comp_meta = bulk_data$meta

comp_clusters = unique(comp_meta$cluster_name)

#remove genes from countmatrix with <0.1 counts/cell in all pseudobulks  --> fikters out very lowly expressed genes

m1 = bulk_data$counts

keep_genes = rownames(m1)[apply(m1, 1, max)>0.1]

comp_counts = m1[keep_genes,]

bulk_data$counts = comp_counts



#######################################
# calculate uncorrected and corrected vst matrix and corresponding gene z-scores
#######################################
# vst = variance stabilising transformation to make sure that data is roughly homoscedastic and variance is similar across low and high expression genes
### extract vst-normalised matrix for combined dataset

dds = DESeqDataSetFromMatrix(comp_counts, colData = comp_meta,
                             design = ~group)

bulk_data$deseq_dataset_groups_clusters_combined = dds

# extract vst norm expression matrix (based on all clusters combined), correct for covariates, calculate gene z-scores

dds = bulk_data$deseq_dataset_groups_clusters_combined

vst_mat = assay(vst(dds))

bulk_data$vst_mat_uncorr = vst_mat

#batch-correct vst matrix for the 5 covariates
#  - categorical covariates (cohort, APOE, CD33, BrainRegion, Sex) via batch =
vst_mat_corr = vst_mat

for (cov1 in covars_corr_cat){
  if (length(unique(comp_meta[[cov1]]))>1){
    vst_mat_corr = limma::removeBatchEffect(vst_mat_corr, batch = comp_meta[[cov1]])
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



#######################################
# plot PCA plot samples (by cluster, based on uncorrected vst matrix)
#######################################
# For each cluster individually and for all clusters combined: Selects top 3000 most variable genes (these drive the most biological variation),
# performs PCA on the vst-normalised expression of these genes, and plots the first 2 principal components with samples colored by group (TREM2 variant and diagnosis groups)
# and labelled by sample ID. This is done separately for the uncorrected vst matrix and the covariate-corrected vst matrix, to assess how well samples cluster by group and whether there are any outliers or batch effects.
meta = bulk_data$meta

samples = unique(bulk_data$gr_tab$sample)
gr = unique(bulk_data$gr_tab$group)

vst_mat = bulk_data$vst_mat_uncorr

#identify top 3000 variable genes
v1 = rowVars(vst_mat)
var_genes = names(v1[order(-v1)][1:3000])


### plot PCA by cluster

pl = list()

cl = "AST_SLC1A2_s4" # changed from HOM...for MIC

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


pdf(file = paste0(sub_out_dir,script_ind, "PCA_plot_samples_by_cluster_uncorr.pdf"),
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


pdf(file = paste0(sub_out_dir,script_ind, "PCA_plot_samples_by_cluster_corr.pdf"),
    width = 7, height = 6)
{
  lapply(pl, function(x){x})
}
dev.off()



#######################################
# Run DESeq2 analysis cluster X vs ref_cluster (common variant only)
#######################################
# compares pseudobulk expression of that cluster against AST_SLC1A2_s0, using only CV samples.
# identifies genes that are differentially expressed between astrocyte subtypes in the absence of TREM2 mutations (regardles of AD vs ctrl)
# NOTE: this baseline comparison does NOT use covariates (design = ~cluster_name),
#       so it is unaffected by the 5-covariate change - kept identical to LD_E02a2.

message("\n\n          *** Running DESeq2 analysis  - ", Sys.time(),"\n\n")

bulk_data$deseq_dataset = list()
bulk_data$deseq_info = list()

### analysis cluster X vs ref_cluster (common variant only)

for (cl in comp_clusters[comp_clusters != ref_cluster]){

  message("\n          *** Running DESeq2 analysis ", cl, " vs ref_cluster - ", Sys.time(),"\n")

  meta_cl = comp_meta[(comp_meta$cluster_name %in% c(cl, ref_cluster) ) & comp_meta$TREM2Variant == "CV",]
  meta_cl$cluster_name = factor(meta_cl$cluster_name, levels = c(ref_cluster, cl))

  counts_cl = comp_counts[, meta_cl$cluster_sample]

  form_full = "~cluster_name"
  form_red = "~1"

  bulk_data$deseq_info[[paste0(cl, "_CV_vs_", ref_cluster)]] = list(form_full = form_full,
                                                                    form_red = form_red)

  #run deseq2
  names(meta_cl)
  dds = DESeqDataSetFromMatrix(counts_cl, colData = meta_cl,
                               design = as.formula(form_full))
  dds = DESeq(dds, test = "LRT", reduced = as.formula(form_red))

  bulk_data$deseq_dataset[[paste0(cl, "_CV_vs_", ref_cluster)]] = dds

}



#######################################
# Run Deseq analysis within clusters CV AD vs Control
#######################################
# for each cluster it tests AD vs ctrl in CV individuals, correcting for the 5 covariates
# (cohort, APOE, CD33, BrainRegion, Sex).
# Q: Which genes change in AD vs Control within each astrocyte cluster, independent of the 5 covariates
cl = "AST_SLC1A2_s3" # changed from HOM_s9 for MIC

for (cl in comp_clusters){

  message("\n          *** Running DESeq2 analysis  CV AD vs Control by cluster - ", cl, " - ", Sys.time(),"\n")

  msg = NULL

  meta_cl = comp_meta[comp_meta$cluster_name == cl & comp_meta$TREM2Variant == "CV",]
  counts_cl = comp_counts[, meta_cl$cluster_sample]

  form_full = "~cohort + APOEgroup + CD33Group + BrainRegion + Sex + NeuropathologicalDiagnosis"
  form_red  = "~cohort + APOEgroup + CD33Group + BrainRegion + Sex"

  #remove spaces from formula to avoid problems with automatic removal of variables
  form_full = str_remove_all(form_full, " ")
  form_red = str_remove_all(form_red, " ")

  #extract variables in model
  model_vars_full = unlist(str_split(str_remove_all(form_full, "~"), "\\+"))
  model_vars_red = unlist(str_split(str_remove_all(form_red, "~"), "\\+"))


  ###check whether variables allow calculating model, remove colinear covariates

  #check whether model variables have multiple levels, if not, drop from model
  for (var1 in model_vars_red){

    meta_cl$model_var = meta_cl[[var1]]

    t2 = meta_cl %>% group_by(model_var) %>% summarise(N_samples = n())

    if (nrow(t2)==1){
      form_full = str_remove_all(form_full, paste0(var1, "\\+"))
      form_red = str_remove_all(form_red, paste0(var1, "\\+"))
      form_red = str_remove_all(form_red, paste0(var1))
      model_vars_red =  model_vars_red[model_vars_red != var1]

      msg = c(msg, paste0("            ! WARNING: Dropping '", var1, "' - All samples have same value \n"))
      message(msg)
    }
  }

  #check whether any reduced model variables are colinear, if so remove first variable
  v1 = model_vars_red

  for (var1 in v1){
    for (var2 in v1[v1 != var1]){

      t2 = meta_cl %>% group_by(pick(all_of(c(var1, var2)))) %>% summarise(N_samples = n())

      if (nrow(t2)<3 & str_detect(form_full, var2)){
        form_full = str_remove_all(form_full, paste0(var1, "\\+"))
        form_red = str_remove_all(form_red, paste0(var1, "\\+"))
        form_red = str_remove_all(form_red, paste0(var1))
        model_vars_red = model_vars_red[model_vars_red != var1]
        v2 = paste0("         WARNING: Colinearity '", var1, "' and '", var2,"' => '", var1, "' dropped")
        message(v2)
        msg = c(msg, v2)
      }

      v1 = v1[v1 != var1]

    }
  }

  #check if covariates are colinear with main predictors of full vs reduced model => if so, omit cluster

  omit_cl = FALSE

  for (var1 in model_vars_full[!(model_vars_full %in% model_vars_red)]){
    for (var2 in model_vars_red){

      t2 = meta_cl %>% group_by(pick(all_of(c(var1, var2)))) %>% summarise(N_samples = n())

      if (nrow(t2)<3){
        omit_cl = TRUE
        v2 = paste0("         WARNING: Colinearity '", var1, "' and '", var2,
                    "' => effect '", var1, "' cannot be determined => cluster ", cl, " omitted")
        message(v2)
        msg = c(msg, v2)
      }
    }
  }

  #if less samples than coefficients (variables+1) omit cluster (model can't be calculated)
  if (nrow(meta_cl) <= length(model_vars_full)+1){
    omit_cl = TRUE
    v1 = paste0("         WARNING: Not enough samples to calculate model => cluster ", cl," ommitted")
    message(v1)
    msg = c(msg, v1)
  }

  #full-rank safeguard: the collinearity checks above only catch single variables
  # with one level or PAIRS of collinear variables. They miss higher-order linear
  # dependencies where 3+ covariates are jointly redundant (e.g. cohort is a linear
  # combination of APOEgroup + BrainRegion + Sex) even though no single pair is
  # collinear. In that case DESeqDataSetFromMatrix() aborts the ENTIRE run with
  # "model matrix is not full rank". Test the design-matrix rank explicitly here and,
  # if rank-deficient, omit the cluster (consistent with the other omit_cl cases).
  if (omit_cl == FALSE){
    form_full = str_remove(form_full, "\\+$")
    form_red  = str_remove(form_red, "\\+$")
    mm_full = model.matrix(as.formula(form_full), data = meta_cl)
    if (qr(mm_full)$rank < ncol(mm_full)){
      omit_cl = TRUE
      v1 = paste0("         WARNING: design matrix not full rank (higher-order covariate collinearity) => cluster ", cl, " omitted")
      message(v1)
      msg = c(msg, v1)
    }
  }

  #if more samples than coefficients (variables+1) run DESeq2 (else model can't be calculated)
  if (omit_cl == FALSE){
    form_full = str_remove(form_full, "\\+$") #remove "+" in end if final variable has been dropped
    form_red = str_remove(form_red, "\\+$")
    dds = DESeqDataSetFromMatrix(counts_cl, colData = meta_cl,
                                 design = as.formula(form_full))
    dds = DESeq(dds, test = "LRT", reduced = as.formula(form_red))

    bulk_data$deseq_dataset[[paste0(cl, "_TREM2_CV_AD_vs_Control")]] = dds
  }

  message(msg)

  bulk_data$deseq_info[[paste0(cl, "_TREM2_CV_AD_vs_Control")]] = list(form_full = form_full,
                                                                    form_red = form_red,
                                                                    msg = msg)
}



#######################################
# Run Deseq analysis within clusters TREM2 R62H vs CV (only AD cases)
#######################################
# for each cluster it tests TREM2 R62H vs CV in AD individuals, correcting for the 5 covariates
# (cohort, APOE, CD33, BrainRegion, Sex).
cl = "AST_GFAP_s2" # changed from IRM_s14 for MIC

for (cl in comp_clusters){

  message("\n          *** Running DESeq2 analysis AD R62H vs CV by cluster - ", cl, " - ", Sys.time(),"\n")

  msg = NULL

  meta_cl = comp_meta[comp_meta$cluster_name == cl &
                        comp_meta$NeuropathologicalDiagnosis == "AD"&
                        comp_meta$TREM2Variant %in% c("CV", "R62H"),]
  meta_cl$TREM2Variant = factor(meta_cl$TREM2Variant, levels = c("CV", "R62H"))
  counts_cl = comp_counts[, meta_cl$cluster_sample]

  form_full = "~cohort + APOEgroup + CD33Group + BrainRegion + Sex + TREM2Variant"
  form_red  = "~cohort + APOEgroup + CD33Group + BrainRegion + Sex"

  #remove spaces from formula to avoid problems with automatic removal of variables
  form_full = str_remove_all(form_full, " ")
  form_red = str_remove_all(form_red, " ")

  #extract variables in model
  model_vars_full = unlist(str_split(str_remove_all(form_full, "~"), "\\+"))
  model_vars_red = unlist(str_split(str_remove_all(form_red, "~"), "\\+"))


  ###check whether variables allow calculating model, remove colinear covariates

  #check whether model variables have multiple levels, if not, drop from model
  for (var1 in model_vars_red){

    meta_cl$model_var = meta_cl[[var1]]

    t2 = meta_cl %>% group_by(model_var) %>% summarise(N_samples = n())

    if (nrow(t2)==1){
      form_full = str_remove_all(form_full, paste0(var1, "\\+"))
      form_red = str_remove_all(form_red, paste0(var1, "\\+"))
      form_red = str_remove_all(form_red, paste0(var1))
      model_vars_red =  model_vars_red[model_vars_red != var1]

      msg = c(msg, paste0("            ! WARNING: Dropping '", var1, "' - All samples have same value \n"))
      message(msg)
    }
  }

  #check whether any reduced model variables are colinear, if so remove first variable
  v1 = model_vars_red

  for (var1 in v1){
    for (var2 in v1[v1 != var1]){

      t2 = meta_cl %>% group_by(pick(all_of(c(var1, var2)))) %>% summarise(N_samples = n())

      if (nrow(t2)<3 & str_detect(form_full, var2)){
        form_full = str_remove_all(form_full, paste0(var1, "\\+"))
        form_red = str_remove_all(form_red, paste0(var1, "\\+"))
        form_red = str_remove_all(form_red, paste0(var1))
        model_vars_red = model_vars_red[model_vars_red != var1]
        v2 = paste0("         WARNING: Colinearity '", var1, "' and '", var2,"' => '", var1, "' dropped")
        message(v2)
        msg = c(msg, v2)
      }

      v1 = v1[v1 != var1]

    }
  }

  #check if covariates are colinear with main predictors of full vs reduced model => if so, omit cluster

  omit_cl = FALSE

  for (var1 in model_vars_full[!(model_vars_full %in% model_vars_red)]){
    for (var2 in model_vars_red){

      t2 = meta_cl %>% group_by(pick(all_of(c(var1, var2)))) %>% summarise(N_samples = n())

      if (nrow(t2)<3){
        omit_cl = TRUE
        v2 = paste0("         WARNING: Colinearity '", var1, "' and '", var2,
                    "' => effect '", var1, "' cannot be determined => cluster ", cl, " omitted")
        message(v2)
        msg = c(msg, v2)
      }
    }
  }

  #if less samples than coefficients (variables+1) omit cluster (model can't be calculated)
  if (nrow(meta_cl) <= length(model_vars_full)+1){
    omit_cl = TRUE
    v1 = paste0("         WARNING: Not enough samples to calculate model => cluster ", cl," ommitted")
    message(v1)
    msg = c(msg, v1)
  }

  #full-rank safeguard: the collinearity checks above only catch single variables
  # with one level or PAIRS of collinear variables. They miss higher-order linear
  # dependencies where 3+ covariates are jointly redundant (e.g. cohort is a linear
  # combination of APOEgroup + BrainRegion + Sex) even though no single pair is
  # collinear. In that case DESeqDataSetFromMatrix() aborts the ENTIRE run with
  # "model matrix is not full rank". Test the design-matrix rank explicitly here and,
  # if rank-deficient, omit the cluster (consistent with the other omit_cl cases).
  if (omit_cl == FALSE){
    form_full = str_remove(form_full, "\\+$")
    form_red  = str_remove(form_red, "\\+$")
    mm_full = model.matrix(as.formula(form_full), data = meta_cl)
    if (qr(mm_full)$rank < ncol(mm_full)){
      omit_cl = TRUE
      v1 = paste0("         WARNING: design matrix not full rank (higher-order covariate collinearity) => cluster ", cl, " omitted")
      message(v1)
      msg = c(msg, v1)
    }
  }

  #if more samples than coefficients (variables+1) run DESeq2 (else model can't be calculated)
  if (omit_cl == FALSE){
    form_full = str_remove(form_full, "\\+$") #remove "+" in end if final variable has been dropped
    form_red = str_remove(form_red, "\\+$")
    dds = DESeqDataSetFromMatrix(counts_cl, colData = meta_cl,
                                  design = as.formula(form_full))
    dds = DESeq(dds, test = "LRT", reduced = as.formula(form_red))

    bulk_data$deseq_dataset[[paste0(cl, "_AD_TREM2_R62H_vs_CV")]] = dds
  }

  message(msg)

  bulk_data$deseq_info[[paste0(cl, "_AD_TREM2_R62H_vs_CV")]] = list(form_full = form_full,
                                                                         form_red = form_red,
                                                                         msg = msg)

}



#######################################
# Run Deseq analysis within clusters TREM2 R47H vs CV (only AD cases)
#######################################
# for each cluster it tests TREM2 R47H vs CV in AD individuals, correcting for the 5 covariates
# (cohort, APOE, CD33, BrainRegion, Sex).
cl = "AST_GFAP_s2" # changend from IRM_s14 for MIC

for (cl in comp_clusters){

  message("\n          *** Running DESeq2 analysis AD R47H vs CV by cluster - ", cl, " - ", Sys.time(),"\n")

  msg = NULL

  meta_cl = comp_meta[comp_meta$cluster_name == cl &
                        comp_meta$NeuropathologicalDiagnosis == "AD"&
                        comp_meta$TREM2Variant %in% c("CV", "R47H"),]
  meta_cl$TREM2Variant = factor(meta_cl$TREM2Variant, levels = c("CV", "R47H"))
  counts_cl = comp_counts[, meta_cl$cluster_sample]

  form_full = "~cohort + APOEgroup + CD33Group + BrainRegion + Sex + TREM2Variant"
  form_red  = "~cohort + APOEgroup + CD33Group + BrainRegion + Sex"

  #remove spaces from formula to avoid problems with automatic removal of variables
  form_full = str_remove_all(form_full, " ")
  form_red = str_remove_all(form_red, " ")

  #extract variables in model
  model_vars_full = unlist(str_split(str_remove_all(form_full, "~"), "\\+"))
  model_vars_red = unlist(str_split(str_remove_all(form_red, "~"), "\\+"))


  ###check whether variables allow calculating model, remove colinear covariates

  #check whether model variables have multiple levels, if not, drop from model
  for (var1 in model_vars_red){

    meta_cl$model_var = meta_cl[[var1]]

    t2 = meta_cl %>% group_by(model_var) %>% summarise(N_samples = n())

    if (nrow(t2)==1){
      form_full = str_remove_all(form_full, paste0(var1, "\\+"))
      form_red = str_remove_all(form_red, paste0(var1, "\\+"))
      form_red = str_remove_all(form_red, paste0(var1))
      model_vars_red =  model_vars_red[model_vars_red != var1]

      msg = c(msg, paste0("            ! WARNING: Dropping '", var1, "' - All samples have same value \n"))
      message(msg)
    }
  }

  #check whether any reduced model variables are colinear, if so remove first variable
  v1 = model_vars_red

  for (var1 in v1){
    for (var2 in v1[v1 != var1]){

      t2 = meta_cl %>% group_by(pick(all_of(c(var1, var2)))) %>% summarise(N_samples = n())

      if (nrow(t2)<3 & str_detect(form_full, var2)){
        form_full = str_remove_all(form_full, paste0(var1, "\\+"))
        form_red = str_remove_all(form_red, paste0(var1, "\\+"))
        form_red = str_remove_all(form_red, paste0(var1))
        model_vars_red = model_vars_red[model_vars_red != var1]
        v2 = paste0("         WARNING: Colinearity '", var1, "' and '", var2,"' => '", var1, "' dropped")
        message(v2)
        msg = c(msg, v2)
      }

      v1 = v1[v1 != var1]

    }
  }

  #check if covariates are colinear with main predictors of full vs reduced model => if so, omit cluster

  omit_cl = FALSE

  for (var1 in model_vars_full[!(model_vars_full %in% model_vars_red)]){
    for (var2 in model_vars_red){

      t2 = meta_cl %>% group_by(pick(all_of(c(var1, var2)))) %>% summarise(N_samples = n())

      if (nrow(t2)<3){
        omit_cl = TRUE
        v2 = paste0("         WARNING: Colinearity '", var1, "' and '", var2,
                    "' => effect '", var1, "' cannot be determined => cluster ", cl, " omitted")
        message(v2)
        msg = c(msg, v2)
      }
    }
  }

  #if less samples than coefficients (variables+1) omit cluster (model can't be calculated)
  if (nrow(meta_cl) <= length(model_vars_full)+1){
    omit_cl = TRUE
    v1 = paste0("         WARNING: Not enough samples to calculate model => cluster ", cl," ommitted")
    message(v1)
    msg = c(msg, v1)
  }

  #full-rank safeguard: the collinearity checks above only catch single variables
  # with one level or PAIRS of collinear variables. They miss higher-order linear
  # dependencies where 3+ covariates are jointly redundant (e.g. cohort is a linear
  # combination of APOEgroup + BrainRegion + Sex) even though no single pair is
  # collinear. In that case DESeqDataSetFromMatrix() aborts the ENTIRE run with
  # "model matrix is not full rank". Test the design-matrix rank explicitly here and,
  # if rank-deficient, omit the cluster (consistent with the other omit_cl cases).
  if (omit_cl == FALSE){
    form_full = str_remove(form_full, "\\+$")
    form_red  = str_remove(form_red, "\\+$")
    mm_full = model.matrix(as.formula(form_full), data = meta_cl)
    if (qr(mm_full)$rank < ncol(mm_full)){
      omit_cl = TRUE
      v1 = paste0("         WARNING: design matrix not full rank (higher-order covariate collinearity) => cluster ", cl, " omitted")
      message(v1)
      msg = c(msg, v1)
    }
  }

  #if more samples than coefficients (variables+1) run DESeq2 (else model can't be calculated)
  if (omit_cl == FALSE){
    form_full = str_remove(form_full, "\\+$") #remove "+" in end if final variable has been dropped
    form_red = str_remove(form_red, "\\+$")
    dds = DESeqDataSetFromMatrix(counts_cl, colData = meta_cl,
                                  design = as.formula(form_full))
    dds = DESeq(dds, test = "LRT", reduced = as.formula(form_red))

    bulk_data$deseq_dataset[[paste0(cl, "_AD_TREM2_R47H_vs_CV")]] = dds
  }

  message(msg)

  bulk_data$deseq_info[[paste0(cl, "_AD_TREM2_R47H_vs_CV")]] = list(form_full = form_full,
                                                                         form_red = form_red,
                                                                         msg = msg)

}




#######################################
# Run Deseq analysis within clusters TREM2 R47H vs R62H (only AD cases)
#######################################
# for each cluster it tests TREM2 R47H vs R62H in AD individuals, correcting for the 5 covariates
# (cohort, APOE, CD33, BrainRegion, Sex).
cl = "AST_GFAP_s2" # changed from IRM_s14 for MIC

for (cl in comp_clusters){

  message("\n          *** Running DESeq2 analysis AD R47H vs R62H by cluster - ", cl, " - ", Sys.time(),"\n")

  msg = NULL

  meta_cl = comp_meta[comp_meta$cluster_name == cl &
                        comp_meta$NeuropathologicalDiagnosis == "AD"&
                        comp_meta$TREM2Variant %in% c("R62H", "R47H"),]
  meta_cl$TREM2Variant = factor(meta_cl$TREM2Variant, levels = c("R62H", "R47H"))
  counts_cl = comp_counts[, meta_cl$cluster_sample]

  form_full = "~cohort + APOEgroup + CD33Group + BrainRegion + Sex + TREM2Variant"
  form_red  = "~cohort + APOEgroup + CD33Group + BrainRegion + Sex"

  #remove spaces from formula to avoid problems with automatic removal of variables
  form_full = str_remove_all(form_full, " ")
  form_red = str_remove_all(form_red, " ")

  #extract variables in model
  model_vars_full = unlist(str_split(str_remove_all(form_full, "~"), "\\+"))
  model_vars_red = unlist(str_split(str_remove_all(form_red, "~"), "\\+"))


  ###check whether variables allow calculating model, remove colinear covariates

  #check whether model variables have multiple levels, if not, drop from model
  for (var1 in model_vars_red){

    meta_cl$model_var = meta_cl[[var1]]

    t2 = meta_cl %>% group_by(model_var) %>% summarise(N_samples = n())

    if (nrow(t2)==1){
      form_full = str_remove_all(form_full, paste0(var1, "\\+"))
      form_red = str_remove_all(form_red, paste0(var1, "\\+"))
      form_red = str_remove_all(form_red, paste0(var1))
      model_vars_red =  model_vars_red[model_vars_red != var1]

      msg = c(msg, paste0("            ! WARNING: Dropping '", var1, "' - All samples have same value \n"))
      message(msg)
    }
  }

  #check whether any reduced model variables are colinear, if so remove first variable
  v1 = model_vars_red

  for (var1 in v1){
    for (var2 in v1[v1 != var1]){

      t2 = meta_cl %>% group_by(pick(all_of(c(var1, var2)))) %>% summarise(N_samples = n())

      if (nrow(t2)<3 & str_detect(form_full, var2)){
        form_full = str_remove_all(form_full, paste0(var1, "\\+"))
        form_red = str_remove_all(form_red, paste0(var1, "\\+"))
        form_red = str_remove_all(form_red, paste0(var1))
        model_vars_red = model_vars_red[model_vars_red != var1]
        v2 = paste0("         WARNING: Colinearity '", var1, "' and '", var2,"' => '", var1, "' dropped")
        message(v2)
        msg = c(msg, v2)
      }

      v1 = v1[v1 != var1]

    }
  }

  #check if covariates are colinear with main predictors of full vs reduced model => if so, omit cluster

  omit_cl = FALSE

  for (var1 in model_vars_full[!(model_vars_full %in% model_vars_red)]){
    for (var2 in model_vars_red){

      t2 = meta_cl %>% group_by(pick(all_of(c(var1, var2)))) %>% summarise(N_samples = n())

      if (nrow(t2)<3){
        omit_cl = TRUE
        v2 = paste0("         WARNING: Colinearity '", var1, "' and '", var2,
                    "' => effect '", var1, "' cannot be determined => cluster ", cl, " omitted")
        message(v2)
        msg = c(msg, v2)
      }
    }
  }

  #if less samples than coefficients (variables+1) omit cluster (model can't be calculated)
  if (nrow(meta_cl) <= length(model_vars_full)+1){
    omit_cl = TRUE
    v1 = paste0("         WARNING: Not enough samples to calculate model => cluster ", cl," ommitted")
    message(v1)
    msg = c(msg, v1)
  }

  #full-rank safeguard: the collinearity checks above only catch single variables
  # with one level or PAIRS of collinear variables. They miss higher-order linear
  # dependencies where 3+ covariates are jointly redundant (e.g. cohort is a linear
  # combination of APOEgroup + BrainRegion + Sex) even though no single pair is
  # collinear. In that case DESeqDataSetFromMatrix() aborts the ENTIRE run with
  # "model matrix is not full rank". Test the design-matrix rank explicitly here and,
  # if rank-deficient, omit the cluster (consistent with the other omit_cl cases).
  if (omit_cl == FALSE){
    form_full = str_remove(form_full, "\\+$")
    form_red  = str_remove(form_red, "\\+$")
    mm_full = model.matrix(as.formula(form_full), data = meta_cl)
    if (qr(mm_full)$rank < ncol(mm_full)){
      omit_cl = TRUE
      v1 = paste0("         WARNING: design matrix not full rank (higher-order covariate collinearity) => cluster ", cl, " omitted")
      message(v1)
      msg = c(msg, v1)
    }
  }

  #if more samples than coefficients (variables+1) run DESeq2 (else model can't be calculated)
  if (omit_cl == FALSE){
    form_full = str_remove(form_full, "\\+$") #remove "+" in end if final variable has been dropped
    form_red = str_remove(form_red, "\\+$")
    dds = DESeqDataSetFromMatrix(counts_cl, colData = meta_cl,
                                 design = as.formula(form_full))
    dds = DESeq(dds, test = "LRT", reduced = as.formula(form_red))

    bulk_data$deseq_dataset[[paste0(cl, "_AD_TREM2_R47H_vs_R62H")]] = dds
  }

  message(msg)

  bulk_data$deseq_info[[paste0(cl, "_AD_TREM2_R47H_vs_R62H")]] = list(form_full = form_full,
                                                                    form_red = form_red,
                                                                    msg = msg)

}




#####################################
# extract DESeq results and DEGs
#####################################
# loose threshold kept identical to LD_E02a2: padj < 0.1, no log2FC cutoff

bulk_data$deseq_results = list()

for (comp1 in names(bulk_data$deseq_dataset)){

  dds = bulk_data$deseq_dataset[[comp1]]

  t0 = as.data.frame(results(dds))
  t1 = cbind(gene = rownames(t0), t0)
  bulk_data$deseq_results[[comp1]] = t1

  t2 = t1[!is.na(t1$padj) & t1$padj<0.1,]
  bulk_data$DEGs[[paste0(comp1, "_up")]] = t2$gene[t2$log2FoldChange>0]
  bulk_data$DEGs[[paste0(comp1, "_down")]] = t2$gene[t2$log2FoldChange<0]

}


#save bulk_dataset with DESeq results

qsave(bulk_data, file = paste0(sub_out_dir,script_ind, "bulk_data.qs"))




#################################################
# summarise all DEGs (all clusters combined and by cluster)
#################################################

v1 = unique(unlist(bulk_data$DEGs))

t1 = tibble(N_DEGs_all = length(v1),
            N_DEGs_all_TF = length(intersect(v1, GOI$TF)),)

write_csv(t1, file = paste0(sub_out_dir,script_ind, "DEGs_all_clusters_combined_N_genes.csv"))



### save table with DEGs by cluster and TFs

l1 = bulk_data$DEGs
l2 = lapply(l1, function(x){x = x[x%in% GOI$TF]})
names(l2) = paste0(names(l1), "_TF")

l3 = c(l1, l2)

m1 = matrix(nrow = max(lengths(l3)), ncol = length(l3))
colnames(m1) = names(l3)

for (i in names(l3)){
  v1 = l3[[i]]
  if(length(v1)>0){
    m1[1:length(v1),i] = v1
  }
}
m1[is.na(m1)] = ""

write_csv(as_tibble(m1), file = paste0(sub_out_dir,script_ind, "DEGs_by_cluster_genes.csv"))

t1 = tibble(gene_set = names(l3), N_genes = lengths(l3))

write_csv(t1, file = paste0(sub_out_dir,script_ind, "DEGs_by_cluster_N.csv"))



###########################################################
# Define the 5 comparison tags (shared by the bar chart, the new summary
# table and the new 2-panel plot below)
###########################################################
# DEG list names follow the pattern: {cluster}_{comparison_tag}_{direction}
# e.g. "AST_SLC1A2_s4_AD_TREM2_R47H_vs_CV_up"

comp_tags = c(
  "CV_vs_AST_SLC1A2_s0"      = "Cluster vs ref (CV only)",
  "TREM2_CV_AD_vs_Control"   = "AD vs Control (CV only)",
  "AD_TREM2_R62H_vs_CV"      = "R62H vs CV (AD only)",
  "AD_TREM2_R47H_vs_CV"      = "R47H vs CV (AD only)",
  "AD_TREM2_R47H_vs_R62H"    = "R47H vs R62H (AD only)"
)

# load cluster order from cluster assignment table (preserves biological ordering)
clust_tab     = read_csv(paste0(main_dir, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv"))
cluster_order = clust_tab$cluster_name

# helper: parse a DEG-list entry name into (cluster, comparison label, direction)
parse_deg_name = function(nm){
  direction    = ifelse(grepl("_up$", nm), "up", "down")
  comp_key_raw = sub("_up$|_down$", "", nm)
  matched_tag = NA; matched_label = NA
  for (tag in names(comp_tags)){
    if (grepl(tag, comp_key_raw, fixed = TRUE)){
      matched_tag = tag; matched_label = comp_tags[[tag]]; break
    }
  }
  if (is.na(matched_tag)) return(NULL)
  cluster = sub(paste0("_", matched_tag, ".*"), "", comp_key_raw)
  list(cluster = cluster, comp_label = matched_label, direction = direction)
}



###########################################################
# NEW: summary table of DEG counts per cluster & comparison, with TF stats
###########################################################
# For each comparison, for each astrocyte subcluster:
#   n_up, n_down, n_total (= n_up + n_down), up:down ratio,
#   n_TF (Fantom5 TFs among the DEGs) and proportion of DEGs that are TFs.
# Plus, per comparison, a POOLED row = union of DEGs across all clusters
#   (like LD_E02a2's unique(unlist(...))), with the same statistics.

message("\n          *** Building DEG summary table (counts + TF stats)... ", Sys.time(), "\n")

# collect gene lists keyed by comparison label -> cluster -> direction
deg_store = list()   # deg_store[[comp_label]][[cluster]][[direction]] = gene vector

for (nm in names(bulk_data$DEGs)){
  genes = bulk_data$DEGs[[nm]]
  info  = parse_deg_name(nm)
  if (is.null(info)) next
  if (!info$cluster %in% cluster_order) next
  deg_store[[info$comp_label]][[info$cluster]][[info$direction]] = genes
}

summary_rows = list()

for (comp_label in names(deg_store)){

  clusters_here = names(deg_store[[comp_label]])

  # ---- per-cluster rows ----
  pooled_up = character(0); pooled_down = character(0)

  for (cl in clusters_here){
    up_genes   = deg_store[[comp_label]][[cl]][["up"]];   if (is.null(up_genes))   up_genes   = character(0)
    down_genes = deg_store[[comp_label]][[cl]][["down"]]; if (is.null(down_genes)) down_genes = character(0)

    all_genes = unique(c(up_genes, down_genes))
    n_up = length(up_genes); n_down = length(down_genes); n_total = length(all_genes)
    n_TF = length(intersect(all_genes, GOI$TF))

    summary_rows[[length(summary_rows)+1]] = tibble(
      comparison    = comp_label,
      cluster       = cl,
      level         = "per_cluster",
      n_up          = n_up,
      n_down        = n_down,
      n_total       = n_total,
      up_down_ratio = ifelse(n_down > 0, n_up / n_down, NA_real_),
      n_TF          = n_TF,
      prop_TF       = ifelse(n_total > 0, n_TF / n_total, NA_real_)
    )

    pooled_up   = c(pooled_up, up_genes)
    pooled_down = c(pooled_down, down_genes)
  }

  # ---- pooled (union across clusters) row ----
  pooled_up   = unique(pooled_up)
  pooled_down = unique(pooled_down)
  pooled_all  = unique(c(pooled_up, pooled_down))
  n_up = length(pooled_up); n_down = length(pooled_down); n_total = length(pooled_all)
  n_TF = length(intersect(pooled_all, GOI$TF))

  summary_rows[[length(summary_rows)+1]] = tibble(
    comparison    = comp_label,
    cluster       = "POOLED_union_all_clusters",
    level         = "pooled_union",
    n_up          = n_up,
    n_down        = n_down,
    n_total       = n_total,
    up_down_ratio = ifelse(n_down > 0, n_up / n_down, NA_real_),
    n_TF          = n_TF,
    prop_TF       = ifelse(n_total > 0, n_TF / n_total, NA_real_)
  )
}

deg_summary_tab = bind_rows(summary_rows) %>%
  mutate(comparison = factor(comparison, levels = unname(comp_tags))) %>%
  arrange(comparison, desc(level == "pooled_union"), cluster)

write_csv(deg_summary_tab,
          file = paste0(sub_out_dir, script_ind, "DEG_summary_counts_TF_per_cluster_and_pooled.csv"))

message("    Saved: ", sub_out_dir, script_ind, "DEG_summary_counts_TF_per_cluster_and_pooled.csv")



###########################################################
# Bar chart: Number of DEGs per cluster per comparison
# Split by direction (up/down) and gene category (TF vs other)
###########################################################
# For each of the 5 DESeq2 comparisons and each astrocyte cluster:
#   - bars above zero  = upregulated DEGs  (solid red = TF, light red = other)
#   - bars below zero  = downregulated DEGs (solid blue = TF, light blue = other)
# Faceted by comparison so you can see which clusters are most strongly affected
# and whether regulatory genes (TFs) are disproportionately involved.

message("\n          *** Plotting DEG count bar chart... ", Sys.time(), "\n")

# ----- parse bulk_data$DEGs into a tidy table -----

deg_rows = list()

for (nm in names(bulk_data$DEGs)){

  genes = bulk_data$DEGs[[nm]]
  if (length(genes) == 0) next

  info = parse_deg_name(nm)
  if (is.null(info)) next
  if (!info$cluster %in% cluster_order) next

  # count TF vs other genes
  n_TF    = sum(genes %in% GOI$TF)
  n_other = length(genes) - n_TF

  deg_rows[[length(deg_rows) + 1]] = tibble(
    cluster    = info$cluster,
    comparison = info$comp_label,
    direction  = info$direction,
    n_TF       = n_TF,
    n_other    = n_other
  )
}

deg_tab = bind_rows(deg_rows)

# pivot to long format (one row per cluster x comparison x direction x gene_cat)
deg_long = deg_tab %>%
  pivot_longer(cols = c(n_TF, n_other),
               names_to  = "gene_cat",
               values_to = "n_genes") %>%
  mutate(
    # downregulated bars go below zero
    n_plot     = ifelse(direction == "up", n_genes, -n_genes),
    cluster    = factor(cluster,    levels = cluster_order),
    comparison = factor(comparison, levels = unname(comp_tags)),
    gene_cat   = factor(gene_cat,   levels = c("n_TF", "n_other")),
    fill_group = factor(paste0(direction, "_", gene_cat),
                        levels = c("up_n_TF", "up_n_other", "down_n_TF", "down_n_other"))
  )

# colour scheme: TF = solid, other = lighter; up = red/salmon, down = blue/lightblue
fill_colors = c(
  "up_n_TF"      = "#C0392B",   # dark red    — TF upregulated
  "up_n_other"   = "#F1948A",   # light red   — other upregulated
  "down_n_TF"    = "#1A5276",   # dark blue   — TF downregulated
  "down_n_other" = "#85C1E9"    # light blue  — other downregulated
)

fill_labels = c(
  "up_n_TF"      = "Up (TF)",
  "up_n_other"   = "Up (other)",
  "down_n_TF"    = "Down (TF)",
  "down_n_other" = "Down (other)"
)

p_deg_bars = ggplot(deg_long,
                    aes(x = cluster, y = n_plot, fill = fill_group)) +
  geom_col(position = "stack", width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey20") +
  scale_fill_manual(values = fill_colors, labels = fill_labels, name = NULL) +
  scale_y_continuous(
    labels = function(x) abs(x),
    name   = "No. of DEGs (padj < 0.1)"
  ) +
  scale_x_discrete(name = NULL) +
  facet_wrap(~ comparison, ncol = 1, scales = "free_y") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "top",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  ) +
  labs(title = "DEGs per cluster and comparison (5-covariate model)",
       caption = "Bars above zero = upregulated; bars below zero = downregulated.\nSolid = transcription factors (Fantom5); lighter = all other genes.")

pdf(file   = paste0(sub_out_dir, script_ind, "DEG_counts_per_cluster_by_comparison.pdf"),
    width  = max(10, length(cluster_order) * 0.55 + 3),
    height = 18)
plot(p_deg_bars)
dev.off()

message("    Saved: ", sub_out_dir, script_ind, "DEG_counts_per_cluster_by_comparison.pdf")



###########################################################
# NEW: 2-panel DEG-count bar plot — AD-vs-Control + R62H-vs-CV stacked
#      with INDEPENDENT (free) y-axes, TFs coloured separately
###########################################################
# AD-vs-Control has far more DEGs than R62H-vs-CV; free y-axes let each panel
# use its own scale so the (smaller) R62H-vs-CV panel is still readable.

message("\n          *** Plotting 2-panel DEG bar chart (AD-vs-Ctrl + R62H-vs-CV, free y)... ", Sys.time(), "\n")

panel_labels = c("AD vs Control (CV only)", "R62H vs CV (AD only)")

deg_long_2panel = deg_long %>%
  filter(comparison %in% panel_labels) %>%
  mutate(comparison = factor(comparison, levels = panel_labels))

p_deg_bars_2panel = ggplot(deg_long_2panel,
                           aes(x = cluster, y = n_plot, fill = fill_group)) +
  geom_col(position = "stack", width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey20") +
  scale_fill_manual(values = fill_colors, labels = fill_labels, name = NULL) +
  scale_y_continuous(
    labels = function(x) abs(x),
    name   = "No. of DEGs (padj < 0.1)"
  ) +
  scale_x_discrete(name = NULL) +
  facet_wrap(~ comparison, ncol = 1, scales = "free_y") +   # free y-axis per panel
  theme_classic(base_size = 11) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold", size = 10),
    legend.position    = "top",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  ) +
  labs(title = "DEGs per cluster: AD-vs-Control vs R62H-vs-CV (5-covariate model)",
       caption = "Independent (free) y-axes per panel.\nBars above zero = upregulated; below zero = downregulated.\nSolid = transcription factors (Fantom5); lighter = all other genes.")

pdf(file   = paste0(sub_out_dir, script_ind, "DEG_counts_2panel_ADvsCtrl_R62HvsCV_freeY.pdf"),
    width  = max(10, length(cluster_order) * 0.55 + 3),
    height = 9)
plot(p_deg_bars_2panel)
dev.off()

message("    Saved: ", sub_out_dir, script_ind, "DEG_counts_2panel_ADvsCtrl_R62HvsCV_freeY.pdf")


#################################################
# plot heatmap top3000 variable genes by cluster_sample
#################################################
# cellwidth and fontsize increased for readable legends/annotation bars;
# PDF width calculated dynamically so the plot is never cut off.

message("\n          *** Plotting top-3000 variable genes heatmap... ", Sys.time(), "\n")

vst_mat = bulk_data$vst_mat

#identify top 3000 variable genes
v1       = rowVars(vst_mat)
pl_genes = names(v1[order(-v1)][1:3000])

pl_mat_X  = bulk_data$gene_Z_scores$clusters_combined[pl_genes,]
lims_X    = 0.3 * c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

# shared cellwidth for dynamic PDF width calculation
cellwidth_val = 5

# Dynamic width: n_cols * cellwidth (pt) / 72 (pt/inch) + buffer for legend + annotation bars
n_cols      = ncol(pl_mat_X)
pdf_width_X = ceiling(n_cols * cellwidth_val / 72) + 6   # +6 inches: legend + dendrogram + margins

pdf(file   = paste0(sub_out_dir, script_ind, "Gene_expr_heatmap_var_genes_top3000.pdf"),
    width  = pdf_width_X, height = 25)

bulkdata_heatmap(
  pl_mat          = pl_mat_X,
  pl_meta         = meta,
  pl_genes        = pl_genes,
  x_col           = "cluster_sample",
  meta_annot_cols = c("TREM2Variant", "NeuropathologicalDiagnosis",
                      "CD33Group", "APOEgroup", "cohort", "Sex", "cluster_name"),
  show_rownames   = FALSE,
  show_colnames   = FALSE,
  cluster_rows    = TRUE, cluster_cols = FALSE,
  color           = viridis(250),
  lims            = lims_X,
  cellwidth       = cellwidth_val,
  cellheight      = 0.5,
  fontsize        = 10,
  title           = "Top 3000 variable genes - Z-score (corrected VST, 5 covariates)"
)

dev.off()

message("    Saved: ", sub_out_dir, script_ind, "Gene_expr_heatmap_var_genes_top3000.pdf")


#################################################
# plot heatmap DEGs combined by cluster_sample
#################################################
# Shows all genes significant (padj < 0.1) in any DESeq2 comparison.
# cluster_rows = FALSE to avoid multi-hour clustering at sub-pixel row height;
# cellheight = NA lets pheatmap auto-scale rows to fill the PDF height.

message("\n          *** Plotting DEGs combined heatmap (all DEGs)... ", Sys.time(), "\n")

pl_genes_DEG = unique(unlist(bulk_data$DEGs))
n_DEGs       = length(pl_genes_DEG)

message("    Total DEGs across all comparisons: ", n_DEGs)

pl_mat_DEG    = bulk_data$gene_Z_scores$clusters_combined[pl_genes_DEG, ]
lims_DEG      = 0.3 * c(-max(abs(pl_mat_DEG), na.rm = TRUE), max(abs(pl_mat_DEG), na.rm = TRUE))
pdf_width_DEG = ceiling(ncol(pl_mat_DEG) * cellwidth_val / 72) + 6

pdf(file   = paste0(sub_out_dir, script_ind, "Gene_expr_heatmap_DEGs_comb.pdf"),
    width  = pdf_width_DEG, height = 40)

bulkdata_heatmap(
  pl_mat          = pl_mat_DEG,
  pl_meta         = meta,
  pl_genes        = pl_genes_DEG,
  x_col           = "cluster_sample",
  meta_annot_cols = c("TREM2Variant", "NeuropathologicalDiagnosis",
                      "CD33Group", "APOEgroup", "cohort", "Sex", "cluster_name"),
  show_rownames   = FALSE,
  show_colnames   = FALSE,
  cluster_rows    = FALSE,   # large gene sets — clustering not useful at sub-pixel row height
  cluster_cols    = FALSE,
  color           = viridis(250),
  lims            = lims_DEG,
  cellwidth       = cellwidth_val,
  cellheight      = NA,      # auto-scale: pheatmap fills the PDF height
  fontsize        = 10,
  title           = paste0("All DEGs combined (n=", n_DEGs, ") - Z-score (corrected VST, 5 covariates)")
)

dev.off()

message("    Saved: ", sub_out_dir, script_ind, "Gene_expr_heatmap_DEGs_comb.pdf")



###########################################################
# plot volcano plots with DEGs highlighted
###########################################################

pl = list()

for (comp in names(bulk_data$deseq_results)){

  t1 = bulk_data$deseq_results[[comp]]
  t1 = t1[!is.na(t1$padj),]
  t1$log10pval = log10(t1$pvalue)
  t1$log10padj = log10(t1$padj)
  t1$DEG = t1$padj <= 0.1
  t1$highlight = t1$DEG

  t1$gene_cat = "Other"
  t1$gene_cat[t1$DEG & t1$log2FoldChange>0] = "up"
  t1$gene_cat[t1$DEG & t1$log2FoldChange<0] = "down"

  t1$plot_label = ""

  t3 = t1[order(t1$log10pval),]
  t4 = t1[order(t1$log2FoldChange),]
  t5 = t1[order(-t1$log2FoldChange),]
  v5 = unique(c(rownames(t3)[1:20], rownames(t4)[1:20], rownames(t5)[1:20]))

  t1$plot_label[rownames(t1) %in% v5] = rownames(t1)[rownames(t1) %in% v5]

  pl[[comp]] = ggplot(t1, aes(x = log2FoldChange, y = -log10pval, color = gene_cat))+
    geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_hline(yintercept = -log10(0.1), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_point(aes(size = highlight), alpha = 0.8)+
    geom_label_repel(aes(label = plot_label), seed = 42, min.segment.length = 0, max.overlaps = Inf,
                     max.time = 5)+
    scale_size_manual(limits = c(FALSE, TRUE), values = c(1, 2))+
    scale_color_manual(limits = c("up", "Other", "down"),
                       values = c("red", "grey40", "blue"))+
    theme_minimal()+
    labs(title = comp)

}


pdf(file = paste0(sub_out_dir,script_ind, "Volcano_plots_DEGs_by_comp_group_cluster.pdf"),
    width = 10, height = 8)
{
  lapply(pl, function(x){x})
}
dev.off()





#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed LD_E02c (5-covariate version) ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")
