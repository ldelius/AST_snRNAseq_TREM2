message("\n\n##########################################################################\n",
        "# Start LD_H02: Variance-driven gene characterisation ", Sys.time(),
        "\n##########################################################################\n\n")

# TREM2Variant is the single primary variable; the other H01 covariates are
# removed from the corrected matrix.
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
library(variancePartition)
library(BiocParallel)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(DOSE)
library(enrichplot)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify output directory
out_dir = paste0(main_dir,"LD_H_VarPartition_output/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

#specify script/output index as prefix for file names
script_ind = "LD_H02_v02_"

### load dataset (LD_H01 output: cleaned meta + varPart_analysis already populated)
bulk_data = qread(file = paste0(out_dir, "LD_H01_v01_bulk_data.qs"))

###define thresholds for minimal explained variance for each gene and primary variable
var_thresholds = c(0.05, 0.10, 0.15, 0.20)

### define formula for covariate correction of the vst matrix (linear mixed model)
# Use random effects for cluster_name and Braak. Two-level categorical and
# numeric covariates are fixed because random effects for two-level factors
# produced singular fits and incomplete removal. H01 uses random effects for
# variance decomposition; this model constructs a covariate-corrected matrix.
form_corr = ~ (1 | cluster_name) + (1 | Braak) +
  cohort + BrainRegion + Sex + CD33Group + APOEgroup +
  PMI_scaled + RNA_counts_scaled + Age_scaled + plaque_dens_scaled +
  pctAT8PositiveArea + pctPHF1PositiveArea + pct4G8PositiveArea

### get Fantom5 transcription-factor panel
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


### create heatmap annotations for Complex Heatmap

create_heatmap_annot = function(annot_df, annot_dim = c("column", "row"),
                                annotation_name_side = c("left", "top")){

  #define colors

  annot_colors = list()

  for (col1 in names(annot_df)){

    v1 = annot_df[[col1]]

    if (is.numeric(v1)){
      v2 = colorRamp2(breaks = seq(from = min(na.omit(v1)), to = max(na.omit(v1)), length.out = 100),
                      colors = viridis(100))
      annot_colors = c(annot_colors, list(v2))

    } else if (is.factor(v1)) {
      v2 = pal(levels(na.omit(v1)))
      names(v2) = levels(na.omit(v1))
      annot_colors = c(annot_colors, list(v2))
    } else {
      stop("   !ERROR: categorical annotations need to be factors")
    }
  }
  names(annot_colors) = names(annot_df)

  # Create the Column Annotation Object
  annot <- HeatmapAnnotation(
    df = annot_df,
    col = annot_colors,
    which = annot_dim,
    annotation_name_side = annotation_name_side[1]
  )

  return(annot)
}


#######################################
# plot variance explained by primary variable per gene (ranked)
#######################################
# Rank genes by the TREM2Variant variance fraction and show extraction thresholds.

t1 = bulk_data$varPart_analysis$varPart

prim_vars = bulk_data$varPart_analysis$primary_var
prim_var  = prim_vars[1]

t2 = t1[t1[[prim_var]] > 0.05, , drop = FALSE]

t3 = data.frame(gene = rownames(t2), var_expl = t2[[prim_var]],
                stringsAsFactors = FALSE)
t3 = t3[order(-t3$var_expl),]
t3$rank = seq_len(nrow(t3))

#label the top genes
t3$label = ""
n_lab = min(30, nrow(t3))
if (n_lab > 0){t3$label[1:n_lab] = t3$gene[1:n_lab]}

p1 = ggplot(t3, aes(x = rank, y = var_expl))+
  geom_point(aes(color = label != ""), alpha = 0.5)+
  geom_hline(yintercept = var_thresholds, color = "red", linetype = "dashed", linewidth = 0.5)+
  geom_text_repel(aes(label = label), color = "red", max.overlaps = 30)+
  scale_color_manual(limits = c(FALSE, TRUE), values = c("grey", "red"))+
  labs(title = paste0("Fraction variance explained by ", prim_var, " per gene (ranked)"),
       x = "Gene rank",
       y = paste0("Fraction explained by ", prim_var))+
  theme_bw()

pdf(file = paste0(out_dir,script_ind, "Fract_variance_explained_by_prim_var.pdf"),
    width = 7, height = 7)
{
  plot(p1)
}
dev.off()


#######################################
# extract genes with >X% of variance explained by the primary variable
#######################################

t1 = bulk_data$varPart_analysis$varPart

for (var_threshold in var_thresholds){

  var_genes1 = rownames(t1)[t1[[prim_var]] >= var_threshold]

  bulk_data$varPart_analysis$var_genes[[paste0(prim_var, "_thr_", var_threshold)]] = var_genes1

}

### save table with var genes and TFs

l1 = bulk_data$varPart_analysis$var_genes
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

write_csv(as_tibble(m1), file = paste0(out_dir,script_ind, "var_genes_most_var_expl_by_prim_var.csv"))

t1 = tibble(gene_set = names(l3), N_genes = lengths(l3))

write_csv(t1, file = paste0(out_dir,script_ind, "var_genes_most_var_expl_by_prim_var_N.csv"))



#######################################
# extract/plot key covariates for all var_gene_sets
#######################################

varPart = bulk_data$varPart_analysis$varPart

vars_full = c(bulk_data$varPart_analysis$candidate_covars,
              bulk_data$varPart_analysis$primary_var)

var_genes_list = bulk_data$varPart_analysis$var_genes[lengths(bulk_data$varPart_analysis$var_genes)>0]

var_genes_var_means_tab = NULL

var_means_by_var_comp_mat = matrix(nrow = length(vars_full), ncol = length(var_genes_list),
                dimnames = list(vars_full, names(var_genes_list)))
pl = list()

for (comp in names(var_genes_list)){

  var_genes = var_genes_list[[comp]]

  m1 = varPart[var_genes, vars_full]

  vp_mean <- sort(colMeans(m1), decreasing = TRUE)

  vp_mean_tab = tibble(covariate = names(vp_mean), Variance_explained_mean = vp_mean)

  var_genes_var_means_tab = rbind(var_genes_var_means_tab, cbind(comp = comp, vp_mean_tab))

  ### add variance means by comp and variable to var_means_by_var_comp_mat matrix
  var_means_by_var_comp_mat[,comp] = vp_mean[match(rownames(var_means_by_var_comp_mat),
                                                   names(vp_mean))]

  ### Plot vp_mean by comp
  pl[[comp]] <- plotVarPart(m1) + ggtitle(paste0(comp, "\n - Variance explained by covariates"))+
    scale_x_discrete(limits = names(vp_mean))

  pl[[paste0(comp, "_max20")]] <- pl[[comp]] + scale_y_continuous(limits = c(0,20)) +
    ggtitle(paste0(comp, "\n - Variance explained by covariates (max 20%)"))

}

write_csv(var_genes_var_means_tab, file = paste0(out_dir, script_ind, "var_genes_var_expl_by_covars_mean_ord.csv"))


pdf(file = paste0(out_dir,script_ind, "var_genes_var_expl_by_covars_mean_ord_by_comp.pdf"),
    width = 10, height = 10)
{
  lapply(pl, function(x){x})
}
dev.off()


pdf(file = paste0(out_dir,script_ind, "var_genes_var_expl_by_covars_mean_ord_all_comps.pdf"),
    width = 7, height = 7)
{
  p1 = pheatmap::pheatmap(var_means_by_var_comp_mat,
                          cluster_rows = TRUE, cluster_cols = TRUE,
                          show_rownames = TRUE, show_colnames = TRUE,
                          color = colorRampPalette(c("white", "red"))(250),
                          breaks = seq(0, max(var_means_by_var_comp_mat), length.out = 251),
                          border_color = NA, cellwidth = 10, cellheight = 10,
                          fontsize = 10,
                          main = "Fraction variance explained")
}
dev.off()



#######################################
# Calculate corrected vst matrix and z-scores with mixed model
#######################################
# correction formula (form_corr) is defined explicitly at the top of the script.

vst_mat = bulk_data$vst_mat_uncorr
meta = as.data.frame(bulk_data$meta)

message("\n   Correction formula:\n   ", deparse1(form_corr), "\n")

#calculate Z-score per gene by pseudobulk (cluster_sample) for uncorrected vst matrix
z_mat = t(apply(vst_mat, 1, scale))
colnames(z_mat) = colnames(vst_mat)
bulk_data$gene_Z_scores_uncorr = z_mat


# Extract covariate-corrected residuals per gene. fitVarPartModel is incompatible
# with lme4 >= 2.0, so fit the same REML model directly with lmer. Failed fits are
# reported and fall back to mean-centring.
param <- SnowParam(6, "SOCK", progressbar = TRUE)

# align meta rows to expression columns, build the response formula (y ~ covariates)
meta = meta[colnames(vst_mat), ]
form_resid = update(form_corr, y ~ .)   # add response: y ~ (1 | cov) + ... + cov

fit_resid_one = function(i, mat, meta, form){
  d = meta
  d$y = mat[i, ]
  tryCatch(
    suppressWarnings(residuals(lme4::lmer(form, data = d, REML = TRUE))),
    error = function(e){ rep(NA_real_, nrow(d)) }
  )
}

residList = bplapply(seq_len(nrow(vst_mat)), fit_resid_one,
                     mat = vst_mat, meta = meta, form = form_resid,
                     BPPARAM = param)

# convert list to matrix
residMatrix = do.call(rbind, residList)
rownames(residMatrix) = rownames(vst_mat)
colnames(residMatrix) = colnames(vst_mat)

# handle any genes whose mixed model failed (fall back to mean-centring, report)
failed = which(rowSums(is.na(residMatrix)) > 0)
if (length(failed) > 0){
  message("\n   !! ", length(failed), " of ", nrow(residMatrix),
          " genes failed the mixed-model fit - mean-centring those instead\n")
  for (i in failed){
    y = vst_mat[i, ]
    residMatrix[i, ] = y - mean(y)
  }
}

bulk_data$vst_mat_corr = residMatrix
bulk_data$vst_mat_corr_formula = form_resid


#calculate Z-score per gene by pseudobulk (cluster_sample) for corrected vst matrix
z_mat = t(apply(residMatrix, 1, scale))
colnames(z_mat) = colnames(residMatrix)
bulk_data$gene_Z_scores_corr = z_mat

qsave(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.qs"))


#################################################
# plot heatmap variable genes combined by cluster_sample (uncorrected)
#################################################

pl_gene_list = bulk_data$varPart_analysis$var_genes

pl = list()

for (comp in names(pl_gene_list)){

  message("   * Plot variable genes ", comp, " - ", Sys.time())

  pl_genes = pl_gene_list[[comp]]

  if (length(pl_genes) < 2){next}

  #define metadata, plot matrix, plot scale

  pl_mat_X = bulk_data$gene_Z_scores_uncorr

  lims_X = 0.1*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

  pl_meta = bulk_data$meta[colnames(pl_mat_X),]

  pl_mat_X = pl_mat_X[pl_genes, rownames(pl_meta)]

  #Create Column Annotation Data (Top Bars)
  annot_vars = c(bulk_data$varPart_analysis$candidate_covars, prim_vars)

  col_anno_df <- as.data.frame(pl_meta[,annot_vars])

  for (col1 in names(col_anno_df)){
    if (!(is.numeric(col_anno_df[[col1]]) | is.factor(col_anno_df[[col1]]))){
      col_anno_df[[col1]] = factor(col_anno_df[[col1]], levels = unique(col_anno_df[[col1]]))}
  }
  rownames(col_anno_df) = pl_meta$cluster_sample

  col_annot = create_heatmap_annot(annot_df = col_anno_df)

  #create heatmap

  pl[[comp]] <- Heatmap(
    pl_mat_X,
    name = "Expr \nZ-score",
    col = colorRamp2(breaks = seq(from = lims_X[1], to = lims_X[2], length.out = 100),
                     colors = viridis(100)),
    cluster_columns = FALSE,
    column_title = paste0(comp, "\n - Expr Z-score"),
    top_annotation = col_annot,
    cluster_rows = TRUE,
    row_title = NULL,
    left_annotation = NULL,
    show_row_names = TRUE,
    show_column_names = FALSE,
    use_raster = FALSE,
    width = 15, height = 15)

}


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmaps_var_genes_uncorr.pdf"),
    width = 15, height = 15)
{
  lapply(pl, function(x){x})
}
dev.off()


#################################################
# plot heatmap variable genes combined by cluster_sample (corrected)
#################################################

pl_gene_list = bulk_data$varPart_analysis$var_genes

pl = list()

for (comp in names(pl_gene_list)){

  message("   * Plot variable genes ", comp, " - ", Sys.time())

  pl_genes = pl_gene_list[[comp]]

  if (length(pl_genes) < 2){next}

  #define metadata, plot matrix, plot scale

  pl_mat_X = bulk_data$gene_Z_scores_corr

  lims_X = 0.1*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

  pl_meta = bulk_data$meta[colnames(pl_mat_X),]

  pl_mat_X = pl_mat_X[pl_genes, rownames(pl_meta)]

  #Create Column Annotation Data (Top Bars)
  annot_vars = c(bulk_data$varPart_analysis$candidate_covars, prim_vars)

  col_anno_df <- as.data.frame(pl_meta[,annot_vars])

  for (col1 in names(col_anno_df)){
    if (!(is.numeric(col_anno_df[[col1]]) | is.factor(col_anno_df[[col1]]))){
      col_anno_df[[col1]] = factor(col_anno_df[[col1]], levels = unique(col_anno_df[[col1]]))}
  }
  rownames(col_anno_df) = pl_meta$cluster_sample

  col_annot = create_heatmap_annot(annot_df = col_anno_df)

  #create heatmap

  pl[[comp]] <- Heatmap(
    pl_mat_X,
    name = "Expr \nZ-score",
    col = colorRamp2(breaks = seq(from = lims_X[1], to = lims_X[2], length.out = 100),
                     colors = viridis(100)),
    cluster_columns = FALSE,
    column_title = paste0(comp, "\n - Z-score var genes (corrected)"),
    top_annotation = col_annot,
    cluster_rows = TRUE,
    row_title = NULL,
    left_annotation = NULL,
    show_row_names = TRUE,
    show_column_names = FALSE,
    use_raster = FALSE,
    width = 15, height = 15)

}


pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmaps_var_genes_corr.pdf"),
    width = 15, height = 15)
{
  lapply(pl, function(x){x})
}
dev.off()




#################################################
### GO-BP over-representation analysis of variable genes by comp
#################################################

GO_list = list()
GO_results_tab = NULL

N_comps = length(bulk_data$varPart_analysis$var_genes)

for (i in 1:N_comps){

  comp = names(bulk_data$varPart_analysis$var_genes)[i]

  message("\n          *** GO analysis var genes by comp ", comp,
          " (", i, " of ",N_comps,  ") - ", Sys.time(),"\n")

  ego = NULL

  go_genes =  bulk_data$varPart_analysis$var_genes[[comp]]

  if (length(go_genes)>2){

    ego = enrichGO(gene         = go_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

    GO_list[[comp]] = ego
  }

  if (!is.null(ego)){

    t1 = ego@result[ego@result$p.adjust<=0.05,]

    if (nrow(t1)>0){
      t2 = cbind(comp = comp, t1)
      GO_results_tab = rbind(GO_results_tab, t2)
    }
  }
}

GO_results_tab = GO_results_tab[GO_results_tab$Count>1,]

write_csv(GO_results_tab, file = paste0(out_dir,script_ind, "GO_results_by_comp.csv"))

bulk_data$varPart_analysis$GO_results$by_comp_GO_list = GO_list
bulk_data$varPart_analysis$GO_results$by_comp_GO_res = GO_results_tab

qsave(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.qs"))


############################################################################
### visualise GO analysis (dotplot top 10 terms by comp, GO term network plot)
############################################################################

t1 = GO_results_tab

t2 = NULL

for (comp in unique(t1$comp)){
  t3 = t1[t1$comp == comp,]
  if (nrow(t3)>10){t3 = t3[1:10,]}
  t2 = rbind(t2, t3)
}

comps = unique(t1$comp)

p1 = ggplot(t2, aes(x = comp, y = Description, size = Count, colour = comp))+
  geom_point()+
  scale_color_manual (limits = comps, values = pal(comps))+
  scale_size_continuous(limits = c(0, max(t2$Count)))+
  scale_x_discrete(limits = comps)+
  scale_y_discrete(limits = unique(t2$Description))+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        axis.text.y  = element_text(size = 12),
        axis.title   = element_text(size = 14),
        plot.title   = element_text(size = 15),
        legend.text  = element_text(size = 11),
        legend.title = element_text(size = 12))+
  labs(title = "Top10 GO terms per comp")

pdf(file = paste0(out_dir,script_ind, "GO_results_by_comp_dotplot_top_terms.pdf"),
    width = 11, height = 12)
plot(p1)
dev.off()



###plot GO term network plot

pl = list()

N_comps = length(GO_list)

for (i in 1:N_comps){

  comp = names(GO_list)[i]

  message("\n          *** GO network var genes by comp ", comp,
          " (", i, " of ",N_comps,  ") - ", Sys.time(),"\n")

  edo = GO_list[[comp]]
  t1 = edo@result
  t1 = t1[t1$p.adjust<0.01,]

  if (nrow(t1)>2){
    edo <- pairwise_termsim(edo)
    p1 = emapplot(edo, showCategory = 100)+labs(title = comp)
    pl[[comp]] <- p1
  }
}

pdf(file = paste0(out_dir,script_ind, "GO_results_by_comp_network_plot.pdf"),
    width = 15, height = 12)
lapply(pl, function(x){x})
dev.off()






#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed LD_H02 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")
