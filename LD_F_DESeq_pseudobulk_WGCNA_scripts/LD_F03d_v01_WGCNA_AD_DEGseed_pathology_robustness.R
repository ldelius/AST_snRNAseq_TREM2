message("\n\n##########################################################################\n",
        "# Start LD_F03d: WGCNA pathology robustness analysis ", Sys.time(),
        "\n##########################################################################\n\n")

library(qs)
library(tidyverse)
library(limma)
library(DESeq2)
library(colorRamps)
library(viridis)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(ggrepel)
library(WGCNA)
library(enrichplot)
library(lme4)
library(lmerTest)


### define directories

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

in_dir       = paste0(main_dir, "LD_F_DESeq_pseudobulk_WGCNA/")
base_out_dir = paste0(main_dir, "LD_F_DESeq_pseudobulk_WGCNA/LD_F03d_v01/")
dir.create(base_out_dir, showWarnings = FALSE, recursive = TRUE)


### load DEseq2 dataset (F02c: AD-only, both variants, 7-covar group-protected corrected matrix)
bulk_data = qread(file = paste0(in_dir, "LD_F02c_v01_bulk_data.qs"))


### get gene of interest gene sets

GOI = list()
t1 = read_csv(paste0(main_dir, "data_TREM2_michael/A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol

t1 = read_csv(paste0(main_dir, "data_TREM2_michael/A_input/GOI_sets_251020.csv"))

for (goi_set in unique(t1$gene_set)){
  GOI[[goi_set]] = t1$gene[t1$gene_set == goi_set]
}


### covariate setup: numeric coercions + Braak conversion
# Done once on bulk_data before the loop; all models share the same base meta.

t1 = bulk_data$meta

t1$PostMortemInterval  = as.numeric(t1$PostMortemInterval)
t1$plaque_dens         = as.numeric(t1$plaque_dens)
t1$Age                 = as.numeric(t1$Age)
t1$pctAT8PositiveArea  = as.numeric(t1$pctAT8PositiveArea)
t1$pctPHF1PositiveArea = as.numeric(t1$pctPHF1PositiveArea)
t1$pct4G8PositiveArea  = as.numeric(t1$pct4G8PositiveArea)
t1$N_cells             = as.numeric(t1$N_cells)

t1$RNA_counts = apply(bulk_data$counts, 2, sum)

Braak_conv = tibble(Braak_stage   = c("0", "I", "II", "III", "IV", "V", "V,VI", "VI"),
                    Braak_numeric = c(0, 1, 2, 3, 4, 5, 5.5, 6))
t1$Braak_numeric = Braak_conv$Braak_numeric[match(t1$Braak, Braak_conv$Braak_stage)]

bulk_data$meta = t1

### cluster names (stable across both models)
cluster_names = unique(bulk_data$meta$cluster_name)


### pathology correction models
# Model A: minimal non-collinear set (Braak dominates pathology variance at 4.8%;
#           plaque_dens excluded due to r = 0.821 collinearity with Braak)
# Model B: full measured pathology profile (all IHC + plaque added on top of A)
path_models = list(
  ModelA = c("Age", "Braak_numeric"),
  ModelB = c("Age", "Braak_numeric", "plaque_dens",
             "pctAT8PositiveArea", "pctPHF1PositiveArea", "pct4G8PositiveArea")
)


###########################################################
# functions
###########################################################

pal = function(v){
  v2 = length(unique(v))
  if (v2 == 2){
    p2 = c("grey20", "dodgerblue")
  } else if (v2 == 3){
    p2 = c("dodgerblue", "grey20", "orange")
  } else if (v2 == 4){
    p2 = c("dodgerblue", "grey20", "orange", "green4")
  } else if (v2 < 6){
    p2 = matlab.like(6)[1:v2]
  } else {
    p2 = matlab.like(v2)
  }
  return(p2)
}


bulkdata_heatmap = function(pl_mat, pl_meta, x_col, pl_genes = NULL,
                            meta_annot_cols = NULL, gene_annot = NULL,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            show_rownames = T, show_colnames = T,
                            color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                            lims = NULL, cellwidth = 15, cellheight = 10,
                            fontsize = 10, title = "Z-score vst-norm gene expression"){

  if (is.null(pl_genes)){pl_genes = rownames(pl_mat)}

  pl_mat  = pl_mat[match(pl_genes, rownames(pl_mat), nomatch = 0),]
  pl_meta = pl_meta[match(pl_meta[[x_col]], colnames(pl_mat), nomatch = 0),]

  if (!is.null(meta_annot_cols)){
    annot_row = NULL
    annot_col = data.frame(row.names = pl_meta[[x_col]])
    for (col1 in meta_annot_cols){
      v1 = pl_meta[match(colnames(pl_mat), pl_meta[[x_col]]),][[col1]]
      v1 = factor(v1, levels = unique(pl_meta[[col1]]))
      annot_col = as.data.frame(cbind(annot_col, v1))
    }
    colnames(annot_col) = meta_annot_cols
    rownames(annot_col) = colnames(pl_mat)
    annot_colors = lapply(meta_annot_cols, function(x){
      v1 = pal(unique(annot_col[[x]]))
      names(v1) = levels(annot_col[[x]])
      return(v1)
    })
    names(annot_colors) = meta_annot_cols
  } else {
    annot_col    = NULL
    annot_row    = NULL
    annot_colors = NULL
  }

  if (!is.null(gene_annot)){
    annot_row = data.frame(row.names = row.names(gene_annot))
    for (col1 in colnames(gene_annot)){
      v1 = gene_annot[[col1]]
      v1 = factor(v1, levels = unique(v1))
      annot_row = as.data.frame(cbind(annot_row, v1))
    }
    colnames(annot_row) = colnames(gene_annot)
    rownames(annot_row) = rownames(gene_annot)
    annot_colors_row = lapply(colnames(gene_annot), function(x){
      v1 = pal(unique(annot_row[[x]]))
      names(v1) = levels(annot_row[[x]])
      return(v1)
    })
    names(annot_colors_row) = colnames(gene_annot)
    annot_colors = c(annot_colors, annot_colors_row)
  }

  if (is.null(lims)){lims = c(-0.7*max(abs(na.omit(pl_mat))), 0.7*max(abs(na.omit(pl_mat))))}

  p1 = pheatmap::pheatmap(pl_mat, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                          color = color,
                          breaks = seq(lims[1], lims[2], length.out = length(color)+1),
                          show_rownames = show_rownames, show_colnames = show_colnames,
                          annotation_col = annot_col, annotation_row = annot_row,
                          annotation_colors = annot_colors,
                          border_color = NA, cellwidth = cellwidth, cellheight = cellheight,
                          fontsize = fontsize, main = title)
  return(p1)
}


create_heatmap_annot = function(annot_df, annot_dim = c("column", "row"),
                                annotation_name_side = c("left", "top")){
  annot_colors = list()
  for (col1 in names(annot_df)){
    v1 = annot_df[[col1]]
    if (is.numeric(v1)){
      v2 = colorRamp2(breaks = seq(from = min(na.omit(v1)), to = max(na.omit(v1)), length.out = 100),
                      colors = viridis(100))
      annot_colors = c(annot_colors, list(v2))
    } else {
      v2 = pal(unique(v1))
      names(v2) = unique(v1)
      annot_colors = c(annot_colors, list(v2))
    }
  }
  names(annot_colors) = names(annot_df)
  annot <- HeatmapAnnotation(
    df = annot_df,
    col = annot_colors,
    which = annot_dim,
    annotation_name_side = annotation_name_side[1]
  )
  return(annot)
}


###########################################################
# Main loop over pathology correction models
###########################################################

for (model_name in names(path_models)){

  message("\n\n##########################################################################\n",
          "# Pathology model: ", model_name, " - ", Sys.time(),
          "\n   Correcting for: ", paste(path_models[[model_name]], collapse = ", "),
          "\n##########################################################################\n\n")

  out_dir    = paste0(base_out_dir, model_name, "/")
  script_ind = paste0("LD_F03d_v01_", model_name, "_")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # fresh copy of bulk_data for this model
  bd = bulk_data


  #######################################
  # Single-pass correction from uncorrected VST matrix
  # Replicates F02c's sequential group-protected removeBatchEffect for the 7 main
  # covariates (same order as F02c's covars_corr), then continues with the
  # model-specific pathology covariates — all in one sequential pass starting
  # from vst_mat_uncorr. This is equivalent to the previous two-step approach
  # (F02c's 7-covar output + F03d pathology correction) because F02c itself uses
  # sequential correction, and sequential removeBatchEffect is order-independent
  # given the same group-protection at each step. The change makes F03d
  # self-contained and removes the dependency on F02c's pre-corrected matrix.
  #######################################

  message("\n   *Single-pass correction: 7 main covariates + pathology (", model_name, ") \n")

  # align meta to uncorrected vst_mat column order
  meta_vst = bd$meta[match(colnames(bd$vst_mat_uncorr), bd$meta$cluster_sample), ]

  # group-protection: preserve TREM2 variant group direction
  design_group = model.matrix(~ 0 + as.factor(meta_vst$group))

  # start from uncorrected matrix
  vst_mat_path = bd$vst_mat_uncorr

  # 7 main covariates — replicating F02c covars_corr order
  # factors use batch=, numerics use covariates= (matching F02c's loop logic)
  main_covars_factor  = c("cohort", "BrainRegion", "APOEgroup", "CD33Group", "Sex")
  main_covars_numeric = c("PostMortemInterval", "log10_nCount_RNA")

  for (cov1 in main_covars_factor){
    if (length(unique(meta_vst[[cov1]])) > 1){
      vst_mat_path = removeBatchEffect(vst_mat_path,
                                       batch  = meta_vst[[cov1]],
                                       design = design_group)
    }
  }
  for (cov1 in main_covars_numeric){
    if (length(unique(meta_vst[[cov1]])) > 1){
      vst_mat_path = removeBatchEffect(vst_mat_path,
                                       covariates = as.numeric(meta_vst[[cov1]]),
                                       design     = design_group)
    }
  }

  # pathology covariates (model-specific); impute NAs with column means
  path_cov_cols = path_models[[model_name]]
  path_cov_mat  = as.matrix(sapply(meta_vst[, path_cov_cols, drop = FALSE], as.numeric))
  rownames(path_cov_mat) = meta_vst$cluster_sample

  for (col_i in colnames(path_cov_mat)){
    na_idx = is.na(path_cov_mat[, col_i])
    if (any(na_idx)){
      message("   Imputing ", sum(na_idx), " NAs in ", col_i, " with column mean")
      path_cov_mat[na_idx, col_i] = mean(path_cov_mat[!na_idx, col_i])
    }
  }

  for (col_i in colnames(path_cov_mat)){
    vst_mat_path = removeBatchEffect(vst_mat_path,
                                     covariates = path_cov_mat[, col_i],
                                     design     = design_group)
  }


  #######################################
  # calculate WGCNA soft-thresholding power
  #######################################

  message("\n\n          *** Calculate soft-thresholding power... ", Sys.time(), "\n\n")

  input_mat = t(vst_mat_path[unique(unlist(bd$DEGs)), ])
  bd$wgcna$input_mat = input_mat

  allowWGCNAThreads()

  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

  sft = pickSoftThreshold(
    input_mat,
    powerVector = powers,
    verbose = 5
  )

  t1 = as_tibble(sft$fitIndices)
  t1$signed_R2 = -sign(t1$slope) * t1$SFT.R.sq

  wgcna_power = t1$Power[t1$signed_R2 > 0.8][1]
  if (wgcna_power < 6){ wgcna_power = 6 }

  bd$wgcna$power_table = t1
  bd$wgcna$power       = wgcna_power

  p1 = ggplot(t1, aes(x = Power, y = signed_R2)) + geom_point() +
    geom_line(color = "grey") + geom_label(aes(label = Power)) +
    geom_hline(yintercept = 0.8, color = "red") +
    geom_vline(xintercept = wgcna_power, color = "red") +
    labs(title = "Scale independence", x = "Soft Threshold (power)",
         y = "Scale Free Topology Model Fit, signed R^2")

  p2 = ggplot(t1, aes(x = Power, y = mean.k.)) + geom_point() +
    geom_line() + geom_label(aes(label = Power)) +
    labs(title = "Mean connectivity", x = "Soft Threshold (power)",
         y = "Mean connectivity")

  pdf(file = paste0(out_dir, script_ind, "Power_thresholding_tests.pdf"),
      width = 9, height = 4)
  { plot(p1 + p2) }
  dev.off()

  qsave(bd, file = paste0(out_dir, script_ind, "bulk_data.qs"))


  #######################################
  # calculate WGCNA with selected soft-thresholding power
  #######################################

  temp_cor <- cor
  cor <- WGCNA::cor

  netwk <- blockwiseModules(input_mat,
                            power             = wgcna_power,
                            networkType       = "signed",
                            deepSplit         = 2,
                            pamRespectsDendro = F,
                            minModuleSize     = 30,
                            maxBlockSize      = 20000,
                            reassignThreshold = 0,
                            mergeCutHeight    = 0.25,
                            saveTOMs          = FALSE,
                            numericLabels     = T,
                            verbose           = 3)

  cor <- temp_cor

  bd$wgcna$network = netwk

  qsave(bd, file = paste0(out_dir, script_ind, "bulk_data.qs"))


  #######################################
  # Basic network characterisation and extraction of modules
  #######################################

  netwk = bd$wgcna$network

  t1 <- data.frame(
    gene          = names(netwk$colors),
    module        = paste0("M", netwk$colors),
    module_number = netwk$colors,
    colors        = labels2colors(netwk$colors)
  )

  mod_gene_tab = t1[order(t1$module_number), ]

  mods = unique(mod_gene_tab$module)

  bd$wgcna$mod_gene_tab = mod_gene_tab

  for (mod1 in mods){
    bd$wgcna$mod_genes[[mod1]] = mod_gene_tab$gene[mod_gene_tab$module == mod1]
  }

  qsave(bd, file = paste0(out_dir, script_ind, "bulk_data.qs"))


  ### save table with module genes and TFs

  l1 = bd$wgcna$mod_genes
  l2 = lapply(l1, function(x){ x = x[x %in% GOI$TF] })
  names(l2) = paste0(names(l1), "_TF")

  l3 = c(l1, l2)

  m1 = matrix(nrow = max(lengths(l3)), ncol = length(l3))
  colnames(m1) = names(l3)

  for (i in names(l3)){
    v1 = l3[[i]]
    if (length(v1) > 0){ m1[1:length(v1), i] = v1 }
  }
  m1[is.na(m1)] = ""

  write_csv(as_tibble(m1), file = paste0(out_dir, script_ind, "Module_genes.csv"))

  t1 = tibble(gene_set = names(l3), N_genes = lengths(l3))
  write_csv(t1, file = paste0(out_dir, script_ind, "Module_genes_N.csv"))


  ### plot network dendrogram

  mergedColors = labels2colors(netwk$colors)

  pdf(file = paste0(out_dir, script_ind, "WGCNA_network_dendrogram.pdf"),
      width = 10, height = 3)
  {
    plotDendroAndColors(
      netwk$dendrograms[[1]],
      mergedColors,
      "Module colors",
      dendroLabels = FALSE,
      hang         = 0.03,
      addGuide     = TRUE,
      guideHang    = 0.05)
  }
  dev.off()


  ### extract module eigengene matrix

  m1 = t(netwk$MEs)
  rownames(m1) = str_replace_all(rownames(m1), "ME", "M")
  m2 = m1[mods, ]

  bd$wgcna$mod_eigengene_mat = m2


  #######################################################
  ### plot module eigengene heatmap by sample and module
  #######################################################

  meta = bd$meta

  pl_mat_X = bd$wgcna$mod_eigengene_mat
  lims_X   = 0.2 * c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

  col_anno_df <- as.data.frame(meta[, c("cluster_name", "APOEgroup", "CD33Group", "BrainRegion",
                                        "plaque_dens", "NeuropathologicalDiagnosis", "TREM2Variant")])
  for (col1 in names(col_anno_df)){
    if (!is.numeric(col_anno_df[[col1]])){
      col_anno_df[[col1]] = factor(col_anno_df[[col1]], levels = unique(col_anno_df[[col1]]))
    }
  }
  rownames(col_anno_df) = meta$cluster_sample

  col_annot = create_heatmap_annot(annot_df = col_anno_df)

  ht_list <- Heatmap(
    pl_mat_X,
    name              = "Z-score",
    col               = colorRamp2(breaks = seq(from = lims_X[1], to = lims_X[2], length.out = 100),
                                   colors = viridis(100)),
    cluster_columns   = FALSE,
    column_title      = "Cluster_samples",
    top_annotation    = col_annot,
    cluster_rows      = TRUE,
    row_title         = "Modules",
    left_annotation   = NULL,
    show_row_names    = TRUE,
    show_column_names = FALSE,
    width = 15, height = 5
  )

  pdf(file = paste0(out_dir, script_ind, "Module_eigengene_heatmap_by_sample.pdf"),
      width = 10, height = 10)
  { draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom") }
  dev.off()


  #######################################################
  ### plot module eigengene heatmap by sample, columns grouped by TREM2 variant
  #######################################################

  trem2_split  = meta$TREM2Variant[match(colnames(pl_mat_X), meta$cluster_sample)]
  trem2_levels = intersect(c("CV", "R47H", "R62H"), unique(as.character(trem2_split)))
  trem2_split  = factor(as.character(trem2_split), levels = trem2_levels)

  ht_list_trem2 <- Heatmap(
    pl_mat_X,
    name                   = "Z-score",
    col                    = colorRamp2(breaks = seq(from = lims_X[1], to = lims_X[2], length.out = 100),
                                        colors = viridis(100)),
    column_split           = trem2_split,
    cluster_columns        = FALSE,
    cluster_column_slices  = FALSE,
    column_title           = "%s",
    top_annotation         = col_annot,
    cluster_rows           = TRUE,
    row_title              = "Modules",
    left_annotation        = NULL,
    show_row_names         = TRUE,
    show_column_names      = FALSE,
    width = 15, height = 5
  )

  pdf(file = paste0(out_dir, script_ind, "Module_eigengene_heatmap_by_TREM2Variant.pdf"),
      width = 10, height = 10)
  { draw(ht_list_trem2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom") }
  dev.off()


  ######################################################################
  ### Exploratory module-trait correlation screen
  ######################################################################

  message("\n\n   *Compute exploratory module-trait correlation screen \n")

  me_mat  = bd$wgcna$mod_eigengene_mat
  meta_mt = bd$meta
  meta_mt$RNA_counts = apply(bd$counts[, meta_mt$cluster_sample], 2, sum)

  meta_mt = meta_mt[match(colnames(me_mat), meta_mt$cluster_sample), ]

  numeric_traits = c("plaque_dens", "pct4G8PositiveArea", "pctAT8PositiveArea",
                     "pctPHF1PositiveArea", "Braak_numeric", "Age",
                     "PostMortemInterval", "N_cells", "RNA_counts")
  trait_df = meta_mt[, numeric_traits, drop = FALSE]

  for (lv in setdiff(levels(meta_mt$TREM2Variant), "CV")){
    dummy = rep(NA_integer_, nrow(meta_mt))
    dummy[meta_mt$TREM2Variant == "CV"] = 0L
    dummy[meta_mt$TREM2Variant == lv]   = 1L
    trait_df[[paste0("TREM2_", lv)]] = dummy
  }

  trait_df[["APOE4_pos"]]  = as.integer(meta_mt$APOEgroup  == "APOE4-pos")
  trait_df[["CD33var"]]    = as.integer(meta_mt$CD33Group   == "CD33var")
  trait_df[["Region_SSC"]] = as.integer(meta_mt$BrainRegion == "SSC")

  cohort_levels = levels(meta_mt$cohort)
  for (lv in cohort_levels[-1]){
    trait_df[[paste0("Cohort_", lv)]] = as.integer(meta_mt$cohort == lv)
  }

  trait_df = trait_df[, sapply(trait_df, function(x){
    x2 = x[!is.na(x)]
    length(x2) > 1 && var(x2) > 0
  }), drop = FALSE]

  trait_mat = as.matrix(trait_df)

  n_mods      = nrow(me_mat)
  n_traits    = ncol(trait_mat)
  mod_names   = rownames(me_mat)
  trait_names = colnames(trait_mat)

  cor_mat  = matrix(NA, nrow = n_mods, ncol = n_traits, dimnames = list(mod_names, trait_names))
  pval_mat = matrix(NA, nrow = n_mods, ncol = n_traits, dimnames = list(mod_names, trait_names))
  n_mat    = matrix(NA, nrow = n_mods, ncol = n_traits, dimnames = list(mod_names, trait_names))

  for (mod1 in mod_names){
    eg = me_mat[mod1, ]
    for (tr1 in trait_names){
      tr   = trait_mat[, tr1]
      ok   = complete.cases(eg, tr)
      n_ok = sum(ok)
      if (n_ok >= 5){
        r = cor(eg[ok], tr[ok], method = "pearson")
        cor_mat[mod1,  tr1] = r
        n_mat[mod1,    tr1] = n_ok
        pval_mat[mod1, tr1] = corPvalueStudent(r, n_ok)
      }
    }
  }

  padj_mat = matrix(p.adjust(as.vector(pval_mat), method = "BH"),
                    nrow = n_mods, ncol = n_traits,
                    dimnames = list(mod_names, trait_names))

  res_long = expand.grid(module = mod_names, trait = trait_names,
                         stringsAsFactors = FALSE) %>%
    mutate(
      r      = as.vector(cor_mat),
      pvalue = as.vector(pval_mat),
      padj   = as.vector(padj_mat),
      n      = as.vector(n_mat)
    )

  write_csv(res_long,
            file = paste0(out_dir, script_ind, "Module_trait_correlation_results.csv"))

  make_cor_col = function(mat){
    mx = max(abs(mat), na.rm = TRUE)
    colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
  }

  sig_label = function(p_raw, p_adj){
    ifelse(is.na(p_adj), "",
      ifelse(p_adj < 0.001, "***",
        ifelse(p_adj < 0.01, "**",
          ifelse(p_adj < 0.05, "*", ""))))
  }

  cell_labels = matrix(sig_label(as.vector(pval_mat), as.vector(padj_mat)),
                       nrow = n_mods, ncol = n_traits,
                       dimnames = list(mod_names, trait_names))

  ht_cor_full = Heatmap(
    cor_mat,
    name              = "Pearson r",
    col               = make_cor_col(cor_mat),
    cell_fun          = function(j, i, x, y, width, height, fill){
      grid.text(cell_labels[i, j], x, y, gp = gpar(fontsize = 8))
    },
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    row_title         = "Modules",
    column_title      = "Exploratory module-trait screen  (BH-FDR: * <0.05, ** <0.01, *** <0.001)",
    show_row_names    = TRUE,
    show_column_names = TRUE,
    na_col            = "grey90",
    column_names_rot  = 45,
    column_names_gp   = gpar(fontsize = 11),
    row_names_gp      = gpar(fontsize = 11),
    column_title_gp   = gpar(fontsize = 10, fontface = "plain")
  )

  pdf(file   = paste0(out_dir, script_ind, "Module_trait_correlation_heatmap.pdf"),
      width  = max(n_traits * 0.5 + 3, 8),
      height = max(n_mods   * 0.4 + 3, 6))
  { draw(ht_cor_full, heatmap_legend_side = "right") }
  dev.off()

  sig_traits = trait_names[apply(padj_mat, 2, function(p) any(p < 0.05, na.rm = TRUE))]

  if (length(sig_traits) >= 2){
    cor_mat_sig  = cor_mat[,     sig_traits, drop = FALSE]
    cell_lbl_sig = cell_labels[, sig_traits, drop = FALSE]

    ht_cor_sig = Heatmap(
      cor_mat_sig,
      name              = "Pearson r",
      col               = make_cor_col(cor_mat_sig),
      cell_fun          = function(j, i, x, y, width, height, fill){
        grid.text(cell_lbl_sig[i, j], x, y, gp = gpar(fontsize = 8))
      },
      cluster_rows      = FALSE,
      cluster_columns   = FALSE,
      row_title         = "Modules",
      column_title      = "Significant traits  (BH-FDR: * <0.05, ** <0.01, *** <0.001)",
      show_row_names    = TRUE,
      show_column_names = TRUE,
      na_col            = "grey90",
      column_names_rot  = 45,
      column_names_gp   = gpar(fontsize = 11),
      row_names_gp      = gpar(fontsize = 11),
      column_title_gp   = gpar(fontsize = 10, fontface = "plain")
    )

    pdf(file   = paste0(out_dir, script_ind, "Module_trait_correlation_heatmap_sig_traits.pdf"),
        width  = max(length(sig_traits) * 0.5 + 3, 6),
        height = max(n_mods * 0.4 + 3, 6))
    { draw(ht_cor_sig, heatmap_legend_side = "right") }
    dev.off()

  } else {
    message("   *Fewer than 2 BH-FDR significant traits found; skipping filtered heatmap")
  }


  ######################################################################
  ### plot module eigengene vs group boxplots by cluster and module
  ######################################################################

  message("\n\n   *Generate module eigengene vs group plots \n")

  m1  = bd$wgcna$mod_eigengene_mat
  meta = bd$meta
  gr   = unique(meta$group)

  pl     = list()
  pw_res = list()
  kw_res = list()

  fmt_p = function(p){
    if (is.na(p)) return("NA")
    if (p < 0.001) return("<0.001")
    formatC(p, digits = 3, format = "f")
  }

  for (mod1 in rownames(m1)){
    for (cl in cluster_names){

      meta_cl = meta[meta$cluster_name == cl, ]
      t1 = cbind(meta_cl, module = mod1, eigengene = m1[mod1, meta_cl$cluster_sample])

      grp_present = intersect(gr, unique(t1$group))

      kw_p = NA_real_
      if (length(grp_present) >= 2){
        grp_ns = table(factor(t1$group, levels = grp_present))
        if (sum(grp_ns >= 2) >= 2){
          kw_p = suppressWarnings(
            kruskal.test(eigengene ~ factor(group, levels = grp_present), data = t1)$p.value
          )
        }
      }
      kw_res[[length(kw_res) + 1]] = data.frame(
        cluster  = cl, module = mod1,
        n_groups = length(grp_present),
        n_total  = nrow(t1),
        kw_pvalue = kw_p,
        stringsAsFactors = FALSE
      )

      pairs   = combn(grp_present, 2, simplify = FALSE)
      pair_df = data.frame(
        cluster = cl, module = mod1,
        g1 = sapply(pairs, `[`, 1), g2 = sapply(pairs, `[`, 2),
        n1 = NA_integer_, n2 = NA_integer_,
        pvalue = NA_real_, stringsAsFactors = FALSE
      )
      for (pi in seq_along(pairs)){
        x = t1$eigengene[t1$group == pairs[[pi]][1]]
        y = t1$eigengene[t1$group == pairs[[pi]][2]]
        pair_df$n1[pi] = length(x)
        pair_df$n2[pi] = length(y)
        if (length(x) >= 2 && length(y) >= 2){
          pair_df$pvalue[pi] = suppressWarnings(
            wilcox.test(x, y, exact = FALSE)$p.value
          )
        }
      }
      pw_res[[length(pw_res) + 1]] = pair_df

      pl[[paste0(mod1, "_", cl)]] = ggplot(t1, aes(x = group, y = eigengene)) +
        geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.6,
                     color = "grey30", width = 0.6) +
        geom_point(aes(fill = group), shape = 21,
                   position = position_jitter(width = 0.15, height = 0),
                   size = 1, stroke = 0.3, color = "grey30") +
        geom_hline(yintercept = 0) +
        scale_x_discrete(limits = gr) +
        scale_color_manual(limits = gr, values = pal(gr)) +
        scale_fill_manual(limits = gr, values = pal(gr)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              legend.position = "none") +
        labs(title = paste0(mod1, " - ", cl, " - Module eigengene by group"))
    }
  }

  kw_long = do.call(rbind, kw_res)
  kw_long$kw_padj = p.adjust(kw_long$kw_pvalue, method = "BH")
  write_csv(kw_long,
            file = paste0(out_dir, script_ind, "Module_eigengene_vs_group_kruskal_wallis.csv"))

  pw_long = do.call(rbind, pw_res)
  pw_long$padj = p.adjust(pw_long$pvalue, method = "BH")
  write_csv(pw_long,
            file = paste0(out_dir, script_ind, "Module_eigengene_vs_group_pairwise_wilcoxon.csv"))

  pdf(file = paste0(out_dir, script_ind, "Module_eigengene_vs_group_scatterplot_by_cluster.pdf"),
      width = 6, height = 3)
  { lapply(pl, function(x){ x }) }
  dev.off()


  ######################################################################
  ### Across-cluster LMM per module (TREM2Variant effect)
  ######################################################################

  # Variant associations are descriptive because modules are seeded on TREM2 DEGs.

  message("\n\n          *** Across-cluster LMM per module (TREM2Variant) - ", Sys.time(), "\n\n")

  me_mat   = bd$wgcna$mod_eigengene_mat
  meta_lmm = bd$meta
  meta_lmm = meta_lmm[match(colnames(me_mat), meta_lmm$cluster_sample), ]

  lmm_omni_list  = list()
  lmm_pairs_list = list()

  for (mod1 in rownames(me_mat)){

    df1 = data.frame(
      eigengene    = as.numeric(me_mat[mod1, ]),
      TREM2Variant = factor(meta_lmm$TREM2Variant),
      cluster_name = factor(meta_lmm$cluster_name),
      sample       = factor(meta_lmm$sample),
      stringsAsFactors = FALSE
    )
    df1 = df1[complete.cases(df1), ]
    df1$TREM2Variant = droplevels(df1$TREM2Variant)
    df1$cluster_name = droplevels(df1$cluster_name)
    df1$sample       = droplevels(df1$sample)

    if (nlevels(df1$TREM2Variant) < 2) next

    fit = tryCatch(
      lmerTest::lmer(eigengene ~ TREM2Variant + (1 | cluster_name) + (1 | sample),
                     data = df1, REML = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    aov1 = tryCatch(anova(fit, type = 2), error = function(e) NULL)
    get_aov = function(col){
      if (!is.null(aov1) && "TREM2Variant" %in% rownames(aov1) && col %in% colnames(aov1)){
        aov1["TREM2Variant", col]
      } else NA_real_
    }

    lmm_omni_list[[mod1]] = data.frame(
      module     = mod1,
      n_obs      = nrow(df1),
      n_clusters = nlevels(df1$cluster_name),
      n_samples  = nlevels(df1$sample),
      F_value    = get_aov("F value"),
      num_df     = get_aov("NumDF"),
      den_df     = get_aov("DenDF"),
      pvalue     = get_aov("Pr(>F)"),
      stringsAsFactors = FALSE
    )

    lv       = levels(df1$TREM2Variant)
    ref      = lv[1]
    fe_names = names(fixef(fit))
    n_fe     = length(fe_names)
    all_pairs = combn(lv, 2, simplify = FALSE)

    pairs_rows = list()
    for (pp in all_pairs){
      g1 = pp[1]; g2 = pp[2]
      L  = rep(0, n_fe)
      if (g1 != ref) L[which(fe_names == paste0("TREM2Variant", g1))] = -1
      if (g2 != ref) L[which(fe_names == paste0("TREM2Variant", g2))] =  1

      ct = tryCatch(
        lmerTest::contest(fit, L = matrix(L, nrow = 1), joint = FALSE),
        error = function(e) NULL
      )
      if (is.null(ct) || nrow(ct) == 0) next

      pairs_rows[[length(pairs_rows) + 1]] = data.frame(
        module   = mod1,
        contrast = paste0(g2, " - ", g1),
        estimate = ct[1, "Estimate"],
        SE       = ct[1, "Std. Error"],
        df       = ct[1, "df"],
        t.ratio  = ct[1, "t value"],
        p.value  = ct[1, "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
    }
    if (length(pairs_rows) > 0){
      lmm_pairs_list[[mod1]] = do.call(rbind, pairs_rows)
    }
  }

  lmm_omni  = do.call(rbind, lmm_omni_list)
  lmm_pairs = do.call(rbind, lmm_pairs_list)

  if (!is.null(lmm_omni) && nrow(lmm_omni) > 0){
    lmm_omni$padj_BH = p.adjust(lmm_omni$pvalue, method = "BH")
    write_csv(lmm_omni,
              file = paste0(out_dir, script_ind, "Module_eigengene_LMM_across_clusters_omnibus.csv"))
  }

  if (!is.null(lmm_pairs) && nrow(lmm_pairs) > 0){
    lmm_pairs$padj_BH = p.adjust(lmm_pairs$p.value, method = "BH")
    write_csv(lmm_pairs,
              file = paste0(out_dir, script_ind, "Module_eigengene_LMM_across_clusters_pairwise.csv"))

    fp = lmm_pairs
    fp$ci_lo     = fp$estimate - 1.96 * fp$SE
    fp$ci_hi     = fp$estimate + 1.96 * fp$SE
    fp$sig_label = ifelse(is.na(fp$padj_BH), "",
                   ifelse(fp$padj_BH < 0.001, "***",
                   ifelse(fp$padj_BH < 0.01,  "**",
                   ifelse(fp$padj_BH < 0.05,  "*", ""))))

    p_fp = ggplot(fp, aes(x = estimate, y = module, color = contrast)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25,
                     position = position_dodge(width = 0.6)) +
      geom_point(position = position_dodge(width = 0.6), size = 2) +
      geom_text(aes(label = sig_label),
                position = position_dodge(width = 0.6),
                hjust = -0.4, size = 4, show.legend = FALSE) +
      scale_y_discrete(limits = rev(rownames(me_mat))) +
      labs(title    = paste0("Across-cluster LMM: TREM2 variant contrasts per module (", model_name, ")"),
           subtitle = "Estimate ± 95% CI; BH-FDR: * <0.05, ** <0.01, *** <0.001",
           x = "Eigengene difference (LMM estimate)", y = "Module") +
      theme_classic()

    pdf(file   = paste0(out_dir, script_ind, "Module_eigengene_LMM_across_clusters_forestplot.pdf"),
        width  = 8,
        height = max(4, length(unique(fp$module)) * 0.4 + 2))
    plot(p_fp)
    dev.off()
  }


  #################################################
  ### GO-BP over-representation analysis of module genes
  #################################################

  GO_list        = list()
  GO_results_tab = NULL

  N_comps = length(bd$wgcna$mod_genes)

  for (i in 1:N_comps){

    comp = names(bd$wgcna$mod_genes)[i]

    message("\n          *** GO analysis: module ", comp,
            " (", i, " of ", N_comps, ") - ", Sys.time(), "\n")

    ego      = NULL
    go_genes = bd$wgcna$mod_genes[[comp]]

    if (length(go_genes) > 2){
      ego = enrichGO(gene          = go_genes,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
      GO_list[[comp]] = ego
    }

    if (!is.null(ego)){
      t1 = ego@result[ego@result$p.adjust <= 0.05, ]
      if (nrow(t1) > 0){
        t2 = cbind(module = comp, t1)
        GO_results_tab = rbind(GO_results_tab, t2)
      }
    }
  }

  GO_results_tab = GO_results_tab[GO_results_tab$Count > 1, ]

  write_csv(GO_results_tab, file = paste0(out_dir, script_ind, "GO_results_by_comp.csv"))

  bd$GO_results$by_comp_GO_list = GO_list
  bd$GO_results$by_comp_GO_res  = GO_results_tab

  qsave(bd, file = paste0(out_dir, script_ind, "bulk_data.qs"))


  ############################################################################
  ### GO dotplot: top 10 terms per module
  ############################################################################

  t1 = GO_results_tab
  t2 = NULL

  for (mod1 in unique(t1$module)){
    t3 = t1[t1$module == mod1, ]
    if (nrow(t3) > 10){ t3 = t3[1:10, ] }
    t2 = rbind(t2, t3)
  }

  mod_colors = unique(mod_gene_tab$colors)

  p1 = ggplot(t2, aes(x = module, y = Description, size = Count, colour = module)) +
    geom_point() +
    scale_color_manual(limits = mods, values = mod_colors) +
    scale_size_continuous(limits = c(0, max(t2$Count))) +
    scale_x_discrete(limits = mods) +
    scale_y_discrete(limits = unique(t2$Description)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste0("Top 10 GO terms per module (", model_name, ")"))

  pdf(file = paste0(out_dir, script_ind, "GO_results_by_module_dotplot_top_terms.pdf"),
      width = 11, height = 12)
  plot(p1)
  dev.off()


  message("\n\n##########################################################################\n",
          "# Completed LD_F03d v01 - Model ", model_name, " - ", Sys.time(),
          "\n##########################################################################\n\n")

} # end model loop


sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed LD_F03d v01 (all models) ", Sys.time(),
        "\n##########################################################################\n\n\n")
