message("\n\n##########################################################################\n",
        "# Start G01a: Merge AST and MIC pseudobulk objects ", Sys.time(),
        "\n##########################################################################\n",
        "\n   Loads the AST pseudobulk object (LD_F03a1_v02_bulk_data.qs) and Michael's\n",
        "   MIC equivalent (F03a1_bulk_data.qs), concatenates raw counts on shared\n",
        "   genes, runs joint VST + sequential limma::removeBatchEffect for covariates\n",
        "   (cohort, APOEgroup, CD33Group, BrainRegion), and computes per-gene Z-scores\n",
        "   across all cluster_samples.\n",
        "\n   Output is used by G03a/b/c pathway-gene heatmap sections to display ligand\n",
        "   expression in MIC clusters and receptor expression in AST clusters in one\n",
        "   heatmap, annotated by TREM2Variant x NeuropathologicalDiagnosis.\n",
        "\n##########################################################################\n\n")


#set environment/load packages
library(qs)
library(tidyverse)
library(DESeq2)
library(limma)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "LD_G01a_"

#specify output directory
out_dir = paste0(main_dir, "LD_G_MIC_AST_communication_analysis_output/")
if (!dir.exists(out_dir)){dir.create(out_dir)}


#covariates to correct vst matrix for (matches LD_E02a2_v02 covariate set)
covars_corr = c("cohort", "APOEgroup", "CD33Group", "BrainRegion")



###########################################################
# 1. load AST and MIC pseudobulk objects
###########################################################

# F03a1 outputs are used (rather than earlier E02a2) so WGCNA module fields are
# available on the AST side if any downstream G03 section needs them. The
# matrices we use (counts, meta) are identical across E02a2/E03a2/E03b2/F02a1/
# F03a1 - they're computed once in E02a2 and carried through.

message("\n\n          *** Load pseudobulk objects... ", Sys.time(), "\n\n")

bulk_ast = qread(file = paste0(main_dir,
                               "LD_F_DESeq_pseudobulk_WGCNA/LD_F03a1_v02_bulk_data.qs"))
bulk_mic = qread(file = paste0(main_dir,
                               "data_TREM2_michael/F_DESeq_pseudobulk_WGCNA/F03a1_bulk_data.qs"))


cat("AST: counts dim =", paste(dim(bulk_ast$counts), collapse = " x "),
    "  meta nrow =", nrow(bulk_ast$meta), "\n")
cat("MIC: counts dim =", paste(dim(bulk_mic$counts), collapse = " x "),
    "  meta nrow =", nrow(bulk_mic$meta), "\n")



###########################################################
# 2. align gene sets across AST and MIC count matrices
###########################################################

# Both pipelines start from the same upstream Seurat (Michael's snRNA-seq), so
# in principle the gene sets should match. Intersecting defensively in case
# upstream filtering (e.g. low-expression filter in E01) differed between
# lineages, and reporting how many genes are dropped on each side.

message("\n\n          *** Align gene sets... ", Sys.time(), "\n\n")


shared_genes = intersect(rownames(bulk_ast$counts), rownames(bulk_mic$counts))

cat("\nGenes in AST counts:    ", nrow(bulk_ast$counts), "\n",
    "Genes in MIC counts:    ", nrow(bulk_mic$counts), "\n",
    "Shared genes:           ", length(shared_genes), "\n",
    "AST-only genes dropped: ", nrow(bulk_ast$counts) - length(shared_genes), "\n",
    "MIC-only genes dropped: ", nrow(bulk_mic$counts) - length(shared_genes), "\n",
    sep = "")


### subset count matrices to shared genes
counts_ast = bulk_ast$counts[shared_genes, , drop = FALSE]
counts_mic = bulk_mic$counts[shared_genes, , drop = FALSE]



###########################################################
# 3. sanity checks before merging
###########################################################

message("\n\n          *** Sanity checks... ", Sys.time(), "\n\n")


### required metadata columns - hard-fail if any missing from either side
covar_cols = c("cluster_sample", "cluster_name", "sample",
               "TREM2Variant", "NeuropathologicalDiagnosis",
               "BrainRegion", "APOEgroup", "CD33Group", "cohort")
 
covar_cols_check = intersect(intersect(covar_cols, colnames(bulk_ast$meta)),
                             colnames(bulk_mic$meta))
 
missing_required = setdiff(covar_cols, covar_cols_check)
if (length(missing_required) > 0){
  stop("Required metadata columns missing from one or both pseudobulk objects: ",
       paste(missing_required, collapse = ", "))
}
 
 
### one row per sample on each side, joined so AST/MIC values sit side by side
# ungroup() defensively: upstream meta tables may carry grouping (e.g. by
# cell_type) which would make distinct(sample, .keep_all = TRUE) return one
# row per sample x group, not one row per sample
t_ast = bulk_ast$meta %>% ungroup() %>%
  distinct(sample, .keep_all = TRUE) %>% select(all_of(covar_cols_check))
t_mic = bulk_mic$meta %>% ungroup() %>%
  distinct(sample, .keep_all = TRUE) %>% select(all_of(covar_cols_check))
 
joined = inner_join(t_ast, t_mic, by = "sample", suffix = c("_ast", "_mic"))
 
cat("\nN samples - AST only:", length(setdiff(t_ast$sample, t_mic$sample)),
    " | MIC only:", length(setdiff(t_mic$sample, t_ast$sample)),
    " | shared:", nrow(joined), "\n")
 
 
### compare each covariate column-pair, NA-safe
# disagree if values differ, or if exactly one side is NA (xor)
# both NA -> agree (mism = NA -> set to FALSE)
# cast to character so factor level-order differences don't trigger spurious
# mismatches between AST and MIC pipelines
mismatched_covars = c()
for (cov in setdiff(covar_cols_check, c("sample", "cluster_sample", "cluster_name"))){
  a = as.character(joined[[paste0(cov, "_ast")]])
  b = as.character(joined[[paste0(cov, "_mic")]])
  
  mism = (a != b) | xor(is.na(a), is.na(b))
  mism[is.na(mism)] = FALSE
  
  if (any(mism)){
    cat("Mismatch in", cov, "for samples:", paste(joined$sample[mism], collapse = ", "), "\n")
    mismatched_covars = c(mismatched_covars, cov)
  }
}
 
if (length(mismatched_covars) > 0){
  stop("Covariate mismatch in: ", paste(mismatched_covars, collapse = ", "),
       ". Resolve before merging.")
} else {
  cat("\n>>> All compared covariates agree across shared samples. <<<\n")
}



###########################################################
# 4. tag each cluster_sample with its source lineage, then merge
###########################################################

# Adding lineage as an explicit column (AST / MIC), so downstream we can filter
# heatmaps by lineage (e.g. show ligand expression in MIC clusters and
# receptor expression in AST clusters).

message("\n\n          *** Concatenate counts and metadata... ", Sys.time(), "\n\n")


meta_ast = bulk_ast$meta[, covar_cols_check, drop = FALSE]
meta_mic = bulk_mic$meta[, covar_cols_check, drop = FALSE]

meta_ast$lineage = "AST"
meta_mic$lineage = "MIC"


### cbind() requires identical row order (genes), rbind() requires identical columns;
### already ensured above by the shared_genes subset and covar_cols_check intersect
counts_comb = cbind(counts_ast, counts_mic)
meta_comb   = rbind(meta_ast, meta_mic)
rownames(meta_comb) = meta_comb$cluster_sample


### enforce column order match between counts and meta
counts_comb = counts_comb[, meta_comb$cluster_sample]
stopifnot(identical(colnames(counts_comb), meta_comb$cluster_sample))
stopifnot(identical(rownames(counts_comb), shared_genes))


### drop cluster_samples with NA in any correction covariate
# DESeq2 cannot handle NAs in design covariates; matches LD_E02a2 behaviour.
keep_cs = complete.cases(meta_comb[, covars_corr])
n_dropped = sum(!keep_cs)
if (n_dropped > 0){
  cat("\nDropping ", n_dropped, " cluster_samples with NA in correction covariates\n", sep = "")
  cat("NA breakdown by covariate:\n")
  for (cv in covars_corr){
    cat("  ", cv, ": ", sum(is.na(meta_comb[[cv]])), " NA\n", sep = "")
  }
}

counts_comb = counts_comb[, keep_cs]
meta_comb   = meta_comb[keep_cs, , drop = FALSE]


### factorise covariates to match levels Michael uses upstream (LD_E02a2 / F02)
meta_comb$cohort                     = factor(meta_comb$cohort,
                                              levels = c("BiogenInitial", "BiogenExtra"))
meta_comb$APOEgroup                  = factor(meta_comb$APOEgroup,
                                              levels = c("APOE4-neg", "APOE4-pos"))
meta_comb$CD33Group                  = factor(meta_comb$CD33Group,
                                              levels = c("CV", "CD33var"))
meta_comb$BrainRegion                = factor(meta_comb$BrainRegion,
                                              levels = c("SSC", "MTG"))
meta_comb$NeuropathologicalDiagnosis = factor(meta_comb$NeuropathologicalDiagnosis,
                                              levels = c("Control", "AD"))
meta_comb$TREM2Variant               = factor(meta_comb$TREM2Variant,
                                              levels = c("CV", "R47H", "R62H"))


cat("\nMerged dataset - cluster_samples per lineage:\n")
print(table(meta_comb$lineage))


### free memory
rm(bulk_ast, bulk_mic, counts_ast, counts_mic, meta_ast, meta_mic); gc()



###########################################################
# 5. joint DESeq2 dataset, VST normalisation
###########################################################

# Fresh VST on the merged raw counts - dispersion is estimated jointly across
# AST and MIC samples, so VST values are on the same scale across lineages
#
# design = ~cluster_name informs vst()'s variance estimation - no DE test is
# run here, the merged pseudobulk is purely for visualisation.
# blind = FALSE matches LD_E02a2 behaviour (uses design for dispersion).

message("\n\n          *** Build joint DESeq2 dataset and run VST... ", Sys.time(), "\n\n")


dds = DESeqDataSetFromMatrix(countData = counts_comb,
                             colData   = meta_comb,
                             design    = ~ cluster_name)

vst_obj = vst(dds, blind = FALSE)
vst_mat = assay(vst_obj)

cat("VST matrix dim: ", paste(dim(vst_mat), collapse = " x "), "\n", sep = "")



###########################################################
# 6. batch-correct VST matrix for covariates
###########################################################

# limma::removeBatchEffect removes covariate-driven variance from the VST matrix
# group = cluster_name preserves cluster-level expression structure during
# correction (the lineage / disease / variant signal we want to keep).
# Note: sequential correction is an approximation - joint correction via the `covariates`
# argument with a model.matrix would be statistically cleaner but would diverge
# from how every other corrected VST in this project was computed.

message("\n\n          *** Batch-correct VST matrix... ", Sys.time(), "\n\n")


vst_mat_corr = vst_mat

for (cov1 in covars_corr){
  if (length(unique(meta_comb[[cov1]])) > 1){
    cat("    Correcting for: ", cov1, "\n", sep = "")
    vst_mat_corr = limma::removeBatchEffect(vst_mat_corr,
                                            batch = meta_comb[[cov1]],
                                            group = meta_comb$cluster_name)
  } else {
    cat("    Skipping (only 1 level): ", cov1, "\n", sep = "")
  }
}



###########################################################
# 7. compute per-gene Z-scores across all cluster_samples
###########################################################

# Z-scores are what gets plotted in the pathway-gene heatmaps - one row per
# gene, scaled across all cluster_samples (AST + MIC together). Computed from
# both uncorrected and corrected VST so heatmaps can show either.
# Constant genes (zero variance) produce NaN rows; pheatmap renders as grey.

message("\n\n          *** Compute per-gene Z-scores... ", Sys.time(), "\n\n")


z_uncorr = t(apply(vst_mat,      1, scale))
z_corr   = t(apply(vst_mat_corr, 1, scale))

colnames(z_uncorr) = colnames(vst_mat)
colnames(z_corr)   = colnames(vst_mat_corr)


###########################################################
# 8. diagnostics: summary tables + library-size histogram
###########################################################

# small CSVs and a histogram summarising the merged dataset

message("\n\n          *** Diagnostic summaries... ", Sys.time(), "\n\n")

### per-cluster summary: how many cluster_samples per cluster, by lineage
t_summary = meta_comb %>%
  dplyr::count(lineage, cluster_name, name = "n_cluster_samples") %>%
  arrange(lineage, desc(n_cluster_samples))
 
cat("\nPer-cluster summary:\n");  print(t_summary, n = Inf)
 
write_csv(t_summary, paste0(out_dir, script_ind, "cluster_sample_summary_by_lineage.csv"))
 
 
### histogram: distribution of per-cluster_sample library sizes by lineage
t_libs = tibble(cluster_sample = colnames(counts_comb),
                lib_size       = colSums(counts_comb),
                lineage        = meta_comb$lineage)
 
pl_hist = ggplot(t_libs, aes(x = lib_size, fill = lineage)) +
  geom_histogram(bins = 40) +
  scale_x_log10() +
  scale_fill_manual(values = c("AST" = "dodgerblue", "MIC" = "orange")) +
  labs(x = "Pseudobulk library size (log10)",
       y = "N cluster_samples",
       title = "Pseudobulk library size distribution by lineage") +
  theme_bw()
 
pdf(file = paste0(out_dir, script_ind, "pseudobulk_lib_size_histogram.pdf"),
    width = 7, height = 4)
print(pl_hist)
dev.off()



###########################################################
# 9. assemble output object and save
###########################################################

# Output structure mirrors LD_F03a1_v02_bulk_data.qs so G03a/b/c can use
# bulk_data$gene_Z_scores$clusters_combined etc. with familiar syntax.
# DESeq2 / WGCNA / GSEA fields from the source AST/MIC objects are NOT carried
# through - this object is for visualisation only; G03 sections that need
# DEG/WGCNA info read those directly from LD_F03a1_v02_bulk_data.qs.

message("\n\n          *** Save merged pseudobulk object... ", Sys.time(), "\n\n")


bulk_data = list(
  counts                = counts_comb,
  meta                  = meta_comb,
  vst_mat_uncorr        = vst_mat,
  vst_mat               = vst_mat_corr,
  gene_Z_scores_uncorr  = list(clusters_combined = z_uncorr),
  gene_Z_scores         = list(clusters_combined = z_corr),
  covars_corr           = covars_corr,
  shared_genes          = shared_genes
)


qsave(bulk_data, file = paste0(out_dir, script_ind, "bulk_data_AST_MIC.qs"))


message("\n\n##########################################################################\n",
        "# Finished G01a ", Sys.time(),
        "\n##########################################################################\n\n")