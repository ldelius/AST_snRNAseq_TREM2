message("\n\n##########################################################################\n",
        "# Start G01b: Merge AST and MIC labelled Seurat objects ", Sys.time(),
        "\n##########################################################################\n",
        "\n   Loads the labelled AST Seurat (LD_B04a_v02) and Michael's labelled MIC Seurat (B04),\n",
        "   into a single object, runs NormalizeData() on the merged RNA assay, and saves \n",
        "   for downstream CellChat analysis (G02).\n",
        "\n##########################################################################\n\n")


#set environment/load packages
library(qs)
library(tidyverse)
library(Seurat)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "LD_G01b_"

#specify output directory
out_dir = paste0(main_dir, "LD_G_MIC_AST_communication_analysis_output/")
if (!dir.exists(out_dir)){dir.create(out_dir)}



###########################################################
# 1. load AST and MIC labelled Seurat objects
###########################################################

message("\n\n          *** Load labelled Seurat objects... ", Sys.time(), "\n\n")

seur_ast = qread(file = paste0(main_dir, "LD_B_AST_analysis_output/LD_B04a_v02_seur.qs"))
seur_mic = qread(file = paste0(main_dir, "data_TREM2_michael/B_load_from_scflow_subcluster/B04_seur.qs"))



###########################################################
# 2. verify RNA data not yet log-normalised in either Seurat
###########################################################

# CellChat expects log-normalised counts in the RNA data layer. Our upstream
# pipelines used SCT (which lives in a separate assay), so the RNA data layer
# should be empty here. We will populate it with NormalizeData() after merging,
# so MIC + AST cells are normalised in one pass.
# Michael used SCT in his pipeline, however, CellChatwebsite recommends NormalizeData()

message("\n\n          *** Check RNA data layer state... ", Sys.time(), "\n\n")


check_rna_data_empty = function(seur, label){
  data_mat = tryCatch(GetAssayData(seur, assay = "RNA", layer = "data"),
                      warning = function(w) NULL,
                      error   = function(e) NULL)
  is_empty = is.null(data_mat) || length(data_mat) == 0 || nrow(data_mat) == 0
  cat(label, ": RNA data layer empty (= not yet normalised)? ", is_empty, "\n", sep = "")
  return(is_empty)
}

if (!check_rna_data_empty(seur_ast, "AST") || !check_rna_data_empty(seur_mic, "MIC")){
  stop("RNA data layer is already populated in one or both Seurats. ",
       "NormalizeData() may have been run upstream - decide before proceeding.")
}



###########################################################
# 3. sanity checks before merging
###########################################################

# Two sanity checks before stitching the two Seurats together:
#  - does each side actually contain only its expected cell type?
#  - do patient-level metadata agree for samples present in both?

message("\n\n          *** Sanity checks... ", Sys.time(), "\n\n")
 
 
### confirm AST contains only Astro, MIC only Micro
if (!all(seur_ast$cluster_celltype == "Astro")){
  stop("AST Seurat contains non-Astro cells.")
}
if (!all(seur_mic$cluster_celltype == "Micro")){
  stop("MIC Seurat contains non-Micro cells.")
}
 
 
### patient-level metadata agreement
covar_cols = c("group", "sample", "SampleID", "StudyID",
               "Age", "Sex",
               "APOE", "APOEgroup", "TREM2Variant",
               "BrainRegion", "BrainBank", "Braak",
               "NeuropathologicalDiagnosis",
               "CD33Group", "cohort", "manifest", "plaque_dens")
 
# only compare columns present in BOTH Seurats
covar_cols_check = intersect(intersect(covar_cols, colnames(seur_ast@meta.data)),
                             colnames(seur_mic@meta.data))
 
cat("\nCovariates compared between AST and MIC: ",
    paste(covar_cols_check, collapse = ", "), "\n")
cat("Not compared (missing in one or both): ",
    paste(setdiff(covar_cols, covar_cols_check), collapse = ", "), "\n")
 
if (!"sample" %in% covar_cols_check){
  stop("'sample' missing from one or both Seurats - cannot align metadata.")
}
 
 
### one row per sample on each side
t_ast = seur_ast@meta.data %>% distinct(sample, .keep_all = TRUE) %>% select(all_of(covar_cols_check))
t_mic = seur_mic@meta.data %>% distinct(sample, .keep_all = TRUE) %>% select(all_of(covar_cols_check))
 
shared_samples = intersect(t_ast$sample, t_mic$sample)
cat("\nN samples - AST only:", length(setdiff(t_ast$sample, t_mic$sample)),
    " | MIC only:", length(setdiff(t_mic$sample, t_ast$sample)),
    " | shared:", length(shared_samples), "\n")
 
 
### compare each covariate row-by-row across shared samples
mismatches = list()
for (s in shared_samples){
  r1 = t_ast[t_ast$sample == s, ]
  r2 = t_mic[t_mic$sample == s, ]
  for (cov in setdiff(covar_cols_check, "sample")){
    if (!identical(r1[[cov]], r2[[cov]])){
      mismatches[[length(mismatches) + 1]] =
        tibble(sample = s, covar = cov, AST = r1[[cov]], MIC = r2[[cov]])
    }
  }
}
 
if (length(mismatches) > 0){
  cat("\n>>> Covariate mismatches: <<<\n")
  print(bind_rows(mismatches))
  stop("Covariate mismatch between AST and MIC. Resolve before merging.")
} else {
  cat("\n>>> All compared covariates agree across shared samples. <<<\n")
}



###########################################################
# 4. tag each cell with its source cell type, then merge
###########################################################

# Adding cell_type_joint as an explicit column with AST and MIC, so that downstream
# we can filter cells by source 

seur_ast$cell_type_joint = "AST"
seur_mic$cell_type_joint = "MIC"


message("\n\n          *** Merge Seurat objects... ", Sys.time(), "\n\n")


# Seurat::merge() concatenates cells. Cell barcodes get prefixed with
# add.cell.ids to avoid collisions. Integration outputs from each side
# (Harmony, joint UMAP) are dropped. CellChat doesn't use these so it doesn't matter here.

seur = merge(seur_ast, y = seur_mic,
             add.cell.ids = c("AST", "MIC"),
             project = "MIC_AST_joint")


# In Seurat v5, merging two Assay5 objects keeps source-specific layers
# (counts.1, counts.2). JoinLayers() collapses them into a single 'counts'
# layer so downstream tools see one consistent layer per type.

seur = JoinLayers(seur, assay = "RNA")


### add cluster_sample column (cluster x patient sample)
seur$cluster_sample = paste0(seur$cluster_name, "_", seur$sample)


cat("\nMerged object - cells per cell type:\n")
print(table(seur$cell_type_joint))
cat("\nRNA layers after JoinLayers:", paste(Layers(seur, assay = "RNA"), collapse = ", "), "\n")


### free memory
rm(seur_ast, seur_mic); gc()



###########################################################
# 5. NormalizeData() on the merged RNA assay
###########################################################

# NormalizeData() is per-cell (Per-cell: each cell's counts get divided by that cell's library size,
# multiplied by 10,000, plus 1, log-transformed.). Running on the merged object gives consistent
# normalisation for AST + MIC in one procedure.

message("\n\n          *** NormalizeData on merged RNA assay... ", Sys.time(), "\n\n")


DefaultAssay(seur) = "RNA"
seur = NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)



###########################################################
# 6. set Idents to cluster_name (CellChat reads from Idents by default)
###########################################################

Idents(seur) = "cluster_name"



###########################################################
# 7. diagnostics for downsampling decision in G02 (cells per cluster_sample ())
###########################################################

# running per cluster summary csv, cap vs totalcell csv, histogram

message("\n\n          *** Cells per cluster_sample diagnostics... ", Sys.time(), "\n\n")


t_csize = seur@meta.data %>%
  count(cell_type_joint, cluster_name, sample, name = "n_cells") %>%
  arrange(desc(n_cells))


### per-cluster summary
t_summary = t_csize %>%
  group_by(cell_type_joint, cluster_name) %>%
  summarise(n_samples    = n(),
            min_cells    = min(n_cells),
            q25_cells    = quantile(n_cells, 0.25),
            median_cells = median(n_cells),
            q75_cells    = quantile(n_cells, 0.75),
            max_cells    = max(n_cells),
            total_cells  = sum(n_cells),
            .groups = "drop") %>%
  arrange(cell_type_joint, cluster_name)


### cap-vs-total-cells trade-off table
caps = c(50, 75, 100, 150, 200, 250, 300, 500)
t_caps = tibble(cap = caps) %>%
  rowwise() %>%
  mutate(total_cells_after_cap = sum(pmin(t_csize$n_cells, cap)),
         n_clusters_capped     = sum(t_csize$n_cells > cap)) %>%
  ungroup()


cat("\nPer-cluster summary:\n");                 print(t_summary, n = Inf)
cat("\nTotal cells after downsampling at candidate caps:\n");  print(t_caps)


### save tables
write_csv(t_summary, paste0(out_dir, script_ind, "cluster_sample_size_summary.csv"))
write_csv(t_csize,   paste0(out_dir, script_ind, "cluster_sample_size_full.csv"))
write_csv(t_caps,    paste0(out_dir, script_ind, "cluster_sample_size_at_caps.csv"))


### histogram, faceted by cluster_name
pl_hist = ggplot(t_csize, aes(x = n_cells, fill = cell_type_joint)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = c(100, 200), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  facet_wrap(~ cluster_name, scales = "free_y") +
  scale_fill_manual(values = c("AST" = "dodgerblue", "MIC" = "orange")) +
  labs(x = "Cells per cluster_sample",
       y = "N samples",
       title = "Distribution of cluster_sample sizes",
       subtitle = "Dashed lines: candidate downsample caps (100, 200)") +
  theme_bw() +
  theme(strip.text = element_text(size = 7))

n_clusters = length(unique(t_csize$cluster_name))
n_cols = ceiling(sqrt(n_clusters))
n_rows = ceiling(n_clusters / n_cols)

pdf(file = paste0(out_dir, script_ind, "cluster_sample_size_histograms.pdf"),
    width = max(8, n_cols * 2), height = max(6, n_rows * 1.8))
print(pl_hist)
dev.off()



###########################################################
# 8. save merged Seurat
###########################################################

message("\n\n          *** Save merged Seurat... ", Sys.time(), "\n\n")


qsave(seur, file = paste0(out_dir, script_ind, "seur_merged.qs"))


message("\n\n##########################################################################\n",
        "# Finished G01b ", Sys.time(),
        "\n##########################################################################\n\n")