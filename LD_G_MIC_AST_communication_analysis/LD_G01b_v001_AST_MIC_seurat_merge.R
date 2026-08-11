message("\n\n##########################################################################\n",
        "# Start LD_G01b: Merge AST and MIC Seurat objects ", Sys.time(),
        "\n##########################################################################\n\n")

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

# For CellChat I want to use log-normalised counts in the RNA data layer as recommended by CellChat website.
# Our upstream pipelines used SCT, so the RNA data layer should be empty here.
# However, double checking here.

message("\n\n          *** Check RNA data layer state... ", Sys.time(), "\n\n")

stopifnot(length(GetAssayData(seur_ast, assay = "RNA", layer = "data")) == 0)
stopifnot(length(GetAssayData(seur_mic, assay = "RNA", layer = "data")) == 0)



###########################################################
# 3. sanity checks before merging
###########################################################

# Two sanity checks before stitching the two Seurats together:
#  - does each side actually contain only its expected cell type?
#  - do patient-level metadata agree for samples present in both?

message("\n\n          *** Sanity checks... ", Sys.time(), "\n\n")
 

### patient-level metadata agreement
covar_cols = c("group", "sample", "SampleID", "StudyID",
               "Age", "Sex",
               "APOE", "APOEgroup", "TREM2Variant",
               "BrainRegion", "BrainBank", "Braak",
               "NeuropathologicalDiagnosis",
               "CD33Group", "cohort", "manifest", "plaque_dens")
 
# only compare columns present in both Seurats
covar_cols_check = intersect(intersect(covar_cols, colnames(seur_ast@meta.data)),
                             colnames(seur_mic@meta.data))
 
cat("\nCovariates compared between AST and MIC: ", paste(covar_cols_check, collapse = ", "), "\n")
cat("Not compared (missing in one or both): ",
    paste(setdiff(covar_cols, covar_cols_check), collapse = ", "), "\n")
 
if (!"sample" %in% covar_cols_check) stop("'sample' missing from one or both Seurats.")
 
 
### one row per sample on each side, joined so AST/MIC values sit side by side
t_ast = seur_ast@meta.data %>% distinct(sample, .keep_all = TRUE) %>% select(all_of(covar_cols_check))
t_mic = seur_mic@meta.data %>% distinct(sample, .keep_all = TRUE) %>% select(all_of(covar_cols_check))
 
joined = inner_join(t_ast, t_mic, by = "sample", suffix = c("_ast", "_mic"))
 
cat("\nN samples - AST only:", length(setdiff(t_ast$sample, t_mic$sample)),
    " | MIC only:", length(setdiff(t_mic$sample, t_ast$sample)),
    " | shared:", nrow(joined), "\n")
 
 
### compare each covariate column-pair, NA-safe
# disagree if values differ, OR if exactly one side is NA (xor)
# both NA -> agree (mism = NA -> set to FALSE)
mismatched_covars = c()
for (cov in setdiff(covar_cols_check, "sample")){
  a = joined[[paste0(cov, "_ast")]]
  b = joined[[paste0(cov, "_mic")]]
  
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
# multiplied by 10,000, plus 1, log-transformed.)

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

# Per-cluster summary + cap-vs-total trade-off table. Used to pick a
# per-cluster_sample cap for G02.
 
message("\n\n          *** Cells per cluster_sample diagnostics... ", Sys.time(), "\n\n")
 
 
t_csize = seur@meta.data %>%
  count(cell_type_joint, cluster_name, sample, name = "n_cells")
 
 
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
 
 
cat("\nPer-cluster summary:\n");                                print(t_summary, n = Inf)
cat("\nTotal cells after downsampling at candidate caps:\n");   print(t_caps)
 
 
### save tables
write_csv(t_summary, paste0(out_dir, script_ind, "cluster_sample_size_summary.csv"))
write_csv(t_caps,    paste0(out_dir, script_ind, "cluster_sample_size_at_caps.csv"))



###########################################################
# 8. save merged Seurat
###########################################################

message("\n\n          *** Save merged Seurat... ", Sys.time(), "\n\n")


qsave(seur, file = paste0(out_dir, script_ind, "seur_merged.qs"))


message("\n\n##########################################################################\n",
        "# Finished G01b ", Sys.time(),
        "\n##########################################################################\n\n")
