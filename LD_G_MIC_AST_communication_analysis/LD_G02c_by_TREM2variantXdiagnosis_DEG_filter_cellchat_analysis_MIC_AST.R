message("\n\n##########################################################################\n",
        "# Start G02c: CellChat by TREM2Variant x NeuropathologicalDiagnosis ", Sys.time(),
        "\n##########################################################################\n",
        "     LD_G02a - ungrouped\n",
        "     LD_G02b - TREM2Variant\n",
        "     LD_G02c - TREM2Variant x NeuropathologicalDiagnosis (this script, 6 splits)\n",
        "\n   All six splits are run, including R47H_Control (N=2 samples).\n",
        "   Splits with very few cells will print warnings during CellChat - kept for\n",
        "   completeness, will be flagged at interpretation time.\n",
        "\n   Parameters: assay='RNA' (LogNormalize from G01b), type='triMean', nboot=100,\n",
        "   filterCommunication(min.cells = 10) post-hoc, cluster_sample cap of 200.\n",
        "\n##########################################################################\n\n")


#set environment/load packages
library(qs)
library(tidyverse)
library(Seurat)
library(CellChat)
library(patchwork)
library(future)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 50 * 1024^3)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

script_ind = "LD_G02c_v001_"

out_dir = paste0(main_dir, "LD_G_MIC_AST_communication_analysis_output/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

in_seur     = paste0(out_dir, "LD_G01b_seur_merged.qs")
in_DEGs_AST = paste0(main_dir, "LD_E_DESeq_pseudobulk/LD_E02a2_v02_bulk_data.qs")
in_DEGs_MIC = paste0(main_dir, "data_TREM2_michael/E_DESeq_pseudobulk/E02a2_bulk_data.qs")


n_workers = 8
future::plan("multisession", workers = n_workers)


cc_type      = "triMean"
cc_nboot     = 100
cc_min_cells = 10
cc_seed      = 42      # passed to computeCommunProb's seed.use (CellChat's built-in mechanism for parallel-safe bootstrap reproducibility)


#cluster_sample cap: keep at most N cells from each cluster_sample combination.
#Necessary because full data (304k cells) caused computeCommunProb to OOM at
#500 GB. Capping per cluster_sample preserves sample-level balance.
cap_cluster_sample = 200



###########################################################
# 1. load merged Seurat, set samples column
###########################################################

message("\n\n          *** Load merged Seurat... ", Sys.time(), "\n\n")


seur = qread(file = in_seur)
DefaultAssay(seur) = "RNA"

seur$samples = seur$sample

seur$cluster_name = droplevels(factor(seur$cluster_name))
Idents(seur) = "cluster_name"


###########################################################
# 1b. cap cells per cluster_sample combination
###########################################################

#Sample-balanced downsampling. Each cluster_sample combination capped at cap_cluster_sample cells.

message("\n\n          *** Cap cells per cluster_sample at ", cap_cluster_sample, "... ", Sys.time(), "\n\n")


set.seed(cc_seed)


t_cs = seur@meta.data %>% as_tibble() %>%
  count(cluster_sample, name = "n_cells") %>%
  mutate(capped = n_cells > cap_cluster_sample,
         n_after = pmin(n_cells, cap_cluster_sample))

cat("cluster_sample summary:\n",
    "  total cluster_samples:        ", nrow(t_cs), "\n",
    "  combos > cap (",  cap_cluster_sample, "):       ", sum(t_cs$capped), "\n",
    "  cells before cap (total):     ", sum(t_cs$n_cells), "\n",
    "  cells after cap (estimated):  ", sum(t_cs$n_after), "\n",
    sep = "")


cell_ids_keep = seur@meta.data %>%
  rownames_to_column("cell_id") %>%
  group_by(cluster_sample) %>%
  slice_sample(n = cap_cluster_sample) %>%
  pull(cell_id)

seur = seur[, cell_ids_keep]

cat("\nCells after cluster_sample cap:", ncol(seur), "\n\n")


seur$cluster_name = droplevels(factor(seur$cluster_name))
Idents(seur) = "cluster_name"

write_csv(t_cs, paste0(out_dir, script_ind, "cluster_sample_cap_summary.csv"))



###########################################################
# 2. cell-count diagnostics per TREM2Variant x Diagnosis split
###########################################################

#Important here: small subgroups (esp. R47H_Control with N=2 samples) may have
#<10 cells in some clusters after the cap, which filterCommunication will then
#drop. Table to cross-reference with G02c output.
 
seur@meta.data %>% as_tibble() %>%
  mutate(group = paste0(TREM2Variant, "_", NeuropathologicalDiagnosis)) %>%
  count(cluster_name, group, name = "n_cells") %>%
  pivot_wider(names_from = group, values_from = n_cells, values_fill = 0) %>%
  write_csv(paste0(out_dir, script_ind, "cells_per_cluster_per_TREM2Variant_x_Diagnosis.csv"))


###########################################################
# 3. load DEGs (AST + MIC), build DEG-filtered CellChatDB
###########################################################

message("\n\n          *** Build DEG-filtered CellChatDB... ", Sys.time(), "\n\n")


DEGs_comb = NULL


### AST DEGs
if (file.exists(in_DEGs_AST)){
  bd_ast = qread(in_DEGs_AST)
  DEGs_comb = unique(c(DEGs_comb, unlist(bd_ast$DEGs)))
  cat("Loaded AST DEGs from:", in_DEGs_AST, "- total unique genes so far:", length(DEGs_comb), "\n")
  rm(bd_ast); gc()
} else {
  stop("AST DEG file not found: ", in_DEGs_AST)
}


### MIC DEGs
if (file.exists(in_DEGs_MIC)){
  bd_mic = qread(in_DEGs_MIC)
  DEGs_comb = unique(c(DEGs_comb, unlist(bd_mic$DEGs)))
  cat("Loaded MIC DEGs from:", in_DEGs_MIC, "- total unique genes after MIC merge:", length(DEGs_comb), "\n")
  rm(bd_mic); gc()
} else {
  stop("MIC DEG file not found: ", in_DEGs_MIC)
}


cat("Total unique DEGs for DB filtering:", length(DEGs_comb), "\n")


### filter CellChatDB.human
#drop = FALSE prevents R collapsing to a vector if exactly one row matches
#(would silently NULL out rownames() and break the downstream filter)
CellChatDB = CellChatDB.human

t1 = CellChatDB
t2 = t1$complex
t3 = apply(t2, c(1,2), `%in%`, DEGs_comb)
complexes_with_DEGs = t3[apply(t3, 1, any), , drop = FALSE]

t2 = t1$cofactor
t3 = apply(t2, c(1,2), `%in%`, DEGs_comb)
cofactors_with_DEGs = t3[apply(t3, 1, any), , drop = FALSE]

t2 = t1$interaction
t3 = t2[t2$ligand    %in% DEGs_comb |
        t2$ligand    %in% rownames(complexes_with_DEGs) |
        t2$receptor  %in% DEGs_comb |
        t2$receptor  %in% rownames(complexes_with_DEGs) |
        t2$agonist   %in% rownames(cofactors_with_DEGs) |
        t2$antagonist %in% rownames(cofactors_with_DEGs), ]

CellChatDB.use = subsetDB(CellChatDB, search = t3$interaction_name, key = "interaction_name")

cat("\nCellChatDB filtering:\n",
    "  Full DB interactions: ", nrow(CellChatDB$interaction), "\n",
    "  After DEG filter:     ", nrow(t3), "\n",
    "  In CellChatDB.use:    ", nrow(CellChatDB.use$interaction), "\n", sep = "")



###########################################################
# 4. helper function
###########################################################

run_cellchat = function(seur_in, run_label, db_use,
                        cc_type, cc_nboot, cc_min_cells, cc_seed,
                        out_dir, script_ind){

  message("\n\n   *** Run CellChat: ", run_label, " - ", Sys.time(), "\n")
  cat("       cells:", ncol(seur_in), "  clusters:", length(unique(Idents(seur_in))),
      "  samples:", length(unique(seur_in$samples)), "\n")

  cellchat = createCellChat(object = seur_in, group.by = "ident", assay = "RNA")
  cellchat@DB = db_use

  cellchat = subsetData(cellchat)
  cellchat = identifyOverExpressedGenes(cellchat)
  cellchat = identifyOverExpressedInteractions(cellchat)

  #seed.use is CellChat's parallel-safe seed (propagates to future workers)
  cellchat = computeCommunProb(cellchat, type = cc_type, nboot = cc_nboot,
                               seed.use = cc_seed)

  cellchat = filterCommunication(cellchat, min.cells = cc_min_cells)
  cellchat = computeCommunProbPathway(cellchat)
  cellchat = aggregateNet(cellchat)

  qsave(cellchat,
        file = paste0(out_dir, script_ind, "cellchat_", run_label, ".qs"))

  message("       Saved: ", out_dir, script_ind, "cellchat_", run_label, ".qs - ", Sys.time())

  return(cellchat)
}



###########################################################
# 5. RUN: split by TREM2Variant x Diagnosis, loop over splits
###########################################################

message("\n\n##########################################################################\n",
        "# RUN: CellChat by TREM2Variant x NeuropathologicalDiagnosis ", Sys.time(),
        "\n##########################################################################\n\n")


stopifnot("TREM2Variant" %in% colnames(seur@meta.data),
          "NeuropathologicalDiagnosis" %in% colnames(seur@meta.data))

seur$group_TREM2_diag = paste0(seur$TREM2Variant, "_", seur$NeuropathologicalDiagnosis)
seur_list_TREM2_diag = SplitObject(seur, split.by = "group_TREM2_diag")

#drop original Seurat now that splits are made - SplitObject duplicates cells
rm(seur); gc()

cat("\nTREM2Variant_x_Diagnosis splits:\n")
for (g in names(seur_list_TREM2_diag)){
  cat("  ", g, ":", ncol(seur_list_TREM2_diag[[g]]), "cells\n")
}


for (g in names(seur_list_TREM2_diag)){

  cc = run_cellchat(
    seur_in      = seur_list_TREM2_diag[[g]],
    run_label    = paste0("TREM2Variant_x_Diagnosis_", g),
    db_use       = CellChatDB.use,
    cc_type      = cc_type,
    cc_nboot     = cc_nboot,
    cc_min_cells = cc_min_cells,
    cc_seed      = cc_seed,
    out_dir      = out_dir,
    script_ind   = script_ind
  )

  rm(cc); gc()
}

rm(seur_list_TREM2_diag); gc()



###########################################################
# 6. companion metadata file
###########################################################

qsave(list(CellChatDB.use = CellChatDB.use,
           DEGs_comb      = DEGs_comb,
           params         = list(type      = cc_type,
                                 nboot     = cc_nboot,
                                 min.cells = cc_min_cells,
                                 seed.use  = cc_seed,
                                 cap_cluster_sample = cap_cluster_sample,
                                 assay     = "RNA",
                                 n_workers = n_workers)),
      file = paste0(out_dir, script_ind, "cellchat_metadata.qs"))



sessionInfo()


message("\n\n##########################################################################\n",
        "# Completed G02c ", Sys.time(),
        "\n##########################################################################\n\n")