# LD_X05b: Verify LD_X05 contrast arithmetic using treatment coding.
# This avoids the degenerate results possible with sccomp string contrasts;
# the smallest unadjusted subtype model is used as a targeted check.

library(tidyverse)
library(qs)
library(Seurat)
library(sccomp)

### paths ---------------------------------------------------------------------
base       = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
b04_path   = file.path(base, "LD_B_AST_analysis_output/LD_B04a_v02_seur.qs")
clust_csv  = file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
group_csv  = file.path(base, "data_TREM2_michael/A_input/group_tab.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
existing_csv = file.path(out_dir, "LD_X05_sccomp_subtypes_unadj.csv")
script_ind = "LD_X05b_v01_"
cores_n    = 8
for (p in c(b04_path, clust_csv, group_csv, existing_csv)) if (!file.exists(p)) stop("Missing input: ", p)

ord = read_csv(clust_csv, show_col_types = FALSE)
subtype_levels = unique(ord$cell_type)
grtab = read_csv(group_csv, show_col_types = FALSE)
gr = unique(grtab$group)

message("Loading B04 object (metadata only)...")
seur = qread(b04_path)
meta = seur@meta.data
if (!"group" %in% names(meta)) meta$group = grtab$group[match(meta$sample, grtab$sample)]
rm(seur); gc()
samp_meta = grtab %>% dplyr::distinct(sample, group)

### per-sample subtype counts (same as LD_X05's counts_tab, subtype level) ----
u = factor(meta$cell_type, levels = subtype_levels); s = factor(meta$sample)
ct = as.data.frame(table(sample = s, cellgroup = u), stringsAsFactors = FALSE)
names(ct)[names(ct) == "Freq"] = "N_cells"
ct = ct %>% dplyr::left_join(samp_meta, by = "sample") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(N_sample = sum(N_cells)) %>%
  dplyr::ungroup()

### treatment-coded model: Control_CV as reference, NO contrasts= -------------
# relevel so Control_CV is first -> its the implicit reference under treatment
# (dummy) coding, so "groupAD_R62H" etc. become direct vs-Control_CV differences
ct$group = droplevels(factor(ct$group, levels = c("Control_CV", setdiff(gr, "Control_CV"))))

message("Fitting treatment-coded model (no contrasts=)...")
res = ct |>
  sccomp_glm(formula_composition = ~ group,
             .sample = sample, .cell_group = cellgroup, .count = N_cells,
             bimodal_mean_variability_association = TRUE, cores = cores_n)

r = as.data.frame(res)
if (!"factor" %in% names(r) && "parameter" %in% names(r)) r$factor = r$parameter
par_col = intersect(c("parameter", "factor"), names(r))[1]
cg_col  = intersect(c("cellgroup", "cell_group"), names(r))[1]
eff_col = intersect(c("c_effect", "effect"), names(r))[1]

verify = r %>%
  dplyr::filter(.data[[par_col]] == "groupAD_R62H") %>%
  dplyr::transmute(cellgroup = .data[[cg_col]],
                   method = "treatment_coding_no_contrasts",
                   effect = .data[[eff_col]], c_FDR = c_FDR)

write_csv(verify, file.path(out_dir, paste0(script_ind, "verify_AD_R62H_vs_CtrlCV_subtypes.csv")))

### compare against LD_X05's existing (contrasts=) result ---------------------
existing = read_csv(existing_csv, show_col_types = FALSE) %>%
  dplyr::filter(comparison == "AD_R62H vs Ctrl_CV") %>%
  dplyr::transmute(cellgroup, method = "cell_means_with_contrasts", effect = NA_real_, c_FDR)

comparison_tab = dplyr::bind_rows(verify, existing) %>% dplyr::arrange(cellgroup, method)
write_csv(comparison_tab, file.path(out_dir, paste0(script_ind, "verification_comparison.csv")))

message("\n=== Verification: AD_R62H vs Control_CV, subtypes, unadjusted ===")
print(as.data.frame(comparison_tab))
message("\nIf c_FDR values per cellgroup are similar in magnitude/ranking between the\n",
       "two methods, LD_X05's contrasts= usage is verified reliable for this pattern.\n",
       "If 'cell_means_with_contrasts' shows FDR ~= 1 / identical estimates across\n",
       "cellgroups while 'treatment_coding_no_contrasts' shows real separation,\n",
       "contrasts= is broken here and LD_X05's results need redoing without it.")
