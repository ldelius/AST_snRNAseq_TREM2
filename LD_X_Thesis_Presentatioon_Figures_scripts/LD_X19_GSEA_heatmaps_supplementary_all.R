# LD_X19: Supplementary GSEA heatmaps with all subclusters, comparisons and
# significant terms. Plot definitions are inherited from LD_X10b.

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")
### inherit results-figure definitions ---------------------------------------
# NB: the eval below also sets base/out_dir/script_ind (to LD_X10b_'s values), so
# this script's own paths and prefix are (re)assigned AFTER it, not before.
x10b = file.path(base, "LD_X_Thesis_Presentatioon_Figures_scripts",
                 "LD_X10b_GSEA_heatmap_combined_prep.R")
if (!file.exists(x10b)) stop("Missing: ", x10b)
src  = readLines(x10b)
stop_at = grep("^ml = create_path_comp_mat_list", src)[1]   # first line that builds/draws
if (is.na(stop_at)) stop("Could not find the build entry point in LD_X10b.")

# This figure has 61 columns against the results figure's ~21, so the inherited
# 8 pt column labels overlap. Patch that one literal into a variable before
# evaluating, rather than duplicating the whole build_heatmap() function here.
COL_LAB_FS = 6      # column (subcluster) label size; LD_X10b uses 8
n_patched = sum(grepl("column_names_gp = gpar(fontsize = 8)", src, fixed = TRUE))
if (n_patched != 1)
  stop("Expected exactly one column_names_gp literal in LD_X10b, found ", n_patched,
       " - check whether that script changed.")
src = sub("column_names_gp = gpar(fontsize = 8)",
          "column_names_gp = gpar(fontsize = COL_LAB_FS)", src, fixed = TRUE)

eval(parse(text = paste(src[seq_len(stop_at - 1)], collapse = "\n")))

# reclaim our own output prefix and paths (the eval set them to LD_X10b's)
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
script_ind = "LD_X19_"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
message("Inherited config from LD_X10b (", stop_at - 1, " lines, nothing drawn).")

### include all comparisons --------------------------------------------------
# The results figure keeps only three contrast blocks. The supplement keeps every
# comparison, so two more blocks are added back and need their own tag/label rules:
#   R47H vs R62H            - the direct variant-vs-variant contrast
#   Subcluster vs SLC1A2_s0 - the cluster-identity contrast (each subcluster's CV
#                             samples against the largest subcluster's)
comp_tag = function(comps) dplyr::case_when(
  grepl("_TREM2_CV_AD_vs_Control$", comps)   ~ "AD vs Control",
  grepl("_AD_TREM2_R62H_vs_CV$",    comps)   ~ "R62H vs CV",
  grepl("_AD_TREM2_R47H_vs_CV$",    comps)   ~ "R47H vs CV",
  grepl("_AD_TREM2_R47H_vs_R62H$",  comps)   ~ "R47H vs R62H",
  grepl("_CV_vs_AST_SLC1A2_s0$",    comps)   ~ "Subcluster vs SLC1A2_s0",
  TRUE ~ "other")
cluster_of = function(comps){
  x = comps
  x = sub("_TREM2_CV_AD_vs_Control$", "", x)
  x = sub("_AD_TREM2_R62H_vs_CV$",    "", x)
  x = sub("_AD_TREM2_R47H_vs_CV$",    "", x)
  x = sub("_AD_TREM2_R47H_vs_R62H$",  "", x)
  x = sub("_CV_vs_AST_SLC1A2_s0$",    "", x)
  x
}
# Subcluster-identity contrast first (far left): it answers a different question
# from the disease/genotype blocks, so it reads as the baseline the rest follow.
tag_levels = c("Subcluster vs SLC1A2_s0",
               "AD vs Control", "R62H vs CV", "R47H vs CV", "R47H vs R62H")

comps_all = unique(gsea_res_tab$comp)
unknown = comps_all[comp_tag(comps_all) == "other"]
if (length(unknown)) stop("Unrecognised comparison naming: ", paste(unknown, collapse = ", "))
# order columns by contrast block, then by cluster name within block
comps_all = comps_all[order(factor(comp_tag(comps_all), levels = tag_levels),
                            cluster_of(comps_all))]

comps_largest = comps_all                      # build_heatmap() checks against this
tag_seq       = factor(comp_tag(comps_all), levels = tag_levels)
col_labels    = cluster_of(comps_all)
message("Comparisons: ", length(comps_all), " (",
        paste(names(table(tag_seq)), table(tag_seq), sep = ": ", collapse = ", "), ")")

### retain all terms ---------------------------------------------------------
sig = gsea_res_tab %>%
  filter(sub_cat == "HALLMARKS", !is.na(padj), padj < SIG_CUT, comp %in% comps_all)
terms_all = unique(sig$pathway)
ungrouped = setdiff(terms_all, unlist(custom_hallmark_groups))

hallmark_drop_terms = character(0)              # nothing dropped in the supplement
custom_hallmark_groups = c(custom_hallmark_groups, list(Other = ungrouped))
custom_group_levels = names(custom_hallmark_groups)
term_to_group = setNames(rep(custom_group_levels, lengths(custom_hallmark_groups)),
                         unlist(custom_hallmark_groups))
# extend the group palette to 9 blocks, same viridis family as the results figure,
# with "Other" in neutral grey so it does not read as another biological block
custom_group_colors = c(
  setNames(scales::viridis_pal(option = "viridis")(length(custom_group_levels) - 1),
           setdiff(custom_group_levels, "Other")),
  Other = "grey75")
message("Hallmark terms: ", length(terms_all), " (", length(ungrouped), " in \"Other\")")

### build and save -----------------------------------------------------------
ml = create_path_comp_mat_list(gsea_res_tab, comps_sel = comps_all, sig_cut = SIG_CUT)

ht_hallmark = build_heatmap(ml[["HALLMARKS"]],     "Hallmark")
ht_green    = build_heatmap(ml[["user_def_sets"]], "Green")

save_heatmap(ht_hallmark, file.path(out_dir, paste0(script_ind, "GSEA_heatmap_ALL_Hallmark_FDR10_allsig")))
save_heatmap(ht_green,    file.path(out_dir, paste0(script_ind, "GSEA_heatmap_ALL_Green_FDR10_allsig")))

message("Done. Hallmark ", nrow(ml[["HALLMARKS"]]), " x ", ncol(ml[["HALLMARKS"]]),
        " | Green ", nrow(ml[["user_def_sets"]]), " x ", ncol(ml[["user_def_sets"]]),
        ". Outputs in: ", out_dir)
