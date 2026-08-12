# LD_X10: re-draw the two largest7 GSEA NES heatmaps (Hallmark + Green) used in
# the thesis Results, for further figure edits.
#
# DATA SOURCE: LD_E03c_GSEA_results.csv (written by LD_E03c_GSEA_5covar.R - the
# fgsea run itself, Hallmark + Green24 astrocyte-state signatures, 5-covariate
# E02c DESeq2 results). This script only re-plots from that already-computed
# table; it does not rerun fgsea/DESeq2/msigdbr. Matrix-building logic
# (create_path_comp_mat_list / hallmark_family ordering) is copied as-is from
# LD_E03c_GSEA_5covar.R; the drawing step is reimplemented with ComplexHeatmap
# instead of pheatmap so comparison groups can get an actual text title above
# each block of columns (pheatmap can only do a colour-coded legend for that).

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

base    = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
gsea_csv = file.path(base, "LD_E_DESeq_pseudobulk/LD_E03c_GSEA_results.csv")
out_dir  = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X10_v02_"
if (!file.exists(gsea_csv)) stop("Missing input: ", gsea_csv)

gsea_res_tab = read_csv(gsea_csv, show_col_types = FALSE)

# comps are named "<cluster>_<comp_tag>". largest7 keeps only comps whose
# subcluster is among the 7 largest by nuclei (same set as LD_E03c).
largest7 = c("AST_SLC1A2_s0", "AST_SLC1A2_s3", "AST_SLC1A2_s4",
             "AST_GFAP_s1", "AST_GFAP_s2", "AST_CHI3L1_s6", "AST_CHI3L1_s9")
all_comps     = unique(gsea_res_tab$comp)
comps_largest = all_comps[vapply(all_comps,
                  function(c1) any(startsWith(c1, paste0(largest7, "_"))), logical(1))]

# drop the "cluster vs AST_SLC1A2_s0" comparisons (cluster-identity contrast,
# not a diagnosis/genotype comparison) - not needed for the presentation
comps_largest = comps_largest[!grepl("_CV_vs_AST_SLC1A2_s0$", comps_largest)]

# contrast-type per comp -> becomes the column-group title drawn above each
# block ("AD vs Control" / "R62H vs CV" / "R47H vs CV" / "R47H vs R62H"), and
# the cluster name alone -> becomes the (short) bottom column label
comp_tag = function(comps) dplyr::case_when(
  grepl("_TREM2_CV_AD_vs_Control$", comps)  ~ "AD vs Control",
  grepl("_AD_TREM2_R62H_vs_CV$",    comps)  ~ "R62H vs CV",
  grepl("_AD_TREM2_R47H_vs_CV$",    comps)  ~ "R47H vs CV",
  grepl("_AD_TREM2_R47H_vs_R62H$",  comps)  ~ "R47H vs R62H",
  TRUE ~ "other")
cluster_of = function(comps){
  x = comps
  x = sub("_TREM2_CV_AD_vs_Control$", "", x)
  x = sub("_AD_TREM2_R62H_vs_CV$",    "", x)
  x = sub("_AD_TREM2_R47H_vs_CV$",    "", x)
  x = sub("_AD_TREM2_R47H_vs_R62H$",  "", x)
  x
}
tag_levels = c("AD vs Control", "R62H vs CV", "R47H vs CV", "R47H vs R62H")
tag_seq    = factor(comp_tag(comps_largest), levels = tag_levels)
col_labels = cluster_of(comps_largest)

SIG_CUT = 0.10   # FDR < 0.1, same threshold as LD_E03c

### create_path_comp_mat_list: pathways x comparisons NES matrix, one per sub_cat
create_path_comp_mat_list = function(gsea_res_tab, comps_sel = unique(gsea_res_tab$comp),
                                     comp_ref = comps_sel, sig_cut = 0.05){
  path_comp_mat_list = list()
  for (cat1 in unique(gsea_res_tab$sub_cat)){
    t2 = gsea_res_tab[gsea_res_tab$sub_cat == cat1 & gsea_res_tab$comp %in% comps_sel &
                          !is.na(gsea_res_tab$padj) & gsea_res_tab$padj < sig_cut,]
    pathways_sel = unique(t2$pathway)
    path_comp_mat = matrix(nrow = length(pathways_sel), ncol = length(comps_sel),
                           dimnames = list(pathways_sel, comps_sel))
    for (comp1 in comps_sel){
      for (path1 in pathways_sel){
        t3 = t2[t2$comp == comp1 & t2$pathway == path1,]
        if (nrow(t3)>0){
          path_comp_mat[path1, comp1] = t2[["NES"]][t2$comp == comp1 & t2$pathway == path1]
        }
      }
    }
    path_comp_mat[is.na(path_comp_mat)] = 0
    m1 = path_comp_mat != 0
    keep_pathways = apply(as.matrix(m1[,colnames(m1) %in% comp_ref]), 1, any)
    m2 = path_comp_mat[keep_pathways,]
    if (is.matrix(m2)){
      if (nrow(m2)>1){path_comp_mat_list[[cat1]] = m2}
    }
  }
  return(path_comp_mat_list)
}

# MSigDB Hallmark process categories (Liberzon et al. 2015, Cell Syst) -> used to
# order Hallmark rows by functional family.
hallmark_family = c(
  HALLMARK_ALLOGRAFT_REJECTION="immune", HALLMARK_COAGULATION="immune",
  HALLMARK_COMPLEMENT="immune", HALLMARK_INTERFERON_ALPHA_RESPONSE="immune",
  HALLMARK_INTERFERON_GAMMA_RESPONSE="immune", HALLMARK_IL6_JAK_STAT3_SIGNALING="immune",
  HALLMARK_INFLAMMATORY_RESPONSE="immune",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB="signaling", HALLMARK_ANDROGEN_RESPONSE="signaling",
  HALLMARK_ESTROGEN_RESPONSE_EARLY="signaling", HALLMARK_ESTROGEN_RESPONSE_LATE="signaling",
  HALLMARK_IL2_STAT5_SIGNALING="signaling", HALLMARK_KRAS_SIGNALING_UP="signaling",
  HALLMARK_KRAS_SIGNALING_DN="signaling", HALLMARK_MTORC1_SIGNALING="signaling",
  HALLMARK_NOTCH_SIGNALING="signaling", HALLMARK_PI3K_AKT_MTOR_SIGNALING="signaling",
  HALLMARK_HEDGEHOG_SIGNALING="signaling", HALLMARK_TGF_BETA_SIGNALING="signaling",
  HALLMARK_WNT_BETA_CATENIN_SIGNALING="signaling",
  HALLMARK_APOPTOSIS="pathway", HALLMARK_HYPOXIA="pathway",
  HALLMARK_PROTEIN_SECRETION="pathway", HALLMARK_UNFOLDED_PROTEIN_RESPONSE="pathway",
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY="pathway",
  HALLMARK_BILE_ACID_METABOLISM="metabolic", HALLMARK_CHOLESTEROL_HOMEOSTASIS="metabolic",
  HALLMARK_FATTY_ACID_METABOLISM="metabolic", HALLMARK_GLYCOLYSIS="metabolic",
  HALLMARK_HEME_METABOLISM="metabolic", HALLMARK_OXIDATIVE_PHOSPHORYLATION="metabolic",
  HALLMARK_XENOBIOTIC_METABOLISM="metabolic",
  HALLMARK_E2F_TARGETS="proliferation", HALLMARK_G2M_CHECKPOINT="proliferation",
  HALLMARK_MYC_TARGETS_V1="proliferation", HALLMARK_MYC_TARGETS_V2="proliferation",
  HALLMARK_P53_PATHWAY="proliferation", HALLMARK_MITOTIC_SPINDLE="proliferation",
  HALLMARK_DNA_REPAIR="DNA_damage", HALLMARK_UV_RESPONSE_DN="DNA_damage",
  HALLMARK_UV_RESPONSE_UP="DNA_damage",
  HALLMARK_ADIPOGENESIS="development", HALLMARK_ANGIOGENESIS="development",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION="development", HALLMARK_MYOGENESIS="development",
  HALLMARK_SPERMATOGENESIS="development", HALLMARK_PANCREAS_BETA_CELLS="development",
  HALLMARK_APICAL_JUNCTION="cellular_component", HALLMARK_APICAL_SURFACE="cellular_component",
  HALLMARK_PEROXISOME="cellular_component")
fam_levels = c("immune", "signaling", "pathway", "metabolic", "proliferation",
               "DNA_damage", "development", "cellular_component", "other")
# same colourblind-safe palette convention as the rest of the thesis (X02 UMAP):
# viridis for >8 categories (9 families here exceeds Okabe-Ito's 8), instead of
# pheatmap's arbitrary auto-generated pastel colours
fam_colors = set_names(scales::viridis_pal(option = "viridis")(length(fam_levels)), fam_levels)

# Hallmark terms dropped for simplicity (not statistically excluded - just
# hidden from this plot; still present in LD_E03c_GSEA_results.csv/other plots)
hallmark_drop_terms = c(
  "HALLMARK_PANCREAS_BETA_CELLS",
  "HALLMARK_ANDROGEN_RESPONSE",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_HEDGEHOG_SIGNALING",
  "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_E2F_TARGETS"
)

# Green24 state display labels - identical text to green_labels in
# LD_X04_B_characterisation_plots.R (order = Ast.1 -> Ast.10)
green_state_labels = c(
  "Green24_Ast.1"  = "Ast1 homeostatic",
  "Green24_Ast.2"  = "Ast2 homeostatic",
  "Green24_Ast.3"  = "Ast3 enhanced mitophagy",
  "Green24_Ast.4"  = "Ast4 reactive",
  "Green24_Ast.5"  = "Ast5 reactive",
  "Green24_Ast.6"  = "Ast6",
  "Green24_Ast.7"  = "Ast7 IFN response",
  "Green24_Ast.8"  = "Ast8 stress response",
  "Green24_Ast.9"  = "Ast9 stress response",
  "Green24_Ast.10" = "Ast10 AD-elevated"
)

# draw one NES heatmap (a single matrix) to PDF and PNG. Columns are split into
# titled blocks by comparison (tag_seq/tag_levels, shared across both plots);
# Hallmark rows are additionally grouped by family with a coloured side-bar
# (row_split only creates the gap/order; its own text titles are suppressed
# via row_title = NULL so the colour + legend does the labelling, as before).
draw_heatmap = function(mat, collection, fbase){
  # column order/labels are fixed (tag_seq/col_labels built from comps_largest,
  # same order used to build `mat`) - never re-clustered
  stopifnot(identical(colnames(mat), comps_largest))
  labels_col = col_labels
  row_labels = rownames(mat)
  left_annot = NULL

  if (collection == "Hallmark"){
    mat = mat[!rownames(mat) %in% hallmark_drop_terms, , drop = FALSE]
    fam = hallmark_family[rownames(mat)]; fam[is.na(fam)] = "other"
    fam = factor(fam, levels = fam_levels)
    row_order = order(fam, rownames(mat))
    mat = mat[row_order, , drop = FALSE]; fam = fam[row_order]
    row_split = fam; cluster_rows = FALSE
    row_labels = sub("^HALLMARK_", "", rownames(mat))   # drop the "HALLMARK_" prefix for display only
    left_annot = rowAnnotation(family = fam, col = list(family = fam_colors),
                               show_annotation_name = FALSE)
  } else {                                                # Green: cluster as original
    row_split = NULL; cluster_rows = TRUE
    row_labels = unname(green_state_labels[rownames(mat)])
  }

  lim = max(abs(mat))
  # same up/down convention as the rest of the thesis (Green24 z-scores, DEG
  # panel, pairwise concordance): blue = down/negative NES, orange = up/positive
  col_fun = colorRamp2(c(-lim, 0, lim), c("#0072B2", "white", "#E69F00"))
  W = ncol(mat) * 0.20 + 5
  H = nrow(mat) * 0.18 + 3

  # ComplexHeatmap's default legend for a colorRamp2 colour function auto-picks
  # its own tick range, which does NOT necessarily match the colour function's
  # actual saturation point (lim) - explicit `at` forces the legend to show the
  # true data range instead of a wider auto-extended one
  legend_at = round(c(-lim, -lim / 2, 0, lim / 2, lim), 1)

  ht = Heatmap(mat, name = "NES", col = col_fun,
              cluster_rows = cluster_rows, cluster_columns = FALSE,
              row_split = row_split, cluster_row_slices = FALSE, row_title = NULL,
              row_gap = unit(1.5, "mm"), left_annotation = left_annot,
              column_split = tag_seq, cluster_column_slices = FALSE,
              column_title_gp = gpar(fontsize = 10, fontface = "bold"),
              column_gap = unit(2, "mm"),
              column_labels = labels_col, column_names_gp = gpar(fontsize = 8),
              row_labels = row_labels, row_names_gp = gpar(fontsize = 8),
              rect_gp = gpar(col = NA), border = TRUE,
              column_title_side = "top", column_names_rot = 45,
              heatmap_legend_param = list(title = "NES", at = legend_at))

  pdf(paste0(fbase, ".pdf"), width = W, height = H); draw(ht); dev.off()
  png(paste0(fbase, ".png"), width = W, height = H, units = "in", res = 300); draw(ht); dev.off()
}

ml = create_path_comp_mat_list(gsea_res_tab, comps_sel = comps_largest, sig_cut = SIG_CUT)

draw_heatmap(ml[["HALLMARKS"]], "Hallmark",
            paste0(out_dir, "/", script_ind, "GSEA_heatmap_largest7_Hallmark_FDR10_allsig"))
draw_heatmap(ml[["user_def_sets"]], "Green",
            paste0(out_dir, "/", script_ind, "GSEA_heatmap_largest7_Green_FDR10_allsig"))

message("Done. Written to ", out_dir)
