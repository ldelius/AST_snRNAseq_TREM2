# LD_X10b: Hallmark and Green24 GSEA heatmaps for the seven largest subclusters.

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)   # viewport/grid.layout/grid.text - used to place the two
                # independent heatmaps side by side in the combined figure

base    = "/Volumes/lvd25/home/AST_scRNAseq_TREM2"
gsea_csv = file.path(base, "LD_E_DESeq_pseudobulk/LD_E03c_GSEA_results.csv")
out_dir  = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X10b_"
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
# not a diagnosis/genotype comparison) - not needed for this figure
comps_largest = comps_largest[!grepl("_CV_vs_AST_SLC1A2_s0$", comps_largest)]

# drop the "R47H vs R62H" comparison block entirely (carried over from PP)
comps_largest = comps_largest[!grepl("_AD_TREM2_R47H_vs_R62H$", comps_largest)]

# contrast-type per comp -> becomes the column-group title drawn above each
# block ("AD vs Control" / "R62H vs CV" / "R47H vs CV"), and the cluster name
# alone -> becomes the (short) bottom column label
comp_tag = function(comps) dplyr::case_when(
  grepl("_TREM2_CV_AD_vs_Control$", comps)  ~ "AD vs Control",
  grepl("_AD_TREM2_R62H_vs_CV$",    comps)  ~ "R62H vs CV",
  grepl("_AD_TREM2_R47H_vs_CV$",    comps)  ~ "R47H vs CV",
  TRUE ~ "other")
cluster_of = function(comps){
  x = comps
  x = sub("_TREM2_CV_AD_vs_Control$", "", x)
  x = sub("_AD_TREM2_R62H_vs_CV$",    "", x)
  x = sub("_AD_TREM2_R47H_vs_CV$",    "", x)
  x
}
tag_levels = c("AD vs Control", "R62H vs CV", "R47H vs CV")   # R47H vs R62H removed
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

# Hallmark terms dropped for simplicity (not statistically excluded - just
# hidden from this plot; still present in LD_E03c_GSEA_results.csv/other
# plots). Thesis base drop list plus the PP version's further slide-deck drops.
hallmark_drop_terms = c(
  "HALLMARK_PANCREAS_BETA_CELLS",
  "HALLMARK_ANDROGEN_RESPONSE",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_HEDGEHOG_SIGNALING",
  "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_PEROXISOME",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_UV_RESPONSE_DN",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_HEME_METABOLISM",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_KRAS_SIGNALING_DN",
  "HALLMARK_XENOBIOTIC_METABOLISM",
  "HALLMARK_TGF_BETA_SIGNALING"
)

# custom Hallmark grouping (hand-picked, replacing the official MSigDB family
# classification) - fixed blocks/order, top to bottom. Must exactly match the
# terms remaining after hallmark_drop_terms.
custom_hallmark_groups = list(
  "Interferon" = c(
    "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  "Other immune" = c(
    "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_COAGULATION", "HALLMARK_COMPLEMENT",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  "Metabolism" = c(
    "HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_CHOLESTEROL_HOMEOSTASIS"),
  "Stress and proteostasis" = c(
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),
  "Translation and growth" = c(
    "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_MTORC1_SIGNALING", "HALLMARK_PI3K_AKT_MTOR_SIGNALING"),
  "Proliferation and damage" = c(
    "HALLMARK_P53_PATHWAY", "HALLMARK_DNA_REPAIR", "HALLMARK_UV_RESPONSE_UP",
    "HALLMARK_APOPTOSIS"),
  "Structural remodelling" = c(
    "HALLMARK_MYOGENESIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_APICAL_JUNCTION", "HALLMARK_ADIPOGENESIS"),
  "Other signalling" = c(
    "HALLMARK_NOTCH_SIGNALING", "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    "HALLMARK_IL2_STAT5_SIGNALING")
)
custom_group_levels = names(custom_hallmark_groups)
term_to_group = setNames(rep(custom_group_levels, lengths(custom_hallmark_groups)),
                         unlist(custom_hallmark_groups))
# ColorBrewer "Dark2" 8-class qualitative palette (colourblind-safe), as in
# the PP version. Its hex values don't overlap with the NES diverging scale's
# blue/orange (#0072B2/#E69F00) at all, so - unlike plain Okabe-Ito, which
# would have reused those exact colours for a different meaning in the same
# figure - there's no same-figure clash here to avoid.
custom_group_colors = c(
  "Interferon" = "#1B9E77", "Other immune" = "#D95F02", "Metabolism" = "#7570B3",
  "Stress and proteostasis" = "#E7298A", "Translation and growth" = "#66A61E",
  "Proliferation and damage" = "#E6AB02", "Structural remodelling" = "#A6761D",
  "Other signalling" = "#666666")

# Green24 state display labels (order = Ast.1 -> Ast.10) - matches the labels
# already fixed in LD_X04c / LD_X04_B_characterisation_plots.R / LD_X10 earlier
green_state_labels = c(
  "Green24_Ast.1"  = "Ast1 homeostatic",
  "Green24_Ast.2"  = "Ast2 homeostatic",
  "Green24_Ast.3"  = "Ast3 enhanced mitophagy",
  "Green24_Ast.4"  = "Ast4 reactive",
  "Green24_Ast.5"  = "Ast5 reactive",
  "Green24_Ast.6"  = "Ast6",
  "Green24_Ast.7"  = "Ast7 Interferon response",
  "Green24_Ast.8"  = "Ast8 stress response",
  "Green24_Ast.9"  = "Ast9 stress response",
  "Green24_Ast.10" = "Ast10 AD-elevated"
)

# draw one NES heatmap (a single matrix) to PDF and PNG. Columns are split into
# titled blocks by comparison (tag_seq/tag_levels, shared across both plots);
# Hallmark rows are additionally grouped by the custom blocks with a coloured
# side-bar (row_split only creates the gap/order; its own text titles are
# suppressed via row_title = NULL so the colour + legend does the labelling).
#
# build_heatmap() returns the Heatmap OBJECT (not drawn yet) plus a suggested
# standalone width/height - split out from drawing so the same object can
# either be saved on its own, or placed into the combined side-by-side figure
# further down. Cell shape: rectangular (narrower than tall), not square -
# columns are fixed-width (comparisons don't vary), so shrinking column width
# is pure width-saving; rows get slightly MORE height to keep it legible.
CELL_W = 0.09   # in, per column - was 0.20, then 0.13; pushed narrower again
CELL_H = 0.20   # in, per row    - was 0.18; a bit taller to stay legible when narrower

build_heatmap = function(mat, collection){
  stopifnot(identical(colnames(mat), comps_largest))
  labels_col = col_labels
  row_labels = rownames(mat)
  left_annot = NULL
  row_rot    = 0   # horizontal in both panels

  if (collection == "Hallmark"){
    mat = mat[!rownames(mat) %in% hallmark_drop_terms, , drop = FALSE]
    if (!setequal(names(term_to_group), rownames(mat)))
      stop("custom_hallmark_groups does not match the term set exactly (missing or extra terms)")
    # exact order as given in custom_hallmark_groups (block by block, term by
    # term within each block) - NOT alphabetical
    row_order = match(unlist(custom_hallmark_groups), rownames(mat))
    mat = mat[row_order, , drop = FALSE]
    grp = factor(term_to_group[rownames(mat)], levels = custom_group_levels)
    row_split = grp; cluster_rows = FALSE
    row_labels = sub("^HALLMARK_", "", rownames(mat))   # drop the "HALLMARK_" prefix for display only
    # nrow = 4: lays the 8 group entries out as 4 rows x 2 columns (filled
    # column-by-column, i.e. first 4 down the left column, next 4 down the
    # right) instead of one long 8-row column - much more compact next to
    # the NES bar in the merged bottom legend row
    left_annot = rowAnnotation(group = grp, col = list(group = custom_group_colors),
                               show_annotation_name = FALSE,
                               annotation_legend_param = list(group = list(nrow = 4)))
  } else {                                                # Green: cluster as original
    row_split = NULL; cluster_rows = TRUE
    row_labels = unname(green_state_labels[rownames(mat)])
  }

  lim = max(abs(mat))
  # same up/down convention as the rest of the thesis (Green24 z-scores, DEG
  # panel, pairwise concordance): blue = down/negative NES, orange = up/positive
  col_fun = colorRamp2(c(-lim, 0, lim), c("#0072B2", "white", "#E69F00"))
  W = ncol(mat) * CELL_W + 5
  H = nrow(mat) * CELL_H + 3

  # ComplexHeatmap's default legend for a colorRamp2 colour function auto-picks
  # its own tick range, which does NOT necessarily match the colour function's
  # actual saturation point (lim) - explicit `at` forces the legend to show the
  # true data range instead of a wider auto-extended one
  legend_at = round(c(-lim, -lim / 2, 0, lim / 2, lim), 1)

  row_fontsize = 9   # same size in both panels (was 8 Hallmark / 13 Green)

  ht = Heatmap(mat, name = "NES", col = col_fun,
              cluster_rows = cluster_rows, cluster_columns = FALSE,
              row_split = row_split, cluster_row_slices = FALSE, row_title = NULL,
              row_gap = unit(1.5, "mm"), left_annotation = left_annot,
              column_split = tag_seq, cluster_column_slices = FALSE,
              column_title_gp = gpar(fontsize = 10, fontface = "bold"),
              column_gap = unit(2, "mm"),
              column_labels = labels_col, column_names_gp = gpar(fontsize = 8),
              row_labels = row_labels, row_names_gp = gpar(fontsize = row_fontsize),
              row_names_rot = row_rot,
              rect_gp = gpar(col = NA), border = TRUE,
              column_title_side = "top", column_names_rot = 45,
              # direction = "horizontal": legend sits at the bottom
              # (heatmap_legend_side="bottom" in save_heatmap) - horizontal
              # keeps the bar compact in HEIGHT (a thin strip) rather than a
              # tall vertical bar eating into the bottom legend row's height
              heatmap_legend_param = list(title = "NES", at = legend_at, direction = "horizontal"))

  list(ht = ht, W = W, H = H)
}

# draws one heatmap to its own PDF/PNG, legends moved BELOW the plot (was: to
# the right) - saves width, which is the whole point here; merge_legend
# combines the NES colourbar and the Hallmark family legend into one row
# instead of two, saving further space for that panel
save_heatmap = function(built, fbase){
  # extra left padding: the leftmost 45deg-rotated column label otherwise gets
  # clipped by the edge of the plotting device (bottom, left, top, right)
  # (order: bottom, left, top, right) - extra left for the rotated column
  # labels; extra right so row labels ending in a full word (e.g. "..._PATHWAY")
  # aren't clipped at the device edge
  draw_padding = unit(c(2, 8, 2, 6), "mm")
  draw_args = list(built$ht, padding = draw_padding, heatmap_legend_side = "bottom",
                   annotation_legend_side = "bottom", merge_legend = TRUE)

  pdf(paste0(fbase, ".pdf"), width = built$W, height = built$H)
  do.call(draw, draw_args); dev.off()
  png(paste0(fbase, ".png"), width = built$W, height = built$H, units = "in", res = 300)
  do.call(draw, draw_args); dev.off()
}

ml = create_path_comp_mat_list(gsea_res_tab, comps_sel = comps_largest, sig_cut = SIG_CUT)

ht_hallmark = build_heatmap(ml[["HALLMARKS"]],     "Hallmark")
ht_green    = build_heatmap(ml[["user_def_sets"]], "Green")

save_heatmap(ht_hallmark, paste0(out_dir, "/", script_ind, "GSEA_heatmap_largest7_Hallmark_FDR10_allsig"))
save_heatmap(ht_green,    paste0(out_dir, "/", script_ind, "GSEA_heatmap_largest7_Green_FDR10_allsig"))

### combined Green and Hallmark heatmaps ------------------------------------
# Absolute viewport dimensions preserve each panel's standalone size.

actual_green_pct = round(100 * ht_green$W / (ht_green$W + ht_hallmark$W), 1)
message("Combined figure: Green/Hallmark width ratio at exact standalone size = ",
       actual_green_pct, "% / ", round(100 - actual_green_pct, 1), "%")

combined_W = ht_green$W + ht_hallmark$W
combined_H = ht_hallmark$H

draw_combined = function(){
  # (order: bottom, left, top, right) - extra left for the rotated column
  # labels; extra right so row labels ending in a full word (e.g. "..._PATHWAY")
  # aren't clipped at the device edge
  draw_padding = unit(c(2, 8, 2, 6), "mm")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2,
      widths  = unit(c(ht_green$W, ht_hallmark$W), "in"),
      heights = unit(c(combined_H - ht_green$H, ht_green$H), "in"))))

  # Green: row 2 (bottom), col 1 - viewport height/width exactly match its
  # own standalone dimensions, so it renders unstretched
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  draw(ht_green$ht, newpage = FALSE, padding = draw_padding,
      heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)
  grid.text("a", x = unit(2, "mm"), y = unit(1, "npc") - unit(2, "mm"),
           just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 16))
  popViewport()

  # Hallmark: rows 1-2, col 2 - full canvas height, its own natural width
  pushViewport(viewport(layout.pos.row = 1:2, layout.pos.col = 2))
  draw(ht_hallmark$ht, newpage = FALSE, padding = draw_padding,
      heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)
  grid.text("b", x = unit(2, "mm"), y = unit(1, "npc") - unit(2, "mm"),
           just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 16))
  popViewport()

  popViewport()
}

pdf(paste0(out_dir, "/", script_ind, "Fig_Green_Hallmark_combined.pdf"), width = combined_W, height = combined_H)
draw_combined(); dev.off()
png(paste0(out_dir, "/", script_ind, "Fig_Green_Hallmark_combined.png"), width = combined_W, height = combined_H, units = "in", res = 300)
draw_combined(); dev.off()

message("Done. Written to ", out_dir)
