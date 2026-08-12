message("\n\n##########################################################################\n",
        "# Start LD_F03e v01: Replot WGCNA Figure 7 ", Sys.time(),
        "\n##########################################################################\n\n")

library(tidyverse)
library(qs)
library(circlize)
library(ComplexHeatmap)
library(viridis)

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

f03c_dir = paste0(main_dir, "LD_F_DESeq_pseudobulk_WGCNA/LD_F03c_v02/")
in_dir   = paste0(main_dir, "LD_F_DESeq_pseudobulk_WGCNA/LD_F03e_v01/")
out_dir  = in_dir
script_ind = "LD_F03e_v01_"

# Modules highlighted in the module-trait panel.
mods_highlight = c("M1", "M5", "M8", "M11", "M12")


### load small inputs only ---------------------------------------------------

replot_bundle = qread(paste0(in_dir, script_ind, "replot_bundle.qs"))
me_mat        = replot_bundle$me_mat
meta          = replot_bundle$meta
mods          = replot_bundle$mods

meta_idx = match(colnames(me_mat), meta$cluster_sample)
if (anyNA(meta_idx)) stop("Metadata are missing for one or more eigengene columns")
meta = meta[meta_idx, ]
if (!identical(as.character(meta$cluster_sample), colnames(me_mat)))
  stop("Eigengene columns and metadata are not aligned")

lmm_path      = read_csv(paste0(in_dir,   script_ind, "Module_eigengene_pathology_LMM.csv"), show_col_types = FALSE)
cor_res       = read_csv(paste0(f03c_dir, "LD_F03c_v02_Module_trait_correlation_results.csv"), show_col_types = FALSE)
# corrected TREM2 LMM (adds (1|BrainBankNetworkIDFormatted): sample = donor x
# brain region, not donor itself, so F03c v02's original (1|sample)-only LMM
# left cross-region correlation within a donor unmodelled - see compute script
trem2_pairs   = read_csv(paste0(in_dir, script_ind, "Module_eigengene_LMM_TREM2_corrected_pairwise.csv"), show_col_types = FALSE)

path_traits = unique(lmm_path$trait)


### helpers -------------------------------------------------------------------

pal = function(v){
  v2 = length(unique(v))
  if (v2 == 2){
    p2 = c("grey20", "dodgerblue")
  } else if (v2 == 3){
    p2 = c("dodgerblue", "grey20", "orange")
  } else if (v2 == 4){
    p2 = c("dodgerblue", "grey20", "orange", "green4")
  } else if (v2 < 6){
    p2 = colorRamps::matlab.like(6)[1:v2]
  } else {
    p2 = colorRamps::matlab.like(v2)
  }
  return(p2)
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
  annot = HeatmapAnnotation(df = annot_df, col = annot_colors, which = annot_dim,
                            annotation_name_side = annotation_name_side[1])
  return(annot)
}

sig_from_padj = function(p){
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", ""))))
}


######################################################################
### 1) Module-trait correlation heatmap: same ComplexHeatmap layout as
### F03c v02's original (cell_fun asterisks, no clustering, same sizing
### formula), just filtered to pathology + TREM2 traits, thesis
### blue/orange scale, and LMM-based (not Pearson) asterisks.
######################################################################

message("\n\n   *Module-trait correlation heatmap (pathology + TREM2 only) \n")

trait_names = c("TREM2_R47H", "TREM2_R62H",
                setdiff(path_traits, c("Age", "PostMortemInterval", "plaque_dens", "pctPHF1PositiveArea")))
trait_names = intersect(trait_names, unique(cor_res$trait))   # keep only traits actually present

mod_names = mods
n_mods    = length(mod_names)
n_traits  = length(trait_names)

cor_mat = matrix(NA_real_, nrow = n_mods, ncol = n_traits,
                 dimnames = list(mod_names, trait_names))
for (i in seq_len(nrow(cor_res))){
  if (cor_res$trait[i] %in% trait_names){
    cor_mat[cor_res$module[i], cor_res$trait[i]] = cor_res$r[i]
  }
}

### re-derive BH-FDR for just the vs-CV pairwise contrasts we actually plot
### here (R47H-CV, R62H-CV). The CSV's own padj_BH was corrected across all
### 3 pairwise contrasts (incl. R62H-R47H, which this heatmap never shows) -
### that's over-conservative for what's being reported on this specific plot,
### so the correction family is narrowed to match exactly the two hypotheses
### this heatmap displays.
trem2_vsCV = trem2_pairs %>%
  filter(contrast %in% c("R47H - CV", "R62H - CV")) %>%
  mutate(padj_BH_vsCV = p.adjust(p.value, method = "BH"))

sig_mat = matrix("", nrow = n_mods, ncol = n_traits, dimnames = list(mod_names, trait_names))
for (mod1 in mod_names){
  # per-variant pairwise contrasts (vs CV), not the shared 2-df omnibus test -
  # otherwise R47H and R62H always show the identical significance, since the
  # omnibus test doesn't distinguish which variant is driving it
  r47h_p = trem2_vsCV$padj_BH_vsCV[trem2_vsCV$module == mod1 & trem2_vsCV$contrast == "R47H - CV"]
  r47h_p = if (length(r47h_p) == 1) r47h_p else NA_real_
  r62h_p = trem2_vsCV$padj_BH_vsCV[trem2_vsCV$module == mod1 & trem2_vsCV$contrast == "R62H - CV"]
  r62h_p = if (length(r62h_p) == 1) r62h_p else NA_real_
  if ("TREM2_R47H" %in% trait_names) sig_mat[mod1, "TREM2_R47H"] = sig_from_padj(r47h_p)
  if ("TREM2_R62H" %in% trait_names) sig_mat[mod1, "TREM2_R62H"] = sig_from_padj(r62h_p)

  for (tr1 in path_traits){
    if (!tr1 %in% trait_names) next
    p1 = lmm_path$padj_BH[lmm_path$module == mod1 & lmm_path$trait == tr1]
    p1 = if (length(p1) == 1) p1 else NA_real_
    sig_mat[mod1, tr1] = sig_from_padj(p1)
  }
}

### significant-results table (padj < 0.05) for citing in the thesis text -
### same LMM p-values driving the asterisks above, just as plain numbers
sig_results_list = list()
for (mod1 in mod_names){
  for (tr1 in trait_names){
    r_val = cor_res$r[cor_res$module == mod1 & cor_res$trait == tr1]
    r_val = if (length(r_val) == 1) r_val else NA_real_

    if (tr1 %in% c("TREM2_R47H", "TREM2_R62H")){
      contrast_lab = if (tr1 == "TREM2_R47H") "R47H - CV" else "R62H - CV"
      row1 = trem2_vsCV[trem2_vsCV$module == mod1 & trem2_vsCV$contrast == contrast_lab, ]
      if (nrow(row1) == 1){
        sig_results_list[[paste(mod1, tr1)]] = data.frame(
          module = mod1, trait = tr1, r = r_val, estimate = row1$estimate,
          n_obs = NA_integer_, pvalue = row1$p.value, padj = row1$padj_BH_vsCV,
          stringsAsFactors = FALSE)
      }
    } else {
      row1 = lmm_path[lmm_path$module == mod1 & lmm_path$trait == tr1, ]
      if (nrow(row1) == 1){
        sig_results_list[[paste(mod1, tr1)]] = data.frame(
          module = mod1, trait = tr1, r = r_val, estimate = row1$estimate,
          n_obs = row1$n_obs, pvalue = row1$pvalue, padj = row1$padj_BH,
          stringsAsFactors = FALSE)
      }
    }
  }
}
sig_results = do.call(rbind, sig_results_list)
sig_results = sig_results[!is.na(sig_results$padj) & sig_results$padj < 0.1, ]
sig_results = sig_results[order(sig_results$padj), ]

write_csv(sig_results, file = paste0(out_dir, script_ind, "Significant_LMM_results_padj_lt_0.1.csv"))
message("\n\n   *LMM results with padj < 0.1 (module-trait heatmap traits only):\n")
print(as.data.frame(sig_results))

# display column labels only (no underscores); matrix dimnames themselves are
# left as-is above since sig_mat/cor_mat indexing above relies on them
trait_display = c(TREM2_R47H = "TREM2 R47H", TREM2_R62H = "TREM2 R62H",
                  plaque_dens = "Plaque density", pct4G8PositiveArea = "4G8+ area",
                  pctAT8PositiveArea = "AT8+ area", pctPHF1PositiveArea = "PHF1+ area",
                  Braak_numeric = "Braak stage", Age = "Age", PostMortemInterval = "PMI")
column_labels_cor = unname(trait_display[trait_names])

# thesis-wide diverging colour convention: blue = negative, orange = positive
mx_cor  = max(abs(cor_mat), na.rm = TRUE)
col_fun_cor = colorRamp2(c(-mx_cor, 0, mx_cor), c("#0072B2", "white", "#E69F00"))

ht_cor = Heatmap(
  cor_mat,
  name              = "Pearson r",
  col               = col_fun_cor,
  cell_fun          = function(j, i, x, y, width, height, fill) {
    grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 8))
  },
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  show_heatmap_legend = FALSE,      # custom legend (with significance key) used instead
  row_title         = "Modules",
  show_row_names    = TRUE,
  show_column_names = TRUE,
  column_labels     = column_labels_cor,
  na_col            = "grey90",
  column_names_rot  = 45,
  column_names_gp   = gpar(fontsize = 11),
  row_names_gp      = gpar(fontsize = 11)
)

# custom legend: colour bar, then the significance key stacked underneath it
# (replaces the column_title text that used to spell this out in words)
lgd_color_cor = Legend(col_fun = col_fun_cor, title = "Pearson r",
                       at = round(c(-mx_cor, -mx_cor / 2, 0, mx_cor / 2, mx_cor), 1))
lgd_sig_cor = Legend(labels = c("*   padj < 0.05", "**  padj < 0.01", "*** padj < 0.001"),
                     title = "Significance (BH-FDR)", type = "grid",
                     legend_gp = gpar(fill = NA, col = NA),
                     grid_height = unit(0, "mm"), grid_width = unit(0, "mm"))
combined_legend_cor = packLegend(lgd_color_cor, lgd_sig_cor, direction = "vertical")

draw_cor_plot = function(newpage = TRUE){
  draw(ht_cor, annotation_legend_list = list(combined_legend_cor), annotation_legend_side = "right",
      newpage = newpage)
  # module-row highlight(s): single box per module, spanning the full plot
  # width, drawn after draw() via decorate_heatmap_body (same technique as
  # the earlier local replot script)
  for (m1 in mods_highlight){
    idx = match(m1, mod_names)
    if (is.na(idx)) next
    y1 = (n_mods - idx) / n_mods
    y2 = (n_mods - idx + 1) / n_mods
    decorate_heatmap_body("Pearson r", {
      grid.rect(x = unit(0.5, "npc"), y = unit((y1 + y2) / 2, "npc"),
               width = unit(1, "npc"), height = unit(y2 - y1, "npc"),
               gp = gpar(fill = NA, col = "black", lwd = 2.5))
    })
  }
}

# width bumped slightly above the plain 2/3 shrink so the rotated leftmost
# column label ("TREM2 R47H") isn't clipped
w_cor_panel = 0.72 * max(n_traits * 0.5 + 3, 8)
h_cor_panel = 0.75 * max(n_mods * 0.4 + 3, 6)

pdf(file  = paste0(out_dir, script_ind, "Module_trait_correlation_pathology_TREM2.pdf"),
    width  = w_cor_panel, height = h_cor_panel)
draw_cor_plot()
dev.off()


######################################################################
### 2) TREM2Variant module eigengene heatmap, top bars reduced to
### cluster_name + TREM2Variant only
######################################################################

message("\n\n   *TREM2Variant eigengene heatmap (2 annotation bars) \n")

pl_mat_X = me_mat
lims_X   = 0.2 * c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

col_anno_df = as.data.frame(meta[, c("cluster_name", "TREM2Variant")])
col_anno_df$cluster_name = gsub("_", " ", col_anno_df$cluster_name)   # display only, no underscores
colnames(col_anno_df) = c("Cluster", "TREM2Variant")   # annotation bar title comes from the column name
for (col1 in names(col_anno_df)){
  col_anno_df[[col1]] = factor(col_anno_df[[col1]], levels = unique(col_anno_df[[col1]]))
}
rownames(col_anno_df) = meta$cluster_sample

col_annot = create_heatmap_annot(annot_df = col_anno_df)

trem2_split  = meta$TREM2Variant[match(colnames(pl_mat_X), meta$cluster_sample)]
trem2_levels = intersect(c("CV", "R47H", "R62H"), unique(as.character(trem2_split)))
trem2_split  = factor(as.character(trem2_split), levels = trem2_levels)

ht_list_trem2 = Heatmap(
  pl_mat_X,
  name = "Z-score",
  col  = colorRamp2(breaks = seq(from = lims_X[1], to = lims_X[2], length.out = 100),
                    colors = viridis(100)),
  column_split          = trem2_split,
  cluster_columns        = FALSE,
  cluster_column_slices  = FALSE,
  column_title           = "%s",
  top_annotation         = col_annot,
  cluster_rows           = TRUE,
  row_title              = "Modules",
  left_annotation         = NULL,
  show_row_names          = TRUE,
  show_column_names       = FALSE,
  width = 15, height = 2.5
)

w_trem2_panel = 10
h_trem2_panel = 5.0

pdf(file = paste0(out_dir, script_ind, "Module_eigengene_heatmap_by_TREM2Variant.pdf"),
    width = w_trem2_panel, height = h_trem2_panel)
{
  draw(ht_list_trem2, heatmap_legend_side = "right", annotation_legend_side = "right",
      merge_legend = TRUE)
}
dev.off()


######################################################################
### 2b) Combined figure: A = TREM2Variant eigengene heatmap, B = module-trait
### correlation heatmap, side by side
######################################################################

message("\n\n   *Combined figure: A (eigengene heatmap) + B (module-trait heatmap) \n")

w_combo = w_trem2_panel + w_cor_panel
h_combo = max(h_trem2_panel, h_cor_panel) + 0.3   # small headroom for the A/B labels

pdf(file = paste0(out_dir, script_ind, "Figure_eigengene_and_module_trait_heatmap_AB.pdf"),
    width = w_combo, height = h_combo)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2,
             widths = unit(c(w_trem2_panel, w_cor_panel), "in"))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_list_trem2, heatmap_legend_side = "right", annotation_legend_side = "right",
    merge_legend = TRUE, newpage = FALSE)
grid.text("a", x = unit(2, "mm"), y = unit(1, "npc") - unit(2, "mm"),
         just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 16))
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw_cor_plot(newpage = FALSE)
grid.text("b", x = unit(2, "mm"), y = unit(1, "npc") - unit(2, "mm"),
         just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 16))
popViewport()

dev.off()
message("\n\nDone. Written to ", out_dir, "\n")
