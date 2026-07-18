# LD_X07: Pseudobulk DESeq2 - DEG threshold sensitivity + pairwise log2FC comparisons.
#
#   PART 1 (decision table): DEG counts per comparison x astrocyte subcluster under
#     four padj x |log2FC| cutoffs, to help choose a DEG threshold.
#       Comparisons: AD vs Ctrl (CV), R62H vs CV (AD), R47H vs CV (AD).
#       Threshold sets: padj < {0.05, 0.1} x |log2FC| > {0, 0.25}.
#       Same DEG definition as LD_E02a2 (LRT padj + raw MLE log2FoldChange); no
#       re-fit, no shrinkage. Reads LD_E02a2_v02_bulk_data.qs ($deseq_results).
#       Scope = PER SUBCLUSTER (the only stored scope for E02a2).
#
#   PART 2 (thesis figures): pairwise log2FC concordance of DESeq2 contrasts, on the
#     POOLED astrocyte pseudobulk (all subclusters, cluster_name as a covariate).
#     Reads the pooled DESeq results checkpointed by LD_E04a ($E04_deseq_res).
#       (2a) MAIN figure - three contrast pairs side by side (E03a2 scatter style,
#            reg_both genes, base covariate model M0):
#              R62H_vs_CV  ~ AD_vs_Ctrl ;  R47H_vs_CV ~ AD_vs_Ctrl ;  R62H_vs_CV ~ R47H_vs_CV
#       (2b) SUPPLEMENTARY figure - the same three pairs as rows x covariate levels as
#            columns (M0 -> +Sex -> +cohort -> +PMI -> +Age -> +Braak), showing the
#            effect survives progressive covariate adjustment.
#
#   NB on thresholds in PART 2: the scatter is restricted to the E03 DEG-union gene
#   universe (E02a2 DEGs, padj < 0.1) and genes are coloured "up"/"down"/"nreg" by
#   NOMINAL p < 0.05 per contrast (classify_reg). This is a descriptive effect-direction
#   concordance plot, NOT an FDR-significant DEG plot - matches E03a2/E04a exactly.
#
#   DATA: LD_E_DESeq_pseudobulk/LD_E02a2_v02_bulk_data.qs        (Part 1; ~3.4 GB)
#         LD_E_DESeq_pseudobulk/LD_E04a_v01_bulk_data.qs         (Part 2; pooled results)
#         LD_E_DESeq_pseudobulk/LD_E02a2_v02_DEGs_by_cluster_genes.csv (Part 2; gene universe)

message("\n\n##########################################################################\n",
        "# Start LD_X07 pseudobulk DESeq2 figures (5-covariate, E02c): ", Sys.time(),
        "\n##########################################################################\n\n")

library(tidyverse)
library(qs)
library(ggrepel)
library(patchwork)

### paths -------------------------------------------------------------------
base       = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
e_out      = file.path(base, "LD_E_DESeq_pseudobulk")
e02_path   = file.path(e_out, "LD_E02a2_v02_bulk_data.qs")
e04_path   = file.path(e_out, "LD_E04a_v01_bulk_data.qs")
deg_csv    = file.path(e_out, "LD_E02a2_v02_DEGs_by_cluster_genes.csv")
clust_csv  = file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X07_"

# E02c (5-covariate) inputs for the DEG-count + summary part (light CSVs; no 3.2 GB load)
e02c_deg_genes_csv = file.path(e_out, "LD_E02c/LD_E02c_v01_DEGs_by_cluster_genes.csv")
tf_csv             = file.path(base, "data_TREM2_michael/A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv")

# Legacy Part 1 (threshold sensitivity, E02a2) and Part 2 (concordance, E04a) still read the
# OLD 3-covariate data, so they are disabled for now. Part 2 will be rebuilt on the 5-covariate
# data next; set RUN_LEGACY = TRUE only to reproduce the old figures.
RUN_LEGACY = FALSE

save_plot = function(p, suffix, w, h) {
  ggsave(file.path(out_dir, paste0(script_ind, suffix, ".pdf")), p, width = w, height = h, useDingbats = FALSE)
  ggsave(file.path(out_dir, paste0(script_ind, suffix, ".png")), p, width = w, height = h, dpi = 300)
}

# cluster order (biological ordering from the cluster-assignment table)
clust_tab     = read_csv(clust_csv, show_col_types = FALSE)
cluster_order = unique(clust_tab$cluster_name)


##########################################################################
# PART A - DEG counts (2-panel, two y-axis versions) + DEG/TF summary CSV
#          5-covariate model (E02c); reads light CSVs, no 3.2 GB object load
##########################################################################

message("\n          *** PART A: DEG counts + TF summary (E02c, 5-covariate)... ", Sys.time(), "\n")

deg_wide = read_csv(e02c_deg_genes_csv, show_col_types = FALSE)
TF       = unique(read_csv(tf_csv, show_col_types = FALSE)$Symbol)

# comparison tags -> readable labels (same convention as E02c); cluster-vs-ref and
# R47H_vs_R62H columns are intentionally not listed, so they are skipped below.
comp_tags = c(
  "TREM2_CV_AD_vs_Control" = "AD vs Control (CV only)",
  "AD_TREM2_R62H_vs_CV"    = "R62H vs CV (AD only)",
  "AD_TREM2_R47H_vs_CV"    = "R47H vs CV (AD only)"
)

# parse a DEG-list column name "{cluster}_{tag}_{up|down}" -> (cluster, label, direction)
parse_deg_col = function(nm){
  if (grepl("_up$", nm))   { direction = "up"   } else
  if (grepl("_down$", nm)) { direction = "down" } else return(NULL)
  key = sub("_up$|_down$", "", nm)
  for (tag in names(comp_tags)){
    if (grepl(tag, key, fixed = TRUE)){
      return(list(cluster = sub(paste0("_", tag, ".*"), "", key),
                  comp_label = comp_tags[[tag]], direction = direction))
    }
  }
  NULL
}

# reconstruct DEG gene lists: store[[comp_label]][[cluster]][[direction]] = gene vector
store = list()
for (nm in colnames(deg_wide)){
  info = parse_deg_col(nm)
  if (is.null(info)) next
  if (!info$cluster %in% cluster_order) next
  g = as.character(deg_wide[[nm]])
  g = g[!is.na(g) & g != ""]
  store[[info$comp_label]][[info$cluster]][[info$direction]] = g
}

# ---- per-cluster TF/other counts (for the 2-panel bar plot) ----
deg_rows = list()
for (comp_label in names(store)){
  for (cl in names(store[[comp_label]])){
    for (dir1 in c("up", "down")){
      g = store[[comp_label]][[cl]][[dir1]]; if (is.null(g)) g = character(0)
      deg_rows[[length(deg_rows)+1]] = tibble(
        cluster = cl, comparison = comp_label, direction = dir1,
        n_TF = sum(g %in% TF), n_other = sum(!(g %in% TF))
      )
    }
  }
}
deg_tab = bind_rows(deg_rows)

# ---- 2-panel data: AD-vs-Control + R62H-vs-CV ----
panel_labels     = c("AD vs Control (CV only)", "R62H vs CV (AD only)")
deg_2p           = deg_tab %>% filter(comparison %in% panel_labels)
clusters_present = cluster_order[cluster_order %in% deg_2p$cluster]

deg_long = deg_2p %>%
  pivot_longer(c(n_TF, n_other), names_to = "gene_cat", values_to = "n_genes") %>%
  mutate(
    n_plot     = ifelse(direction == "up", n_genes, -n_genes),
    cluster    = factor(cluster,    levels = clusters_present),
    comparison = factor(comparison, levels = panel_labels),
    fill_group = factor(paste0(direction, "_", gene_cat),
                        levels = c("up_n_TF", "up_n_other", "down_n_TF", "down_n_other"))
  )

fill_colors = c("up_n_TF" = "#C0392B", "up_n_other" = "#F1948A",
                "down_n_TF" = "#1A5276", "down_n_other" = "#85C1E9")
fill_labels = c("up_n_TF" = "Up (TF)", "up_n_other" = "Up (other)",
                "down_n_TF" = "Down (TF)", "down_n_other" = "Down (other)")

deg_theme = theme_classic(base_size = 11) +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
        strip.background   = element_rect(fill = "grey92", color = NA),
        strip.text         = element_text(face = "bold", size = 10),
        legend.position    = "top",
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3))

deg_base = function() {
  list(geom_col(position = "stack", width = 0.75),
       geom_hline(yintercept = 0, linewidth = 0.4, color = "grey20"),
       scale_fill_manual(values = fill_colors, labels = fill_labels, name = NULL),
       scale_y_continuous(labels = function(x) abs(x), name = "No. of DEGs (padj < 0.1)"),
       scale_x_discrete(name = NULL, drop = FALSE), deg_theme)
}

pdf_w = max(10, length(clusters_present) * 0.55 + 3)

# version 1: DAMPENED-proportional panel heights. Each panel keeps its own readable scale,
# but the height ratio = (DEG-extent ratio)^prop_damp, so fully proportional (prop_damp = 1)
# is softened. prop_damp = 0 -> equal panels; 0.5 -> square-root (gentle); 1 -> fully proportional.
prop_damp = 0.5

ext = deg_long %>%
  group_by(comparison, cluster) %>%
  summarise(up_ext =  sum(n_plot[n_plot > 0]),
            dn_ext = -sum(n_plot[n_plot < 0]), .groups = "drop") %>%
  group_by(comparison) %>%
  summarise(extent = max(up_ext) + max(dn_ext), .groups = "drop")
ad_ext  = ext$extent[ext$comparison == panel_labels[1]]
r62_ext = ext$extent[ext$comparison == panel_labels[2]]
hratio  = (ad_ext / r62_ext) ^ prop_damp

mk_panel = function(lab, show_x){
  p = ggplot(filter(deg_long, comparison == lab),
             aes(cluster, n_plot, fill = fill_group)) + deg_base() + facet_wrap(~ comparison)
  if (!show_x) p = p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}
p_prop = (mk_panel(panel_labels[1], FALSE) / mk_panel(panel_labels[2], TRUE)) +
  plot_layout(heights = c(hratio, 1), guides = "collect") +
  plot_annotation(
    title   = "DEGs per cluster: AD-vs-Control vs R62H-vs-CV (5-covariate)",
    caption = "Panel heights dampened-proportional to DEG number (prop_damp = 0.5).\nUp = above zero, down = below zero. Solid = TF (Fantom5), lighter = other genes.") &
  theme(legend.position = "top")
save_plot(p_prop, "DEG_counts_2panel_proportional", pdf_w, 9)

# version 2: SHARED fixed y-axis (same scale both panels; R62H bars appear small)
p_shared = ggplot(deg_long, aes(cluster, n_plot, fill = fill_group)) + deg_base() +
  facet_wrap(~ comparison, ncol = 1, scales = "fixed") +
  labs(title = "DEGs per cluster: AD-vs-Control vs R62H-vs-CV (5-covariate, shared scale)",
       caption = "Shared y-axis: R62H-vs-CV bars are small relative to AD-vs-Control.\nUp = above zero, down = below zero. Solid = TF (Fantom5), lighter = other genes.")
save_plot(p_shared, "DEG_counts_2panel_sharedscale", pdf_w, 9)

message("    Saved 2-panel DEG-count figures (proportional + shared scale)")

# ---- DEG/TF summary CSV: Overall + per contrast (pooled = union across subclusters) ----
# pooled up/down are unions across subclusters; a gene can be up in one subcluster and
# down in another, so n_up + n_down may exceed n_DEG (which is the union of all DEGs).
pooled_stats = function(comp_label){
  cls  = names(store[[comp_label]])
  up   = unique(unlist(lapply(cls, function(c) store[[comp_label]][[c]][["up"]])))
  down = unique(unlist(lapply(cls, function(c) store[[comp_label]][[c]][["down"]])))
  list(up = up[!is.na(up)], down = down[!is.na(down)])
}
summ_row = function(level, up, down){
  up = unique(up); down = unique(down); all = unique(c(up, down))
  n_tf = length(intersect(all, TF))
  tibble(level = level, n_DEG = length(all), n_up = length(up), n_down = length(down),
         up_down_ratio = ifelse(length(down) > 0, length(up) / length(down), NA_real_),
         n_TF = n_tf, pct_TF = ifelse(length(all) > 0, 100 * n_tf / length(all), NA_real_))
}

csv_contrasts = c("AD vs Control (CV only)" = "AD_vs_Control",
                  "R62H vs CV (AD only)"    = "R62H_vs_CV",
                  "R47H vs CV (AD only)"    = "R47H_vs_CV")
rows = list(); ov_up = character(0); ov_down = character(0)
for (lab in names(csv_contrasts)){
  st = pooled_stats(lab)
  ov_up = c(ov_up, st$up); ov_down = c(ov_down, st$down)
  rows[[length(rows)+1]] = summ_row(csv_contrasts[[lab]], st$up, st$down)
}
deg_tf_summary = bind_rows(summ_row("Overall", ov_up, ov_down), bind_rows(rows))

write_csv(deg_tf_summary, file.path(out_dir, paste0(script_ind, "DEG_TF_summary.csv")))
message("    Saved: ", file.path(out_dir, paste0(script_ind, "DEG_TF_summary.csv")))
message("\n##### DEG / TF summary (5-covariate, pooled across subclusters) #####\n")
print(as.data.frame(deg_tf_summary), row.names = FALSE)


##########################################################################
# PART B - pairwise log2FC concordance scatters (5-covariate)
#   B1 pooled MAIN (3 pairs, M0)           from E04c $E04_deseq_res
#   B2 pooled incremental GRID (M0 .. +Braak) from E04c
#   B3 per-subcluster MULTIPAGE (3 pages)  from E02c $deseq_results (M0 / 5-cov)
# reg_both genes (padj < 0.1 in BOTH contrasts), classified concordant/discordant,
# with a fitted regression line and Pearson r.
##########################################################################

message("\n          *** PART B: log2FC concordance scatters... ", Sys.time(), "\n")

e04c_path = file.path(e_out, "LD_E04c_bulk_data.qs")
e02c_path = file.path(e_out, "LD_E02c/LD_E02c_v01_bulk_data.qs")
DOT       = 0.6   # point size (was 2 in the old version)

# three scatter pairs (y vs x); contrast names as stored in E04c
pairs = list(
  list(y = "R62H_vs_CV", x = "CV_AD_vs_Control", ylab = "R62H vs CV", xlab = "AD vs Control"),
  list(y = "R47H_vs_CV", x = "CV_AD_vs_Control", ylab = "R47H vs CV", xlab = "AD vs Control"),
  list(y = "R62H_vs_CV", x = "R47H_vs_CV",       ylab = "R62H vs CV", xlab = "R47H vs CV")
)
# E02c per-subcluster comparison tags (key = "<cluster>_<tag>")
e02c_tag = c(CV_AD_vs_Control = "TREM2_CV_AD_vs_Control",
             R62H_vs_CV       = "AD_TREM2_R62H_vs_CV",
             R47H_vs_CV       = "AD_TREM2_R47H_vs_CV")
# original E03a2 colours: concordant up=red, concordant down=blue, discordant=magenta3
reg_levels = c("down_down", "up_up", "down_up", "up_down",
               "down_nreg", "up_nreg", "nreg_down", "nreg_up", "nreg_nreg")
reg_cols = c(down_down = "blue", up_up = "red", down_up = "magenta3", up_down = "magenta3",
             down_nreg = "grey40", up_nreg = "grey40", nreg_down = "grey40",
             nreg_up = "grey40", nreg_nreg = "grey")
covar_model  = "cohort + APOE + CD33 + BrainRegion + Sex"

norm_res = function(res){
  res = as.data.frame(res)
  g = if ("gene" %in% names(res)) as.character(res$gene) else rownames(res)
  tibble(gene = g, log2FC = res$log2FoldChange, padj = res$padj, pval = res$pvalue)
}
# reg_both table: genes significant in BOTH contrasts, classified into reg_group by
# sign exactly as E03a2 so the colours match the original. Significance criterion is
# parameterised: pooled scatters use FDR (padj < 0.1); the per-subcluster trend-check
# uses nominal p < 0.05 (as the original), because single subclusters lack the power
# for genes to reach FDR < 0.1 in both contrasts.
build_pair = function(res_y, res_x, sig_col = "padj", sig_cut = 0.1){
  empty = tibble(gene = character(), log2FC_y = numeric(), log2FC_x = numeric(),
                 reg_group = factor(character(), levels = reg_levels))
  if (is.null(res_y) || is.null(res_x)) return(empty)
  d = inner_join(norm_res(res_y), norm_res(res_x), by = "gene", suffix = c("_y", "_x"))
  d$sig_y = d[[paste0(sig_col, "_y")]]
  d$sig_x = d[[paste0(sig_col, "_x")]]
  d = d[!is.na(d$log2FC_y) & !is.na(d$log2FC_x) & !is.na(d$sig_y) & !is.na(d$sig_x), ]
  d = d[d$sig_y < sig_cut & d$sig_x < sig_cut, ]
  if (nrow(d) == 0) return(empty)
  reg_y = ifelse(d$log2FC_y > 0, "up", "down")
  reg_x = ifelse(d$log2FC_x > 0, "up", "down")
  d$reg_group = factor(paste0(reg_y, "_", reg_x), levels = reg_levels)
  d
}
safe_cor = function(x, y) if (length(x) >= 3) suppressWarnings(cor(x, y)) else NA_real_

# shared scatter layers, matching the original E03a2 aesthetic (colours, dashed
# +/-log2(1.2) guides, solid 0 lines, grey lm line, theme_minimal). Only agreed
# changes: smaller points (DOT) and a single overall lm line (matches the single r).
scatter_base = function(dot = DOT){
  list(geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2),
       geom_hline(yintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2),
       geom_vline(xintercept = 0, linewidth = 0.3, color = "grey30"),
       geom_hline(yintercept = 0, linewidth = 0.3, color = "grey30"),
       geom_smooth(method = "lm", formula = y ~ x, color = "grey30", linewidth = 0.5, se = FALSE),
       geom_point(aes(color = reg_group), size = dot, alpha = 0.8),
       scale_color_manual(limits = reg_levels, values = reg_cols, name = "regulation"),
       theme_minimal(base_size = 10))
}

### B1 + B2: pooled, from E04c ------------------------------------------------
e04 = qread(e04c_path)
res_pooled = function(level, contrast) e04$E04_deseq_res[[paste("pooled", level, contrast, sep = "|")]]

# B1: main figure, M0, three pairs side by side
b1 = lapply(pairs, function(p){
  d = build_pair(res_pooled("M0_base", p$y), res_pooled("M0_base", p$x))
  lim = if (nrow(d) > 0) max(abs(c(d$log2FC_x, d$log2FC_y)), na.rm = TRUE) else 1
  ggplot(d, aes(log2FC_x, log2FC_y)) + scatter_base() +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +   # symmetric axes, as original
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, size = 3,
             label = sprintf("r = %.2f, n = %d", safe_cor(d$log2FC_x, d$log2FC_y), nrow(d))) +
    labs(x = paste0("log2FC (", p$xlab, ")"), y = paste0("log2FC (", p$ylab, ")"),
         title = paste0(p$ylab, " vs ", p$xlab))
})
p_main = wrap_plots(b1, nrow = 1, guides = "collect") +
  plot_annotation(
    title   = "Pairwise log2FC concordance (pooled astrocytes, 5-covariate)",
    caption = paste0("Pooled DESeq2; covariates: ", covar_model,
                     " + subcluster (+ contrast variable). reg_both genes (padj < 0.1 in both).")) &
  theme(legend.position = "bottom")
save_plot(p_main, "log2FC_concordance_pooled_main", 14, 5)

# B2: incremental-covariate grid (rows = pairs, cols = levels)
levels_ord = c("M0_base", "M1_Age", "M2_Age_PMI", "M3_Age_PMI_Braak")
level_lab  = c(M0_base = "M0", M1_Age = "+Age", M2_Age_PMI = "+PMI", M3_Age_PMI_Braak = "+Braak")
pair_labs  = vapply(pairs, function(p) paste0(p$ylab, " ~ ", p$xlab), character(1))
grid_df = bind_rows(lapply(pairs, function(p){
  bind_rows(lapply(levels_ord, function(L){
    d = build_pair(res_pooled(L, p$y), res_pooled(L, p$x))
    if (nrow(d) == 0) return(NULL)
    d$pair = paste0(p$ylab, " ~ ", p$xlab); d$level = L; d
  }))
}))
grid_df$level = factor(level_lab[grid_df$level], levels = unname(level_lab))
grid_df$pair  = factor(grid_df$pair, levels = pair_labs)
r_grid = grid_df %>% group_by(pair, level) %>%
  summarise(r = safe_cor(log2FC_x, log2FC_y), .groups = "drop")
p_grid = ggplot(grid_df, aes(log2FC_x, log2FC_y)) + scatter_base(dot = 0.35) +
  geom_text(data = r_grid, aes(x = -Inf, y = Inf, label = sprintf("r=%.2f", r)),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.4, size = 2.7) +
  facet_grid(pair ~ level) +
  labs(x = "log2FC (x-axis contrast)", y = "log2FC (y-axis contrast)",
       title = "log2FC concordance under incremental covariate adjustment (pooled)",
       caption = paste0("M0 = ", covar_model,
                        " + subcluster; columns add one covariate cumulatively. reg_both genes (padj < 0.1 in both).")) +
  theme(legend.position = "bottom")
save_plot(p_grid, "log2FC_concordance_pooled_covariate_grid", 15, 9)
rm(e04); gc()

### B3: per-subcluster multipage, from E02c ----------------------------------
message("\n          *** PART B3: per-subcluster scatters (loading E02c ~3.2 GB)... ", Sys.time(), "\n")
e02  = qread(e02c_path)
dres = e02$deseq_results

pdf(file.path(out_dir, paste0(script_ind, "log2FC_concordance_per_subcluster.pdf")),
    width = 12, height = 10)
for (p in pairs){
  cls = cluster_order[ vapply(cluster_order, function(cl)
    !is.null(dres[[paste0(cl, "_", e02c_tag[p$y])]]) &&
    !is.null(dres[[paste0(cl, "_", e02c_tag[p$x])]]), logical(1)) ]
  per = bind_rows(lapply(cls, function(cl){
    # per-subcluster trend-check uses nominal p < 0.05 (as the original); FDR<0.1
    # leaves too few genes per single subcluster to see the trend.
    d = build_pair(dres[[paste0(cl, "_", e02c_tag[p$y])]], dres[[paste0(cl, "_", e02c_tag[p$x])]],
                   sig_col = "pval", sig_cut = 0.05)
    if (nrow(d) == 0) return(NULL)
    d$cluster = cl; d
  }))
  if (nrow(per) == 0) next
  per$cluster = factor(per$cluster, levels = cluster_order)
  r_cl = per %>% group_by(cluster) %>% summarise(r = safe_cor(log2FC_x, log2FC_y), .groups = "drop")
  pg = ggplot(per, aes(log2FC_x, log2FC_y)) + scatter_base() +
    geom_text(data = r_cl, aes(x = -Inf, y = Inf, label = sprintf("r=%.2f", r)),
              inherit.aes = FALSE, hjust = -0.1, vjust = 1.4, size = 2.7) +
    facet_wrap(~ cluster, scales = "free") +
    labs(x = paste0("log2FC (", p$xlab, ")"), y = paste0("log2FC (", p$ylab, ")"),
         title = paste0(p$ylab, " vs ", p$xlab, " - per subcluster (5-covariate)"),
         caption = "Per-subcluster DESeq2 (E02c). reg_both genes (nominal p < 0.05 in both contrasts).") +
    theme(legend.position = "bottom")
  print(pg)
}
dev.off()
message("    Saved per-subcluster multipage PDF (one page per contrast pair)")


##########################################################################
# PART C - volcano plots (same aesthetic as the original E02c volcano, with
#   slightly smaller dots). Two multipage PDFs, one volcano per page:
#   C1 per subcluster x contrast (incl. cluster-vs-ref), from E02c
#   C2 pooled, 4 contrasts at the 5-covariate model (M0), from E04c
##########################################################################
message("\n          *** PART C: volcano plots... ", Sys.time(), "\n")

VDOT_OTHER = 0.7   # non-DEG point size (original = 1)
VDOT_DEG   = 1.4   # DEG point size      (original = 2)

# one volcano for a single results table; replicates the original E02c volcano
volcano_plot = function(res, title){
  t1 = as.data.frame(res)
  t1$gene = if ("gene" %in% names(t1)) as.character(t1$gene) else rownames(t1)
  t1 = t1[!is.na(t1$padj), ]
  if (nrow(t1) == 0) return(NULL)
  t1$neglog10p = -log10(t1$pvalue)
  t1$DEG = t1$padj <= 0.1
  t1$gene_cat = "Other"
  t1$gene_cat[t1$DEG & t1$log2FoldChange > 0] = "up"
  t1$gene_cat[t1$DEG & t1$log2FoldChange < 0] = "down"
  # label top 20 by p-value, top 20 by +log2FC, top 20 by -log2FC (as the original)
  lab = unique(c(t1$gene[order(t1$pvalue)][1:20],
                 t1$gene[order(-t1$log2FoldChange)][1:20],
                 t1$gene[order(t1$log2FoldChange)][1:20]))
  lab = lab[!is.na(lab)]
  t1$plot_label = ifelse(t1$gene %in% lab, t1$gene, "")
  ggplot(t1, aes(x = log2FoldChange, y = neglog10p, color = gene_cat)) +
    geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2) +
    geom_hline(yintercept = -log10(0.1), linewidth = 0.3, color = "grey30", linetype = 2) +
    geom_point(aes(size = DEG), alpha = 0.8) +
    geom_label_repel(aes(label = plot_label), seed = 42, min.segment.length = 0,
                     max.overlaps = Inf, max.time = 5) +
    scale_size_manual(limits = c(FALSE, TRUE), values = c(VDOT_OTHER, VDOT_DEG)) +
    scale_color_manual(limits = c("up", "Other", "down"), values = c("red", "grey40", "blue")) +
    theme_minimal() +
    labs(title = title, x = "log2 fold change", y = "-log10(p)")
}

# C1: per subcluster x contrast (reuse the E02c object already loaded)
pdf(file.path(out_dir, paste0(script_ind, "Volcano_per_subcluster.pdf")), width = 10, height = 8)
for (comp in names(dres)){
  v = volcano_plot(dres[[comp]], comp)
  if (!is.null(v)) print(v)
}
dev.off()
message("    Saved per-subcluster volcano PDF (", length(dres), " comparisons)")
rm(e02, dres); gc()

# C2: pooled, 4 contrasts at the 5-covariate model (M0); reload the small E04c object
e04v = qread(e04c_path)
pooled_contrasts = c(CV_AD_vs_Control = "AD vs Control", R62H_vs_CV = "R62H vs CV",
                     R47H_vs_CV = "R47H vs CV", R47H_vs_R62H = "R47H vs R62H")
pdf(file.path(out_dir, paste0(script_ind, "Volcano_pooled.pdf")), width = 10, height = 8)
for (cn in names(pooled_contrasts)){
  v = volcano_plot(e04v$E04_deseq_res[[paste("pooled", "M0_base", cn, sep = "|")]],
                   paste0("Pooled: ", pooled_contrasts[[cn]]))
  if (!is.null(v)) print(v)
}
dev.off()
message("    Saved pooled volcano PDF (", length(pooled_contrasts), " contrasts)")
rm(e04v); gc()


if (RUN_LEGACY) {
##########################################################################
# PART 1 - DEG threshold-sensitivity table (per subcluster, from E02a2)
##########################################################################

message("\n          *** PART 1: loading E02a2 DESeq results... ", Sys.time(), "\n")
bulk_e02      = qread(e02_path)
deseq_results = bulk_e02$deseq_results
rm(bulk_e02); gc()                                   # drop the heavy object (carries dds)
stopifnot(length(deseq_results) > 0)

# the three comparisons, with the exact tag used in the stored result names
# ( result keys follow "{cluster}_{tag}" )
comp_tags = c(
  "TREM2_CV_AD_vs_Control" = "AD vs Ctrl (CV)",
  "AD_TREM2_R62H_vs_CV"    = "R62H vs CV (AD)",
  "AD_TREM2_R47H_vs_CV"    = "R47H vs CV (AD)"
)

# threshold sets
sets = tibble(
  set_label = c("padj<0.05 & |log2FC|>0",
                "padj<0.05 & |log2FC|>0.25",
                "padj<0.1 & |log2FC|>0",
                "padj<0.1 & |log2FC|>0.25"),
  padj_cut  = c(0.05, 0.05, 0.10, 0.10),
  lfc_cut   = c(0,    0.25, 0,    0.25)
)

# for one results table + one threshold set: n_DEG / n_up / n_down / ratio
count_set = function(res, padj_cut, lfc_cut) {
  res = res[!is.na(res$padj), ]
  sig = res$padj < padj_cut
  n_up   = sum(sig & res$log2FoldChange >  lfc_cut)
  n_down = sum(sig & res$log2FoldChange < -lfc_cut)
  tibble(n_DEG = n_up + n_down,
         n_up = n_up, n_down = n_down,
         up_down_ratio = if (n_down > 0) round(n_up / n_down, 2) else NA_real_)
}

message("\n          *** PART 1: tabulating threshold sensitivity... ", Sys.time(), "\n")
rows = list()
for (tag in names(comp_tags)) {
  for (cl in cluster_order) {
    key = paste0(cl, "_", tag)
    if (is.null(deseq_results[[key]])) next          # comparison absent for this cluster
    res = deseq_results[[key]]
    for (i in seq_len(nrow(sets))) {
      cnt = count_set(res, sets$padj_cut[i], sets$lfc_cut[i])
      rows[[length(rows) + 1]] = tibble(
        comparison = comp_tags[[tag]],
        cluster    = cl,
        set_label  = sets$set_label[i],
        padj_cut   = sets$padj_cut[i],
        lfc_cut    = sets$lfc_cut[i]
      ) %>% bind_cols(cnt)
    }
  }
}

tab = bind_rows(rows) %>%
  mutate(comparison = factor(comparison, levels = unname(comp_tags)),
         cluster    = factor(cluster, levels = cluster_order),
         set_label  = factor(set_label, levels = sets$set_label)) %>%
  arrange(comparison, cluster, set_label)

write_csv(tab, file.path(out_dir, paste0(script_ind, "DEG_threshold_sensitivity.csv")))
message("      wrote per-subcluster table: ", nrow(tab), " rows")

# union across subclusters (NOT a pooled model): unique genes passing a set in ANY subcluster
union_rows = list()
for (tag in names(comp_tags)) {
  keys = paste0(cluster_order, "_", tag)
  keys = keys[keys %in% names(deseq_results)]
  for (i in seq_len(nrow(sets))) {
    up_g = character(0); dn_g = character(0)
    for (key in keys) {
      res = deseq_results[[key]]
      res = res[!is.na(res$padj) & res$padj < sets$padj_cut[i], ]
      up_g = union(up_g, res$gene[res$log2FoldChange >  sets$lfc_cut[i]])
      dn_g = union(dn_g, res$gene[res$log2FoldChange < -sets$lfc_cut[i]])
    }
    union_rows[[length(union_rows) + 1]] = tibble(
      comparison    = comp_tags[[tag]],
      set_label     = sets$set_label[i],
      padj_cut      = sets$padj_cut[i],
      lfc_cut       = sets$lfc_cut[i],
      n_unique_DEG  = length(union(up_g, dn_g)),
      n_unique_up   = length(up_g),
      n_unique_down = length(dn_g)
    )
  }
}
union_tab = bind_rows(union_rows) %>%
  mutate(comparison = factor(comparison, levels = unname(comp_tags)),
         set_label  = factor(set_label, levels = sets$set_label)) %>%
  arrange(comparison, set_label)
write_csv(union_tab, file.path(out_dir, paste0(script_ind, "DEG_threshold_sensitivity_union_across_subclusters.csv")))
message("      wrote union-across-subclusters table: ", nrow(union_tab), " rows")

# figure: tile heatmap of n_DEG (cluster x set, faceted by comparison)
message("\n          *** PART 1: plotting threshold-sensitivity tiles... ", Sys.time(), "\n")
p_tile = ggplot(tab, aes(set_label, cluster, fill = n_DEG)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = n_DEG), size = 2.4) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "n DEG") +
  facet_wrap(~ comparison, nrow = 1) +
  scale_y_discrete(limits = rev(levels(tab$cluster))) +
  labs(title = "DEG threshold sensitivity (per subcluster)",
       x = NULL, y = NULL,
       caption = "Same DEG definition as LD_E02a2 (LRT padj + raw log2FC); per-subcluster pseudobulk DESeq2") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5, face = "bold"),
        plot.caption = element_text(size = 7),
        panel.grid = element_blank())
save_plot(p_tile, "DEG_threshold_sensitivity_tiles", 11, 7)

rm(deseq_results); gc()


##########################################################################
# PART 2 - pairwise log2FC concordance (POOLED astrocytes, from E04a)
##########################################################################

message("\n          *** PART 2: loading E04a pooled DESeq results... ", Sys.time(), "\n")
bulk_e04  = qread(e04_path)
deseq_res = bulk_e04$E04_deseq_res                    # keys: "scope|level|contrast"
rm(bulk_e04); gc()
stopifnot(length(deseq_res) > 0)

# gene universe = exactly the DEG-union E03 plotted (E02a2 DEGs, padj < 0.1)
deg_tab      = read_csv(deg_csv, show_col_types = FALSE)
deg_universe = unique(unlist(deg_tab, use.names = FALSE))
deg_universe = deg_universe[!is.na(deg_universe) & deg_universe != ""]
message("      gene universe (E03 DEG-union): ", length(deg_universe), " genes")

# covariate levels present for the pooled scope, in incremental order
keys_split  = strsplit(names(deseq_res), "|", fixed = TRUE)
pooled_lv   = unique(vapply(keys_split, function(x) if (x[1] == "pooled") x[2] else NA_character_, character(1)))
pooled_lv   = pooled_lv[!is.na(pooled_lv)]
lv_order    = c("M0_base", "M1_Sex", "M2_Sex_cohort", "M3_Sex_cohort_PMI", "M4_..Age", "M5_..Braak")
covar_levels = lv_order[lv_order %in% pooled_lv]
main_level   = "M0_base"
# short, human-readable column labels for the covariate ladder
lv_labels = c(M0_base = "base", M1_Sex = "+Sex", M2_Sex_cohort = "+cohort",
              M3_Sex_cohort_PMI = "+PMI", M4_..Age = "+Age", M5_..Braak = "+Braak")

# the three contrast pairs (comp1 = y-axis, comp_ref = x-axis), in the requested order
pairs = list(
  "R62H_vs_CV ~ AD_vs_Ctrl" = c("R62H_vs_CV", "CV_AD_vs_Control"),
  "R47H_vs_CV ~ AD_vs_Ctrl" = c("R47H_vs_CV", "CV_AD_vs_Control"),
  "R62H_vs_CV ~ R47H_vs_CV" = c("R62H_vs_CV", "R47H_vs_CV")
)

# regulation classification + colours (verbatim from E03a2/E04a; nominal p < 0.05)
classify_reg = function(res){
  reg = rep("nreg", nrow(res))
  reg[!is.na(res$pvalue) & res$pvalue < 0.05 & res$log2FoldChange > 0] = "up"
  reg[!is.na(res$pvalue) & res$pvalue < 0.05 & res$log2FoldChange < 0] = "down"
  reg
}
reg_levels = c("down_down", "up_up", "down_up", "up_down",
               "down_nreg", "up_nreg", "nreg_down", "nreg_up", "nreg_nreg")
reg_cols = c("down_down" = "blue", "up_up" = "red", "down_up" = "magenta3",
             "up_down" = "magenta3", "down_nreg" = "grey40", "up_nreg" = "grey40",
             "nreg_down" = "grey40", "nreg_up" = "grey40", "nreg_nreg" = "grey")

# combined table for one (scope, pair) over the requested covariate levels
build_pair_tab = function(scope, pair, levels){
  comp1 = pair[1]; comp_ref = pair[2]
  out = NULL
  for (level in levels){
    r1 = deseq_res[[paste(scope, level, comp1,   sep = "|")]]
    r2 = deseq_res[[paste(scope, level, comp_ref, sep = "|")]]
    if (is.null(r1) || is.null(r2)) next
    g = intersect(intersect(r1$gene, r2$gene), deg_universe)   # E03 DEG-union restriction
    r1 = r1[match(g, r1$gene), ]; r2 = r2[match(g, r2$gene), ]
    t = tibble(level = level, gene = g,
               log2FC = r1$log2FoldChange, pval = r1$pvalue,
               log2FC_ref = r2$log2FoldChange, pval_ref = r2$pvalue)
    t = t[!is.na(t$log2FC) & !is.na(t$log2FC_ref) & !is.na(t$pval) & !is.na(t$pval_ref), ]
    if (nrow(t) == 0) next
    t$reg       = classify_reg(tibble(log2FoldChange = t$log2FC,     pvalue = t$pval))
    t$reg_ref   = classify_reg(tibble(log2FoldChange = t$log2FC_ref, pvalue = t$pval_ref))
    t$reg_group = paste0(t$reg, "_", t$reg_ref)
    out = rbind(out, t)
  }
  out
}

### ---- 2a: MAIN figure - three pairs side by side (E03a2 style, size = 2) ----
message("\n          *** PART 2a: main pairwise log2FC scatter (pooled, ", main_level, ")... ", Sys.time(), "\n")

# one E03a2-style scatter for a single pair (reg_both genes), per-panel symmetric axes
make_main_scatter = function(pair_name, pair){
  tab1 = build_pair_tab("pooled", pair, main_level)
  if (is.null(tab1) || nrow(tab1) == 0) return(NULL)
  d = tab1[tab1$reg != "nreg" & tab1$reg_ref != "nreg", ]   # reg_both
  if (nrow(d) < 3) return(NULL)
  # label up to 10 genes furthest from origin per regulation group
  d$label = ""
  for (rg in reg_levels){
    idx = which(d$reg_group == rg)
    if (!length(idx)) next
    dist = sqrt(d$log2FC[idx]^2 + d$log2FC_ref[idx]^2)
    top  = idx[order(-dist)][seq_len(min(10, length(idx)))]
    d$label[top] = d$gene[top]
  }
  d = d[order(match(d$reg_group, reg_levels)), ]           # concordant genes to foreground
  lim = max(abs(c(d$log2FC, d$log2FC_ref)), na.rm = TRUE)
  ggplot(d, aes(x = log2FC_ref, y = log2FC, color = reg_group)) +
    geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2) +
    geom_hline(yintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey30") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey30") +
    geom_smooth(color = "grey30", method = "lm", formula = y ~ x) +
    geom_point(alpha = 0.8, size = 2) +
    geom_label_repel(aes(label = label), seed = 42, size = 2.4,
                     min.segment.length = 0.2, max.overlaps = Inf, max.time = 5) +
    scale_color_manual(limits = reg_levels, values = reg_cols, name = "regulation\n(nominal p<0.05)") +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    theme_minimal() +
    labs(title = pair_name,
         x = paste0("log2FC  ", pair[2]),
         y = paste0("log2FC  ", pair[1]))
}

main_plots = Filter(Negate(is.null), imap(pairs, ~ make_main_scatter(.y, .x)))
p_main = wrap_plots(main_plots, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Pairwise log2FC concordance - pooled astrocytes (base model M0)",
                  caption = "reg_both genes (nominal p<0.05 in both contrasts); E03 DEG-union universe (E02a2 padj<0.1). Grey dashed = +/-log2(1.2).")
save_plot(p_main, "log2FC_pairwise_pooled_main", 16, 6)

### ---- 2b: SUPPLEMENTARY figure - pairs (rows) x covariate levels (cols) ----
message("\n          *** PART 2b: covariate-robustness grid (pooled)... ", Sys.time(), "\n")

# combined reg_both table across all pairs x all covariate levels
grid_d = NULL
for (pn in names(pairs)){
  tb = build_pair_tab("pooled", pairs[[pn]], covar_levels)
  if (is.null(tb) || nrow(tb) == 0) next
  tb = tb[tb$reg != "nreg" & tb$reg_ref != "nreg", ]      # reg_both
  if (nrow(tb) == 0) next
  tb$pair = pn
  grid_d  = rbind(grid_d, tb)
}
grid_d$pair  = factor(grid_d$pair,  levels = names(pairs))
grid_d$level = factor(grid_d$level, levels = covar_levels)

# label up to 6 genes furthest from origin per (pair, level) facet
grid_d$label = ""
for (pn in levels(grid_d$pair)) for (lv in levels(grid_d$level)){
  idx = which(grid_d$pair == pn & grid_d$level == lv)
  if (!length(idx)) next
  dist = sqrt(grid_d$log2FC[idx]^2 + grid_d$log2FC_ref[idx]^2)
  top  = idx[order(-dist)][seq_len(min(6, length(idx)))]
  grid_d$label[top] = grid_d$gene[top]
}
grid_d = grid_d[order(match(grid_d$reg_group, reg_levels)), ]

# per-facet slope / r / n annotation (base-R split; version-safe)
grid_stat = do.call(rbind, lapply(split(grid_d, list(grid_d$pair, grid_d$level), drop = TRUE), function(x){
  tibble(pair = x$pair[1], level = x$level[1], n = nrow(x),
         slope = if (nrow(x) >= 3) unname(coef(lm(log2FC ~ log2FC_ref, x))[2]) else NA_real_,
         r     = if (nrow(x) >= 3) suppressWarnings(cor(x$log2FC_ref, x$log2FC)) else NA_real_)
}))
grid_stat$pair  = factor(grid_stat$pair,  levels = names(pairs))
grid_stat$level = factor(grid_stat$level, levels = covar_levels)
grid_stat$lab   = paste0("n=", grid_stat$n, "\nslope=", round(grid_stat$slope, 2), "\nr=", round(grid_stat$r, 2))

lim = max(abs(c(grid_d$log2FC, grid_d$log2FC_ref)), na.rm = TRUE)

p_grid = ggplot(grid_d, aes(x = log2FC_ref, y = log2FC, colour = reg_group)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_smooth(colour = "grey30", method = "lm", formula = y ~ x, se = TRUE) +
  geom_point(alpha = 0.8, size = 2) +
  geom_label_repel(aes(label = label), seed = 42, size = 2.2,
                   min.segment.length = 0.2, max.overlaps = Inf, max.time = 3) +
  geom_text(data = grid_stat, aes(x = -lim, y = lim, label = lab),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.6, colour = "grey20") +
  scale_colour_manual(limits = reg_levels, values = reg_cols, name = "regulation\n(nominal p<0.05)") +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  facet_grid(pair ~ level, labeller = labeller(level = lv_labels)) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 0)) +
  labs(title = "log2FC concordance is robust to progressive covariate adjustment (pooled astrocytes)",
       x = "log2FC (x-axis contrast)", y = "log2FC (y-axis contrast)",
       caption = "reg_both genes (nominal p<0.05 in both contrasts). Columns add one covariate left to right (base = APOE+CD33+Region+cluster).")

save_plot(p_grid, "log2FC_pairwise_pooled_covariate_grid", 17, 9)

} # end if (RUN_LEGACY)

message("\n\n##########################################################################\n",
        "# Completed LD_X07 ", Sys.time(),
        "\n##########################################################################\n\n")

sessionInfo()
