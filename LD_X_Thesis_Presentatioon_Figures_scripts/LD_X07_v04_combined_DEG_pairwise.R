# LD_X07_v04: Combined DEG-count and pooled log2FC concordance figure.

library(tidyverse)
library(qs)
library(patchwork)

# Accept either project-directory spelling.
base       = c("/rds/general/user/lvd25/home/AST_snRNAseq_TREM2",
               "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2")
base       = base[dir.exists(base)]
stopifnot("no project directory found at either sn/sc path" = length(base) >= 1)
base       = base[1]
e_out      = file.path(base, "LD_E_DESeq_pseudobulk")
e04c_path  = file.path(e_out, "LD_E04c_bulk_data.qs")
e02c_deg_genes_csv = file.path(e_out, "LD_E02c/LD_E02c_v01_DEGs_by_cluster_genes.csv")
clust_csv  = file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
script_ind = "LD_X07_v04_"
for (p in c(e04c_path, e02c_deg_genes_csv, clust_csv)) if (!file.exists(p)) stop("Missing input: ", p)

clust_tab     = read_csv(clust_csv, show_col_types = FALSE)
cluster_order = unique(clust_tab$cluster_name)

### Panel A: DEG counts ------------------------------------------------------

deg_wide = read_csv(e02c_deg_genes_csv, show_col_types = FALSE)
# Transcription factors are included in the up/down totals.

comp_tags = c(
  "TREM2_CV_AD_vs_Control" = "AD vs Control (CV only)",
  "AD_TREM2_R62H_vs_CV"    = "R62H vs CV (AD only)",
  "AD_TREM2_R47H_vs_CV"    = "R47H vs CV (AD only)"
)

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

store = list()
for (nm in colnames(deg_wide)){
  info = parse_deg_col(nm)
  if (is.null(info)) next
  if (!info$cluster %in% cluster_order) next
  g = as.character(deg_wide[[nm]])
  g = g[!is.na(g) & g != ""]
  store[[info$comp_label]][[info$cluster]][[info$direction]] = g
}

deg_rows = list()
for (comp_label in names(store)){
  for (cl in names(store[[comp_label]])){
    for (dir1 in c("up", "down")){
      g = store[[comp_label]][[cl]][[dir1]]; if (is.null(g)) g = character(0)
      deg_rows[[length(deg_rows)+1]] = tibble(
        cluster = cl, comparison = comp_label, direction = dir1, n_genes = length(g)
      )
    }
  }
}
deg_tab = bind_rows(deg_rows)

panel_labels     = c("AD vs Control (CV only)", "R62H vs CV (AD only)")
deg_2p           = deg_tab %>% filter(comparison %in% panel_labels)
clusters_present = cluster_order[cluster_order %in% deg_2p$cluster]

deg_long = deg_2p %>%
  mutate(
    n_plot     = ifelse(direction == "up", n_genes, -n_genes),
    cluster    = factor(cluster,    levels = clusters_present),
    comparison = factor(comparison, levels = panel_labels),
    direction  = factor(direction, levels = c("up", "down"))
  )

# up/down colours matched to the rest of the thesis's up/down convention
# (orange = up, blue = down - same pair as the Green24 diverging scale)
fill_colors = c(up = "#E69F00", down = "#0072B2")
fill_labels = c(up = "Up", down = "Down")

deg_theme_compact = theme_classic(base_size = 8) +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y        = element_text(size = 6),
        axis.title         = element_text(size = 9.5),
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold", size = 8.5),
        legend.position    = "bottom",
        legend.text        = element_text(size = 6.5),
        legend.title        = element_blank(),
        legend.key.size    = unit(3, "mm"),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        plot.margin        = margin(t = 2, r = 4, b = 2, l = 2))

# each comparison keeps its OWN y-axis (so R62H-vs-CV's real effects stay
# visible - a shared/fixed axis crushed them to nothing against AD-vs-Control's
# ~10x larger range) AND equal panel width (no size differentiation between
# the two - the actual axis tick values already show the scale difference to
# anyone reading them).
p_deg = ggplot(deg_long, aes(cluster, n_plot, fill = direction)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey20") +
  scale_fill_manual(values = fill_colors, labels = fill_labels, name = NULL) +
  scale_y_continuous(labels = function(x) abs(x), name = "No. DEGs") +
  # drop = TRUE (default): each facet only shows the clusters it actually has
  # data for - AD-vs-Control has 13 subclusters, R62H-vs-CV only has 12
  # (no DESeq2 result exists for AST_GFAP_s12's R62H contrast); previously
  # drop = FALSE forced R62H's panel to reserve an empty slot for it too.
  scale_x_discrete(name = NULL) +
  facet_wrap(~ comparison, nrow = 1, scales = "free") +
  deg_theme_compact +
  labs(tag = "a")

### Panel B: pooled log2FC pairwise concordance ------------------------------

DOT = 0.6
# axis labels state the subset each contrast is restricted to: "AD vs Control"
# is fit within CV genotype only; "R62H/R47H vs CV" is fit within AD diagnosis
# only - same "(CV only)"/"(AD only)" convention already used in comp_tags above
# x_pos/y_pos = the group name that a POSITIVE value on that axis means
# ("R62H_vs_CV" positive = higher in R62H); used both in the axis label and to
# build the specific quadrant labels below (e.g. "AD up / R62H down"), instead
# of a generic "discordant" that doesn't say what's actually discordant.
pairs = list(
  list(y = "R62H_vs_CV", x = "CV_AD_vs_Control",
       ylab = "R62H vs CV (log₂FC, AD donors)", xlab = "AD vs Control (log₂FC, CV donors)",
       y_pos = "R62H", x_pos = "AD"),
  list(y = "R47H_vs_CV", x = "CV_AD_vs_Control",
       ylab = "R47H vs CV (log₂FC, AD donors)", xlab = "AD vs Control (log₂FC, CV donors)",
       y_pos = "R47H", x_pos = "AD"),
  list(y = "R62H_vs_CV", x = "R47H_vs_CV",
       ylab = "R62H vs CV (log₂FC, AD donors)", xlab = "R47H vs CV (log₂FC, AD donors)",
       y_pos = "R62H", x_pos = "R47H")
)

reg_levels = c("down_down", "up_up", "down_up", "up_down")
# same up/down palette as the rest of the thesis (Green24 z-scores, DEG panel
# above): blue = down, orange = up. Discordant gets a third, unclaimed
# Okabe-Ito colour (reddish-purple/magenta) rather than base-R "magenta3".
reg_cols   = c(down_down = "#0072B2", up_up = "#E69F00", down_up = "#CC79A7", up_down = "#CC79A7")

norm_res = function(res){
  res = as.data.frame(res)
  g = if ("gene" %in% names(res)) as.character(res$gene) else rownames(res)
  tibble(gene = g, log2FC = res$log2FoldChange, padj = res$padj, pval = res$pvalue)
}
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

scatter_base = function(dot = DOT){
  list(geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2),
       geom_hline(yintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2),
       geom_vline(xintercept = 0, linewidth = 0.3, color = "grey30"),
       geom_hline(yintercept = 0, linewidth = 0.3, color = "grey30"),
       geom_smooth(method = "lm", formula = y ~ x, color = "grey30", linewidth = 0.5, se = FALSE),
       geom_point(aes(color = reg_group), size = dot, alpha = 0.8),
       scale_color_manual(limits = reg_levels, values = reg_cols,
                          breaks = c("down_down", "up_up", "down_up"),
                          labels = c("Down in both", "Up in both", "Discordant"),
                          name = NULL),
       theme_minimal(base_size = 10))
}

e04 = qread(e04c_path)
res_pooled = function(level, contrast) e04$E04_deseq_res[[paste("pooled", level, contrast, sep = "|")]]

# manual tags (not the automatic a/b/c/d sequence): only the first (leftmost)
# panel gets "b"; the middle panel gets no tag at all; the third panel (a
# genuinely different comparison, R62H vs R47H) is "c"
pair_tags = c("b", "", "c")

# build all three pairs' data first, then use ONE shared lim (max across all
# three) for every panel's axes - previously each panel computed its own lim,
# so the three panels weren't on the same scale and weren't visually comparable
pair_data = lapply(pairs, function(p) build_pair(res_pooled("M0_base", p$y), res_pooled("M0_base", p$x)))
lim = max(vapply(pair_data, function(d)
  if (nrow(d) > 0) max(abs(c(d$log2FC_x, d$log2FC_y)), na.rm = TRUE) else 1, numeric(1)))

quad_size = 3.4
quad_pad  = 0.97   # fraction of `lim` - labels sit right at the plot edge;
                   # per-corner hjust/vjust (below) keeps the TEXT growing
                   # inward from that anchor point so it doesn't get clipped

b1 = lapply(seq_along(pairs), function(i){
  p = pairs[[i]]; d = pair_data[[i]]
  qx = lim * quad_pad; qy = lim * quad_pad
  # specific quadrant labels (not a generic "discordant") - name the actual
  # group and direction on each axis with an arrow, e.g. top-right = "AD↑ / R62H↑".
  # hjust/vjust anchor each label to the CORNER closest to it and grow inward
  # (e.g. top-right: hjust=1 grows left, vjust=1 grows down) so long text stays
  # inside the panel instead of overflowing past the plot edge.
  quad_labels = tibble::tibble(
    x = c( qx,  qx, -qx, -qx),
    y = c( qy, -qy,  qy, -qy),
    hjust = c(1, 1, 0, 0),
    vjust = c(1, 0, 1, 0),
    label = c(paste0(p$x_pos, "↑ / ", p$y_pos, "↑"),
             paste0(p$x_pos, "↑ / ", p$y_pos, "↓"),
             paste0(p$x_pos, "↓ / ", p$y_pos, "↑"),
             paste0(p$x_pos, "↓ / ", p$y_pos, "↓")))
  ggplot(d, aes(log2FC_x, log2FC_y)) + scatter_base() +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    geom_text(data = quad_labels, aes(x, y, label = label, hjust = hjust, vjust = vjust),
             inherit.aes = FALSE, size = quad_size, colour = "grey25", fontface = "italic") +
    labs(x = p$xlab, y = p$ylab, tag = pair_tags[i])
})

p_pairwise = wrap_plots(b1, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom") &
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3)))

### combine panels -----------------------------------------------------------

# tags are set manually per-panel (mk_deg_panel = "a" above; b1 panels =
# "b"/""/"c" above), so no tag_levels here - that would try to auto-number
# on top of them
fig = (p_deg / p_pairwise) +
  patchwork::plot_layout(heights = c(0.28, 0.72)) &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(out_dir, paste0(script_ind, "Fig_DEG_pairwise_combined.pdf")),
       fig, width = 15, height = 8, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "Fig_DEG_pairwise_combined.png")),
       fig, width = 15, height = 8, dpi = 300)

message("Done. Combined DEG + pairwise figure written to ", out_dir)
