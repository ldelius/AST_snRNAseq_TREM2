# LD_X03a: Round-1 marker dot plot showing clusters removed before round 2.
# A cached summary avoids reloading the large Seurat object during plotting.

library(tidyverse)

### paths -------------------------------------------------------------------
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")

in_dir     = file.path(base, "data_TREM2_michael/A_input")
b02_path   = file.path(base, "LD_B_AST_analysis_output/LD_B02a_seur.qs")
marker_csv = file.path(in_dir, "cell_type_markers_241219_w_astr_subtype_markers.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X03a_"
cache_path = file.path(out_dir, paste0(script_ind, "dotplot_data_res1.5.csv"))
message("Using base: ", base)

### what to show ------------------------------------------------------------
res_col  = "SCT_snn_res.1.5"                              # round-1 clustering resolution
# clusters removed in LD_B03a line 35 (expressed non-astrocyte markers)
removed  = c(14, 15, 17, 21, 26, 30, 31, 33, 38, 40)

# marker gene -> cell type, taken from the same source table and in the same order
# the original DotPlot used (level %in% cell_types / neuronal_lineage).
markers = read_csv(marker_csv, show_col_types = FALSE) %>%
  filter(level %in% c("cell_types", "neuronal_lineage")) %>%
  select(gene, cell_type) %>%
  distinct(gene, .keep_all = TRUE)

# short, readable panel labels for the marker groups
group_labs = c(AST_RG = "Astro/RG", AST_oRG = "oRG", AST = "Astro", AST_mat = "Astro mature",
               NEU = "Neuron", ODC_OPC = "OPC", ODC = "Oligo", MIC = "Micro",
               END = "Endo", PER = "Peri", PROL = "Prolif", IPC = "IPC",
               NEU_Cor_exc = "Exc neuron", NEU_Cor_inh = "Inh neuron")

### stage 1: cache build (HPC only, runs once) ------------------------------
if (!file.exists(cache_path)) {
  message("Cache missing - building it from ", basename(b02_path),
          " (heavy step, HPC; runs once).")
  if (!file.exists(b02_path)) stop("Missing input: ", b02_path)
  library(qs); library(Seurat)

  seur = qread(b02_path)
  DefaultAssay(seur) = "SCT"
  if (!res_col %in% names(seur@meta.data))
    stop("Column ", res_col, " not in the object metadata.")

  feats = intersect(markers$gene, rownames(seur))   # same intersect as the original
  # DotPlot is used purely as the summariser: $data holds avg.exp.scaled (z-scored
  # per gene across clusters) and pct.exp, which is all the plot needs.
  dp = DotPlot(seur, features = feats, group.by = res_col, scale.by = "size")$data

  dp %>%
    transmute(gene = as.character(features.plot), cluster = as.character(id),
              avg_exp_scaled = avg.exp.scaled, pct_exp = pct.exp) %>%
    write_csv(cache_path)
  message("Wrote cache: ", cache_path, " (", nrow(dp), " rows)")
  rm(seur); gc()
} else {
  message("Using existing cache: ", cache_path)
}

### stage 2: plot (laptop-friendly) -----------------------------------------
dp = read_csv(cache_path, show_col_types = FALSE) %>%
  left_join(markers, by = "gene") %>%
  mutate(
    gene    = factor(gene, levels = markers$gene[markers$gene %in% gene]),
    group   = factor(recode(cell_type, !!!group_labs),
                     levels = unname(group_labs[unique(markers$cell_type)])),
    # descending level order puts cluster 0 at the top; do NOT reverse the scale
    # afterwards, or the factor positions and the highlight rectangles desynchronise
    cluster = factor(as.integer(cluster), levels = rev(sort(unique(as.integer(cluster))))),
    removed = as.integer(as.character(cluster)) %in% removed
  )

n_clust = nlevels(dp$cluster)

# black outline spanning the full row, for each cluster dropped in round 2
# one rectangle per removed cluster per facet: with zero panel spacing (set in the
# theme) these butt together and read as a single box spanning the whole row
box_df = dp %>% distinct(cluster, group, removed) %>% filter(removed) %>%
  mutate(y = as.integer(cluster))

# Colour: diverging blue -> grey -> orange, centred on 0, because DotPlot's colour
# is a z-score across clusters (negative = below the gene's mean). Same scale and
# hex values as the z-scaled module-score heatmap in LD_X04c, and the same
# thesis-wide down/up convention (#0072B2 low, #E69F00 high).
p = ggplot(dp, aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct_exp, colour = avg_exp_scaled)) +
  geom_rect(data = box_df, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
            fill = NA, colour = "black", linewidth = 0.5) +
  facet_grid(~ group, scales = "free_x", space = "free_x", switch = "x") +
  scale_colour_gradient2(low = "#0072B2", mid = "grey90", high = "#E69F00",
                         midpoint = 0, name = "Scaled mean\nexpression") +
  scale_size_continuous(range = c(0, 4.5), name = "% cells\nexpressing") +
  labs(x = NULL, y = "Cluster (resolution 1.5)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                       size = 8, face = "italic"),
        axis.text.y     = element_text(size = 8),
        panel.grid      = element_line(linewidth = 0.2, colour = "grey92"),
        panel.spacing.x = unit(0, "pt"),   # so the row highlight boxes join up
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text.x    = element_text(size = 7.5, margin = margin(2, 1, 2, 1)),
        strip.placement = "outside",
        legend.key.size = unit(4, "mm"),
        legend.title    = element_text(size = 8),
        legend.text     = element_text(size = 7))

### save --------------------------------------------------------------------
W = 11; H = 7
ggsave(file.path(out_dir, paste0(script_ind, "Supp_Fig2a_marker_dotplot_res1.5.pdf")),
       p, width = W, height = H)
ggsave(file.path(out_dir, paste0(script_ind, "Supp_Fig2a_marker_dotplot_res1.5.png")),
       p, width = W, height = H, dpi = 300)

message("Done. ", n_clust, " clusters, ", nlevels(dp$gene), " genes; ",
        length(removed), " clusters outlined as removed. Outputs in: ", out_dir)
