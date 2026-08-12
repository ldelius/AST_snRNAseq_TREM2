# LD_X12: Supplementary marker plots and clustering-resolution UMAPs.
# Cached summaries avoid reloading the large Seurat objects during plotting.

library(tidyverse)
library(patchwork)

### paths -------------------------------------------------------------------
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")

in_dir     = file.path(base, "data_TREM2_michael/A_input")
b_out      = file.path(base, "LD_B_AST_analysis_output")
b02_path   = file.path(b_out, "LD_B02a_seur.qs")
b04_path   = file.path(b_out, "LD_B04a_v02_seur.qs")
clust_csv  = file.path(b_out, "LD_B03a_cluster_assignment.csv")   # subcluster labels/order
marker_csv = file.path(in_dir, "cell_type_markers_241219_w_astr_subtype_markers.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X12_"

# panel a cache keeps the LD_X03a_ name: that script already builds it, so a run
# submitted from there is reused here rather than repeated.
cache_a = file.path(out_dir, "LD_X03a_dotplot_data_res1.5.csv")
cache_b = file.path(out_dir, paste0(script_ind, "umap_resolutions.rds"))
cache_c = file.path(out_dir, paste0(script_ind, "dotplot_subtype_markers.csv"))
message("Using base: ", base)

### what to show ------------------------------------------------------------
res_col_a   = "SCT_snn_res.1.5"                            # panel a: round-1 resolution
removed     = c(14, 15, 17, 21, 26, 30, 31, 33, 38, 40)    # dropped in LD_B03a:35
res_to_plot = c("0.2", "0.3", "0.5")                       # panel b: resolutions shown
res_chosen  = "0.3"                                        # marked "(selected)"
key_genes   = c("GFAP", "SLC1A2", "CHI3L1")                # panel c: bold (family names)

mk_tab = read_csv(marker_csv, show_col_types = FALSE)

# panel a genes: broad cell-type markers, grouped by their cell_type column
markers_a = mk_tab %>% filter(level %in% c("cell_types", "neuronal_lineage")) %>%
  select(gene, group_raw = cell_type) %>% distinct(gene, .keep_all = TRUE)
group_labs_a = c(AST_RG = "Astro/RG", AST_oRG = "oRG", AST = "Astro", AST_mat = "Astro mature",
                 NEU = "Neuron", ODC_OPC = "OPC", ODC = "Oligo", MIC = "Micro",
                 END = "Endo", PER = "Peri", PROL = "Prolif", IPC = "IPC",
                 NEU_Cor_exc = "Exc neuron", NEU_Cor_inh = "Inh neuron")

# panel c genes: astrocyte-subtype markers, grouped by their subtype column
markers_c = mk_tab %>% filter(level == "Astrocyte_subtypes") %>%
  select(gene, group_raw = subtype) %>% distinct(gene, .keep_all = TRUE)
group_labs_c = c(A1 = "A1", A2 = "A2", pan_react = "Pan-reactive", DAA = "DAA",
                 Ast1_2_homeostatic = "Ast1/2 homeostatic", Ast3_enh_mitophagy = "Ast3 mitophagy",
                 Ast4_react = "Ast4 reactive", Ast5_react = "Ast5 reactive", Ast6 = "Ast6",
                 Ast7_IFN_resp = "Ast7 IFN", Ast8_stress_resp = "Ast8 stress",
                 Ast9_stress_resp = "Ast9 stress", Ast10_stress_resp = "Ast10 stress")

### stage 1: cache builds (HPC only, each runs once) ------------------------
if (!file.exists(cache_a)) {
  message("Panel a cache missing - building from ", basename(b02_path), " (heavy, HPC).")
  if (!file.exists(b02_path)) stop("Missing input: ", b02_path)
  library(qs); library(Seurat)
  seur = qread(b02_path); DefaultAssay(seur) = "SCT"
  if (!res_col_a %in% names(seur@meta.data)) stop("Column ", res_col_a, " not found.")
  # DotPlot used purely as the summariser: $data has avg.exp.scaled and pct.exp
  DotPlot(seur, features = intersect(markers_a$gene, rownames(seur)),
          group.by = res_col_a, scale.by = "size")$data %>%
    transmute(gene = as.character(features.plot), cluster = as.character(id),
              avg_exp_scaled = avg.exp.scaled, pct_exp = pct.exp) %>%
    write_csv(cache_a)
  message("  wrote ", basename(cache_a)); rm(seur); gc()
} else message("Panel a: using cache ", basename(cache_a))

if (!file.exists(cache_b) || !file.exists(cache_c)) {
  message("Panel b/c cache missing - building from ", basename(b04_path),
          " (heavy, HPC; one load serves both).")
  if (!file.exists(b04_path)) stop("Missing input: ", b04_path)
  library(qs); library(Seurat)
  s = qread(b04_path)

  if (!file.exists(cache_b)) {                       # panel b: UMAP + resolutions
    emb = Embeddings(s, "umap")[, 1:2]
    d = data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2])
    res_cols = paste0("SCT_snn_res.", res_to_plot)
    miss = res_cols[!res_cols %in% colnames(s@meta.data)]
    if (length(miss)) stop("Missing clustering columns: ", paste(miss, collapse = ", "),
                           "\nAvailable: ", paste(grep("SCT_snn_res", colnames(s@meta.data),
                                                       value = TRUE), collapse = ", "))
    for (rc in res_cols) d[[rc]] = as.character(s@meta.data[[rc]])
    saveRDS(d, cache_b)   # RDS: compact, and readable with base R on the laptop
    message("  wrote ", basename(cache_b), " (", nrow(d), " nuclei)")
  }

  if (!file.exists(cache_c)) {                       # panel c: subtype-marker dot plot
    DefaultAssay(s) = "SCT"
    if (!"cluster_name" %in% colnames(s@meta.data))
      stop("Column 'cluster_name' not found in the B04 object.")
    DotPlot(s, features = intersect(markers_c$gene, rownames(s)),
            group.by = "cluster_name", scale.by = "size")$data %>%
      transmute(gene = as.character(features.plot), cluster_name = as.character(id),
                avg_exp_scaled = avg.exp.scaled, pct_exp = pct.exp) %>%
      write_csv(cache_c)
    message("  wrote ", basename(cache_c))
  }
  rm(s); gc()
} else message("Panel b/c: using caches ", basename(cache_b), " + ", basename(cache_c))

### shared style ------------------------------------------------------------
okabe_ito = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
              "#0072B2", "#D55E00", "#CC79A7", "#000000")
# Cluster identity in panel b uses turbo, not the thesis default viridis. Viridis is
# perceptually UNIFORM - excellent for continuous quantities, but for ~20 discrete
# clusters it reads as one purple-to-green wash and neighbouring clusters are hard to
# tell apart. Turbo spans a much wider hue range, so each cluster is distinct. Losing
# viridis' colourblind-safety is acceptable here specifically because colour is not
# the only channel: every cluster is also labelled with its number at its centroid.
# begin/end trim turbo's very dark navy and dark red ends, which would otherwise
# render the largest cluster near-black and swallow its centroid label.
make_pal = function(levels, option = "viridis", begin = 0, end = 1) {
  n = length(levels)
  cols = if (option == "okabe" && n <= length(okabe_ito)) okabe_ito[seq_len(n)]
         else scales::viridis_pal(option = option, begin = begin, end = end)(n)
  set_names(cols, levels)
}

theme_umap = function() theme_classic(base_size = 11) +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line  = element_line(linewidth = 0.3),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "none")

# shared dot-plot scales, so panels a and c encode expression identically
dot_scales = function(size_range)
  list(scale_colour_gradient2(low = "#0072B2", mid = "grey90", high = "#E69F00",
                              midpoint = 0, name = "Scaled mean\nexpression"),
       scale_size_continuous(range = size_range, name = "% cells\nexpressing"))

theme_dot = function(base = 11) theme_bw(base_size = base) +
  theme(axis.text.y      = element_text(size = 8),
        panel.grid       = element_line(linewidth = 0.2, colour = "grey92"),
        legend.key.size  = unit(4, "mm"),
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 7))

PT = 0.15; AL = 0.6
use_rast = requireNamespace("ggrastr", quietly = TRUE)
point_layer = function() {
  if (use_rast) ggrastr::rasterise(geom_point(size = PT, alpha = AL, stroke = 0), dpi = 300)
  else          geom_point(size = PT, alpha = AL, stroke = 0)
}

### panel a: marker dot plot ------------------------------------------------
dp_a = read_csv(cache_a, show_col_types = FALSE) %>%
  left_join(markers_a, by = "gene") %>%
  mutate(
    gene  = factor(gene, levels = markers_a$gene[markers_a$gene %in% gene]),
    group = factor(recode(group_raw, !!!group_labs_a),
                   levels = unname(group_labs_a[unique(markers_a$group_raw)])),
    # descending levels put cluster 0 on top; do NOT also reverse the scale, or
    # the factor positions and the highlight rectangles desynchronise
    cluster = factor(as.integer(cluster), levels = rev(sort(unique(as.integer(cluster))))),
    removed = as.integer(as.character(cluster)) %in% removed)

# Removed clusters are identified by bold y-axis labels.
faces_a = ifelse(as.integer(levels(dp_a$cluster)) %in% removed, "bold", "plain")

pa = ggplot(dp_a, aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct_exp, colour = avg_exp_scaled)) +
  facet_grid(~ group, scales = "free_x", space = "free_x", switch = "x") +
  dot_scales(c(0, 4.5)) +
  labs(x = NULL, y = "Cluster (resolution 1.5)") +
  theme_dot() +
  theme(axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                        size = 8, face = "italic"),
        axis.text.y      = element_text(size = 8, face = faces_a),
        panel.spacing.x  = unit(0, "pt"),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text.x     = element_text(size = 7.5, margin = margin(2, 1, 2, 1)),
        strip.placement  = "outside")

### panel b: resolution UMAPs ------------------------------------------------
df_umap = readRDS(cache_b)

mk_res = function(df, rc, title) {
  cl_lv = as.character(sort(as.integer(unique(df[[rc]]))))
  df[[rc]] = factor(df[[rc]], levels = cl_lv)
  set.seed(1234); d = df[sample(nrow(df)), ]
  cent = d %>% group_by(.cl = .data[[rc]]) %>%
    summarise(x = median(UMAP_1), y = median(UMAP_2), .groups = "drop")
  ggplot(d, aes(UMAP_1, UMAP_2, colour = .data[[rc]])) +
    point_layer() +
    # labels on a translucent white plate so they stay readable over any cluster
    # colour (plain black text disappears on the darker ones)
    geom_label(data = cent, aes(x, y, label = .cl), inherit.aes = FALSE,
               size = 3, fontface = "bold", colour = "black",
               fill = alpha("white", 0.65), label.size = 0,
               label.padding = unit(0.6, "mm")) +
    scale_colour_manual(values = make_pal(cl_lv, option = "turbo",
                                          begin = 0.12, end = 0.93)) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_umap()
}

panels_b = lapply(res_to_plot, function(r)
  mk_res(df_umap, paste0("SCT_snn_res.", r),
         paste0("Resolution ", r, if (r == res_chosen) " (selected)" else "")))

### panel c: subtype-marker dot plot -----------------------------------------
# subcluster order (SLC1A2 -> GFAP -> CHI3L1 families) comes from the manual
# assignment table, so the y axis matches every other subcluster figure
clust_ord = read_csv(clust_csv, show_col_types = FALSE)$cluster_name

dp_c = read_csv(cache_c, show_col_types = FALSE) %>%
  left_join(markers_c, by = "gene") %>%
  mutate(gene         = factor(gene, levels = markers_c$gene[markers_c$gene %in% gene]),
         group        = factor(recode(group_raw, !!!group_labs_c),
                               levels = unname(group_labs_c[unique(markers_c$group_raw)])),
         cluster_name = factor(cluster_name, levels = rev(clust_ord)))

# No subtype-group annotation here (no strips, bands or separators): with only 18
# rows the dots have room, and the grouping is already conveyed by the gene order.
# Genes stay ordered by subtype so related markers remain adjacent.
gene_lv = levels(dp_c$gene)
n_cl    = nlevels(dp_c$cluster_name)
faces   = ifelse(gene_lv %in% key_genes, "bold", "italic")   # the three family genes

pc = ggplot(dp_c, aes(x = gene, y = cluster_name)) +
  geom_point(aes(size = pct_exp, colour = avg_exp_scaled)) +
  dot_scales(c(0, 3.5)) +
  labs(x = NULL, y = "Astrocyte subcluster") +
  theme_dot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   size = 5.5, face = faces))

### assemble -----------------------------------------------------------------
pa            = pa + labs(tag = "a")
panels_b[[1]] = panels_b[[1]] + labs(tag = "b")
pb            = wrap_plots(panels_b, nrow = 1)
pc_tag        = pc + labs(tag = "c")

# panel a gets the most height: 43 cluster rows, which otherwise pack tight enough
# that neighbouring dots overlap. c needs least (18 rows).
fig = wrap_plots(list(pa, pb, pc_tag), ncol = 1, heights = c(2.6, 1.0, 1.05)) &
  theme(plot.tag = element_text(face = "bold", size = 14))

# c has 111 genes on the x axis (the source plot was 25 in wide), so the whole
# figure is landscape/wide; a and b are stretched to match.
W = 19; H = 21
ggsave(file.path(out_dir, paste0(script_ind, "Supp_Fig2.pdf")), fig,
       width = W, height = H, useDingbats = FALSE, limitsize = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "Supp_Fig2.png")), fig,
       width = W, height = H, dpi = 300, limitsize = FALSE)

# individual panels too, for flexible placement in the thesis
ggsave(file.path(out_dir, paste0(script_ind, "Supp_Fig2a_marker_dotplot.pdf")),
       pa, width = 11, height = 10)
ggsave(file.path(out_dir, paste0(script_ind, "Supp_Fig2b_resolution_UMAPs.pdf")),
       pb, width = 12, height = 4.2)
ggsave(file.path(out_dir, paste0(script_ind, "Supp_Fig2c_subtype_marker_dotplot.pdf")),
       pc, width = 19, height = 5.2, limitsize = FALSE)

message("Done. a: ", nlevels(dp_a$cluster), " clusters x ", nlevels(dp_a$gene), " genes (",
        length(removed), " outlined) | b: ", nrow(df_umap), " nuclei, res ",
        paste(res_to_plot, collapse = "/"), " | c: ", n_cl, " subclusters x ",
        nlevels(dp_c$gene), " genes. Outputs in: ", out_dir)
