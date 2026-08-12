# LD_X02: Overview UMAPs for the whole dataset and astrocyte families/subclusters.
# The whole-dataset panel uses the B01 LIGER embedding; astrocyte panels use B04
# Harmony. Large objects are loaded sequentially to limit peak memory.

library(tidyverse)
library(qs)
library(Seurat)
library(patchwork)

### paths (all on RDS) ------------------------------------------------------
base     = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
b01_path = file.path(base, "data_TREM2_michael/B_load_from_scflow_subcluster/B01_seur.qs")
b04_path = file.path(base, "LD_B_AST_analysis_output/LD_B04a_v02_seur.qs")
clust_csv= file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
out_dir  = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X02_v03_"
for (p in c(b01_path, b04_path, clust_csv))
  if (!file.exists(p)) stop("Missing input: ", p)
message("base: ", base, "\nout:  ", out_dir)

### shared plotting style ----------------------------------------------------
PT  = 0.20   # point size (same in all panels)
AL  = 0.80   # point alpha (raised from 0.55: at low alpha, heavy overplotting on
             # white left dense cluster cores looking washed out / desaturated)
use_rast  = requireNamespace("ggrastr", quietly = TRUE)   # keeps PDF small if available
use_repel = requireNamespace("ggrepel", quietly = TRUE)   # nicer label placement if available

point_layer = function() {
  if (use_rast) ggrastr::rasterise(geom_point(size = PT, alpha = AL, stroke = 0), dpi = 300)
  else          geom_point(size = PT, alpha = AL, stroke = 0)
}

# place each category's name at the centre (median position) of its points,
# so the cluster is labelled directly on the UMAP instead of (only) in a legend.
label_layer = function(df, col, emphasise = NULL, emph_size = 5.2) {
  cen = df %>%
    group_by(.lab = .data[[col]]) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop") %>%
    mutate(.size = if (is.null(emphasise)) 3.3 else if_else(.lab %in% emphasise, emph_size, 3.3))
  geom = if (use_repel)
    ggrepel::geom_text_repel(data = cen, aes(UMAP_1, UMAP_2, label = .lab, size = .size),
                             inherit.aes = FALSE, fontface = "bold",
                             colour = "black", bg.color = "white", bg.r = 0.12,
                             min.segment.length = 0, segment.size = 0.2, max.overlaps = Inf,
                             show.legend = FALSE)
  else
    geom_text(data = cen, aes(UMAP_1, UMAP_2, label = .lab, size = .size),
              inherit.aes = FALSE, fontface = "bold", colour = "black",
              show.legend = FALSE)
  list(geom, scale_size_identity())
}

theme_umap = function() theme_classic(base_size = 11) +
  theme(axis.text   = element_blank(),
        axis.ticks  = element_blank(),
        axis.line   = element_line(linewidth = 0.3),
        plot.title  = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "bottom",
        legend.title    = element_blank(),
        legend.text     = element_text(size = 7),
        legend.key.size = unit(3, "mm"))

# Use Okabe-Ito for up to eight categories and viridis for larger sets.
okabe_ito = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
              "#0072B2", "#D55E00", "#CC79A7", "#000000")
make_pal = function(levels, cb = FALSE) {
  n = length(levels)
  cols = if (n <= length(okabe_ito)) okabe_ito[seq_len(n)]
         else                        scales::viridis_pal(option = "viridis")(n)
  set_names(cols, levels)
}

# fixed (name-keyed, not position-keyed) palette for the three astrocyte families,
# so the same three hex values can be reused verbatim in any other figure that
# colours by family, keeping the palette consistent across the thesis.
astro_family_pal = c(AST_SLC1A2 = "#E69F00", AST_GFAP = "#0072B2", AST_CHI3L1 = "#009E73")

mk_umap = function(df, col, title, leg_ncol = 4, label = FALSE, show_legend = TRUE, cb = FALSE, emphasise = NULL, pal_override = NULL) {
  lv = levels(df[[col]]); if (is.null(lv)) lv = sort(unique(df[[col]]))
  cols = if (!is.null(pal_override)) pal_override[lv] else make_pal(lv, cb = cb)
  p = ggplot(df, aes(UMAP_1, UMAP_2, colour = .data[[col]])) +
    point_layer() +
    scale_colour_manual(values = cols, drop = FALSE) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_umap()
  if (label) p = p + label_layer(df, col, emphasise = emphasise)
  if (show_legend)
    p = p + guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = leg_ncol))
  else
    p = p + theme(legend.position = "none")
  p
}

### helper: pull a 2-D embedding aligned to metadata rows -------------------
get_emb = function(seur, name) {
  red = seur@reductions[[name]]
  if (is.null(red)) stop("Reduction '", name, "' not found. Have: ",
                         paste(names(seur@reductions), collapse = ", "))
  m = if (inherits(red, "DimReduc")) Embeddings(red) else as.matrix(as.data.frame(red))
  cells = rownames(seur@meta.data)
  if (!is.null(rownames(m)) && all(cells %in% rownames(m))) m = m[cells, , drop = FALSE]
  data.frame(UMAP_1 = as.numeric(m[, 1]), UMAP_2 = as.numeric(m[, 2]))
}

### (a) whole dataset, cell type --------------------------------------------
message("Loading B01 (whole dataset)...")
s = qread(b01_path)
df_all = cbind(get_emb(s, "UMAP_Liger"),
               cell_type = as.character(s$cluster_celltype))
df_all$cell_type = dplyr::recode(df_all$cell_type,
                                  "Astro" = "Astrocytes",
                                  "Micro" = "Microglia")  # spell out for figure; rest kept as abbreviations
df_all$cell_type = factor(df_all$cell_type, levels = sort(unique(df_all$cell_type)))
rm(s); gc()
message("  panel (a) total nuclei: ", nrow(df_all))

### (b,c) astrocytes, cell_type & cluster_name ------------------------------
message("Loading B04 (astrocytes)...")
ord = read_csv(clust_csv, show_col_types = FALSE)   # gives a sensible legend order
s = qread(b04_path)
df_ast = cbind(get_emb(s, "umap"),
               cell_type    = as.character(s$cell_type),
               cluster_name = as.character(s$cluster_name))
rm(s); gc()
# order factors so the legend is grouped by astrocyte subtype family
df_ast$cell_type    = factor(df_ast$cell_type,    levels = unique(ord$cell_type))
df_ast$cluster_name = factor(df_ast$cluster_name, levels = unique(ord$cluster_name))
message("  panels (b,c) astrocyte nuclei: ", nrow(df_ast))
message("  astrocyte subclusters (cluster_name): ", length(unique(as.character(df_ast$cluster_name))),
        "  | subtypes (cell_type): ", length(unique(as.character(df_ast$cell_type))))

### assemble one-row figure --------------------------------------------------
pa = mk_umap(df_all, "cell_type",    "All cell types",        label = TRUE, show_legend = FALSE,
             emphasise = "Astrocytes")
pb = mk_umap(df_ast, "cell_type",    "Astrocyte families",   label = TRUE, show_legend = FALSE, pal_override = astro_family_pal)
pc = mk_umap(df_ast, "cluster_name", "Astrocyte subclusters", label = TRUE, show_legend = FALSE)

fig = (pa | pb | pc) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(out_dir, paste0(script_ind, "Fig_overview_UMAPs.pdf")),
       fig, width = 16, height = 6, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "Fig_overview_UMAPs.png")),
       fig, width = 16, height = 6, dpi = 300)

message("Done. Figure written to ", out_dir)
