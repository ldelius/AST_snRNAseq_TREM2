# LD_X03: Resolution-selection and integration-QC UMAPs.
# Uses B03 round-2 clusters stored in the B04 object; panel a is assembled separately.

library(tidyverse)
library(qs)
library(Seurat)
library(patchwork)

### paths (all on RDS) ------------------------------------------------------
base     = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
b04_path = file.path(base, "LD_B_AST_analysis_output/LD_B04a_v02_seur.qs")
out_dir  = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X03_"
if (!file.exists(b04_path)) stop("Missing input: ", b04_path)

### what to plot ------------------------------------------------------------
res_to_plot = c("0.2", "0.3", "0.5")   # resolutions (cols SCT_snn_res.<x>)
res_chosen  = "0.3"                     # marked "(selected)" in the title
qc_var      = "cohort"                  # metadata column for the QC panel

### style -------------------------------------------------------------------
PT = 0.15; AL = 0.6
use_rast = requireNamespace("ggrastr", quietly = TRUE)
point_layer = function() {
  if (use_rast) ggrastr::rasterise(geom_point(size = PT, alpha = AL, stroke = 0), dpi = 300)
  else          geom_point(size = PT, alpha = AL, stroke = 0)
}
theme_umap = function() theme_classic(base_size = 12) +
  theme(axis.text   = element_blank(),
        axis.ticks  = element_blank(),
        axis.line   = element_line(linewidth = 0.3),
        plot.title  = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "bottom",
        legend.title    = element_blank(),
        legend.text     = element_text(size = 12),
        legend.key.size = unit(5, "mm"))
make_pal = function(levels) set_names(scales::hue_pal(l = 55, c = 110)(length(levels)), levels)

### load object and retain plotting data ------------------------------------
message("Loading B04 astrocyte object...")
s = qread(b04_path)

emb = Embeddings(s, "umap")[, 1:2]
df  = data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2])

res_cols = paste0("SCT_snn_res.", res_to_plot)
have = res_cols %in% colnames(s@meta.data)
if (!all(have)) stop("Missing clustering columns: ",
                     paste(res_cols[!have], collapse = ", "),
                     "\nAvailable: ", paste(grep("SCT_snn_res", colnames(s@meta.data), value = TRUE), collapse = ", "))
for (rc in res_cols) df[[rc]] = as.character(s@meta.data[[rc]])

if (!qc_var %in% colnames(s@meta.data))
  stop("QC column '", qc_var, "' not found. Available: ",
       paste(colnames(s@meta.data), collapse = ", "))
df[[qc_var]] = as.character(s@meta.data[[qc_var]])
rm(s); gc()
message("  astrocyte nuclei: ", nrow(df))

### panel builders ----------------------------------------------------------
# resolution panel: coloured by cluster, no legend, cluster numbers at centroids
mk_res = function(df, rc, title) {
  cl_lv = as.character(sort(as.integer(unique(df[[rc]]))))
  df[[rc]] = factor(df[[rc]], levels = cl_lv)
  set.seed(1234); d = df[sample(nrow(df)), ]
  cent = d %>% group_by(.cl = .data[[rc]]) %>%
    summarise(x = median(UMAP_1), y = median(UMAP_2), .groups = "drop")
  ggplot(d, aes(UMAP_1, UMAP_2, colour = .data[[rc]])) +
    point_layer() +
    geom_text(data = cent, aes(x, y, label = .cl), inherit.aes = FALSE,
              size = 4.8, fontface = "bold") +
    scale_colour_manual(values = make_pal(cl_lv)) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_umap() + theme(legend.position = "none")
}

# QC panel: coloured by cohort. The larger cohort would otherwise overplot (hide)
# the smaller one, making the panel look single-coloured. So we downsample each
# cohort to equal size, then shuffle, so genuine mixing is visible and fair.
mk_qc = function(df, var, title) {
  lv = sort(unique(df[[var]]))
  df[[var]] = factor(df[[var]], levels = lv)
  set.seed(1234)
  n_min = min(table(df[[var]]))
  idx = unlist(lapply(lv, function(g) {
    ii = which(df[[var]] == g); sample(ii, min(length(ii), n_min))
  }))
  d = df[idx, ]; d = d[sample(nrow(d)), ]
  message("  QC panel: downsampled to ", n_min, " nuclei per '", var, "' group")
  ggplot(d, aes(UMAP_1, UMAP_2, colour = .data[[var]])) +
    point_layer() +
    scale_colour_manual(values = make_pal(lv)) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_umap() +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

### assemble ----------------------------------------------------------------
panels = lapply(seq_along(res_to_plot), function(i) {
  r = res_to_plot[i]
  ttl = paste0("Resolution ", r, if (r == res_chosen) " (selected)" else "")
  mk_res(df, paste0("SCT_snn_res.", r), ttl)
})
qc = mk_qc(df, qc_var, "Integration QC (by cohort)")

# manual panel letters: the three resolution UMAPs are collectively panel "b"
# (tag on the first only); the QC UMAP is panel "c". Panel "a" (dot plot) is
# added separately, so it is not produced here.
panels[[1]] = panels[[1]] + labs(tag = "b")
qc          = qc          + labs(tag = "c")

fig = wrap_plots(c(panels, list(qc)), nrow = 1) &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(out_dir, paste0(script_ind, "Supp_resolution_QC_UMAPs.pdf")),
       fig, width = 20, height = 5.5, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "Supp_resolution_QC_UMAPs.png")),
       fig, width = 20, height = 5.5, dpi = 300)

message("Done. Figure written to ", out_dir)
