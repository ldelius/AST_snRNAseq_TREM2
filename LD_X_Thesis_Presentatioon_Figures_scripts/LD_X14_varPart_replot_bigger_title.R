# LD_X14: Replot the LD_H01 variance-partition result with adjusted text sizing.

library(tidyverse)

### paths -------------------------------------------------------------------
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")

h_out   = file.path(base, "LD_H_VarPartition_output")
out_dir = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X14_"

src        = "LD_H01_v02"                     # which variancePartition run to redraw
by_gene    = file.path(h_out, paste0(src, "_varPart_variance_expl_by_gene.csv"))
mean_ord   = file.path(h_out, paste0(src, "_varPart_variance_expl_mean_ord.csv"))
for (f in c(by_gene, mean_ord)) if (!file.exists(f)) stop("Missing input: ", f)
message("Using base: ", base, " | redrawing ", src)

### plot settings ------------------------------------------------------------
# The thesis caption supplies the title.
TITLE      = NULL
TITLE_SIZE = 20
AXIS_SIZE  = 12

### load ---------------------------------------------------------------------
vp  = read_csv(by_gene, show_col_types = FALSE)
ord = read_csv(mean_ord, show_col_types = FALSE)$covariate   # source column order
if (!setequal(ord, names(vp)))
  stop("Covariates in the two CSVs do not match:\n  by_gene: ", paste(names(vp), collapse = ", "),
       "\n  mean_ord: ", paste(ord, collapse = ", "))
vp = vp[, ord]        # preserve source ordering
message("  ", nrow(vp), " genes x ", ncol(vp), " covariates")

### draw ---------------------------------------------------------------------
big_title = theme(plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
                  axis.text.x = element_text(size = AXIS_SIZE),
                  axis.text.y = element_text(size = AXIS_SIZE),
                  axis.title.y = element_text(size = AXIS_SIZE + 1))

if (requireNamespace("variancePartition", quietly = TRUE)) {
  message("  Drawing with variancePartition::plotVarPart()")
  # plotVarPart() sets its own title; NULL removes it rather than leaving the default
  p = variancePartition::plotVarPart(as.data.frame(vp)) + ggtitle(TITLE) + big_title
} else {
  message("  variancePartition not installed here - using the ggplot fallback redraw.")
  long = vp %>%
    pivot_longer(everything(), names_to = "covariate", values_to = "frac") %>%
    mutate(covariate = factor(covariate, levels = ord),
           pct = 100 * frac)
  p = ggplot(long, aes(covariate, pct, fill = covariate)) +
    geom_violin(scale = "width", linewidth = 0.2, colour = "grey30") +
    geom_boxplot(width = 0.07, outlier.size = 0.25, fill = "white", linewidth = 0.2) +
    scale_fill_manual(values = scales::hue_pal()(length(ord)), guide = "none") +
    scale_y_continuous(limits = c(0, 100)) +
    labs(title = TITLE, x = NULL, y = "Variance explained (%)") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    big_title
}

### save ---------------------------------------------------------------------
W = 10; H = 7
ggsave(file.path(out_dir, paste0(script_ind, "varPart_variance_explained.pdf")), p,
       width = W, height = H)
ggsave(file.path(out_dir, paste0(script_ind, "varPart_variance_explained.png")), p,
       width = W, height = H, dpi = 300)

# Maximum-20% zoom.
p20 = p + coord_cartesian(ylim = c(0, 20)) +
  ggtitle(if (is.null(TITLE)) NULL else paste0(TITLE, " (max 20%)"))
ggsave(file.path(out_dir, paste0(script_ind, "varPart_variance_explained_max20.pdf")), p20,
       width = W, height = H)

message("Done. Outputs in: ", out_dir)
