# LD_X07_v02: Replot pooled log2FC concordance from the E04c checkpoint.
# Opposite-direction categories share one legend entry.

library(tidyverse)
library(qs)
library(patchwork)

base       = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
e_out      = file.path(base, "LD_E_DESeq_pseudobulk")
e04c_path  = file.path(e_out, "LD_E04c_bulk_data.qs")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
script_ind = "LD_X07_v02_"
if (!file.exists(e04c_path)) stop("Missing input: ", e04c_path)

DOT = 0.6

pairs = list(
  list(y = "R62H_vs_CV", x = "CV_AD_vs_Control", ylab = "R62H vs CV", xlab = "AD vs Control"),
  list(y = "R47H_vs_CV", x = "CV_AD_vs_Control", ylab = "R47H vs CV", xlab = "AD vs Control"),
  list(y = "R62H_vs_CV", x = "R47H_vs_CV",       ylab = "R62H vs CV", xlab = "R47H vs CV")
)

reg_levels = c("down_down", "up_up", "down_up", "up_down")   # the only 4 combos that ever occur
reg_cols   = c(down_down = "blue", up_up = "red", down_up = "magenta3", up_down = "magenta3")

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

b1 = lapply(pairs, function(p){
  d = build_pair(res_pooled("M0_base", p$y), res_pooled("M0_base", p$x))
  lim = if (nrow(d) > 0) max(abs(c(d$log2FC_x, d$log2FC_y)), na.rm = TRUE) else 1
  ggplot(d, aes(log2FC_x, log2FC_y)) + scatter_base() +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, size = 3,
             label = sprintf("r = %.2f, n = %d", safe_cor(d$log2FC_x, d$log2FC_y), nrow(d))) +
    labs(x = paste0("log2FC (", p$xlab, ")"), y = paste0("log2FC (", p$ylab, ")"),
         title = paste0(p$ylab, " vs ", p$xlab))
})

p_main = wrap_plots(b1, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom") &
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3)))

ggsave(file.path(out_dir, paste0(script_ind, "log2FC_concordance_pooled_main.pdf")),
       p_main, width = 14, height = 5, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "log2FC_concordance_pooled_main.png")),
       p_main, width = 14, height = 5, dpi = 300)

message("Done. Written to ", out_dir)
