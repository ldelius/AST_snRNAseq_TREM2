message("\n\n##########################################################################\n",
        "# Start LD_X09: Hallmark interferon alpha/gamma dotplots ", Sys.time(),
        "\n##########################################################################\n\n")

# R47H_vs_R62H uses R62H as reference, so positive log2FC indicates higher
# expression in R47H. Colour and size scales are shared within each scope but
# not across scopes because per-subcluster values have a wider range.
library(tidyverse)
library(qs)
library(msigdbr)


### directories / index ------------------------------------------------------

base       = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
setwd(base)

script_ind = "LD_X09_"
e_out      = file.path(base, "LD_E_DESeq_pseudobulk")
e02c_path  = file.path(e_out, "LD_E02c/LD_E02c_v01_bulk_data.qs")
e04c_path  = file.path(e_out, "LD_E04c_bulk_data.qs")
clust_csv  = file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


### config --------------------------------------------------------------------

clust_tab     = read_csv(clust_csv, show_col_types = FALSE)
cluster_order = unique(clust_tab$cluster_name)   # SLC1A2 -> GFAP -> CHI3L1, matches X07/X08

# per-subcluster comparison tags (key = "<cluster>_<tag>"), matching X07/X08
e02c_tag = c(R62H_vs_CV   = "AD_TREM2_R62H_vs_CV",
             R47H_vs_CV   = "AD_TREM2_R47H_vs_CV",
             R47H_vs_R62H = "AD_TREM2_R47H_vs_R62H")

contrasts = names(e02c_tag)   # same 3 names used for the pooled scope

SIG_CUT = 0.05
SIZE_CAP = 10


### load DESeq2 results (both scopes) -----------------------------------------

e02  = qread(e02c_path)
dres = e02$deseq_results
rm(e02); gc()
message("   loaded ", length(dres), " E02c per-subcluster DESeq result tables")

e04           = qread(e04c_path)
e04_deseq_res = e04$E04_deseq_res
rm(e04); gc()
res_pooled = function(contrast) e04_deseq_res[[paste("pooled", "M0_base", contrast, sep = "|")]]


### Hallmark interferon gene sets ---------------------------------------------

hallmark = msigdbr(species = "Homo sapiens", collection = "H")
alpha_set = unique(hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE"])
gamma_set = unique(hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE"])

set_alpha      = alpha_set
set_all        = union(alpha_set, gamma_set)
set_gamma_only = setdiff(gamma_set, alpha_set)

message("   HALLMARK_INTERFERON_ALPHA_RESPONSE: ", length(alpha_set), " genes")
message("   HALLMARK_INTERFERON_GAMMA_RESPONSE: ", length(gamma_set), " genes")
message("   intersect(alpha, gamma): ", length(intersect(alpha_set, gamma_set)), " genes")
message("   set_alpha: ", length(set_alpha), " | set_all: ", length(set_all),
        " | set_gamma_only: ", length(set_gamma_only))

gene_sets = list(alpha = set_alpha, all = set_all, gamma_only = set_gamma_only)


### restrict gene sets to genes actually tested, per scope --------------------

detected_subcluster = rownames(dres[[1]])
detected_pooled      = rownames(res_pooled("R62H_vs_CV"))

restrict_sets = function(sets, detected, scope_label){
  lapply(names(sets), function(nm){
    g = sets[[nm]]
    g_keep = intersect(g, detected)
    message("   [", scope_label, "] set '", nm, "': ", length(g_keep), " / ", length(g),
            " genes detected (", length(g) - length(g_keep), " dropped, not tested)")
    g_keep
  }) %>% setNames(names(sets))
}

sets_subcluster = restrict_sets(gene_sets, detected_subcluster, "subcluster")
sets_pooled     = restrict_sets(gene_sets, detected_pooled, "pooled")


### build long tables ----------------------------------------------------------

# one row per (cluster, contrast, gene) for the per-subcluster scope
build_long_subcluster = function(genes){
  out = list()
  for (cl in cluster_order){
    for (cn in contrasts){
      r = dres[[paste0(cl, "_", e02c_tag[cn])]]
      if (is.null(r)) next
      r = r[r$gene %in% genes, ]
      if (nrow(r) == 0) next
      out[[length(out) + 1]] = tibble(
        cluster = cl, contrast = cn, gene = r$gene,
        log2FoldChange = r$log2FoldChange, pvalue = r$pvalue, padj = r$padj)
    }
  }
  bind_rows(out)
}

# one row per (contrast, gene) for the pooled scope
build_long_pooled = function(genes){
  out = list()
  for (cn in contrasts){
    r = res_pooled(cn)
    if (is.null(r)) next
    r = as.data.frame(r)
    g = if ("gene" %in% names(r)) as.character(r$gene) else rownames(r)
    r = r[g %in% genes, ]
    if (nrow(r) == 0) next
    out[[length(out) + 1]] = tibble(
      contrast = cn, gene = if ("gene" %in% names(r)) r$gene else rownames(r),
      log2FoldChange = r$log2FoldChange, pvalue = r$pvalue, padj = r$padj)
  }
  bind_rows(out)
}

finalize_long = function(d){
  d = d[!is.na(d$log2FoldChange) & !is.na(d$padj), ]
  d$negLogP = pmin(-log10(d$padj), SIZE_CAP)
  d$sig     = d$padj < SIG_CUT
  d
}

long_subcluster = lapply(sets_subcluster, build_long_subcluster) %>% lapply(finalize_long)
long_pooled     = lapply(sets_pooled,     build_long_pooled)     %>% lapply(finalize_long)


### gene ordering (fixed once per scope, from R62H_vs_CV over set_all) --------

# per-subcluster: median log2FoldChange across clusters
ord_sub_tab = long_subcluster$all[long_subcluster$all$contrast == "R62H_vs_CV", ]
ord_sub = ord_sub_tab %>% group_by(gene) %>% summarise(m = median(log2FoldChange, na.rm = TRUE)) %>%
  arrange(desc(m))
gene_order_subcluster = ord_sub$gene

# pooled: single log2FoldChange value per gene
ord_pool_tab = long_pooled$all[long_pooled$all$contrast == "R62H_vs_CV", ]
ord_pool = ord_pool_tab %>% arrange(desc(log2FoldChange))
gene_order_pooled = ord_pool$gene

apply_gene_factor = function(d, gene_order){
  d$gene = factor(d$gene, levels = rev(gene_order[gene_order %in% unique(d$gene)]))
  d
}
long_subcluster = lapply(long_subcluster, apply_gene_factor, gene_order = gene_order_subcluster)
long_pooled     = lapply(long_pooled,     apply_gene_factor, gene_order = gene_order_pooled)


### shared colour/size scale, computed once per scope (over set_all) ---------

scale_from = function(d_all){
  lim       = max(abs(d_all$log2FoldChange), na.rm = TRUE)
  size_rng  = range(d_all$negLogP, na.rm = TRUE)
  list(lim = lim, size_rng = size_rng)
}
scale_subcluster = scale_from(long_subcluster$all)
scale_pooled     = scale_from(long_pooled$all)

message("   [subcluster] colour limit +/-", round(scale_subcluster$lim, 2),
        " | size range ", paste(round(scale_subcluster$size_rng, 2), collapse = "-"))
message("   [pooled] colour limit +/-", round(scale_pooled$lim, 2),
        " | size range ", paste(round(scale_pooled$size_rng, 2), collapse = "-"))


### plotting -------------------------------------------------------------------

contrast_lab = c(R62H_vs_CV = "R62H vs CV", R47H_vs_CV = "R47H vs CV", R47H_vs_R62H = "R47H vs R62H")

make_plot = function(d, scale_info, x_var, facet, title){
  sig_d = d[d$sig, ]
  p = ggplot(d, aes(x = .data[[x_var]], y = gene, colour = log2FoldChange, size = negLogP)) +
    geom_point() +
    { if (nrow(sig_d) > 0) geom_point(data = sig_d, shape = 1, colour = "black", stroke = 0.6) } +
    scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0,
                           limits = c(-scale_info$lim, scale_info$lim), name = "log2FC") +
    scale_size_continuous(limits = scale_info$size_rng, range = c(0.5, 5),
                          name = paste0("-log10(padj)\n(capped at ", SIZE_CAP, ")")) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = NULL, y = NULL)
  if (facet) p = p + facet_wrap(~ contrast, nrow = 1,
                                labeller = labeller(contrast = contrast_lab))
  p
}

set_title = c(alpha = "HALLMARK_INTERFERON_ALPHA_RESPONSE", all = "Alpha + Gamma union",
             gamma_only = "Gamma-only (not in alpha set)")

plots_subcluster = list()
for (nm in names(long_subcluster)){
  d = long_subcluster[[nm]]
  plots_subcluster[[nm]] = make_plot(d, scale_subcluster, x_var = "cluster", facet = TRUE,
                                     title = paste0(set_title[[nm]], " - per subcluster"))
}

plots_pooled = list()
for (nm in names(long_pooled)){
  d = long_pooled[[nm]]
  d$contrast = factor(contrast_lab[as.character(d$contrast)], levels = unname(contrast_lab))
  plots_pooled[[nm]] = make_plot(d, scale_pooled, x_var = "contrast", facet = FALSE,
                                 title = paste0(set_title[[nm]], " - pooled"))
}

# dynamic page height from gene count (per-scope, since set sizes differ);
# width fixed by x-axis cardinality (many clusters x3 facets vs 3 contrasts)
height_for = function(d) max(6, 0.16 * length(unique(d$gene)) + 2)

pdf_names = c(alpha = "IFN_alpha_dotplot", all = "IFN_all_dotplot", gamma_only = "IFN_gamma_only_dotplot")

for (nm in names(plots_subcluster)){
  h = height_for(long_subcluster[[nm]])
  pdf(file.path(out_dir, paste0(script_ind, pdf_names[[nm]], "_subcluster.pdf")), width = 12, height = h)
  print(plots_subcluster[[nm]])
  dev.off()
}
for (nm in names(plots_pooled)){
  h = height_for(long_pooled[[nm]])
  pdf(file.path(out_dir, paste0(script_ind, pdf_names[[nm]], "_pooled.pdf")), width = 6, height = h)
  print(plots_pooled[[nm]])
  dev.off()
}
message("    Saved 6 dotplot PDFs (3 gene sets x 2 scopes)")


### combined long CSV + qsave of plot objects ----------------------------------

csv_tab = bind_rows(
  bind_rows(lapply(names(long_subcluster), function(nm)
    long_subcluster[[nm]] %>% mutate(set = nm, scope = "subcluster", gene = as.character(gene)))),
  bind_rows(lapply(names(long_pooled), function(nm)
    long_pooled[[nm]] %>% mutate(set = nm, scope = "pooled", gene = as.character(gene), cluster = NA_character_)))
)
write_csv(csv_tab, file.path(out_dir, paste0(script_ind, "IFN_dotplot_data.csv")))
message("    Saved combined long CSV: ", nrow(csv_tab), " rows")

all_plots = c(setNames(plots_subcluster, paste0(names(plots_subcluster), "_subcluster")),
             setNames(plots_pooled, paste0(names(plots_pooled), "_pooled")))
qsave(all_plots, file.path(out_dir, paste0(script_ind, "IFN_dotplots.qs")))
message("    Saved qs bundle of all 6 plot objects")


sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed LD_X09 ", Sys.time(),
        "\n##########################################################################\n\n\n")
