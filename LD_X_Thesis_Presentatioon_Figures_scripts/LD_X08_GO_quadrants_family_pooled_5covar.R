message("\n\n##########################################################################\n",
        "# Start LD_X08: Pooled astrocyte-family GO analysis ", Sys.time(),
        "\n##########################################################################\n\n")

library(tidyverse)
library(qs)
library(clusterProfiler)
library(org.Hs.eg.db)


### directories / index

base       = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
setwd(base)

script_ind = "LD_X08_"
e_out      = file.path(base, "LD_E_DESeq_pseudobulk")
e02c_path  = file.path(e_out, "LD_E02c/LD_E02c_v01_bulk_data.qs")
e04c_path  = file.path(e_out, "LD_E04c_bulk_data.qs")
deg_csv    = file.path(e_out, "LD_E02c/LD_E02c_v01_DEGs_by_cluster_genes.csv")
clust_csv  = file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


### config

clust_tab = read_csv(clust_csv, show_col_types = FALSE)
families  = c(SLC1A2 = "AST_SLC1A2", GFAP = "AST_GFAP", CHI3L1 = "AST_CHI3L1")
family_clusters = lapply(families, function(ct) clust_tab$cluster_name[clust_tab$cell_type == ct])
message("   family subcluster counts: ",
        paste(names(family_clusters), lengths(family_clusters), sep = " = ", collapse = ", "))

# E02c per-subcluster comparison tags (key = "<cluster>_<tag>"), matching LD_X07
e02c_tag = c(CV_AD_vs_Control = "TREM2_CV_AD_vs_Control",
             R62H_vs_CV       = "AD_TREM2_R62H_vs_CV",
             R47H_vs_CV       = "AD_TREM2_R47H_vs_CV")

# pairs: c(comp1 = y-axis, comp_ref = x-axis), matching LD_X07 Part B orientation
plot_pairs = list(
  P1_R62H_vs_CV__vs__AD_vs_Control = c("R62H_vs_CV", "CV_AD_vs_Control"),
  P2_R47H_vs_CV__vs__AD_vs_Control = c("R47H_vs_CV", "CV_AD_vs_Control"),
  P3_R62H_vs_CV__vs__R47H_vs_CV    = c("R62H_vs_CV", "R47H_vs_CV")
)

quadrants  = c("up_up", "down_down", "up_down", "down_up")
quad_label = c(up_up = "up_up (shared up)", down_down = "down_down (shared down)",
               up_down = "up_down (divergent)", down_up = "down_up (divergent)")

pair_label = c(P1_R62H_vs_CV__vs__AD_vs_Control = "R62H vs CV\nvs AD vs Ctrl",
               P2_R47H_vs_CV__vs__AD_vs_Control = "R47H vs CV\nvs AD vs Ctrl",
               P3_R62H_vs_CV__vs__R47H_vs_CV    = "R62H vs CV\nvs R47H vs CV")


### load E02c per-subcluster DESeq results (basic 5-covariate model, no ladder)

e02  = qread(e02c_path)
dres = e02$deseq_results
rm(e02); gc()
message("   loaded ", length(dres), " E02c per-subcluster DESeq result tables")

### load E04c pooled DESeq results (single model across ALL astrocyte subclusters,
### cluster_name as covariate - the true "pooled" scope used in LD_X07 Part B1 /
### LD_E04b's original "pooled" scope; NOT a union of per-subcluster results)

e04           = qread(e04c_path)
e04_deseq_res = e04$E04_deseq_res
rm(e04); gc()
res_pooled = function(level, contrast) e04_deseq_res[[paste("pooled", level, contrast, sep = "|")]]


### gene universes
# DEG-union (E02c padj < 0.1, same restriction LD_X07 uses) for the foreground;
# expressed/tested genes (all genes across all loaded result tables) as the ORA
# background universe - matches LD_E04b conventions exactly.

deg_tab      = read_csv(deg_csv, show_col_types = FALSE)
deg_universe = unique(unlist(deg_tab, use.names = FALSE))
deg_universe = deg_universe[!is.na(deg_universe) & deg_universe != ""]

gene_universe = unique(unlist(lapply(dres, function(x) if (!is.null(x)) x$gene)))
message("   DEG-union (foreground pool): ", length(deg_universe),
        " | expressed background: ", length(gene_universe), " genes")


### helpers

classify_reg = function(lfc, p){
  r = rep("nreg", length(p))
  r[!is.na(p) & p < 0.05 & lfc > 0] = "up"
  r[!is.na(p) & p < 0.05 & lfc < 0] = "down"
  r
}

# reg_both quadrant gene list for one subcluster/pair (DEG-union restricted,
# nominal p < 0.05 per contrast) - identical logic to LD_E04b's get_quadrant_genes
get_quadrant_genes_cluster = function(cl, pair){
  r1 = dres[[paste0(cl, "_", e02c_tag[pair[1]])]]
  r2 = dres[[paste0(cl, "_", e02c_tag[pair[2]])]]
  if (is.null(r1) || is.null(r2)) return(NULL)
  g = intersect(intersect(r1$gene, r2$gene), deg_universe)
  if (length(g) == 0) return(NULL)
  r1 = r1[match(g, r1$gene), ]; r2 = r2[match(g, r2$gene), ]
  d = tibble(gene = g,
             lfc = r1$log2FoldChange, p = r1$pvalue,
             lfc_ref = r2$log2FoldChange, p_ref = r2$pvalue)
  d = d[!is.na(d$lfc) & !is.na(d$lfc_ref) & !is.na(d$p) & !is.na(d$p_ref), ]
  d$reg     = classify_reg(d$lfc, d$p)
  d$reg_ref = classify_reg(d$lfc_ref, d$p_ref)
  d = d[d$reg != "nreg" & d$reg_ref != "nreg", ]
  d$quadrant = paste0(d$reg, "_", d$reg_ref)
  split(d$gene, factor(d$quadrant, levels = quadrants))
}

# pool (union) quadrant gene sets across all subclusters in one family
get_quadrant_genes_family = function(family_name, pair){
  cls = family_clusters[[family_name]]
  pooled = setNames(vector("list", length(quadrants)), quadrants)
  for (cl in cls){
    ql = get_quadrant_genes_cluster(cl, pair)
    if (is.null(ql)) next
    for (q in quadrants) pooled[[q]] = union(pooled[[q]], ql[[q]])
  }
  pooled
}

# reg_both quadrant gene list from the single ALL-astrocyte pooled DESeq2 model
# (M0_base, LD_E04c) - no per-subcluster union needed, this is already one fit
get_quadrant_genes_allpooled = function(pair){
  r1 = res_pooled("M0_base", pair[1])
  r2 = res_pooled("M0_base", pair[2])
  if (is.null(r1) || is.null(r2)) return(setNames(vector("list", length(quadrants)), quadrants))
  g = intersect(intersect(r1$gene, r2$gene), deg_universe)
  if (length(g) == 0) return(setNames(vector("list", length(quadrants)), quadrants))
  r1 = r1[match(g, r1$gene), ]; r2 = r2[match(g, r2$gene), ]
  d = tibble(gene = g,
             lfc = r1$log2FoldChange, p = r1$pvalue,
             lfc_ref = r2$log2FoldChange, p_ref = r2$pvalue)
  d = d[!is.na(d$lfc) & !is.na(d$lfc_ref) & !is.na(d$p) & !is.na(d$p_ref), ]
  d$reg     = classify_reg(d$lfc, d$p)
  d$reg_ref = classify_reg(d$lfc_ref, d$p_ref)
  d = d[d$reg != "nreg" & d$reg_ref != "nreg", ]
  d$quadrant = paste0(d$reg, "_", d$reg_ref)
  split(d$gene, factor(d$quadrant, levels = quadrants))
}

# dispatch: "ALL_pooled" = single all-astrocyte model; otherwise per-family union
get_quadrant_genes_scope = function(scope, pair){
  if (scope == "ALL_pooled") get_quadrant_genes_allpooled(pair) else get_quadrant_genes_family(scope, pair)
}

# GO-BP over-representation on one gene set (matches LD_E04b params)
run_go = function(genes){
  if (length(genes) < 3) return(NULL)
  ego = tryCatch(
    enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
             ont = "BP", pAdjustMethod = "BH", universe = gene_universe,
             pvalueCutoff = 0.01, qvalueCutoff = 0.05),
    error = function(e) NULL)
  if (is.null(ego) || is.null(ego@result)) return(NULL)
  res = ego@result[ego@result$p.adjust <= 0.05 & ego@result$Count > 1, ]
  if (nrow(res) == 0) return(NULL)
  res
}


### run GO over all scope (ALL_pooled + 3 families) x pair x quadrant

scopes = c("ALL_pooled", names(families))

GO_all      = NULL
gene_counts = NULL

for (scope in scopes){
  for (pn in names(plot_pairs)){

    ql = get_quadrant_genes_scope(scope, plot_pairs[[pn]])

    for (q in quadrants){
      genes = ql[[q]]
      gene_counts = rbind(gene_counts, tibble(
        scope = scope, pair = pn, quadrant = q, n_genes = length(genes)))

      message("   GO: ", scope, " | ", pn, " | ", q,
              " (", length(genes), " pooled genes) - ", Sys.time())

      res = run_go(genes)
      if (!is.null(res)){
        GO_all = rbind(GO_all, cbind(scope = scope, pair = pn, quadrant = q, res))
      }
    }
  }
}

write_csv(gene_counts, file = file.path(out_dir, paste0(script_ind, "GO_quadrants_family_pooled_gene_counts.csv")))
if (!is.null(GO_all)) write_csv(GO_all, file = file.path(out_dir, paste0(script_ind, "GO_quadrants_family_pooled_results.csv")))


### dotplot: top 10 GO terms per quadrant, three pairs side by side, one page per scope
# NB: quadrant is mapped to a discrete x-axis; ggplot's default scale_x_discrete()
# has drop = TRUE, which silently removes any quadrant with zero significant terms
# from the axis (e.g. up_up often has 0 hits - see gene_counts vs GO_all). Setting
# drop = FALSE forces all 4 quadrants to always appear, empty or not, so a missing
# category reads as "0 significant terms", not as if that quadrant was never tested.

make_go_dotplot_scope = function(scope, title){
  if (is.null(GO_all)) return(NULL)
  t = GO_all[GO_all$scope == scope, ]
  if (nrow(t) == 0) return(NULL)

  # top 10 terms per (pair, quadrant) by adjusted p
  t = t %>% group_by(pair, quadrant) %>%
    slice_min(p.adjust, n = 10, with_ties = FALSE) %>% ungroup()

  t$quadrant = factor(t$quadrant, levels = quadrants, labels = quad_label[quadrants])
  t$pair     = factor(t$pair, levels = names(plot_pairs), labels = pair_label[names(plot_pairs)])

  ggplot(t, aes(x = quadrant, y = Description, size = Count, colour = -log10(p.adjust))) +
    geom_point() +
    facet_wrap(~ pair, nrow = 1, scales = "free_y") +
    scale_x_discrete(drop = FALSE, limits = quad_label[quadrants]) +
    scale_colour_viridis_c() +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(title, "  |  GO-BP per quadrant"),
         x = "quadrant", y = "GO biological process",
         colour = "-log10 padj", size = "gene count")
}

scope_title = c(ALL_pooled = "All astrocytes pooled (single model, cluster_name covariate)",
                setNames(paste0(names(families), "  (pooled union across ",
                                lengths(family_clusters), " subclusters)"), names(families)))

pdf(file = file.path(out_dir, paste0(script_ind, "GO_quadrants_family_pooled.pdf")), width = 16, height = 9)
for (scope in scopes){
  p = make_go_dotplot_scope(scope, scope_title[[scope]])
  if (!is.null(p)) plot(p)
}
dev.off()

message("    Saved GO-quadrant PDF (", length(scopes), " scopes: ALL_pooled + ", length(names(families)), " families)")


sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed LD_X08 ", Sys.time(),
        "\n##########################################################################\n\n\n")
