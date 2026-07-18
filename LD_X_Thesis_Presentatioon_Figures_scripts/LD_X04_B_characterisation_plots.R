# LD_X04: Astrocyte cluster-characterisation plots for the thesis (B scripts).
#   Umbrella script - add further characterisation plots as new sections below.
#
#   Plot 1 (MAIN): Green et al. (2024) astrocyte-state module-score dot plot,
#                  with renamed state labels.
#
# DATA: LD_B04a_v02_seur.qs (cleaned round-2 astrocytes; B03 embedding/clusters).
#   Module scores are RECOMPUTED here because the saved object does not store them
#   (B04a's only qsave precedes AddModuleScore). Run on the HPC R (qs + Seurat).

library(tidyverse)
library(qs)
library(Seurat)

### paths (all on RDS) ------------------------------------------------------
base       = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
b04_path   = file.path(base, "LD_B_AST_analysis_output/LD_B04a_v02_seur.qs")
green_csv  = file.path(base, "data_TREM2_michael/A_input/Green24_S2_subpopulation_markers.csv")
clust_csv  = file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X04_"
for (p in c(b04_path, green_csv, clust_csv))
  if (!file.exists(p)) stop("Missing input: ", p)

### load object + cluster order ---------------------------------------------
message("Loading B04 astrocyte object...")
seur = qread(b04_path)
DefaultAssay(seur) = "SCT"
ord = read_csv(clust_csv, show_col_types = FALSE)
cluster_names = unique(ord$cluster_name)   # subtype-grouped order (SLC1A2, GFAP, CHI3L1)

################################################################################
# Plot 1 (MAIN): Green et al. (2024) astrocyte-state module-score dot plot
################################################################################
# Build per-state signatures (same filter as B04a), score each with AddModuleScore,
# then plot module score per astrocyte subcluster.

green = read_csv(green_csv, show_col_types = FALSE)
ast_states = paste0("Ast.", 1:10)          # numeric order Ast.1 ... Ast.10

sig = lapply(ast_states, function(s) {
  g = green$gene[green$state == s &
                 green$avg_log2FC > log2(1.2) &
                 green$p_val_adj < 0.05]
  intersect(unique(g), rownames(seur))     # keep only genes present in the object
})
names(sig) = ast_states
message("Green24 astrocyte signature sizes (genes per state):")
print(lengths(sig))

set.seed(1234)
seur = AddModuleScore(seur, features = sig, name = "greenAst")
# AddModuleScore names columns greenAst1..greenAst10 in the order of `sig`
score_cols = paste0("greenAst", seq_along(ast_states))   # greenAst1 == Ast.1, ...

# display labels for the 10 Green astrocyte states (order = Ast.1 -> Ast.10)
green_labels = c(
  "AST1_homeostatic",
  "AST2_homeostatic",
  "Ast3_enh_mitophagy",
  "Ast4_reactive",
  "Ast5_reactive",
  "Ast6",
  "Ast7_IFN_response",
  "Ast8_stress response",
  "Ast9_stress response",
  "Ast10_disease_driving")
names(green_labels) = score_cols

p_green = DotPlot(seur, features = score_cols, group.by = "cluster_name",
                  scale.by = "size") +
  scale_x_discrete(labels = green_labels) +
  scale_y_discrete(limits = rev(cluster_names)) +
  labs(x = "Green et al. (2024) astrocyte state",
       y = "Astrocyte subcluster",
       title = "Module-score similarity to Green et al. astrocyte states") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5, face = "bold", size = 13),
        panel.grid  = element_line(linewidth = 0.2))

ggsave(file.path(out_dir, paste0(script_ind, "Green24_module_score_dotplot.pdf")),
       p_green, width = 8, height = 7, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "Green24_module_score_dotplot.png")),
       p_green, width = 8, height = 7, dpi = 300)
message("Plot 1 (Green24 dot plot) written.")

################################################################################
# Table: per-subcluster abundance (% of its subtype, and % of all astrocytes)
################################################################################
# "overcluster" = transcriptomic subtype (cell_type: AST_SLC1A2 / GFAP / CHI3L1).
abund = seur@meta.data %>%
  dplyr::count(cell_type, cluster_name, name = "n_cells") %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(pct_of_subtype = 100 * n_cells / sum(n_cells)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pct_of_total = 100 * n_cells / sum(n_cells)) %>%
  dplyr::arrange(cell_type, dplyr::desc(n_cells))

write_csv(abund, file.path(out_dir, paste0(script_ind, "cluster_abundance_pct.csv")))
message("Per-subcluster abundance (n, % of subtype, % of total):")
print(as.data.frame(abund), digits = 3)

################################################################################
# Plot: GO (BP) over-representation - two largest subclusters per subtype, AND
#       all subclusters that have significant enrichment
################################################################################
library(clusterProfiler)
library(org.Hs.eg.db)

# pick the two most abundant subclusters in each subtype (from the table above)
sel = abund %>% dplyr::group_by(cell_type) %>%
  dplyr::slice_max(n_cells, n = 2, with_ties = FALSE) %>% dplyr::ungroup()
sel_clusters = sel$cluster_name[order(match(sel$cluster_name, cluster_names))]  # SLC1A2 > GFAP > CHI3L1 order
message("GO: two largest subclusters per subtype: ", paste(sel_clusters, collapse = ", "))

# de novo marker genes per subcluster (same definition as B04a; SCT assay)
seur = PrepSCTFindMarkers(seur)            # required before marker testing on SCT
Idents(seur) = "cluster_name"
mk = FindAllMarkers(seur, only.pos = TRUE) %>%   # Wilcoxon rank-sum (default)
  dplyr::filter(p_val_adj < 0.001, avg_log2FC > 0.5)

gene_lists = split(mk$gene, mk$cluster)                    # ALL subclusters with markers
gene_lists = lapply(gene_lists, unique)
gene_lists = gene_lists[!vapply(gene_lists, is.null, logical(1))]
gene_lists = gene_lists[order(match(names(gene_lists), cluster_names))]  # biology-driven order
message("Marker genes per subcluster:"); print(lengths(gene_lists))

# GO BP over-representation, all subclusters (clusters with no enriched term are
# simply absent from the result)
cc = compareCluster(geneClusters = gene_lists, fun = "enrichGO",
                    OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                    universe = rownames(seur),
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

write_csv(as.data.frame(cc), file.path(out_dir, paste0(script_ind, "GO_top_terms_table.csv")))

# --- build the dot plot manually: lets us (i) keep clusters with NO significant
#     terms as empty columns (e.g. s0), and (ii) control the term (y) order ------
res = as.data.frame(cc)
res$GeneRatio_num = vapply(strsplit(res$GeneRatio, "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))

# make a top-N GO dot plot over a given set of subclusters (x-axis) and save it.
# Subclusters in `clusters` with no enriched term appear as empty columns; to
# exclude term-less subclusters, just pass only the clusters that have terms.
make_go_dotplot = function(res, clusters, n_terms, suffix, title) {
  res = res %>% dplyr::filter(Cluster %in% clusters)
  if (nrow(res) == 0) { message("No GO terms for ", suffix, " - skipped"); return(invisible()) }
  topn = res %>% dplyr::group_by(Cluster) %>%
    dplyr::slice_min(p.adjust, n = n_terms, with_ties = FALSE) %>% dplyr::ungroup()
  terms_show = unique(topn$Description)

  # plotting data = every significant (subcluster, term) pair among the shown terms
  pdat = res %>% dplyr::filter(Description %in% terms_show)
  pdat$Cluster = factor(pdat$Cluster, levels = clusters)

  # y-axis order: group each term by the subtype of its most significant subcluster,
  # ordered SLC1A2 -> GFAP -> CHI3L1 so that CHI3L1 terms end up at the TOP of the axis
  term_ord = topn %>% dplyr::group_by(Description) %>%
    dplyr::slice_min(p.adjust, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>%
    dplyr::mutate(subtype = sub("_s[0-9]+$", "", Cluster),
                  srank   = dplyr::recode(subtype, AST_SLC1A2 = 1, AST_GFAP = 2, AST_CHI3L1 = 3)) %>%
    dplyr::arrange(srank, p.adjust)
  pdat$Description = factor(pdat$Description, levels = term_ord$Description)

  p = ggplot(pdat, aes(x = Cluster, y = Description,
                       size = Count, colour = -log10(p.adjust))) +   # match B04a: size = gene count, colour = -log10(adj p)
    geom_point() +
    # add = 0.55 trims the empty margin on either side so columns sit closer together
    scale_x_discrete(limits = clusters, drop = FALSE, expand = expansion(add = 0.55)) +
    scale_size_continuous(name = "Gene count", range = c(1, 6)) +
    scale_colour_gradient(name = "-log10(adj. p)", low = "blue", high = "red") +   # red = more significant, as in B04a
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11),
          plot.title  = element_text(hjust = 0.5, face = "bold", size = 13,
                                     margin = margin(b = 6)),
          plot.margin = margin(t = 10, r = 8, b = 4, l = 6))  # top space keeps the title from being clipped

  h = max(6, 0.30 * length(terms_show) + 2.4)  # height scales with number of terms (extra for title)
  # per-column width pulls the subcluster columns together; the +3 base reserves
  # room for the (now larger) GO-term labels on the y-axis
  w = max(6, 0.52 * length(clusters) + 3)      # width scales with number of clusters
  ggsave(file.path(out_dir, paste0(script_ind, suffix, ".pdf")),
         p, width = w, height = h, useDingbats = FALSE)
  ggsave(file.path(out_dir, paste0(script_ind, suffix, ".png")),
         p, width = w, height = h, dpi = 300)
  message("GO dot plot written: ", suffix, " (", length(clusters), " clusters, ", length(terms_show), " terms)")
}

# two largest subclusters per subtype (s0 etc. shown as empty if no terms)
make_go_dotplot(res, sel_clusters, 5,  "GO_top5_main_subclusters",
                "Top-5 GO BP terms: two largest subclusters per subtype")
make_go_dotplot(res, sel_clusters, 10, "GO_top10_main_subclusters",
                "Top-10 GO BP terms: two largest subclusters per subtype")

# ALL subclusters that have significant GO terms (term-less ones excluded)
sig_clusters = cluster_names[cluster_names %in% unique(as.character(res$Cluster))]
make_go_dotplot(res, sig_clusters, 5, "GO_top5_all_sig_subclusters",
                "Top-5 GO BP terms: subclusters with significant enrichment")

################################################################################
# Plot N: <add next characterisation plot here>
################################################################################

message("Done. Outputs in ", out_dir)
