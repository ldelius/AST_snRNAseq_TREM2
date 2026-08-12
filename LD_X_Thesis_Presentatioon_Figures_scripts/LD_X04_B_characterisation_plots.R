# LD_X04: Astrocyte-state module scores and subcluster GO enrichment.
# Module scores are recomputed because they are not stored in the B04a object.

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
script_ind = "LD_X04_v03_"
for (p in c(b04_path, green_csv, clust_csv))
  if (!file.exists(p)) stop("Missing input: ", p)

### load object + cluster order ---------------------------------------------
message("Loading B04 astrocyte object...")
seur = qread(b04_path)
DefaultAssay(seur) = "SCT"
ord = read_csv(clust_csv, show_col_types = FALSE)
cluster_names = unique(ord$cluster_name)   # subtype-grouped order (SLC1A2, GFAP, CHI3L1)

### Green et al. astrocyte-state module scores -------------------------------

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
  "Ast1 homeostatic",
  "Ast2 homeostatic",
  "Ast3 enhanced mitophagy",
  "Ast4 reactive",
  "Ast5 reactive",
  "Ast6",
  "Ast7 IFN response",
  "Ast8 stress response",
  "Ast9 stress response",
  "Ast10 AD-elevated")
names(green_labels) = score_cols

# draw a rectangle around a block of the dot grid, given the group of x labels
# (Green state display names) and y labels (cluster_name values) it should span.
# x_order/y_order must be the same vectors used in scale_x_discrete()/scale_y_discrete()
# above, since box position is just the matched index range on the discrete axes.
x_order = unname(green_labels[score_cols])
y_order = rev(cluster_names)
highlight_box = function(x_names, y_names, colour = "black", linewidth = 0.9) {
  xi = match(x_names, x_order); yi = match(y_names, y_order)
  if (anyNA(xi) || anyNA(yi)) stop("highlight_box: name not found in axis order")
  annotate("rect", xmin = min(xi) - 0.5, xmax = max(xi) + 0.5,
                    ymin = min(yi) - 0.5, ymax = max(yi) + 0.5,
                    fill = NA, colour = colour, linewidth = linewidth)
}

slc1a2_clusters  = cluster_names[grepl("^AST_SLC1A2", cluster_names)]
gfap_clusters    = cluster_names[grepl("^AST_GFAP",   cluster_names)]
chi3l1_clusters  = cluster_names[grepl("^AST_CHI3L1", cluster_names)]

p_green = DotPlot(seur, features = score_cols, group.by = "cluster_name",
                  scale.by = "size") +
  scale_x_discrete(labels = green_labels) +
  scale_y_discrete(limits = rev(cluster_names)) +
  highlight_box(c("Ast1 homeostatic", "Ast2 homeostatic"), slc1a2_clusters) +
  highlight_box("Ast4 reactive", gfap_clusters) +
  highlight_box("Ast5 reactive", chi3l1_clusters) +
  highlight_box("Ast3 enhanced mitophagy", c("AST_CHI3L1_s9", "AST_CHI3L1_s16"),
                colour = "grey50") +
  highlight_box("Ast7 IFN response", chi3l1_clusters, colour = "grey50") +
  highlight_box(c("Ast8 stress response", "Ast9 stress response", "Ast10 AD-elevated"),
                "AST_CHI3L1_s9", colour = "grey50") +
  labs(x = "Green et al. (2024) astrocyte state",
       y = "Astrocyte subcluster",
       colour = "Mean module score\n(z-scaled across subclusters)",
       size   = "% cells, score > 0") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_line(linewidth = 0.2))

ggsave(file.path(out_dir, paste0(script_ind, "Green24_module_score_dotplot.pdf")),
       p_green, width = 8, height = 7, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "Green24_module_score_dotplot.png")),
       p_green, width = 8, height = 7, dpi = 300)
message("Plot 1 (Green24 dot plot) written.")

### per-subcluster abundance -------------------------------------------------
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

### GO BP over-representation ------------------------------------------------
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
make_go_dotplot = function(res, clusters, n_terms, suffix, title, exclude_terms = NULL,
                           family_sep = FALSE, term_order = NULL) {
  res = res %>% dplyr::filter(Cluster %in% clusters)
  if (nrow(res) == 0) { message("No GO terms for ", suffix, " - skipped"); return(invisible()) }
  topn = res %>% dplyr::group_by(Cluster) %>%
    dplyr::slice_min(p.adjust, n = n_terms, with_ties = FALSE) %>% dplyr::ungroup()
  terms_show = unique(topn$Description)

  # plotting data = every significant (subcluster, term) pair among the shown terms
  pdat = res %>% dplyr::filter(Description %in% terms_show)
  if (!is.null(exclude_terms)) {
    # Remove selected generic terms without reranking and omit empty subclusters.
    pdat = pdat %>% dplyr::filter(!Description %in% exclude_terms)
    clusters = clusters[clusters %in% unique(as.character(pdat$Cluster))]
    terms_show = terms_show[!terms_show %in% exclude_terms]
  }
  pdat$Cluster = factor(pdat$Cluster, levels = clusters)

  if (!is.null(term_order)) {
    # manual curation order: term_order is given top-to-bottom as it should read on
    # the plot; ggplot's discrete y-axis puts factor level 1 at the BOTTOM, so reverse it
    if (!setequal(term_order, terms_show))
      stop("term_order does not match the term set exactly (missing or extra terms)")
    pdat$Description = factor(pdat$Description, levels = rev(term_order))
  } else {
    # y-axis order: group each term by the subtype of its most significant subcluster,
    # ordered SLC1A2 -> GFAP -> CHI3L1 so that CHI3L1 terms end up at the TOP of the axis
    term_ord = topn %>% dplyr::group_by(Description) %>%
      dplyr::slice_min(p.adjust, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>%
      dplyr::mutate(subtype = sub("_s[0-9]+$", "", Cluster),
                    srank   = dplyr::recode(subtype, AST_SLC1A2 = 1, AST_GFAP = 2, AST_CHI3L1 = 3)) %>%
      dplyr::arrange(srank, p.adjust)
    pdat$Description = factor(pdat$Description, levels = term_ord$Description)
  }

  p = ggplot(pdat, aes(x = Cluster, y = Description,
                       size = Count, colour = -log10(p.adjust))) +   # match B04a: size = gene count, colour = -log10(adj p)
    geom_point() +
    # add = 0.55 trims the empty margin on either side so columns sit closer together
    scale_x_discrete(limits = clusters, drop = FALSE, expand = expansion(add = 0.55)) +
    scale_size_continuous(name = "Gene count", range = c(1, 6)) +
    # grey -> black: GO enrichment here is magnitude-only (all gene lists are
    # only.pos = TRUE markers, so there is no "down" direction to represent) -
    # an achromatic ramp avoids the blue/orange hue family used for the genuinely
    # signed Green24 z-scores elsewhere in this figure group, so darker here can't
    # be misread as "toward the negative/blue end" of that unrelated scale
    scale_colour_gradient(name = "-log10(adj. p)", low = "grey80", high = "black") +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 13)

  if (family_sep) {
    # dashed vertical line at every point where the subtype prefix (SLC1A2/GFAP/CHI3L1)
    # changes between adjacent columns, so the three cluster families read as visual groups
    subtype    = sub("_s[0-9]+$", "", clusters)
    boundaries = which(diff(match(subtype, unique(subtype))) != 0) + 0.5
    p = p + geom_vline(xintercept = boundaries, linetype = "dashed",
                       colour = "grey40", linewidth = 0.4)
  }

  p = p +
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

# Curated top-five terms after removing selected generic terms without reranking.
curated_exclude_terms = c(
  "muscle contraction",
  "muscle system process",
  "stem cell differentiation",
  "regulation of plasma membrane organization",
  "regulation of receptor recycling",
  "regulation of cell-substrate adhesion",
  "axonogenesis",
  "regulation of axon extension involved in axon guidance",
  "regulation of axonogenesis",
  "regulation of timing of cell differentiation",
  "regulation of neurogenesis",
  "potassium ion transport",
  "regulation of developmental growth",
  "regulation of nervous system development",
  "extracellular structure organization",
  "external encapsulating structure organization",
  "regulation of trans-synaptic signaling",
  "cellular response to type II interferon",
  "negative regulation of lipase activity"
)
# manual row order, top-to-bottom: CHI3L1-specific, then shared-immune, then
# non-discriminative, then the SLC1A2/neuronal-associated block at the bottom
curated_term_order = c(
  # CHI3L1-specific
  "response to type II interferon",
  "ERK1 and ERK2 cascade",
  "ribosome biogenesis",
  "response to oxidative stress",
  "negative regulation of catalytic activity",
  "cellular response to chemical stress",
  # shared immune
  "cytokine-mediated signaling pathway",
  "wound healing",
  "positive regulation of cytokine production",
  "response to interleukin-1",
  "chemotaxis",
  # non-discriminative
  "cell-substrate adhesion",
  "extracellular matrix organization",
  "morphogenesis of a branching structure",
  "synapse organization",
  # bottom block
  "regulation of postsynaptic membrane potential",
  "ionotropic glutamate receptor signaling pathway",
  "axon development",
  "proteoglycan biosynthetic process",
  "potassium ion transmembrane transport",
  "regulation of neuron projection development",
  "regulation of membrane potential",
  "modulation of chemical synaptic transmission"
)
make_go_dotplot(res, sig_clusters, 5, "GO_top5_all_sig_subclusters_curated", title = NULL,
                exclude_terms = curated_exclude_terms, family_sep = TRUE,
                term_order = curated_term_order)

message("Done. Outputs in ", out_dir)
