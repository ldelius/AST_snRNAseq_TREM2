# LD_X04c: sensitivity check for the Green et al. (2024) astrocyte-state mapping
# (LD_X04 Plot 1) - does the state <-> subcluster ranking survive removing the
# three subtype-defining genes (SLC1A2, GFAP, CHI3L1) from the state signatures?
#
# Deliberately kept separate from LD_X04_B_characterisation_plots.R: that script's
# FindAllMarkers/GO section is what was hitting the OnDemand session memory limit,
# and this check doesn't need any of it - just the B04 object + AddModuleScore.
#
# For each of the 10 Green states:
#   - "full"  score = AddModuleScore using the state's full gene list
#   - "loo"   score = AddModuleScore using the state's gene list with
#                     SLC1A2/GFAP/CHI3L1 removed (leave-one-[gene]-out)
#   mean score per subcluster -> z-scaled across subclusters (within that state,
#   full and loo separately) -> compare: same top-ranked subcluster? correlated
#   ranking overall (Pearson r on the z-profiles)?
#
# DATA: LD_B04a_v02_seur.qs. Run on the HPC R (qs + Seurat).

library(tidyverse)
library(qs)
library(Seurat)
library(patchwork)

### paths (all on RDS) --------------------------------------------------------
base       = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
b04_path   = file.path(base, "LD_B_AST_analysis_output/LD_B04a_v02_seur.qs")
green_csv  = file.path(base, "data_TREM2_michael/A_input/Green24_S2_subpopulation_markers.csv")
clust_csv  = file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X04c_v01_"
for (p in c(b04_path, green_csv, clust_csv))
  if (!file.exists(p)) stop("Missing input: ", p)

### load object + cluster order ------------------------------------------------
message("Loading B04 astrocyte object...")
seur = qread(b04_path)
DefaultAssay(seur) = "SCT"
ord = read_csv(clust_csv, show_col_types = FALSE)
cluster_names = unique(ord$cluster_name)

### build full vs LOO (SLC1A2/GFAP/CHI3L1-excluded) state signatures -----------
# same filter as LD_X04 Plot 1, so the "full" scores here reproduce that plot.
green = read_csv(green_csv, show_col_types = FALSE)
ast_states = paste0("Ast.", 1:10)

sig = lapply(ast_states, function(s) {
  g = green$gene[green$state == s &
                 green$avg_log2FC > log2(1.2) &
                 green$p_val_adj < 0.05]
  intersect(unique(g), rownames(seur))
})
names(sig) = ast_states

subtype_marker_genes = c("SLC1A2", "GFAP", "CHI3L1")
sig_loo = lapply(sig, setdiff, y = subtype_marker_genes)
names(sig_loo) = ast_states

n_removed = mapply(function(a, b) length(a) - length(b), sig, sig_loo)
message("Subtype-marker genes present (and removed) per state (of ", length(subtype_marker_genes), " possible):")
print(n_removed)

# display labels for the 10 Green astrocyte states (order = Ast.1 -> Ast.10),
# same as LD_X04 Plot 1
green_labels = c(
  "Ast1 homeostatic", "Ast2 homeostatic", "Ast3 enhanced mitophagy", "Ast4 reactive",
  "Ast5 reactive", "Ast6", "Ast7 IFN response", "Ast8 stress response",
  "Ast9 stress response", "Ast10 AD-elevated")
names(green_labels) = ast_states

### score with both signature sets ---------------------------------------------
set.seed(1234)
seur = AddModuleScore(seur, features = sig,     name = "Green24_")   # -> Green24_1 .. Green24_10
set.seed(1234)
seur = AddModuleScore(seur, features = sig_loo, name = "GreenLOO_")  # -> GreenLOO_1 .. GreenLOO_10

### dot plots: full vs LOO, same style as LD_X04 Plot 1 -------------------------
# DotPlot just reads already-computed metadata columns, so this is cheap - fine
# to keep in this lightweight script.
score_cols_full = paste0("Green24_",  seq_along(ast_states))
score_cols_loo  = paste0("GreenLOO_", seq_along(ast_states))

slc1a2_clusters = cluster_names[grepl("^AST_SLC1A2", cluster_names)]
gfap_clusters   = cluster_names[grepl("^AST_GFAP",   cluster_names)]
chi3l1_clusters = cluster_names[grepl("^AST_CHI3L1", cluster_names)]
y_order = rev(cluster_names)

# builds a highlight_box() closure bound to a given set of score columns/labels,
# so the same helper works for both the full and LOO dot plots.
make_highlight_box = function(score_cols, labels) {
  x_order = unname(labels[score_cols])
  function(x_names, y_names, colour = "black", linewidth = 0.9) {
    xi = match(x_names, x_order); yi = match(y_names, y_order)
    if (anyNA(xi) || anyNA(yi)) stop("highlight_box: name not found in axis order")
    annotate("rect", xmin = min(xi) - 0.5, xmax = max(xi) + 0.5,
                      ymin = min(yi) - 0.5, ymax = max(yi) + 0.5,
                      fill = NA, colour = colour, linewidth = linewidth)
  }
}

make_green_dotplot = function(score_cols, labels) {
  hb = make_highlight_box(score_cols, labels)
  DotPlot(seur, features = score_cols, group.by = "cluster_name", scale.by = "size") +
    scale_x_discrete(labels = labels[score_cols]) +
    scale_y_discrete(limits = rev(cluster_names)) +
    hb(c("Ast1 homeostatic", "Ast2 homeostatic"), slc1a2_clusters) +
    hb("Ast4 reactive", gfap_clusters) +
    hb("Ast5 reactive", chi3l1_clusters) +
    # diverging, colourblind-safe orange-blue scale (Wong 2011, Nat Methods - the
    # recommended CVD-safe alternative to red-green diverging scales), centred on
    # 0 since this is a genuinely signed z-score (above/below average for that
    # state), not a one-directional magnitude. Reuses the exact Okabe-Ito hex
    # values already used for AST_GFAP (#0072B2) / AST_SLC1A2 (#E69F00) in X02.
    scale_colour_gradient2(low = "#0072B2", mid = "grey90", high = "#E69F00",
                          midpoint = 0, name = "Mean module score\n(z-scaled across subclusters)") +
    labs(x = "Green et al. (2024) astrocyte state", y = NULL,
         size = "% cells, score > 0") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_line(linewidth = 0.2))
}

# score_cols index (1..10) maps 1:1 onto ast_states/green_labels (Ast.1..Ast.10)
labels_full = setNames(green_labels, score_cols_full)
labels_loo  = setNames(green_labels, score_cols_loo)

p_full = make_green_dotplot(score_cols_full, labels_full)
p_loo  = make_green_dotplot(score_cols_loo,  labels_loo)

ggsave(file.path(out_dir, paste0(script_ind, "Green24_dotplot_full.pdf")),
       p_full, width = 8, height = 7, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "Green24_dotplot_full.png")),
       p_full, width = 8, height = 7, dpi = 300)
ggsave(file.path(out_dir, paste0(script_ind, "Green24_dotplot_LOO_excl_subtype_markers.pdf")),
       p_loo, width = 8, height = 7, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "Green24_dotplot_LOO_excl_subtype_markers.png")),
       p_loo, width = 8, height = 7, dpi = 300)
message("Dot plots written (full + LOO).")

################################################################################
# Combined figure: Green24 module-score plot (A) + curated GO term plot (B)
################################################################################
# Built directly here (not via a separate saveRDS-cached combining script - that
# approach silently serialised the whole Seurat object via ggplot's captured
# plot_env, which is what caused the earlier multi-GB / stuck-saving problem).
# The GO panel is rebuilt from the already-written GO_top_terms_table.csv (from
# LD_X04_B_characterisation_plots.R's compareCluster/enrichGO result) - cheap,
# CSV-only, no Seurat/enrichGO rerun. `cluster_names`, `out_dir`, `base`,
# `script_ind`, and `p_full` are all already available from earlier in this script.
go_csv = file.path(base, "LD_X_Thesis_Presentation_output/LD_X04_v02_GO_top_terms_table.csv")

if (!file.exists(go_csv)) {
  message("Skipping combined Green+GO figure - missing: ", go_csv)
} else {
  res = read_csv(go_csv, show_col_types = FALSE)
  sig_clusters = cluster_names[cluster_names %in% unique(as.character(res$Cluster))]

  # same exclusion list / curated row order as LD_X04_v02_replot_curated_GO.R
  curated_exclude_terms = c(
    "muscle contraction", "muscle system process", "stem cell differentiation",
    "regulation of plasma membrane organization", "regulation of receptor recycling",
    "regulation of cell-substrate adhesion", "axonogenesis",
    "regulation of axon extension involved in axon guidance", "regulation of axonogenesis",
    "regulation of timing of cell differentiation", "regulation of neurogenesis",
    "potassium ion transport", "regulation of developmental growth",
    "regulation of nervous system development", "extracellular structure organization",
    "external encapsulating structure organization", "regulation of trans-synaptic signaling",
    "cellular response to type II interferon", "negative regulation of lipase activity",
    "cellular response to chemical stress", "regulation of postsynaptic membrane potential")

  curated_term_order = c(
    "ribosome biogenesis", "response to oxidative stress", "response to type II interferon",
    "ERK1 and ERK2 cascade", "negative regulation of catalytic activity",
    "cytokine-mediated signaling pathway", "positive regulation of cytokine production",
    "response to interleukin-1", "wound healing", "chemotaxis",
    "cell-substrate adhesion", "extracellular matrix organization",
    "proteoglycan biosynthetic process", "morphogenesis of a branching structure",
    "synapse organization", "regulation of neuron projection development",
    "axon development", "potassium ion transmembrane transport",
    "ionotropic glutamate receptor signaling pathway", "regulation of membrane potential",
    "modulation of chemical synaptic transmission")
  n_terms = 5
  term_hline_after = 5   # dashed line after the 5th term (end of CHI3L1-only block)

  topn = res %>% dplyr::filter(Cluster %in% sig_clusters) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::slice_min(p.adjust, n = n_terms, with_ties = FALSE) %>% dplyr::ungroup()
  terms_show = unique(topn$Description)

  pdat_go = res %>% dplyr::filter(Cluster %in% sig_clusters, Description %in% terms_show,
                                  !Description %in% curated_exclude_terms)
  go_clusters = sig_clusters[sig_clusters %in% unique(as.character(pdat_go$Cluster))]
  terms_show  = terms_show[!terms_show %in% curated_exclude_terms]

  if (!setequal(curated_term_order, terms_show))
    stop("curated_term_order does not match the term set exactly (missing or extra terms)")

  pdat_go$Cluster     = factor(pdat_go$Cluster, levels = go_clusters)
  pdat_go$Description = factor(pdat_go$Description, levels = rev(curated_term_order))

  p_go = ggplot(pdat_go, aes(x = Cluster, y = Description,
                             size = Count, colour = -log10(p.adjust))) +
    geom_point() +
    scale_x_discrete(limits = go_clusters, drop = FALSE, expand = expansion(add = 0.55)) +
    scale_size_continuous(name = "Gene count", range = c(1, 6)) +
    scale_colour_gradient(name = "-log10(adj. p)", low = "grey80", high = "black") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 13)

  go_subtype   = sub("_s[0-9]+$", "", go_clusters)
  v_boundaries = which(diff(match(go_subtype, unique(go_subtype))) != 0) + 0.5
  p_go = p_go + geom_vline(xintercept = v_boundaries, linetype = "dashed",
                           colour = "grey40", linewidth = 0.4)

  h_boundary = length(curated_term_order) - term_hline_after + 0.5
  p_go = p_go + geom_hline(yintercept = h_boundary, linetype = "dashed",
                           colour = "grey40", linewidth = 0.4)

  p_go = p_go +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 12.5),
          plot.margin = margin(t = 10, r = 8, b = 4, l = 6))

  # side by side: Green (A) = 45% width, GO (B) = 55%; heights deliberately
  # NOT matched between the two panels. Legend TEXT (not key size) shrunk here,
  # combined-figure-only - the legend keys are already slim, it's the text that
  # was padding out the legend's width and stealing space from the panels.
  fig = (p_full | p_go) +
    patchwork::plot_layout(widths = c(0.45, 0.55)) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag    = element_text(face = "bold", size = 14),
          legend.text  = element_text(size = 7),
          legend.title = element_text(size = 8))

  ggsave(file.path(out_dir, paste0(script_ind, "Fig_Green_GO_combined.pdf")),
        fig, width = 16, height = 9, useDingbats = FALSE)
  ggsave(file.path(out_dir, paste0(script_ind, "Fig_Green_GO_combined.png")),
        fig, width = 16, height = 9, dpi = 300)
  message("Combined Green+GO figure written.")
}

### full vs LOO comparison, per state -------------------------------------------
md = seur@meta.data

get_z = function(pat) {
  md %>%
    dplyr::select(cluster_name, dplyr::matches(pat)) %>%
    tidyr::pivot_longer(-cluster_name, names_to = "state", values_to = "score") %>%
    dplyr::mutate(state = sub(pat, "", state)) %>%
    dplyr::group_by(state, cluster_name) %>%
    dplyr::summarise(m = mean(score), .groups = "drop") %>%
    dplyr::group_by(state) %>%
    dplyr::mutate(z = as.numeric(scale(m))) %>%
    dplyr::ungroup()
}

full = get_z("^Green24_")
loo  = get_z("^GreenLOO_")

# per state: top-ranked subcluster, full vs LOO signature
top_by_state = full %>%
  dplyr::group_by(state) %>% dplyr::slice_max(z, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>% dplyr::select(state, full_top = cluster_name) %>%
  dplyr::left_join(
    loo %>% dplyr::group_by(state) %>% dplyr::slice_max(z, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>% dplyr::select(state, loo_top = cluster_name),
    by = "state") %>%
  dplyr::mutate(state_ix  = as.integer(state),
                state_lab = unname(green_labels[paste0("Ast.", state_ix)]),
                genes_removed = unname(n_removed[paste0("Ast.", state_ix)]),
                top_unchanged = as.character(full_top) == as.character(loo_top)) %>%
  dplyr::arrange(state_ix) %>%
  dplyr::select(state_lab, genes_removed, full_top, loo_top, top_unchanged)

# per state: correlation of the full z-profile against the LOO z-profile across
# all subclusters (not just the top one) - how much does the WHOLE ranking move
rank_corr_by_state = full %>% dplyr::rename(z_full = z) %>%
  dplyr::select(state, cluster_name, z_full) %>%
  dplyr::left_join(loo %>% dplyr::select(state, cluster_name, z_loo = z),
                   by = c("state", "cluster_name")) %>%
  dplyr::group_by(state) %>%
  dplyr::summarise(r = as.numeric(cor(z_full, z_loo)), .groups = "drop") %>%
  dplyr::mutate(state_ix  = as.integer(state),
                state_lab = unname(green_labels[paste0("Ast.", state_ix)]),
                genes_removed = unname(n_removed[paste0("Ast.", state_ix)])) %>%
  dplyr::arrange(state_ix) %>%
  dplyr::select(state_lab, genes_removed, r)

# per state: FULL ranking of every subcluster (not just the winner), full vs LOO,
# so you can see how far any given subcluster moved, not just whether the top
# pick changed. Rank 1 = highest z (top-scoring) within that state.
full_ranking_by_state = full %>% dplyr::rename(z_full = z) %>%
  dplyr::select(state, cluster_name, z_full) %>%
  dplyr::left_join(loo %>% dplyr::select(state, cluster_name, z_loo = z),
                   by = c("state", "cluster_name")) %>%
  dplyr::group_by(state) %>%
  dplyr::mutate(full_rank = dplyr::min_rank(dplyr::desc(z_full)),
                loo_rank  = dplyr::min_rank(dplyr::desc(z_loo))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(state_ix  = as.integer(state),
                state_lab = unname(green_labels[paste0("Ast.", state_ix)]),
                genes_removed = unname(n_removed[paste0("Ast.", state_ix)]),
                rank_shift = full_rank - loo_rank) %>%
  dplyr::arrange(state_ix, full_rank) %>%
  dplyr::select(state_lab, genes_removed, cluster_name,
               z_full, full_rank, z_loo, loo_rank, rank_shift)

message("Top-ranked subcluster per state, full signature vs marker-excluded (LOO):")
print(as.data.frame(top_by_state))
message("Correlation of full vs LOO z-scaled subcluster profile, per state:")
print(as.data.frame(rank_corr_by_state))
message("Full subcluster ranking per state, full vs LOO (rank 1 = top-scoring):")
print(as.data.frame(full_ranking_by_state))

write_csv(top_by_state,      file.path(out_dir, paste0(script_ind, "LOO_sensitivity_top_subcluster_by_state.csv")))
write_csv(rank_corr_by_state, file.path(out_dir, paste0(script_ind, "LOO_sensitivity_rank_correlation_by_state.csv")))
write_csv(full_ranking_by_state, file.path(out_dir, paste0(script_ind, "LOO_sensitivity_full_ranking_by_state.csv")))

### small diagnostic plot: correlation per state --------------------------------
p_corr = ggplot(rank_corr_by_state, aes(x = reorder(state_lab, r), y = r)) +
  geom_col(fill = "grey40") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  coord_flip(ylim = c(min(0, min(rank_corr_by_state$r)), 1)) +
  labs(x = NULL, y = "Pearson r (full vs marker-excluded subcluster z-profile)",
       title = "Green24 state mapping: sensitivity to removing SLC1A2/GFAP/CHI3L1") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, paste0(script_ind, "LOO_sensitivity_correlation_barplot.pdf")),
       p_corr, width = 7, height = 5, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "LOO_sensitivity_correlation_barplot.png")),
       p_corr, width = 7, height = 5, dpi = 300)

message("Done. Outputs in ", out_dir)
