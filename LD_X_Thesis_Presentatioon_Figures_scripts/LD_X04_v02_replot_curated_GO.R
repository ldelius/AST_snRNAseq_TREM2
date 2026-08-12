# LD_X04_v02: Replot curated subcluster GO enrichment from saved results.

library(tidyverse)

base       = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
go_csv     = file.path(base, "LD_X_Thesis_Presentation_output/LD_X04_v02_GO_top_terms_table.csv")
clust_csv  = file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
script_ind = "LD_X04_v02_"
for (p in c(go_csv, clust_csv)) if (!file.exists(p)) stop("Missing input: ", p)

res  = read_csv(go_csv, show_col_types = FALSE)
ord  = read_csv(clust_csv, show_col_types = FALSE)
cluster_names = unique(ord$cluster_name)
sig_clusters  = cluster_names[cluster_names %in% unique(as.character(res$Cluster))]

# same exclusion list as before (unchanged)
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
  "negative regulation of lipase activity",
  "cellular response to chemical stress",
  "regulation of postsynaptic membrane potential"
)

# new manual row order, top-to-bottom:
#   - CHI3L1-only block, reordered (rows 1-5)
#   - horizontal dashed line after row 5
#   - shared-immune block, "wound healing" moved 2 down
#   - non-discriminative block, "proteoglycan biosynthetic process" moved up next
#     to "extracellular matrix organization"
#   - bottom block, with the neuronal/synaptic trio promoted to the top of the block
curated_term_order = c(
  # CHI3L1-specific
  "ribosome biogenesis",
  "response to oxidative stress",
  "response to type II interferon",
  "ERK1 and ERK2 cascade",
  "negative regulation of catalytic activity",
  # shared immune
  "cytokine-mediated signaling pathway",
  "positive regulation of cytokine production",
  "response to interleukin-1",
  "wound healing",
  "chemotaxis",
  # non-discriminative
  "cell-substrate adhesion",
  "extracellular matrix organization",
  "proteoglycan biosynthetic process",
  "morphogenesis of a branching structure",
  "synapse organization",
  # bottom block (promoted trio first)
  "regulation of neuron projection development",
  "axon development",
  "potassium ion transmembrane transport",
  "ionotropic glutamate receptor signaling pathway",
  "regulation of membrane potential",
  "modulation of chemical synaptic transmission"
)
n_terms = 5
term_hline_after = 5   # dashed horizontal line after the 5th term (end of CHI3L1-only block)

topn = res %>% dplyr::filter(Cluster %in% sig_clusters) %>%
  dplyr::group_by(Cluster) %>%
  dplyr::slice_min(p.adjust, n = n_terms, with_ties = FALSE) %>% dplyr::ungroup()
terms_show = unique(topn$Description)

pdat = res %>% dplyr::filter(Cluster %in% sig_clusters, Description %in% terms_show,
                             !Description %in% curated_exclude_terms)
clusters   = sig_clusters[sig_clusters %in% unique(as.character(pdat$Cluster))]
terms_show = terms_show[!terms_show %in% curated_exclude_terms]

if (!setequal(curated_term_order, terms_show))
  stop("curated_term_order does not match the term set exactly (missing or extra terms)")

pdat$Cluster     = factor(pdat$Cluster, levels = clusters)
pdat$Description = factor(pdat$Description, levels = rev(curated_term_order))

p = ggplot(pdat, aes(x = Cluster, y = Description,
                     size = Count, colour = -log10(p.adjust))) +
  geom_point() +
  scale_x_discrete(limits = clusters, drop = FALSE, expand = expansion(add = 0.55)) +
  scale_size_continuous(name = "Gene count", range = c(1, 6)) +
  # grey -> black: magnitude-only (all markers are only.pos = TRUE, no "down"
  # direction), kept off the blue/orange hue family used for the signed Green24
  # z-scores elsewhere in this figure group - see LD_X04_B_characterisation_plots.R
  scale_colour_gradient(name = "-log10(adj. p)", low = "grey80", high = "black") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 13)

# vertical dashed lines: cluster-family boundaries (SLC1A2/GFAP/CHI3L1)
subtype    = sub("_s[0-9]+$", "", clusters)
v_boundaries = which(diff(match(subtype, unique(subtype))) != 0) + 0.5
p = p + geom_vline(xintercept = v_boundaries, linetype = "dashed",
                   colour = "grey40", linewidth = 0.4)

# horizontal dashed line: separates the CHI3L1-only term block from the rest
h_boundary = length(curated_term_order) - term_hline_after + 0.5
p = p + geom_hline(yintercept = h_boundary, linetype = "dashed",
                   colour = "grey40", linewidth = 0.4)

p = p +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        plot.margin = margin(t = 10, r = 8, b = 4, l = 6))

h = max(6, 0.30 * length(terms_show) + 2.4)
w = max(6, 0.52 * length(clusters) + 3)
ggsave(file.path(out_dir, paste0(script_ind, "GO_top5_all_sig_subclusters_curated.pdf")),
       p, width = w, height = h, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "GO_top5_all_sig_subclusters_curated.png")),
       p, width = w, height = h, dpi = 300)

message("Done. Curated GO plot rewritten (", length(clusters), " clusters, ", length(terms_show), " terms).")
