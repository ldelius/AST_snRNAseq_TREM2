# LD_X15: Replot uncurated GO results using the curated figure's colour scheme.

library(tidyverse)

### paths -------------------------------------------------------------------
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")

go_csv     = file.path(base, "LD_X_Thesis_Presentation_output/LD_X04_v02_GO_top_terms_table.csv")
clust_csv  = file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
script_ind = "LD_X15_"
for (p in c(go_csv, clust_csv)) if (!file.exists(p)) stop("Missing input: ", p)
message("Using base: ", base)

### data ---------------------------------------------------------------------
res = read_csv(go_csv, show_col_types = FALSE)
ord = read_csv(clust_csv, show_col_types = FALSE)
cluster_names = unique(ord$cluster_name)
sig_clusters  = cluster_names[cluster_names %in% unique(as.character(res$Cluster))]

n_terms = 5

# top n terms per subcluster by adjusted p, then show every significant
# (subcluster, term) pair among those terms - identical selection to the original
topn = res %>% dplyr::filter(Cluster %in% sig_clusters) %>%
  dplyr::group_by(Cluster) %>%
  dplyr::slice_min(p.adjust, n = n_terms, with_ties = FALSE) %>% dplyr::ungroup()
terms_show = unique(topn$Description)

pdat = res %>% dplyr::filter(Cluster %in% sig_clusters, Description %in% terms_show)
clusters = sig_clusters[sig_clusters %in% unique(as.character(pdat$Cluster))]
pdat$Cluster = factor(pdat$Cluster, levels = clusters)

# y-axis order (same rule as the un-curated original): group each term by the
# subtype of its most significant subcluster, SLC1A2 -> GFAP -> CHI3L1, so CHI3L1
# terms end up at the TOP of the axis
term_ord = topn %>% dplyr::group_by(Description) %>%
  dplyr::slice_min(p.adjust, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>%
  dplyr::mutate(subtype = sub("_s[0-9]+$", "", Cluster),
                srank   = dplyr::recode(subtype, AST_SLC1A2 = 1, AST_GFAP = 2, AST_CHI3L1 = 3)) %>%
  dplyr::arrange(srank, p.adjust)
pdat$Description = factor(pdat$Description, levels = term_ord$Description)

### plot ---------------------------------------------------------------------
p = ggplot(pdat, aes(x = Cluster, y = Description,
                     size = Count, colour = -log10(p.adjust))) +
  geom_point() +
  scale_x_discrete(limits = clusters, drop = FALSE, expand = expansion(add = 0.55)) +
  scale_size_continuous(name = "Gene count", range = c(1, 6)) +
  # THE CHANGE: grey -> black, matching the curated main-results figure. Magnitude
  # only (all markers are only.pos = TRUE, so there is no "down" direction), and
  # kept off the blue/orange hue family used for the signed Green24 z-scores.
  scale_colour_gradient(name = "-log10(adj. p)", low = "grey80", high = "black") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 13)

# vertical dashed lines at the cluster-family boundaries (SLC1A2/GFAP/CHI3L1),
# as in the curated figure. No horizontal separator: that marked the curated
# CHI3L1-only block, which does not exist in the un-curated term set.
subtype      = sub("_s[0-9]+$", "", clusters)
v_boundaries = which(diff(match(subtype, unique(subtype))) != 0) + 0.5
p = p + geom_vline(xintercept = v_boundaries, linetype = "dashed",
                   colour = "grey40", linewidth = 0.4)

p = p +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        plot.margin = margin(t = 10, r = 8, b = 4, l = 6))

### save ---------------------------------------------------------------------
h = max(6, 0.30 * nlevels(pdat$Description) + 2.4)
w = max(6, 0.52 * length(clusters) + 3)
ggsave(file.path(out_dir, paste0(script_ind, "GO_top5_all_sig_subclusters.pdf")),
       p, width = w, height = h, useDingbats = FALSE, limitsize = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "GO_top5_all_sig_subclusters.png")),
       p, width = w, height = h, dpi = 300, limitsize = FALSE)

message("Done. Un-curated GO plot: ", length(clusters), " subclusters, ",
        nlevels(pdat$Description), " terms (", round(w, 1), " x ", round(h, 1), " in). ",
        "Outputs in: ", out_dir)
