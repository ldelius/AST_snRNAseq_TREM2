# LD_X22: supplementary figure - GO biological processes per log2FC quadrant,
# ALL astrocytes pooled (not per family).
#
# DATA: LD_X08_GO_quadrants_family_pooled_results.csv, written by
# LD_X08_GO_quadrants_family_pooled_5covar.R. That script computes four scopes -
# ALL_pooled, SLC1A2, GFAP, CHI3L1 - but only plots the three family ones. This
# takes the ALL_pooled scope, which is the true pooled model from LD_E04c
# (single DESeq2 fit across all subclusters with cluster_name as a covariate,
# 5-covariate base), i.e. the same scope as the LD_X07 concordance panel.
#
# Nothing is recomputed: no enrichGO, no DESeq2, no .qs. Reads one CSV.
#
# COLOURS: grey -> black for -log10(adj. p), matching the curated GO figure
# (LD_X04_v02 / LD_X15). GO over-representation here is magnitude-only - the
# direction is already encoded by the quadrant - so an achromatic ramp is used
# and the blue/orange hue family is left for genuinely signed quantities. This is
# the only change from LD_X08's own styling.
#
# All four quadrants are always drawn. up_up has genes in every pair (362 / 479 /
# 1044) but yields no enriched terms, so it appears as an empty column: that is a
# result, not a missing category.

library(tidyverse)

### paths -------------------------------------------------------------------
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")
out_dir  = file.path(base, "LD_X_Thesis_Presentation_output")
go_csv   = file.path(out_dir, "LD_X08_GO_quadrants_family_pooled_results.csv")
script_ind = "LD_X22_"
if (!file.exists(go_csv)) stop("Missing input: ", go_csv)
message("Using base: ", base)

### knobs --------------------------------------------------------------------
TOP_N   = 10     # GO terms per quadrant per pair
ROW_IN  = 0.30   # inches per term row (LAYOUT "row"); per term column if "rotated"
TERM_PT = 11     # GO term label size

# "row"      ORIGINAL formatting, as LD_X08 / LD_E04b: quadrant on x, GO terms on
#            y, the three contrast pairs side by side in one row. Landscape.
# "rotated"  the same data turned 90 degrees - GO terms along x, quadrant on y -
#            with the three pairs stacked underneath each other. Portrait.
LAYOUT = "row"

### labels -------------------------------------------------------------------
# short labels: with 11 pt term text the panels are narrow, and long quadrant
# labels overlap while long strip titles get truncated
quad_label = c(up_up     = "Shared\nup",
               down_down = "Shared\ndown",
               up_down   = "Up |\ndown",
               down_up   = "Down |\nup")
pair_label = c(P1_R62H_vs_CV__vs__AD_vs_Control = "R62H vs CV | AD vs Ctrl",
               P2_R47H_vs_CV__vs__AD_vs_Control = "R47H vs CV | AD vs Ctrl",
               P3_R62H_vs_CV__vs__R47H_vs_CV    = "R62H vs CV | R47H vs CV")

### data ---------------------------------------------------------------------
go = read_csv(go_csv, show_col_types = FALSE) %>% filter(scope == "ALL_pooled")
if (nrow(go) == 0) stop("No ALL_pooled rows in ", basename(go_csv))

t = go %>%
  group_by(pair, quadrant) %>%
  slice_min(p.adjust, n = TOP_N, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(pair     = factor(pair, levels = names(pair_label), labels = pair_label),
         quadrant = factor(quadrant, levels = names(quad_label), labels = quad_label))

# term order: within each pair, group by the quadrant the term is most
# significant in, then by significance - so each panel reads as blocks rather
# than an arbitrary alphabetical list
term_ord = t %>% group_by(pair, Description) %>%
  slice_min(p.adjust, n = 1, with_ties = FALSE) %>% ungroup() %>%
  arrange(pair, quadrant, desc(p.adjust)) %>%
  distinct(pair, Description)
t = t %>% mutate(Description = factor(Description, levels = unique(term_ord$Description)))

### plot ---------------------------------------------------------------------
# shared encodings: size = gene count, colour = -log10(adj. p) on the thesis
# achromatic ramp; only the axes and facet direction differ between layouts
common = list(
  scale_colour_gradient(name = "-log10(adj. p)", low = "grey80", high = "black"),
  scale_size_continuous(name = "Gene count", range = c(1, 6)),
  labs(x = NULL, y = NULL),
  theme_bw(base_size = 11),
  theme(panel.grid       = element_line(linewidth = 0.2, colour = "grey92"),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text       = element_text(size = 10, face = "bold"),
        legend.key.size  = unit(4, "mm"),
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 7)))

n_terms = t %>% distinct(pair, Description) %>% nrow()
n_pairs = nlevels(t$pair)

if (LAYOUT == "row") {
  # ORIGINAL formatting: quadrant on x, terms on y, pairs side by side
  p = ggplot(t, aes(x = quadrant, y = Description, size = Count, colour = -log10(p.adjust))) +
    geom_point() +
    facet_wrap(~ pair, nrow = 1, scales = "free_y") +
    scale_x_discrete(drop = FALSE) +   # keep empty quadrants (up_up) as columns
    common +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = TERM_PT))
  # height from the busiest panel; width fixed - three label columns need room
  H = max(6, ROW_IN * max(table(t$pair[!duplicated(paste(t$pair, t$Description))])) + 2.0)
  W = 20
} else if (LAYOUT == "rotated") {
  # turned 90 degrees: terms along x, quadrant on y, pairs stacked
  p = ggplot(t, aes(x = Description, y = quadrant, size = Count, colour = -log10(p.adjust))) +
    geom_point() +
    facet_wrap(~ pair, ncol = 1, scales = "free_x") +
    scale_y_discrete(drop = FALSE) +
    common +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = TERM_PT - 2),
          axis.text.y = element_text(size = 10))
  W = max(8.3, 0.16 * n_terms + 3)
  H = 3.4 * n_pairs
} else stop("LAYOUT must be 'row' or 'rotated'")

### save ---------------------------------------------------------------------
suffix = paste0("top", TOP_N, "_", LAYOUT)
ggsave(file.path(out_dir, paste0(script_ind, "GO_quadrants_pooled_", suffix, ".pdf")),
       p, width = W, height = H, limitsize = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "GO_quadrants_pooled_", suffix, ".png")),
       p, width = W, height = H, dpi = 300, limitsize = FALSE)

write_csv(t %>% select(pair, quadrant, Description, Count, p.adjust, geneID),
          file.path(out_dir, paste0(script_ind, "GO_quadrants_pooled_", suffix, "_terms.csv")))

message("Done. ", n_terms, " term entries across ", n_pairs, " pairs, layout ", LAYOUT,
        ", ", round(W, 1), " x ", round(H, 1), " in. Outputs in: ", out_dir)
