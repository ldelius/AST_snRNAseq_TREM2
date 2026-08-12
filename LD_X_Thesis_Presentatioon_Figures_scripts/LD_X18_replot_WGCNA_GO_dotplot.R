# LD_X18: re-draw the LD_F03c_v02 WGCNA "top GO terms per module" dot plot with a
# taller canvas and larger row labels.
#
# NOTHING IS RECOMPUTED. WGCNA and the GO enrichment already ran in
# LD_F03c_v02_WGCNA_AD_both_variants_7covar_corr_DEGseed.R; this reads the GO
# result table it wrote and re-plots it.
#
# TWO KNOBS:
#   TOP_N      terms per module (10 = as the original; 5 = fits an A4 page)
#   ROW_IN     inches of height per y-axis row; label size scales with it
# The original was 11 x 12 in for 82 terms, i.e. 0.15 in per row, which is why
# the row labels had to be small.
#
# MODULE COLOURS: the original coloured the dots by raw WGCNA colour name
# (turquoise/blue/brown/...), which lives only inside the 856 MB bulk_data.qs and
# does not fit the thesis palette. Replaced with viridis (thesis convention for
# categorical identity), so this script needs neither qs nor the big object and
# runs on a laptop in seconds. See the module-colours section below.

library(tidyverse)

### paths -------------------------------------------------------------------
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")

f_dir      = file.path(base, "LD_F_DESeq_pseudobulk_WGCNA/LD_F03c_v02")
go_csv     = file.path(f_dir, "LD_F03c_v02_GO_results_by_comp.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X18_"
if (!file.exists(go_csv)) stop("Missing input: ", go_csv)

### knobs --------------------------------------------------------------------
TOP_N  = 10      # terms per module (original: 10)
ROW_IN = 0.22    # inches per row (original: ~0.15)
LAB_PT = 9       # y-axis label size (original: theme default, ~8.8 at base 11)

### data ---------------------------------------------------------------------
# same selection as the original: Count > 1, then the first TOP_N rows per module
# (the GO table is already ordered by significance within module)
go = read_csv(go_csv, show_col_types = FALSE) %>% filter(Count > 1)

# Replicate the original's row order exactly: it loops over the modules in the
# order they appear in the GO table (NOT alphabetically - dplyr::group_by would
# give M0, M1, M11, M12, ... M2, which reshuffles the y axis) and takes the first
# TOP_N rows of each.
t2 = do.call(rbind, lapply(unique(go$module), function(m) {
  d = go[go$module == m, ]
  if (nrow(d) > TOP_N) d[seq_len(TOP_N), ] else d
}))

mods = unique(t2$module)   # x-axis order: as in the GO table, as the original

### module colours -----------------------------------------------------------
# The original used the raw WGCNA colour names (turquoise/blue/brown/...), which
# are arbitrary labels and clash with the thesis palette. Replaced here with the
# thesis convention for categorical identity: viridis beyond 8 categories, as in
# LD_X02 (cluster identity) and LD_X10b (Hallmark groups). Viridis is also kept
# clear of the blue/orange used for signed quantities elsewhere.
# M0 stays grey: in WGCNA, module 0 is the unassigned ("grey") set, not a real
# module, so it should read as background rather than as one colour among many.
mods_real = setdiff(mods, "M0")
pal = set_names(scales::viridis_pal(option = "viridis", begin = 0.05, end = 0.92)(length(mods_real)),
                mods_real)
if ("M0" %in% mods) pal = c(M0 = "grey70", pal)
pal = pal[mods]

# NOT reversed: ggplot puts the first level of a discrete scale at the BOTTOM, and
# the original passed limits = unique(t2$Description) - so M0's terms sit at the
# bottom of the axis and the last module's at the top.
t2 = t2 %>% mutate(module = factor(module, levels = mods),
                   Description = factor(Description, levels = unique(Description)))

### plot ---------------------------------------------------------------------
p = ggplot(t2, aes(x = module, y = Description, size = Count, colour = module)) +
  geom_point() +
  scale_colour_manual(limits = mods, values = pal, guide = "none") +
  scale_size_continuous(limits = c(0, max(t2$Count))) +
  scale_x_discrete(limits = mods) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = LAB_PT),
        panel.grid  = element_line(linewidth = 0.2, colour = "grey92")) +
  labs(x = NULL, y = NULL)

### save ---------------------------------------------------------------------
n_terms = nlevels(t2$Description)
H = max(6, ROW_IN * n_terms + 1.6)
W = 11
suffix = paste0("top", TOP_N)
ggsave(file.path(out_dir, paste0(script_ind, "WGCNA_GO_dotplot_", suffix, ".pdf")),
       p, width = W, height = H, limitsize = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "WGCNA_GO_dotplot_", suffix, ".png")),
       p, width = W, height = H, dpi = 300, limitsize = FALSE)

message("Done. ", n_terms, " terms x ", length(mods), " modules, ",
        W, " x ", round(H, 1), " in (", ROW_IN, " in/row, labels ", LAB_PT, " pt).",
        if (H > 11.7) "  NB: taller than an A4 page." else "",
        " Outputs in: ", out_dir)
