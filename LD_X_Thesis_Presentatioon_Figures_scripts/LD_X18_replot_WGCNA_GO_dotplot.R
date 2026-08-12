# LD_X18: Replot WGCNA GO enrichment with readable row spacing.

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

### settings -----------------------------------------------------------------
TOP_N  = 10      # terms per module
ROW_IN = 0.22    # inches per row
LAB_PT = 9       # y-axis label size

### data ---------------------------------------------------------------------
# The GO table is ordered by significance within module.
go = read_csv(go_csv, show_col_types = FALSE) %>% filter(Count > 1)

# Preserve input row order by looping over modules in the
# order they appear in the GO table (NOT alphabetically - dplyr::group_by would
# give M0, M1, M11, M12, ... M2, which reshuffles the y axis) and takes the first
# TOP_N rows of each.
t2 = do.call(rbind, lapply(unique(go$module), function(m) {
  d = go[go$module == m, ]
  if (nrow(d) > TOP_N) d[seq_len(TOP_N), ] else d
}))

mods = unique(t2$module)   # input order

### module colours -----------------------------------------------------------
# Use viridis for categorical module identity, separate from the blue/orange
# scale used for signed quantities.
# M0 stays grey: in WGCNA, module 0 is the unassigned ("grey") set, not a real
# module, so it should read as background rather than as one colour among many.
mods_real = setdiff(mods, "M0")
pal = set_names(scales::viridis_pal(option = "viridis", begin = 0.05, end = 0.92)(length(mods_real)),
                mods_real)
if ("M0" %in% mods) pal = c(M0 = "grey70", pal)
pal = pal[mods]

# Keep the first input term at the bottom of the discrete axis.
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

message("Done. Outputs in: ", out_dir)
