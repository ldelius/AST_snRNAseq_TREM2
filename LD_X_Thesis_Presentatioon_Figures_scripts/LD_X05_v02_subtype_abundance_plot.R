# LD_X05_v02: Descriptive subtype abundance plot from LD_X05 outputs.
# A shared y-axis preserves differences in abundance among subtypes.
# No brackets are drawn because no subtype-level contrasts were significant.

library(tidyverse)

base       = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
csv_path   = file.path(base, "LD_X_Thesis_Presentation_output/LD_X05_abundance_by_subtype.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
script_ind = "LD_X05_v02_"
if (!file.exists(csv_path)) stop("Missing input: ", csv_path)

ct = read_csv(csv_path, show_col_types = FALSE)

gr             = c("Control_CV", "Control_R47H", "Control_R62H", "AD_CV", "AD_R47H", "AD_R62H")
subtype_levels = c("AST_SLC1A2", "AST_GFAP", "AST_CHI3L1")   # same order as the original figure
group_labels   = set_names(gsub("_", " ", gr), gr)           # display only - underlying factor levels unchanged

# fill = family identity (same hue as its facet's colour in the X02 UMAP:
# SLC1A2 = orange, GFAP = blue, CHI3L1 = teal) x diagnosis (light = Control,
# dark = AD). Genotype (CV/R47H/R62H) is distinguished by x-axis position/label
# only, not colour. AD colours are darker than the plain UMAP hex (rather than
# reusing it as-is) for stronger visual contrast with Control - GFAP's base hue
# is already dark enough that it's kept as-is.
ct = ct %>%
  dplyr::mutate(cellgroup = factor(cellgroup, levels = subtype_levels),
                group     = factor(group,     levels = gr),
                diagnosis = ifelse(grepl("^Control", group), "Control", "AD"),
                fill_grp  = paste(as.character(cellgroup), diagnosis))

fill_cols = c(
  "AST_SLC1A2 Control" = "#F4D48C", "AST_SLC1A2 AD" = "#BC7C00",
  "AST_GFAP Control"   = "#8CC0DC", "AST_GFAP AD"   = "#0072B2",
  "AST_CHI3L1 Control" = "#8CD3C0", "AST_CHI3L1 AD" = "#008960")

p = ggplot(ct, aes(group, fract_sample)) +
  geom_boxplot(aes(fill = fill_grp), outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.6) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_x_discrete(labels = group_labels) +
  facet_wrap(~ cellgroup, ncol = 3) +   # no `scales =` -> shared/fixed y-axis across facets
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = NULL, y = "Fraction of astrocyte nuclei (per sample)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 10.5),
        axis.text.y     = element_text(size = 10.5),
        axis.title      = element_text(size = 12),
        strip.background = element_blank(),
        strip.text      = element_text(size = 13, face = "bold"))

ggsave(file.path(out_dir, paste0(script_ind, "sig_subtype.pdf")), p, width = 9, height = 4, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "sig_subtype.png")), p, width = 9, height = 4, dpi = 300)

### median subtype-fraction stats, for stating in the thesis text -----------
# overall (pooled across all groups/samples) - the headline number, e.g.
# "SLC1A2 astrocytes comprised the majority of the population (median ~70%)"
median_overall = ct %>%
  dplyr::group_by(cellgroup) %>%
  dplyr::summarise(median_pct = 100 * median(fract_sample), n_samples = dplyr::n(), .groups = "drop")

# Group-specific medians.
median_by_group = ct %>%
  dplyr::group_by(cellgroup, group) %>%
  dplyr::summarise(median_pct = 100 * median(fract_sample), n_samples = dplyr::n(), .groups = "drop")

write_csv(median_overall,  file.path(out_dir, paste0(script_ind, "subtype_median_fraction_overall.csv")))
write_csv(median_by_group, file.path(out_dir, paste0(script_ind, "subtype_median_fraction_by_group.csv")))

message("Median subtype fraction (overall, pooled across groups):")
print(as.data.frame(median_overall))
message("Median subtype fraction, by group:")
print(as.data.frame(median_by_group))

message("Done. Written to ", out_dir)
