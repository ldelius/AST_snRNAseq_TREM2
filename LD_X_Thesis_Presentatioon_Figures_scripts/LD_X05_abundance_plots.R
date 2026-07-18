# LD_X05: Astrocyte abundance (fraction of sample, by group) - all abundance plots.
#   For each grouping (all 18 subclusters / two-largest-per-subtype / subtypes /
#   subtypes-by-region) this makes:
#     (A) a descriptive dodged-bar overview (mean +- SD + per-sample points), and
#     (B) a bracket boxplot annotated with sccomp differential-abundance contrasts.
#
# Statistics: sccomp_glm(~ 0 + group, contrasts = ...) - Bayesian composition
#   model (sccomp 1.6, as B04). Contrasts:
#     - every group vs Control_CV
#     - AD risk variants vs AD_CV
#     - AD risk variants vs each other
#   Only SIGNIFICANT contrasts (FDR < 0.05) are drawn as brackets (labelled by
#   FDR); the full contrast results are written to CSV ("excel list").
#
# DATA: LD_B04a_v02_seur.qs (metadata only). NB: sccomp fits a Stan model per
#   grouping -> heavy/slow (~30-60 min, >=64 GB, 8 cores).

library(tidyverse)
library(qs)
library(Seurat)
library(sccomp)

### paths -------------------------------------------------------------------
base       = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
b04_path   = file.path(base, "LD_B_AST_analysis_output/LD_B04a_v02_seur.qs")
clust_csv  = file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
group_csv  = file.path(base, "data_TREM2_michael/A_input/group_tab.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X05_"
cores_n    = 8
ref_group  = "Control_CV"
fdr_cut    = 0.05

ord = read_csv(clust_csv, show_col_types = FALSE)
cluster_names  = unique(ord$cluster_name)
subtype_levels = unique(ord$cell_type)
grtab = read_csv(group_csv, show_col_types = FALSE)
gr = unique(grtab$group)
group_cols = set_names(scales::hue_pal()(length(gr)), gr)

comparisons = tribble(
  ~name,                  ~g1,            ~g2,
  "Ctrl_R47H vs Ctrl_CV", "Control_R47H", "Control_CV",
  "Ctrl_R62H vs Ctrl_CV", "Control_R62H", "Control_CV",
  "AD_CV vs Ctrl_CV",     "AD_CV",        "Control_CV",
  "AD_R47H vs Ctrl_CV",   "AD_R47H",      "Control_CV",
  "AD_R62H vs Ctrl_CV",   "AD_R62H",      "Control_CV",
  "AD_R47H vs AD_CV",     "AD_R47H",      "AD_CV",
  "AD_R62H vs AD_CV",     "AD_R62H",      "AD_CV",
  "AD_R62H vs AD_R47H",   "AD_R62H",      "AD_R47H")

# contrasts use cell-means coding (~ 0 + group): every pairwise comparison is a
# difference of group coefficients, e.g. "groupAD_R62H - groupAD_CV".

### load metadata -----------------------------------------------------------
message("Loading B04 object (metadata only)...")
seur = qread(b04_path)
meta = seur@meta.data
if (!"group" %in% names(meta))       meta$group       = grtab$group[match(meta$sample, grtab$sample)]
if (!"BrainRegion" %in% names(meta)) meta$BrainRegion = grtab$BrainRegion[match(meta$sample, grtab$sample)]
rm(seur); gc()
# sample-level covariates (from the cohort table) for the adjusted model
covars = c("Sex", "BrainRegion", "APOEgroup", "CD33Group", "cohort")
samp_meta = grtab %>%
  dplyr::distinct(sample, group, Sex, BrainRegion, APOEgroup, CD33Group, cohort)

### per-sample counts + fractions (zeros included) --------------------------
counts_tab = function(meta, unit_col, unit_levels, samp_meta) {
  u = factor(meta[[unit_col]], levels = unit_levels); s = factor(meta$sample)
  tab = as.data.frame(table(sample = s, cellgroup = u), stringsAsFactors = FALSE)
  names(tab)[names(tab) == "Freq"] = "N_cells"
  tab %>% dplyr::left_join(samp_meta, by = "sample") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(N_sample = sum(N_cells),
                  fract_sample = ifelse(N_sample > 0, N_cells / N_sample, 0)) %>%
    dplyr::ungroup()
}

### sccomp contrasts -> tidy (cellgroup, comparison, c_FDR) ------------------
run_contrasts = function(ct, covars = character(0)) {
  present = unique(ct$group)
  ct$group = droplevels(factor(ct$group, levels = gr))
  comps = comparisons %>% dplyr::filter(g1 %in% present, g2 %in% present)
  if (nrow(comps) == 0) return(NULL)
  cstr = setNames(paste0("group", comps$g1, " - group", comps$g2), comps$name)
  # keep only covariates that exist and vary here (e.g. BrainRegion is constant
  # within a single-region model and must be dropped); drop NA-covariate samples
  covars = covars[vapply(covars, function(v) v %in% names(ct) && dplyr::n_distinct(ct[[v]]) > 1, logical(1))]
  if (length(covars)) ct = ct[stats::complete.cases(ct[, covars, drop = FALSE]), ]
  form = stats::as.formula(paste("~ 0 + group",
                                 if (length(covars)) paste("+", paste(covars, collapse = " + ")) else ""))
  res = ct |>
    sccomp_glm(formula_composition = form,
               .sample = sample, .cell_group = cellgroup, .count = N_cells,
               contrasts = cstr,
               bimodal_mean_variability_association = TRUE, cores = cores_n)
  r = as.data.frame(res)
  if (!"factor" %in% names(r) && "parameter" %in% names(r)) r$factor = r$parameter
  par_col = intersect(c("parameter", "factor"), names(r))[1]
  cg_col  = intersect(c("cellgroup", "cell_group"), names(r))[1]
  r %>% dplyr::transmute(cellgroup = .data[[cg_col]],
                         comparison = .data[[par_col]], c_FDR = c_FDR) %>%
    dplyr::mutate(comparison = ifelse(comparison %in% comps$name, comparison,
                                      names(cstr)[match(comparison, cstr)]))
}

### (A) descriptive dodged-bar overview -------------------------------------
plot_dodge = function(ct, unit_order, title, by_region = FALSE) {
  ct = ct %>% dplyr::mutate(unit = factor(cellgroup, levels = unit_order), group = factor(group, levels = gr))
  gvars = if (by_region) c("unit", "group", "BrainRegion") else c("unit", "group")
  summ = ct %>% dplyr::group_by(dplyr::across(dplyr::all_of(gvars))) %>%
    dplyr::summarise(mean_fract = mean(fract_sample), sd_fract = sd(fract_sample), .groups = "drop")
  p = ggplot() +
    geom_col(data = summ, aes(unit, mean_fract, color = group, group = group),
             fill = "grey90", position = position_dodge(0.7), width = 0.7, linewidth = 0.3) +
    geom_errorbar(data = summ, aes(unit, ymin = mean_fract - sd_fract, ymax = mean_fract + sd_fract,
                                   color = group, group = group),
                  position = position_dodge(0.7), width = 0.3, linewidth = 0.2) +
    geom_point(data = ct, aes(unit, fract_sample, color = group, group = group),
               position = position_dodge(0.7), size = 0.5, stroke = 0.3) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = group_cols, name = "Group") +
    labs(title = title, x = NULL, y = "Fraction of sample",
         caption = "descriptive composition (mean ± SD, points = samples; no statistical model)") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
  if (by_region) p = p + facet_wrap(~ BrainRegion)
  p
}

### (B) bracket boxplot (significant contrasts only) ------------------------
bracket_boxplot = function(ct, unit_order, sig, title, ncol = NULL, by_region = FALSE, model = "") {
  ct = ct %>% dplyr::mutate(cellgroup = factor(cellgroup, levels = unit_order), group = factor(group, levels = gr))
  gkeys = if (by_region) c("cellgroup", "BrainRegion") else "cellgroup"
  ymax = ct %>% dplyr::group_by(dplyr::across(dplyr::all_of(gkeys))) %>%
    dplyr::summarise(ymax = max(fract_sample, na.rm = TRUE), .groups = "drop")
  brk = NULL
  if (!is.null(sig) && nrow(sig) > 0) {
    brk = sig %>% dplyr::filter(!is.na(c_FDR), c_FDR < fdr_cut) %>%   # significant only
      dplyr::left_join(comparisons, by = c("comparison" = "name")) %>%
      dplyr::filter(!is.na(g1), !is.na(g2)) %>%
      dplyr::mutate(cellgroup = factor(cellgroup, levels = unit_order),
                    x = match(g1, gr), xend = match(g2, gr),
                    label = sprintf("%.2g", c_FDR)) %>%
      dplyr::filter(!is.na(cellgroup)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(gkeys))) %>%
      dplyr::arrange(abs(xend - x), .by_group = TRUE) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>% dplyr::ungroup() %>%
      dplyr::left_join(ymax, by = gkeys) %>%
      dplyr::mutate(y = ymax * (1 + 0.12 * rank))
  }
  p = ggplot(ct, aes(group, fract_sample)) +
    geom_boxplot(aes(fill = group), outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = group_cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
    labs(title = title, x = NULL, y = "Fraction of sample",
         caption = paste0("model: ", model, "   |   brackets: significant sccomp contrasts (FDR < 0.05)")) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.caption = element_text(size = 7)) +
    (if (by_region) facet_grid(cellgroup ~ BrainRegion, scales = "free_y")
     else facet_wrap(~ cellgroup, scales = "free_y", ncol = ncol))
  if (!is.null(brk) && nrow(brk) > 0)
    p = p +
      geom_segment(data = brk, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.25) +
      geom_text(data = brk, aes(x = (x + xend) / 2, y = y, label = label), inherit.aes = FALSE, vjust = -0.2, size = 2)
  p
}

### build tables ------------------------------------------------------------
ct_cluster = counts_tab(meta, "cluster_name", cluster_names,  samp_meta)
ct_subtype = counts_tab(meta, "cell_type",    subtype_levels, samp_meta)
write_csv(ct_cluster, file.path(out_dir, paste0(script_ind, "abundance_by_subcluster.csv")))
write_csv(ct_subtype, file.path(out_dir, paste0(script_ind, "abundance_by_subtype.csv")))

sel = meta %>% dplyr::count(cell_type, cluster_name, name = "n") %>%
  dplyr::group_by(cell_type) %>% dplyr::slice_max(n, n = 2, with_ties = FALSE) %>% dplyr::ungroup()
sel_clusters = sel$cluster_name[order(match(sel$cluster_name, cluster_names))]

### sccomp (unadjusted + covariate-adjusted) --------------------------------
safe = function(expr) tryCatch(expr, error = function(e) { message("sccomp failed: ", conditionMessage(e)); NULL })
run_region = function(ct, cv = character(0)) dplyr::bind_rows(lapply(c("MTG", "SSC"), function(rg) {
  r = safe(run_contrasts(dplyr::filter(ct, BrainRegion == rg), cv))
  if (!is.null(r)) r$BrainRegion = rg
  r
}))

message("sccomp subclusters (unadjusted)...");        sig_cluster_un = safe(run_contrasts(ct_cluster))
message("sccomp subclusters (adjusted)...");          sig_cluster_aj = safe(run_contrasts(ct_cluster, covars))
message("sccomp subtypes (unadjusted)...");           sig_subtype_un = safe(run_contrasts(ct_subtype))
message("sccomp subtypes (adjusted)...");             sig_subtype_aj = safe(run_contrasts(ct_subtype, covars))
message("sccomp subtypes by region (unadjusted)..."); sig_region_un  = run_region(ct_subtype)
message("sccomp subtypes by region (adjusted)...");   sig_region_aj  = run_region(ct_subtype, covars)

wr = function(d, n) if (!is.null(d) && nrow(d) > 0) write_csv(d, file.path(out_dir, paste0(script_ind, n, ".csv")))
wr(sig_cluster_un, "sccomp_subclusters_unadj"); wr(sig_cluster_aj, "sccomp_subclusters_adj")
wr(sig_subtype_un, "sccomp_subtypes_unadj");    wr(sig_subtype_aj, "sccomp_subtypes_adj")
wr(sig_region_un,  "sccomp_subtypes_by_region_unadj"); wr(sig_region_aj, "sccomp_subtypes_by_region_adj")

### figures -----------------------------------------------------------------
save_plot = function(p, suffix, w, h) {
  ggsave(file.path(out_dir, paste0(script_ind, suffix, ".pdf")), p, width = w, height = h, useDingbats = FALSE)
  ggsave(file.path(out_dir, paste0(script_ind, suffix, ".png")), p, width = w, height = h, dpi = 300)
}
# NULL-safe subset of a (possibly failed -> NULL) sccomp result table
sig_filt = function(sig, keep) if (is.null(sig)) NULL else dplyr::filter(sig, cellgroup %in% keep)

# (A) descriptive overviews
save_plot(plot_dodge(ct_cluster, cluster_names, "Subcluster abundance (all)"), "overview_all_subclusters", 12, 4)
save_plot(plot_dodge(dplyr::filter(ct_cluster, cellgroup %in% sel_clusters), sel_clusters,
                     "Subcluster abundance (two largest per subtype)"), "overview_two_largest", 7, 4)
save_plot(plot_dodge(ct_subtype, subtype_levels, "Subtype abundance (pooled)"), "overview_subtype", 5, 4)
save_plot(plot_dodge(ct_subtype, subtype_levels, "Subtype abundance by region", by_region = TRUE), "overview_subtype_by_region", 8, 4)

# (B) bracket boxplots with significant sccomp contrasts: unadjusted vs adjusted
# model formula strings for the plot captions
f_un     = "~ 0 + group"
f_aj     = paste("~ 0 + group +", paste(covars, collapse = " + "))
f_aj_reg = paste("~ 0 + group +", paste(setdiff(covars, "BrainRegion"), collapse = " + "))  # region constant within facet
ct_two = dplyr::filter(ct_cluster, cellgroup %in% sel_clusters)
# subtypes
save_plot(bracket_boxplot(ct_subtype, subtype_levels, sig_subtype_un, "Subtype abundance (unadjusted)", ncol = 3, model = f_un), "sig_subtype_unadj", 9, 4)
save_plot(bracket_boxplot(ct_subtype, subtype_levels, sig_subtype_aj, "Subtype abundance (covariate-adjusted)", ncol = 3, model = f_aj), "sig_subtype_adj", 9, 4)
# two largest per subtype
save_plot(bracket_boxplot(ct_two, sel_clusters, sig_filt(sig_cluster_un, sel_clusters), "Two largest per subtype (unadjusted)", ncol = 3, model = f_un), "sig_two_largest_unadj", 10, 7)
save_plot(bracket_boxplot(ct_two, sel_clusters, sig_filt(sig_cluster_aj, sel_clusters), "Two largest per subtype (covariate-adjusted)", ncol = 3, model = f_aj), "sig_two_largest_adj", 10, 7)
# all subclusters
save_plot(bracket_boxplot(ct_cluster, cluster_names, sig_cluster_un, "All subclusters (unadjusted)", ncol = 6, model = f_un), "sig_all_subclusters_unadj", 16, 12)
save_plot(bracket_boxplot(ct_cluster, cluster_names, sig_cluster_aj, "All subclusters (covariate-adjusted)", ncol = 6, model = f_aj), "sig_all_subclusters_adj", 16, 12)
# subtypes by region
save_plot(bracket_boxplot(ct_subtype, subtype_levels, sig_region_un, "Subtype by region (unadjusted)", by_region = TRUE, model = f_un), "sig_subtype_by_region_unadj", 8, 7)
save_plot(bracket_boxplot(ct_subtype, subtype_levels, sig_region_aj, "Subtype by region (covariate-adjusted)", by_region = TRUE, model = f_aj_reg), "sig_subtype_by_region_adj", 8, 7)

message("Done. Abundance overviews, sccomp bracket boxplots and contrast tables written to ", out_dir)
