# LD_X06: Astrocyte clusters in relation to pathology.
#   PART 1  Distribution of pathology measures across donor groups, faceted by
#           brain region, with a covariate-adjusted linear model + contrast tests.
#   PART 2  AD-only test of whether subcluster ABUNDANCE changes with pathology
#           burden and whether the relationship DIFFERS BY TREM2 VARIANT
#           (sccomp composition model with a pathology x variant interaction).
#
# Pathology measures: amyloid plaque density (TotalDensity), and tau % positive
#   area for AT8 (pSer202/Thr205) and PHF1 (pSer396/404).
#
# CAVEATS (built into comments below):
#   - R47H has ~8 AD donors -> its slopes / R47H-CV contrasts are exploratory.
#   - TREM2 variant and pathology are correlated (R62H carries higher burden), so
#     variants span different pathology ranges; slope comparisons assume a common
#     linear relationship (some extrapolation).
#   - Pathology is per-sample (donor x region); donors contribute 2 samples, and
#     sccomp 1.6 has no donor random effect (mild pseudoreplication).
#   - geom_smooth lines are UNADJUSTED descriptive views; the formal test is the
#     adjusted sccomp interaction, whose slopes need not match the raw lines.
#   - sccomp FDR is per model; the 3 pathology measures are separate analyses.
#
# DATA: LD_B04a_v02_seur.qs (metadata for cluster counts) + samplesheet + plaque CSV.

library(tidyverse)
library(qs)
library(Seurat)
library(sccomp)

### paths -------------------------------------------------------------------
base       = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
b04_path   = file.path(base, "LD_B_AST_analysis_output/LD_B04a_v02_seur.qs")
clust_csv  = file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
group_csv  = file.path(base, "data_TREM2_michael/A_input/group_tab.csv")
plaque_csv = file.path(base, "data_TREM2_michael/A_input/TREM2_plaque_data_Sam.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X06_"
cores_n    = 8

ord = read_csv(clust_csv, show_col_types = FALSE)
cluster_names  = unique(ord$cluster_name)
subtype_levels = unique(ord$cell_type)
grtab = read_csv(group_csv, show_col_types = FALSE)

# pathology measures (column, label)
path_measures = tribble(
  ~col,                    ~label,
  "plaque_dens",           "Amyloid plaque density",
  "pctAT8PositiveArea",    "AT8 (% area, pTau)",
  "pctPHF1PositiveArea",   "PHF1 (% area, pTau)")

variant_cols = c(CV = "dodgerblue", R62H = "orange", R47H = "grey20")
covars_adj   = c("Sex", "BrainRegion", "APOEgroup", "CD33Group", "cohort")

### sample-level table: covariates + pathology ------------------------------
plaque = read_csv(plaque_csv, show_col_types = FALSE) %>%
  dplyr::select(BrainBankNetworkIDFormatted, BrainRegion, plaque_dens = TotalDensity)

samp = grtab %>%
  dplyr::distinct(sample, group, NeuropathologicalDiagnosis, TREM2Variant, BrainRegion,
                  BrainBankNetworkIDFormatted, Sex, APOEgroup, CD33Group, cohort,
                  pctAT8PositiveArea, pctPHF1PositiveArea) %>%
  # join amyloid plaque density by donor AND region (region-specific)
  dplyr::left_join(plaque, by = c("BrainBankNetworkIDFormatted", "BrainRegion")) %>%
  dplyr::mutate(across(c(plaque_dens, pctAT8PositiveArea, pctPHF1PositiveArea), as.numeric),
                TREM2Variant = factor(TREM2Variant, levels = c("CV", "R62H", "R47H")))

################################################################################
# PART 1: pathology distribution across groups (faceted by region) + adjusted LM
################################################################################
message("\n#### PART 1: pathology distribution + adjusted linear model ", Sys.time(), "\n")

groups4 = c("Control_CV", "AD_CV", "AD_R62H", "AD_R47H")
pcomparisons = tribble(
  ~name,                ~g1,        ~g2,
  "AD_CV vs Ctrl_CV",   "AD_CV",     "Control_CV",
  "AD_R62H vs Ctrl_CV", "AD_R62H",   "Control_CV",
  "AD_R47H vs Ctrl_CV", "AD_R47H",   "Control_CV",
  "AD_R62H vs AD_CV",   "AD_R62H",   "AD_CV",
  "AD_R47H vs AD_CV",   "AD_R47H",   "AD_CV",
  "AD_R62H vs AD_R47H", "AD_R62H",   "AD_R47H")

# adjusted pairwise group contrast from an lm (treatment coding, ref = groups4[1])
lm_contrast = function(fit, g1, g2, ref) {
  b = coef(fit); V = vcov(fit)
  L = setNames(rep(0, length(b)), names(b))
  n1 = paste0("group", g1); n2 = paste0("group", g2)
  if (g1 != ref) { if (!n1 %in% names(b)) return(NULL); L[n1] = L[n1] + 1 }
  if (g2 != ref) { if (!n2 %in% names(b)) return(NULL); L[n2] = L[n2] - 1 }
  est = sum(L * b); se = sqrt(as.numeric(t(L) %*% V %*% L))
  tibble(estimate = est, se = se, t = est / se,
         p = 2 * pt(-abs(est / se), fit$df.residual))
}

p1_path = samp %>% dplyr::filter(group %in% groups4) %>%
  dplyr::mutate(group = factor(group, levels = groups4))

# adjusted model + contrasts, per measure x region (region constant within facet)
p1_stats = list()
for (m in path_measures$col) for (rg in c("MTG", "SSC")) {
  d = p1_path %>% dplyr::filter(BrainRegion == rg) %>%
    dplyr::mutate(value = .data[[m]]) %>%
    dplyr::filter(!is.na(value), !is.na(Sex), !is.na(APOEgroup), !is.na(CD33Group), !is.na(cohort))
  d$group = droplevels(d$group)
  if (dplyr::n_distinct(d$group) < 2) next
  fit = tryCatch(lm(value ~ group + Sex + APOEgroup + CD33Group + cohort, data = d),
                 error = function(e) NULL)
  if (is.null(fit)) next
  rows = pcomparisons %>% dplyr::filter(g1 %in% levels(d$group), g2 %in% levels(d$group))
  res = purrr::pmap_dfr(list(rows$name, rows$g1, rows$g2), function(nm, a, bb) {
    r = lm_contrast(fit, a, bb, levels(d$group)[1]); if (is.null(r)) NULL else dplyr::mutate(r, comparison = nm)
  })
  if (nrow(res)) { res$measure = m; res$BrainRegion = rg; p1_stats[[paste(m, rg)]] = res }
}
p1_stats = dplyr::bind_rows(p1_stats)
if (nrow(p1_stats))
  p1_stats = p1_stats %>% dplyr::group_by(measure, BrainRegion) %>%
    dplyr::mutate(FDR = p.adjust(p, "BH")) %>% dplyr::ungroup()
write_csv(p1_stats, file.path(out_dir, paste0(script_ind, "pathology_group_LM_contrasts.csv")))

# plot: facet measure x region; significant adjusted contrasts as brackets
plot_path_dist = function() {
  long = p1_path %>%
    tidyr::pivot_longer(dplyr::all_of(path_measures$col), names_to = "measure", values_to = "value") %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::left_join(path_measures, by = c("measure" = "col")) %>%
    dplyr::mutate(label = factor(label, levels = path_measures$label))
  ymax = long %>% dplyr::group_by(label, BrainRegion) %>%
    dplyr::summarise(ymax = max(value), .groups = "drop")
  brk = NULL
  if (nrow(p1_stats) > 0) {
    brk = p1_stats %>% dplyr::filter(FDR < 0.05) %>%
      dplyr::left_join(pcomparisons, by = c("comparison" = "name")) %>%
      dplyr::left_join(path_measures, by = c("measure" = "col")) %>%
      dplyr::mutate(label = factor(label, levels = path_measures$label),
                    x = match(g1, groups4), xend = match(g2, groups4),
                    txt = sprintf("%.2g", FDR)) %>%
      dplyr::group_by(label, BrainRegion) %>% dplyr::arrange(abs(xend - x), .by_group = TRUE) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>% dplyr::ungroup() %>%
      dplyr::left_join(ymax, by = c("label", "BrainRegion")) %>%
      dplyr::mutate(y = ymax * (1 + 0.12 * rank))
  }
  p = ggplot(long, aes(factor(group, levels = groups4), value)) +
    geom_boxplot(aes(fill = group), outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.6) +
    facet_grid(label ~ BrainRegion, scales = "free_y") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(title = "Pathology across donor groups (adjusted LM contrasts)",
         x = NULL, y = NULL, caption = "brackets: FDR<0.05 from covariate-adjusted linear model") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.caption = element_text(size = 7))
  if (!is.null(brk) && nrow(brk) > 0)
    p = p + geom_segment(data = brk, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.25) +
            geom_text(data = brk, aes(x = (x + xend) / 2, y = y, label = txt), inherit.aes = FALSE, vjust = -0.2, size = 2)
  p
}
ggsave(file.path(out_dir, paste0(script_ind, "pathology_distribution.pdf")), plot_path_dist(), width = 8, height = 8, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(script_ind, "pathology_distribution.png")), plot_path_dist(), width = 8, height = 8, dpi = 300)
message("Part 1 written.")

################################################################################
# PART 2: AD-only abundance ~ pathology x TREM2 variant (sccomp)
################################################################################
# Estimates are read DIRECTLY from model coefficients (NO contrast arithmetic,
# which sccomp 1.6 evaluates incorrectly -> degenerate values):
#   per-variant slopes  <- ~ 0 + TREM2Variant + TREM2Variant:path_z + covars
#   slope differences   <- ~ path_z * TREM2Variant + covars  (interaction terms)
message("\n#### PART 2: AD-only sccomp pathology x variant ", Sys.time(), "\n")

# per-sample per-cluster and per-subtype counts (need cluster_name/cell_type per cell)
message("Loading B04 object (metadata only)...")
seur = qread(b04_path); meta = seur@meta.data; rm(seur); gc()
counts_sub = meta %>% dplyr::count(sample, cluster_name, name = "N_cells") %>%
  tidyr::complete(sample, cluster_name = cluster_names, fill = list(N_cells = 0)) %>%
  dplyr::rename(cg = cluster_name)
counts_typ = meta %>% dplyr::count(sample, cell_type, name = "N_cells") %>%
  tidyr::complete(sample, cell_type = subtype_levels, fill = list(N_cells = 0)) %>%
  dplyr::rename(cg = cell_type)

# two most abundant subclusters per subtype
sel = meta %>% dplyr::count(cell_type, cluster_name, name = "n") %>%
  dplyr::group_by(cell_type) %>% dplyr::slice_max(n, n = 2, with_ties = FALSE) %>% dplyr::ungroup()
sel_clusters = sel$cluster_name[order(match(sel$cluster_name, cluster_names))]

ad = samp %>% dplyr::filter(NeuropathologicalDiagnosis == "AD")
ad_sub = counts_sub %>% dplyr::inner_join(ad, by = "sample")
ad_typ = counts_typ %>% dplyr::inner_join(ad, by = "sample")

# drop NA pathology/covariate rows, z-score pathology within the AD subset
prep = function(ct, pcol) {
  ct = ct %>% dplyr::filter(!is.na(.data[[pcol]]))
  for (v in covars_adj) ct = ct %>% dplyr::filter(!is.na(.data[[v]]))
  ct$path_z = as.numeric(scale(ct[[pcol]]))
  ct$TREM2Variant = droplevels(ct$TREM2Variant)
  ct
}
sccomp_cols = function(r) {
  if (!"factor" %in% names(r) && "parameter" %in% names(r)) r$factor = r$parameter
  list(par = intersect(c("parameter", "factor"), names(r))[1],
       cg  = intersect(c("cg", "cell_group", "cellgroup"), names(r))[1],
       keep = intersect(c("c_lower", "c_effect", "c_upper", "c_FDR"), names(r)))
}

# model A: per-variant absolute slopes (one coefficient per variant) -> forests
fit_slopes = function(ct, pcol) {
  ct = prep(ct, pcol)
  form = stats::as.formula(paste("~ 0 + TREM2Variant + TREM2Variant:path_z +", paste(covars_adj, collapse = " + ")))
  r = as.data.frame(ct |> sccomp_glm(formula_composition = form, .sample = sample,
        .cell_group = cg, .count = N_cells, bimodal_mean_variability_association = TRUE, cores = cores_n))
  cc = sccomp_cols(r)
  r %>% dplyr::filter(grepl("path_z", .data[[cc$par]]), grepl("TREM2Variant", .data[[cc$par]])) %>%
    dplyr::transmute(cellgroup = .data[[cc$cg]],
                     TREM2Variant = stringr::str_extract(.data[[cc$par]], "CV|R62H|R47H"),
                     dplyr::across(dplyr::all_of(cc$keep)))
}
# model B: slope differences vs CV (interaction terms) + CV slope -> exemplars/CSV
fit_diffs = function(ct, pcol) {
  ct = prep(ct, pcol)
  form = stats::as.formula(paste("~ path_z * TREM2Variant +", paste(covars_adj, collapse = " + ")))
  r = as.data.frame(ct |> sccomp_glm(formula_composition = form, .sample = sample,
        .cell_group = cg, .count = N_cells, bimodal_mean_variability_association = TRUE, cores = cores_n))
  cc = sccomp_cols(r)
  r %>% dplyr::filter(grepl("path_z", .data[[cc$par]])) %>%
    dplyr::mutate(variant = stringr::str_extract(.data[[cc$par]], "R62H|R47H"),
                  term = ifelse(is.na(variant), "slope_CV", paste0("diff_", variant, "_CV"))) %>%
    dplyr::transmute(cellgroup = .data[[cc$cg]], term, dplyr::across(dplyr::all_of(cc$keep)))
}

safe = function(expr) tryCatch(expr, error = function(e) { message("sccomp failed: ", conditionMessage(e)); NULL })
runset = function(fn, ct) dplyr::bind_rows(lapply(seq_len(nrow(path_measures)), function(i) {
  m = path_measures$col[i]; message("sccomp: ", m, " ...")
  r = safe(fn(ct, m)); if (!is.null(r)) r$measure = m; r
}))

slopes_sub = runset(fit_slopes, ad_sub)   # per-variant slopes, subclusters
slopes_typ = runset(fit_slopes, ad_typ)   # per-variant slopes, subtypes
diffs_sub  = runset(fit_diffs,  ad_sub)   # slope differences (interaction), subclusters
write_csv(slopes_sub, file.path(out_dir, paste0(script_ind, "slopes_subcluster.csv")))
write_csv(slopes_typ, file.path(out_dir, paste0(script_ind, "slopes_subtype.csv")))
write_csv(diffs_sub,  file.path(out_dir, paste0(script_ind, "slope_differences_subcluster.csv")))

### forest plots (per-variant slopes) ---------------------------------------
make_forest = function(sl, units, title, suffix, w, h) {
  if (is.null(sl) || nrow(sl) == 0) return(invisible())
  sl = sl %>% dplyr::left_join(path_measures, by = c("measure" = "col")) %>%
    dplyr::mutate(cellgroup = factor(cellgroup, levels = rev(units)),
                  TREM2Variant = factor(TREM2Variant, levels = c("CV", "R62H", "R47H")),
                  sig = !is.na(c_FDR) & c_FDR < 0.05,
                  label = factor(label, levels = path_measures$label))
  p = ggplot(sl, aes(c_effect, cellgroup, color = TREM2Variant)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_pointrange(aes(xmin = c_lower, xmax = c_upper, size = sig, alpha = sig),
                    position = position_dodge(0.6), fatten = 1) +
    facet_wrap(~ label, scales = "free_x") +
    scale_color_manual(values = variant_cols) +
    scale_size_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.9), guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 0.55, `TRUE` = 1), guide = "none") +
    labs(title = title, x = "slope (logit composition per SD pathology)", y = NULL,
         caption = "thick/solid = FDR<0.05; R47H exploratory (small n)") +
    theme_bw(base_size = 9) + theme(plot.caption = element_text(size = 7))
  ggsave(file.path(out_dir, paste0(script_ind, suffix, ".pdf")), p, width = w, height = h, useDingbats = FALSE)
  ggsave(file.path(out_dir, paste0(script_ind, suffix, ".png")), p, width = w, height = h, dpi = 300)
}
make_forest(slopes_sub, cluster_names, "Abundance-pathology slope per variant (all subclusters, AD)", "forest_all_subclusters", 12, 9)
make_forest(dplyr::filter(slopes_sub, cellgroup %in% sel_clusters), sel_clusters,
            "Abundance-pathology slope per variant (two largest per subtype, AD)", "forest_two_largest", 11, 4)
make_forest(slopes_typ, subtype_levels, "Abundance-pathology slope per variant (subtypes, AD)", "forest_subtype", 11, 3.5)

### regression exemplars: clusters with notable R62H-CV slope difference -----
if (!is.null(diffs_sub) && nrow(diffs_sub) > 0) {
  frac = ad_sub %>% dplyr::group_by(sample) %>% dplyr::mutate(fract = N_cells / sum(N_cells)) %>% dplyr::ungroup()
  for (i in seq_len(nrow(path_measures))) {
    m = path_measures$col[i]; lab = path_measures$label[i]
    hits = diffs_sub %>% dplyr::filter(measure == m, term == "diff_R62H_CV", !is.na(c_FDR), c_FDR < 0.1) %>% dplyr::pull(cellgroup)
    if (length(hits) == 0) next
    d = frac %>% dplyr::filter(cg %in% hits, !is.na(.data[[m]]))
    ex = ggplot(d, aes(.data[[m]], fract, color = TREM2Variant)) +
      geom_point(size = 0.6, alpha = 0.6) + geom_smooth(method = "lm", se = TRUE, linewidth = 0.5) +
      facet_wrap(~ cg, scales = "free_y") + scale_color_manual(values = variant_cols) +
      labs(title = paste0("Abundance vs ", lab, " (descriptive, unadjusted)"), x = lab, y = "Fraction of sample",
           caption = "R62H-CV interaction FDR<0.1; lines unadjusted - formal test is the sccomp interaction") +
      theme_bw(base_size = 9) + theme(plot.caption = element_text(size = 7))
    ggsave(file.path(out_dir, paste0(script_ind, "regression_", m, ".pdf")), ex, width = 9, height = 6, useDingbats = FALSE)
    ggsave(file.path(out_dir, paste0(script_ind, "regression_", m, ".png")), ex, width = 9, height = 6, dpi = 300)
  }
}

message("Done. Pathology analysis outputs written to ", out_dir)
