message("\n\n##########################################################################\n",
        "# Start LD_E04a: Incremental covariate adjustment ", Sys.time(),
        "\n##########################################################################\n\n")

# Covariates are added incrementally because the per-cluster, AD-only subsets
# cannot support the full model without collinearity or unidentifiable terms.
# Technical covariates are added before Braak because pathology correlates with
# TREM2Variant and may remove part of the effect of interest.
# Pooled results are exploratory because donors contribute several cluster
# pseudobulks; cluster_name absorbs baseline between-cluster differences.
library(qs)
library(tidyverse)
library(DESeq2)
library(colorRamps)
library(ggrepel)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "LD_E04a_v01_"

#specify output directory
out_dir = paste0(main_dir, "LD_E_DESeq_pseudobulk/")


### load pseudobulk dataset

bulk_data = qread(file = paste0(out_dir, "LD_E01_v02_bulk_data.qs"))

ref_cluster = "AST_SLC1A2_s0"


###########################################################
# prepare metadata: factor categoricals, build numeric covariates
###########################################################

meta = as.data.frame(bulk_data$meta)
rownames(meta) = meta$cluster_sample

# primary categorical variables / covariates already used in E02a2
meta$NeuropathologicalDiagnosis = factor(meta$NeuropathologicalDiagnosis, levels = c("Control", "AD"))
meta$APOEgroup   = factor(meta$APOEgroup,   levels = c("APOE4-neg", "APOE4-pos"))
meta$CD33Group   = factor(meta$CD33Group,   levels = c("CV", "CD33var"))
meta$BrainRegion = factor(meta$BrainRegion, levels = c("SSC", "MTG"))
meta$cohort      = factor(meta$cohort)
meta$Sex         = factor(meta$Sex)
meta$cluster_name = factor(meta$cluster_name)   # used as a covariate in the pooled scope
# TREM2Variant levels set per-contrast below (reference differs per contrast)

# numeric covariates (stored as character in the AST metadata) -> coerce + scale
meta$Age               = as.numeric(meta$Age)
meta$PostMortemInterval = as.numeric(meta$PostMortemInterval)

# Braak is recorded as Roman numerals (0, I-VI) with occasional ranges ("V, VI").
# Map to a numeric severity score (mean of the listed stages) -> 1 df covariate.
braak_map = c("0" = 0, "I" = 1, "II" = 2, "III" = 3, "IV" = 4, "V" = 5, "VI" = 6)
parse_braak = function(x){
  x = trimws(as.character(x))
  vapply(x, function(v){
    if (is.na(v) || v == "" || toupper(v) == "NA") return(NA_real_)
    parts = trimws(unlist(strsplit(v, "[,/;]")))
    vals  = braak_map[parts]
    if (all(is.na(vals))) return(NA_real_)
    mean(vals, na.rm = TRUE)
  }, numeric(1))
}
meta$Braak_raw = parse_braak(meta$Braak)

# scaled (centred + unit-variance) continuous covariates
meta$PMI_scaled = as.numeric(scale(meta$PostMortemInterval))
meta$Age_scaled = as.numeric(scale(meta$Age))
meta$Braak_num  = as.numeric(scale(meta$Braak_raw))

bulk_data$meta = meta


### apply the E02a2 base filter so M0 reproduces E03
# Missing values in added covariates are handled within each adjustment level.

e02_model_vars = c("cohort", "APOEgroup", "CD33Group", "BrainRegion",
                   "NeuropathologicalDiagnosis", "TREM2Variant")

for (v in e02_model_vars){
  meta = meta[!is.na(meta[[v]]), ]
}
for (v in e02_model_vars){
  meta$.mv = meta[[v]]
  t2 = meta %>% group_by(cluster_name, .mv) %>% summarise(N = n(), .groups = "drop")
  t3 = t2[t2$N >= 2, ] %>% group_by(cluster_name) %>% summarise(N_levels = n(), .groups = "drop")
  meta = meta[meta$cluster_name %in% t3$cluster_name[t3$N_levels > 1], ]
}
meta$.mv = NULL
meta = droplevels(meta)

bulk_data$meta = meta


### gene universe = exactly the DEG-union E03 plotted
# E03 restricts the scatter to unique(unlist(bulk_data$DEGs)) (union of DEGs,
# padj<0.1, across ALL E02a2 comparisons). That same gene list was written by
# E02a2 to DEGs_by_cluster_genes.csv (the TF columns are a subset, so the union
# of all cells == unique(unlist($DEGs))). Reading the 1.8 MB CSV avoids loading
# the 3.4 GB E02a2 .qs object just to extract a list of gene names.

deg_tab = read_csv(paste0(out_dir, "LD_E02a2_v02_DEGs_by_cluster_genes.csv"),
                   show_col_types = FALSE)
deg_universe = unique(unlist(deg_tab, use.names = FALSE))
deg_universe = deg_universe[!is.na(deg_universe) & deg_universe != ""]
message("   Gene universe (E03 DEG-union): ", length(deg_universe), " genes")


### filter count matrix (drop genes with <0.1 counts/cell in all pseudobulks; as E02a2)

m1 = bulk_data$counts
m1 = m1[, meta$cluster_sample]
keep_genes = rownames(m1)[apply(m1, 1, max) > 0.1]
comp_counts = m1[keep_genes, ]

comp_clusters = unique(as.character(meta$cluster_name))


###########################################################
# define covariate levels (cumulative; "clean -> severity")
###########################################################

base_covars = c("APOEgroup", "CD33Group", "BrainRegion")

covar_sets = list(
  M0_base            = base_covars,
  M1_Sex             = c(base_covars, "Sex"),
  M2_Sex_cohort      = c(base_covars, "Sex", "cohort"),
  M3_Sex_cohort_PMI  = c(base_covars, "Sex", "cohort", "PMI_scaled"),
  M4_..Age           = c(base_covars, "Sex", "cohort", "PMI_scaled", "Age_scaled"),
  M5_..Braak         = c(base_covars, "Sex", "cohort", "PMI_scaled", "Age_scaled", "Braak_num")
)
covar_levels = names(covar_sets)


###########################################################
# define the 4 DESeq2 contrasts feeding the requested plots
###########################################################
# Each contrast = a subset of the data + a tested variable (and its reference
# level). Reference levels preserve the E02a2 log2FC direction.

contrast_defs = list(
  CV_AD_vs_Control = list(
    subset = function(m) m[m$TREM2Variant == "CV", ],
    tested = "NeuropathologicalDiagnosis",
    set_levels = function(m){ m$NeuropathologicalDiagnosis = factor(m$NeuropathologicalDiagnosis, levels = c("Control", "AD")); m }
  ),
  R62H_vs_CV = list(
    subset = function(m) m[m$NeuropathologicalDiagnosis == "AD" & m$TREM2Variant %in% c("CV", "R62H"), ],
    tested = "TREM2Variant",
    set_levels = function(m){ m$TREM2Variant = factor(m$TREM2Variant, levels = c("CV", "R62H")); m }
  ),
  R47H_vs_CV = list(
    subset = function(m) m[m$NeuropathologicalDiagnosis == "AD" & m$TREM2Variant %in% c("CV", "R47H"), ],
    tested = "TREM2Variant",
    set_levels = function(m){ m$TREM2Variant = factor(m$TREM2Variant, levels = c("CV", "R47H")); m }
  ),
  R47H_vs_R62H = list(
    subset = function(m) m[m$NeuropathologicalDiagnosis == "AD" & m$TREM2Variant %in% c("R62H", "R47H"), ],
    tested = "TREM2Variant",
    set_levels = function(m){ m$TREM2Variant = factor(m$TREM2Variant, levels = c("R62H", "R47H")); m }
  )
)

# scatter pairs: c(comp1 = y-axis, comp_ref = x-axis)  (matches E03a2 orientation)
plot_pairs = list(
  "P1_R62H_vs_CV__vs__R47H_vs_CV"        = c("R62H_vs_CV",   "R47H_vs_CV"),
  "P2_R47H_vs_CV__vs__AD_vs_Control"     = c("R47H_vs_CV",   "CV_AD_vs_Control"),
  "P3_R62H_vs_CV__vs__AD_vs_Control"     = c("R62H_vs_CV",   "CV_AD_vs_Control"),
  "P4_R47H_vs_R62H__vs__AD_vs_Control"   = c("R47H_vs_R62H", "CV_AD_vs_Control")
)


###########################################################
# helper: guarded DESeq2 LRT run (refactored from E02a2)
###########################################################
# Drops single-level / collinear covariates, omits the run if a covariate is
# collinear with the tested variable or there are too few samples. Returns the
# results table (or NULL) plus the final formulas and a log message.

run_deseq_guarded = function(meta_sub, counts_sub, covars, tested, add_cluster = FALSE){

  msg = character(0)

  # add cluster_name as covariate for the pooled scope
  red_vars = covars
  if (add_cluster) red_vars = c(red_vars, "cluster_name")

  # keep only samples with complete data in all modelled variables
  model_vars = c(red_vars, tested)
  meta_sub = meta_sub[complete.cases(meta_sub[, model_vars, drop = FALSE]), ]
  meta_sub = droplevels(meta_sub)
  if (nrow(meta_sub) < 3) return(list(res = NULL, n = nrow(meta_sub),
                                      form_full = NA, form_red = NA,
                                      msg = "too few samples after NA removal"))

  # drop reduced-model covariates with a single observed level
  for (v in red_vars){
    if (!is.numeric(meta_sub[[v]]) && length(unique(meta_sub[[v]])) < 2){
      red_vars = red_vars[red_vars != v]
      msg = c(msg, paste0("dropped '", v, "' (single level)"))
    }
  }

  # drop collinear pairs of categorical reduced covariates (keep the second)
  cat_vars = red_vars[!vapply(red_vars, function(v) is.numeric(meta_sub[[v]]), logical(1))]
  done = character(0)
  for (v1 in cat_vars){
    for (v2 in setdiff(cat_vars, c(v1, done))){
      if (!(v1 %in% red_vars) || !(v2 %in% red_vars)) next
      tab = meta_sub %>% group_by(pick(all_of(c(v1, v2)))) %>% summarise(n = n(), .groups = "drop")
      if (nrow(tab) < 3){
        red_vars = red_vars[red_vars != v1]
        msg = c(msg, paste0("dropped '", v1, "' (collinear with '", v2, "')"))
      }
    }
    done = c(done, v1)
  }

  # omit run if tested variable is collinear with any remaining covariate
  omit = FALSE
  for (v in red_vars){
    tab = meta_sub %>% group_by(pick(all_of(c(tested, v)))) %>% summarise(n = n(), .groups = "drop")
    if (nrow(tab) < 3 && !is.numeric(meta_sub[[v]])){
      omit = TRUE
      msg = c(msg, paste0("OMIT: tested '", tested, "' collinear with '", v, "'"))
    }
  }
  if (omit) return(list(res = NULL, n = nrow(meta_sub), form_full = NA, form_red = NA,
                        msg = paste(msg, collapse = "; ")))

  form_red  = paste0("~", paste(red_vars, collapse = "+"))
  form_full = paste0(form_red, "+", tested)

  # omit if there are not more samples than model coefficients (model unidentifiable).
  # uses the actual design-matrix rank (accounts for multi-level factors e.g. cluster_name)
  n_coef = tryCatch(ncol(model.matrix(as.formula(form_full), data = meta_sub)),
                    error = function(e) Inf)
  if (nrow(meta_sub) <= n_coef){
    msg = c(msg, paste0("OMIT: too few samples (", nrow(meta_sub), ") for ", n_coef, " coefficients"))
    return(list(res = NULL, n = nrow(meta_sub), form_full = form_full, form_red = form_red,
                msg = paste(msg, collapse = "; ")))
  }

  counts_sub = counts_sub[, meta_sub$cluster_sample]

  dds = DESeqDataSetFromMatrix(counts_sub, colData = meta_sub, design = as.formula(form_full))
  dds = DESeq(dds, test = "LRT", reduced = as.formula(form_red))

  res = as.data.frame(results(dds))
  res = cbind(gene = rownames(res), res)

  list(res = res, n = nrow(meta_sub), form_full = form_full, form_red = form_red,
       msg = paste(msg, collapse = "; "))
}


###########################################################
# run all contrasts x covariate levels x scopes (parallel + resumable)
###########################################################
# Independent LRTs are checkpointed for resumable parallel execution.
# Delete the checkpoint directory if model or filtering settings change.

scopes = c(as.character(comp_clusters), "pooled")

ckpt_dir = paste0(out_dir, script_ind, "ckpt/")
if (!dir.exists(ckpt_dir)) dir.create(ckpt_dir, recursive = TRUE)

# filesystem-safe checkpoint path for one job
ckpt_file = function(scope, level, cn){
  key = gsub("[^A-Za-z0-9_.-]", "-", paste(scope, level, cn, sep = "__"))
  paste0(ckpt_dir, key, ".rds")
}

# worker: compute one job and checkpoint it (skips if already on disk)
run_job = function(job){
  f = ckpt_file(job$scope, job$level, job$contrast)
  if (file.exists(f)) return(invisible(NULL))

  pooled = job$scope == "pooled"
  meta_scope = if (pooled) meta else meta[meta$cluster_name == job$scope, ]
  cdef = contrast_defs[[job$contrast]]
  meta_sub = cdef$set_levels(cdef$subset(meta_scope))

  run = tryCatch(
    run_deseq_guarded(meta_sub, comp_counts, covar_sets[[job$level]], cdef$tested,
                      add_cluster = pooled),
    error = function(e) list(res = NULL, n = nrow(meta_sub), form_full = NA,
                             form_red = NA, msg = paste0("ERROR: ", conditionMessage(e)))
  )
  out = list(scope = job$scope, level = job$level, contrast = job$contrast,
             res = run$res, n = run$n, form_full = run$form_full,
             form_red = run$form_red, msg = run$msg)
  # write atomically (temp then rename) so a kill mid-write can't corrupt a ckpt
  tmp = paste0(f, ".tmp"); saveRDS(out, tmp); file.rename(tmp, f)
  message("   done: ", job$scope, " | ", job$level, " | ", job$contrast,
          " (n=", run$n, ") - ", Sys.time())
  invisible(NULL)
}

# full job list, then keep only the ones not yet checkpointed
jobs_all = list()
for (scope in scopes) for (level in covar_levels) for (cn in names(contrast_defs)){
  jobs_all[[length(jobs_all) + 1]] = list(scope = scope, level = level, contrast = cn)
}
jobs = Filter(function(j) !file.exists(ckpt_file(j$scope, j$level, j$contrast)), jobs_all)
message("   ", length(jobs), " of ", length(jobs_all), " DESeq jobs to run (rest already checkpointed)")

# cores from the PBS allocation (fallback to detectCores-1)
n_cores = suppressWarnings(as.integer(Sys.getenv("NCPUS")))
if (is.na(n_cores) || n_cores < 1) n_cores = max(1, parallel::detectCores() - 1)
message("   running on ", n_cores, " cores - ", Sys.time())

if (length(jobs) > 0){
  invisible(parallel::mclapply(jobs, run_job, mc.cores = n_cores, mc.preschedule = FALSE))
}

# collect all checkpoints into deseq_res + model_log
deseq_res = list()
model_log = NULL
for (scope in scopes) for (level in covar_levels) for (cn in names(contrast_defs)){
  f = ckpt_file(scope, level, cn)
  if (!file.exists(f)) next
  o = readRDS(f)
  deseq_res[[paste(scope, level, cn, sep = "|")]] = o$res
  model_log = rbind(model_log, tibble(
    scope = scope, level = level, contrast = cn, n_samples = o$n,
    form_full = o$form_full, form_red = o$form_red, note = o$msg))
}

write_csv(model_log, file = paste0(out_dir, script_ind, "deseq_model_log.csv"))

bulk_data$E04_deseq_res = deseq_res
bulk_data$E04_model_log = model_log
qsave(bulk_data, file = paste0(out_dir, script_ind, "bulk_data.qs"))


###########################################################
# build per-pair log2FC tables (reg classification as in E03a2: pval<0.05)
###########################################################

# classify regulation direction from uncorrected p-value (matches E03a2)
classify_reg = function(res){
  reg = rep("nreg", nrow(res))
  reg[!is.na(res$pvalue) & res$pvalue < 0.05 & res$log2FoldChange > 0] = "up"
  reg[!is.na(res$pvalue) & res$pvalue < 0.05 & res$log2FoldChange < 0] = "down"
  reg
}

reg_levels = c("down_down", "up_up", "down_up", "up_down",
               "down_nreg", "up_nreg", "nreg_down", "nreg_up", "nreg_nreg")
reg_cols = c("down_down" = "blue", "up_up" = "red", "down_up" = "magenta3",
             "up_down" = "magenta3", "down_nreg" = "grey40", "up_nreg" = "grey40",
             "nreg_down" = "grey40", "nreg_up" = "grey40", "nreg_nreg" = "grey")

# assemble a combined table per (scope, pair): one row per gene per covariate level
build_pair_tab = function(scope, pair){
  comp1   = pair[1]   # y-axis
  comp_ref = pair[2]  # x-axis
  out = NULL
  for (level in covar_levels){
    r1 = deseq_res[[paste(scope, level, comp1,   sep = "|")]]
    r2 = deseq_res[[paste(scope, level, comp_ref, sep = "|")]]
    if (is.null(r1) || is.null(r2)) next
    g = intersect(intersect(r1$gene, r2$gene), deg_universe)   # E03 DEG-union restriction
    r1 = r1[match(g, r1$gene), ]; r2 = r2[match(g, r2$gene), ]
    t = tibble(level = level, gene = g,
               log2FC = r1$log2FoldChange, pval = r1$pvalue,
               log2FC_ref = r2$log2FoldChange, pval_ref = r2$pvalue)
    t = t[!is.na(t$log2FC) & !is.na(t$log2FC_ref) & !is.na(t$pval) & !is.na(t$pval_ref), ]
    if (nrow(t) == 0) next
    t$reg      = classify_reg(tibble(log2FoldChange = t$log2FC,     pvalue = t$pval))
    t$reg_ref  = classify_reg(tibble(log2FoldChange = t$log2FC_ref, pvalue = t$pval_ref))
    t$reg_group = paste0(t$reg, "_", t$reg_ref)
    out = rbind(out, t)
  }
  out
}


###########################################################
# robustness summary: slope / r / concordance vs covariate level
###########################################################

summary_tab = NULL

for (scope in scopes){
  for (pn in names(plot_pairs)){
    tab = build_pair_tab(scope, plot_pairs[[pn]])
    if (is.null(tab)) next
    for (level in covar_levels){
      d = tab[tab$level == level, ]
      if (nrow(d) < 3) next
      # genes regulated at p < 0.05 in both contrasts
      db = d[d$reg != "nreg" & d$reg_ref != "nreg", ]
      slope_fun = function(x) if (nrow(x) >= 3) unname(coef(lm(log2FC ~ log2FC_ref, x))[2]) else NA_real_
      r_fun     = function(x) if (nrow(x) >= 3) suppressWarnings(cor(x$log2FC_ref, x$log2FC)) else NA_real_
      conc_fun  = function(x) if (nrow(x) >= 1) mean(sign(x$log2FC) == sign(x$log2FC_ref)) else NA_real_
      summary_tab = rbind(summary_tab, tibble(
        scope = scope, pair = pn, level = level,
        n_shared = nrow(d), r_shared = r_fun(d), slope_shared = slope_fun(d),
        n_reg_both = nrow(db), r_reg_both = r_fun(db), slope_reg_both = slope_fun(db),
        frac_concordant_reg_both = conc_fun(db)
      ))
    }
  }
}

summary_tab$level = factor(summary_tab$level, levels = covar_levels)
write_csv(summary_tab, file = paste0(out_dir, script_ind, "effect_robustness_summary.csv"))


### robustness plot: slope & Pearson r across covariate levels (reg_both genes)

sl = summary_tab %>%
  select(scope, pair, level, slope_reg_both, r_reg_both) %>%
  pivot_longer(c(slope_reg_both, r_reg_both), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, slope_reg_both = "regression slope", r_reg_both = "Pearson r"))

p_rob = ggplot(sl, aes(x = level, y = value, group = metric, colour = metric)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line() + geom_point() +
  facet_grid(scope ~ pair) +
  scale_colour_manual(values = c("regression slope" = "dodgerblue", "Pearson r" = "orange")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Robustness of log2FC correlation to covariate adjustment (reg_both genes)",
       x = "covariate adjustment level", y = "slope / Pearson r")

pdf(file = paste0(out_dir, script_ind, "effect_robustness_slope_r_vs_covar_level.pdf"),
    width = 16, height = 10)
plot(p_rob)
dev.off()


###########################################################
# scatter plots: reg_both genes, faceted by covariate level
###########################################################
# One page per (scope x pair). Facets = covariate levels (M0..M5) so any
# weakening / sign change of the slope is visible across the row. Fixed axis
# limits across facets so magnitude changes are comparable.

make_facet_plot = function(tab, title){
  d = tab[tab$reg != "nreg" & tab$reg_ref != "nreg", ]   # reg_both
  if (is.null(d) || nrow(d) < 3) return(NULL)
  d$level = factor(d$level, levels = covar_levels)

  # label up to 6 genes furthest from origin per facet
  d$label = ""
  for (lv in unique(d$level)){
    idx = which(d$level == lv)
    dist = sqrt(d$log2FC[idx]^2 + d$log2FC_ref[idx]^2)
    top = idx[order(-dist)][seq_len(min(6, length(idx)))]
    d$label[top] = d$gene[top]
  }
  d = d[order(match(d$reg_group, reg_levels)), ]

  # per-facet slope/r annotation
  stat = do.call(rbind, lapply(split(d, d$level), function(x){
    tibble(level = x$level[1], n = nrow(x),
           slope = if (nrow(x) >= 3) unname(coef(lm(log2FC ~ log2FC_ref, x))[2]) else NA_real_,
           r = if (nrow(x) >= 3) suppressWarnings(cor(x$log2FC_ref, x$log2FC)) else NA_real_)
  }))
  stat$level = factor(stat$level, levels = covar_levels)
  stat$lab = paste0("n=", stat$n,
                    "\nslope=", round(stat$slope, 2),
                    "\nr=", round(stat$r, 2))

  lim = max(abs(c(d$log2FC, d$log2FC_ref)), na.rm = TRUE)

  ggplot(d, aes(x = log2FC_ref, y = log2FC, colour = reg_group)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
    geom_smooth(colour = "grey30", method = "lm", formula = y ~ x, se = TRUE) +
    geom_point(alpha = 0.8, size = 1.6) +
    geom_label_repel(aes(label = label), seed = 42, size = 2.4,
                     min.segment.length = 0.2, max.overlaps = Inf, max.time = 3) +
    geom_text(data = stat, aes(x = -lim, y = lim, label = lab),
              inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.8, colour = "grey20") +
    scale_colour_manual(limits = reg_levels, values = reg_cols) +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    facet_wrap(~ level, nrow = 1) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(title = title,
         x = "log2FC (x-axis contrast)", y = "log2FC (y-axis contrast)")
}

### per-cluster PDF

pl = list()
for (pn in names(plot_pairs)){
  for (cl in as.character(comp_clusters)){
    tab = build_pair_tab(cl, plot_pairs[[pn]])
    if (is.null(tab)) next
    ttl = paste0(cl, "  |  ", pn, "   (y = ", plot_pairs[[pn]][1],
                 " , x = ", plot_pairs[[pn]][2], ")")
    p = make_facet_plot(tab, ttl)
    if (!is.null(p)) pl[[paste0(pn, "__", cl)]] = p
  }
}

pdf(file = paste0(out_dir, script_ind, "Log2FC_corr_reg_both_per_cluster.pdf"),
    width = 18, height = 5)
for (p in pl) plot(p)
dev.off()


### pooled PDF

pl = list()
for (pn in names(plot_pairs)){
  tab = build_pair_tab("pooled", plot_pairs[[pn]])
  if (is.null(tab)) next
  ttl = paste0("POOLED (all clusters)  |  ", pn, "   (y = ", plot_pairs[[pn]][1],
               " , x = ", plot_pairs[[pn]][2], ")")
  p = make_facet_plot(tab, ttl)
  if (!is.null(p)) pl[[pn]] = p
}

pdf(file = paste0(out_dir, script_ind, "Log2FC_corr_reg_both_pooled.pdf"),
    width = 18, height = 5)
for (p in pl) plot(p)
dev.off()


#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed LD_E04a ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")
