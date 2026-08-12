message("\n\n##########################################################################\n",
        "# Start LD_E04c: Incremental covariate adjustment (5-covariate base) ", Sys.time(),
        "\n##########################################################################\n\n")

# Age and PMI are added before Braak because Braak is a severity proxy correlated
# with TREM2 genotype and may constitute over-adjustment. The pooled model is
# exploratory because DESeq2 treats a donor's cluster pseudobulks independently.

library(qs)
library(tidyverse)
library(DESeq2)
library(colorRamps)
library(ggrepel)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "LD_E04c_"

#specify output directory
out_dir = paste0(main_dir, "LD_E_DESeq_pseudobulk/")


### load pseudobulk dataset (same object that produced Michael's plot)

bulk_data = qread(file = paste0(out_dir, "LD_E01_v02_bulk_data.qs"))

ref_cluster = "AST_SLC1A2_s0"


### prepare metadata --------------------------------------------------------

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


### apply the SAME sample/cluster filtering as LD_E02a2 (so M0 reproduces E03)
# E02a2 drops samples with NA in any of its model_vars up front, then keeps only
# clusters with >=2 samples in >=2 levels for each model_var (E02a2 lines 50-78).
# This defines the shared base sample/cluster set; the extra covariates' NAs are
# handled per-level inside run_deseq_guarded (a no-op at M0).

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


### gene universe = the DEG-union of the 5-covariate E02c results
# Union of DEGs (padj<0.1) across all E02c comparisons, written by E02c to
# DEGs_by_cluster_genes.csv (the TF columns are a subset, so the union of all
# cells == unique(unlist($DEGs))). Reading the CSV avoids loading the 3.2 GB
# E02c .qs object just to extract a list of gene names.

deg_tab = read_csv(paste0(out_dir, "LD_E02c/LD_E02c_v01_DEGs_by_cluster_genes.csv"),
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


### cumulative covariate levels ---------------------------------------------

# baseline = the agreed 5-covariate model (E02c); increments added cumulatively
# in "clean -> severity" order: Age and PMI (complete, demographic/technical) first,
# Braak (severity proxy, correlates with TREM2 genotype) last as the over-correction test.
base_covars = c("cohort", "APOEgroup", "CD33Group", "BrainRegion", "Sex")

covar_sets = list(
  M0_base          = base_covars,
  M1_Age           = c(base_covars, "Age_scaled"),
  M2_Age_PMI       = c(base_covars, "Age_scaled", "PMI_scaled"),
  M3_Age_PMI_Braak = c(base_covars, "Age_scaled", "PMI_scaled", "Braak_num")
)
covar_levels = names(covar_sets)

# display labels for the incremental columns (M0 model goes in the figure caption)
covar_level_labels = c(M0_base = "M0", M1_Age = "+Age",
                       M2_Age_PMI = "+PMI", M3_Age_PMI_Braak = "+Braak")


### DESeq2 contrasts ---------------------------------------------------------
# Each contrast = a subset of the data + a tested variable (and its reference
# level). Reference levels chosen so log2FC signs match the original E02a2 plot.

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


### guarded DESeq2 LRT -------------------------------------------------------
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


### run contrasts across covariate levels and scopes ------------------------
# Each (scope, level, contrast) is an independent DESeq2 LRT. There are several
# hundred and each takes ~30-45s, so single-threaded this is ~4-5h. We therefore
# (a) run them in PARALLEL across the allocated cores and (b) CHECKPOINT every
# result to its own .rds file. If the job is killed (e.g. walltime), just
# resubmit: finished runs are detected on disk and skipped, so it resumes where
# it stopped. (If you change the MODELS - covar_sets/contrasts/filtering - delete
# the ckpt/ directory first so stale results are not reused.)

# Refit M0 in every scope so the adjustment ladder uses one gene universe.
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

### robustness summary ------------------------------------------------------
# Written as a small CSV so the supplementary table can be built on a laptop
# without loading the .qs (which needs qs + this environment).
#
# TWO gene sets per cell, because they answer different questions:
#   r_shared   - all genes of the E02c DEG universe tested in BOTH contrasts
#   r_sig_both - genes significant in BOTH contrasts at padj < 0.1. This is the
#                set the LD_X07 concordance panel plots and labels with r, so
#                these are the values a reader can check against the figure.

norm_res = function(res){
  res = as.data.frame(res)
  g = if ("gene" %in% names(res)) as.character(res$gene) else rownames(res)
  tibble(gene = g, log2FC = res$log2FoldChange, padj = res$padj)
}

# the three pairs shown in the LD_X07 concordance figure (y vs x)
summary_pairs = list(
  P_R62H_vs_AD  = c("R62H_vs_CV", "CV_AD_vs_Control"),
  P_R47H_vs_AD  = c("R47H_vs_CV", "CV_AD_vs_Control"),
  P_R62H_vs_R47H = c("R62H_vs_CV", "R47H_vs_CV")
)
SIG_CUT = 0.1

r_of = function(x, y) if (length(x) >= 3) suppressWarnings(cor(x, y)) else NA_real_

robust_tab = NULL
for (scope in scopes) for (pn in names(summary_pairs)) for (level in covar_levels){
  pr = summary_pairs[[pn]]
  r1 = deseq_res[[paste(scope, level, pr[1], sep = "|")]]
  r2 = deseq_res[[paste(scope, level, pr[2], sep = "|")]]
  if (is.null(r1) || is.null(r2)) next
  d = dplyr::inner_join(norm_res(r1), norm_res(r2), by = "gene", suffix = c("_y", "_x"))
  d = d[d$gene %in% deg_universe, ]                       # same DEG universe as everywhere else
  d = d[!is.na(d$log2FC_y) & !is.na(d$log2FC_x), ]
  db = d[!is.na(d$padj_y) & !is.na(d$padj_x) & d$padj_y < SIG_CUT & d$padj_x < SIG_CUT, ]
  robust_tab = rbind(robust_tab, tibble(
    scope = scope, pair = pn, level = level,
    n_shared = nrow(d),  r_shared   = r_of(d$log2FC_x,  d$log2FC_y),
    n_sig_both = nrow(db), r_sig_both = r_of(db$log2FC_x, db$log2FC_y)))
}
if (!is.null(robust_tab)) {
  robust_tab$level = factor(robust_tab$level, levels = covar_levels)
  write_csv(robust_tab, file = paste0(out_dir, script_ind, "effect_robustness_summary.csv"))
  message("   wrote effect_robustness_summary.csv (", nrow(robust_tab), " rows)")
}

message("\n\n##########################################################################\n",
        "# Completed LD_E04c ", Sys.time(),
        "\n##########################################################################\n\n")

sessionInfo()
