# LD_X16: Supplementary tables for sccomp differential-abundance analyses.

library(tidyverse)

### paths -------------------------------------------------------------------
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
script_ind = "LD_X16_"
FDR_CUT    = 0.05

# WHICH LD_X05 RESULT SET TO USE. Pinned deliberately, NOT auto-latest: the v05
# rerun (which adds effect sizes) should not silently replace the tables already
# built from v03. Set to "v05" once that rerun has finished and been checked, or
# to "latest" to always take the highest version present.
RESULT_SET = "v03"

hits = list.files(out_dir, pattern = "^LD_X05_v[0-9]+_sccomp_.*\\.csv$",
                  recursive = TRUE, full.names = TRUE)   # top level or per-run subfolder
if (length(hits) == 0) stop("No LD_X05 sccomp result CSVs found under ", out_dir)
vnum = as.integer(sub(".*LD_X05_v([0-9]+)_sccomp_.*", "\\1", basename(hits)))

version = if (identical(RESULT_SET, "latest")) sprintf("v%02d", max(vnum)) else RESULT_SET
keep    = grepl(paste0("LD_X05_", version, "_sccomp_"), basename(hits), fixed = TRUE)
if (!any(keep))
  stop("Result set '", version, "' not found. Available: ",
       paste(sort(unique(sprintf("v%02d", vnum))), collapse = ", "))
res_dir = unique(dirname(hits[keep]))[1]
pre     = paste0("LD_X05_", version, "_sccomp_")
message("Using base: ", base, " | sccomp result set: ", version, " (pinned) in ", res_dir)

### the model inventory ------------------------------------------------------
# One row per sccomp model. `file` is the CSV suffix written by LD_X05.
covars_samp  = "Sex + BrainRegion + APOEgroup + CD33Group + cohort"
covars_donor = "Sex + APOEgroup + CD33Group + cohort"

models = tribble(
  ~file,                          ~Unit,   ~Grouping,        ~Adjustment,  ~Covariates,
  "subclusters_unadj",            "Sample", "18 subclusters", "Unadjusted", "",
  "subclusters_adj",              "Sample", "18 subclusters", "Adjusted",   covars_samp,
  "subtypes_unadj",               "Sample", "3 subtypes",     "Unadjusted", "",
  "subtypes_adj",                 "Sample", "3 subtypes",     "Adjusted",   covars_samp,
  "subclusters_DONORLEVEL_unadj", "Donor",  "18 subclusters", "Unadjusted", "",
  "subclusters_DONORLEVEL_adj",   "Donor",  "18 subclusters", "Adjusted",   covars_donor,
  "subtypes_DONORLEVEL_unadj",    "Donor",  "3 subtypes",     "Unadjusted", "",
  "subtypes_DONORLEVEL_adj",      "Donor",  "3 subtypes",     "Adjusted",   covars_donor
)
# The per-region models (subtypes fitted separately in MTG and SSC) are NOT
# reported: excluded from the thesis. That leaves a clean 2 x 2 x 2 design -
# unit (sample/donor) x grouping (subcluster/subtype) x adjustment - so every
# model tests the same number of contrasts per cell group.

### read every result set ----------------------------------------------------
read_one = function(f) {
  p = file.path(res_dir, paste0(pre, f, ".csv"))
  if (!file.exists(p)) { message("  missing (skipped): ", basename(p)); return(NULL) }
  read_csv(p, show_col_types = FALSE) %>% mutate(file = f)
}
res = map(models$file, read_one) %>% compact() %>% bind_rows()
if (nrow(res) == 0) stop("No sccomp result CSVs found with prefix ", pre)

has_effect = all(c("c_effect", "c_lower", "c_upper") %in% names(res))
message("  effect sizes available: ", has_effect,
        if (!has_effect) "  (FDR-only table; switch RESULT_SET to the v05 rerun to add them)" else "")

### Table A: model inventory -------------------------------------------------
counts = res %>% group_by(file) %>%
  summarise(n_tests = n(), n_sig = sum(c_FDR < FDR_CUT, na.rm = TRUE), .groups = "drop")

tabA = models %>% left_join(counts, by = "file") %>%
  mutate(Formula = ifelse(nzchar(Covariates),
                          paste0("~ 0 + group + ", Covariates), "~ 0 + group"),
         Tests = n_tests, Significant = n_sig) %>%
  select(Unit, `Cell grouping` = Grouping, Adjustment, Formula, Tests, Significant)

### Table B: significant contrasts ------------------------------------------
tabB = res %>%
  filter(c_FDR < FDR_CUT) %>%
  left_join(models %>% select(file, Unit, Grouping, Adjustment), by = "file") %>%
  arrange(Unit, Grouping, Adjustment, c_FDR) %>%
  mutate(FDR = signif(c_FDR, 2))

if (has_effect) {
  tabB = tabB %>%
    mutate(Effect = sprintf("%.2f", c_effect),
           `95% CrI` = sprintf("%.2f to %.2f", c_lower, c_upper)) %>%
    select(Unit, `Cell grouping` = Grouping, Adjustment,
           `Cell group` = cellgroup, Contrast = comparison, Effect, `95% CrI`, FDR)
} else {
  tabB = tabB %>%
    select(Unit, `Cell grouping` = Grouping, Adjustment,
           `Cell group` = cellgroup, Contrast = comparison, FDR)
}

### save CSVs ----------------------------------------------------------------
write_csv(tabA, file.path(out_dir, paste0(script_ind, "TableA_sccomp_model_inventory.csv")))
write_csv(tabB, file.path(out_dir, paste0(script_ind, "TableB_sccomp_significant.csv")))
# full result set, all models pooled, for the appendix
write_csv(res %>% left_join(models %>% select(file, Unit, Grouping, Adjustment), by = "file"),
          file.path(out_dir, paste0(script_ind, "TableC_sccomp_all_results.csv")))

### render -------------------------------------------------------------------
capA = paste0("Supplementary Table 6. Differential-abundance models fitted with sccomp, ",
              "crossing unit of analysis (sample or donor), cell grouping (18 subclusters ",
              "or 3 subtype families) and covariate adjustment. Each model was fitted with ",
              "cell-means coding and the same eight group contrasts, giving 18 x 8 = 144 ",
              "tests per subcluster model and 3 x 8 = 24 per subtype model. Donor-level ",
              "models pool each donor's counts across their samples, and omit brain region ",
              "as a covariate because it is not defined once regions are pooled.")
capB = paste0("Supplementary Table 7. All contrasts reaching FDR < ", FDR_CUT,
              " across the models in Supplementary Table 6.",
              if (has_effect) " Effect is the sccomp composition estimate on the log-ratio scale, with its 95% credible interval."
              else " Effect sizes were not retained in this run and are therefore not shown.")

if (!requireNamespace("flextable", quietly = TRUE)) {
  message("flextable not installed - CSVs written only.")
} else {
  library(flextable)
  sty = function(df, cap, small = FALSE) {
    flextable(df) %>%
      bold(part = "header") %>% italic(part = "header") %>%
      border_remove() %>%
      hline_top(part = "header", border = fp_border_default(width = 1.5)) %>%
      hline_bottom(part = "header", border = fp_border_default(width = 1.5)) %>%
      hline_bottom(part = "body", border = fp_border_default(width = 1.5)) %>%
      fontsize(size = if (small) 8 else 9, part = "all") %>%
      set_caption(cap) %>% autofit()
  }
  ftA = sty(tabA, capA)
  ftB = sty(tabB, capB, small = TRUE)
  ok = tryCatch({ save_as_docx(ftA, ftB,
        path = file.path(out_dir, paste0(script_ind, "sccomp_supplementary_tables.docx"))); TRUE },
        error = function(e) { message("  FAILED docx: ", conditionMessage(e)); FALSE })
  if (ok) message("  wrote sccomp_supplementary_tables.docx")
  try(save_as_image(ftA, path = file.path(out_dir, paste0(script_ind, "TableA.png")), res = 300), silent = TRUE)
  try(save_as_image(ftB, path = file.path(out_dir, paste0(script_ind, "TableB.png")), res = 300), silent = TRUE)
}

message("Done. ", nrow(tabA), " models, ", nrow(tabB), " significant contrasts (FDR < ",
        FDR_CUT, "), ", nrow(res), " tests total. Outputs in: ", out_dir)
