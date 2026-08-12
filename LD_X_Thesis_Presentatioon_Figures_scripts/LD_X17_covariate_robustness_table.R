# LD_X17: Log2FC concordance under cumulative covariate adjustment.

library(tidyverse)

### paths -------------------------------------------------------------------
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")
e_out      = file.path(base, "LD_E_DESeq_pseudobulk")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
clust_csv  = file.path(base, "LD_B_AST_analysis_output/LD_B03a_cluster_assignment.csv")
script_ind = "LD_X17_"

# SWITCH TO "E04c" once LD_E04c has been rerun with per-subcluster scopes.
SOURCE = "E04c"

# WHICH GENES ENTER EACH CORRELATION. Only applies to SOURCE = "E04c", and needs
# LD_X17b_recompute_robustness_measures.R to have been run (it derives all three
# from LD_E04c's checkpoints without refitting anything).
#   "nominal"  nominal p < 0.05 in BOTH contrasts. Fills every cell (median ~170
#              genes) and still selects for genes with an effect in both. This is
#              the rule the Results text was written from.
#   "sig_both" padj < 0.1 in both contrasts - matches the LD_X07 figure exactly,
#              but at subcluster level the median is 0 genes, so only 39 of 156
#              cells are computable. Fine for the pooled row, not for a table.
#   "shared"   no per-pair filter: every gene of the DEG universe present in both.
#              Complete, unbiased, but correlations are much weaker (|r| ~0.1-0.3)
#              because they include genes with no effect in either contrast.
R_MEASURE = "nominal"

### source-specific configuration -------------------------------------------
# Each source has its own covariate-ladder names, r column and pair names; the
# rest of the script is identical for both.
cfg = list(
  E04c = list(
    # the ALLMEASURES file (from LD_X17b) carries all three rules; fall back to
    # LD_E04c's own summary, which has r_shared / r_sig_both only
    csv    = { f = file.path(e_out, "LD_E04c_effect_robustness_summary_ALLMEASURES.csv")
               if (file.exists(f)) f else file.path(e_out, "LD_E04c_effect_robustness_summary.csv") },
    r_col  = switch(R_MEASURE, nominal = "r_nom_both", sig_both = "r_sig_both", shared = "r_shared",
                    stop("R_MEASURE must be 'nominal', 'sig_both' or 'shared'")),
    n_col  = switch(R_MEASURE, nominal = "n_nom_both", sig_both = "n_sig_both", shared = "n_shared"),
    levels = c(M0_base = "Base (5 covariates)", M1_Age = "+ Age",
               M2_Age_PMI = "+ Age + PMI", M3_Age_PMI_Braak = "+ Age + PMI + Braak"),
    pairs  = c(P_R62H_vs_AD  = "R62H vs CV  |  AD vs Control",
               P_R47H_vs_AD  = "R47H vs CV  |  AD vs Control",
               P_R62H_vs_R47H = "R62H vs CV  |  R47H vs CV"),
    note   = switch(R_MEASURE,
                    nominal  = "genes with nominal p < 0.05 in both contrasts",
                    sig_both = "genes with padj < 0.1 in both contrasts",
                    shared   = "all genes of the DEG union present in both contrasts")),
  E04a = list(
    csv    = file.path(e_out, "LD_E04a_v01_effect_robustness_summary.csv"),
    r_col  = "r_reg_both", n_col = "n_reg_both",
    # E04a builds up from a 3-covariate base, so the 5-covariate model is its M2;
    # it also adds PMI BEFORE Age (E04c adds Age first). Same sets at base and end.
    levels = c(M2_Sex_cohort = "Base (5 covariates)", M3_Sex_cohort_PMI = "+ PMI",
               `M4_..Age` = "+ PMI + Age", `M5_..Braak` = "+ PMI + Age + Braak"),
    pairs  = c(P3_R62H_vs_CV__vs__AD_vs_Control  = "R62H vs CV  |  AD vs Control",
               P2_R47H_vs_CV__vs__AD_vs_Control  = "R47H vs CV  |  AD vs Control",
               P1_R62H_vs_CV__vs__R47H_vs_CV     = "R62H vs CV  |  R47H vs CV"),
    note   = "genes regulated (p < 0.05) in both contrasts")
)[[SOURCE]]

if (!file.exists(cfg$csv))
  stop("Source '", SOURCE, "' not available: ", cfg$csv,
       if (SOURCE == "E04c") "\n  -> rerun LD_E04c with per-subcluster scopes first." else "")
message("Using base: ", base, " | source: ", SOURCE, " (pinned)")

### data ---------------------------------------------------------------------
ord = read_csv(clust_csv, show_col_types = FALSE)
scope_order = c("pooled", unique(ord$cluster_name))     # pooled first, then families in order

raw = read_csv(cfg$csv, show_col_types = FALSE)
dat = raw %>%
  filter(level %in% names(cfg$levels), pair %in% names(cfg$pairs)) %>%
  mutate(Model    = factor(unname(cfg$levels[as.character(level)]), levels = unname(cfg$levels)),
         Contrast = factor(unname(cfg$pairs[as.character(pair)]),   levels = unname(cfg$pairs)),
         Scope    = factor(scope, levels = scope_order),
         r = .data[[cfg$r_col]], n = .data[[cfg$n_col]]) %>%
  filter(!is.na(Scope)) %>%
  arrange(Scope, Contrast, Model)

lab_scope = function(x) ifelse(x == "pooled", "All astrocytes (pooled)", as.character(x))

# wide: one row per scope, one column per (contrast, model)
wide_for = function(contrasts_keep, scopes_keep = levels(dat$Scope)) {
  dat %>%
    filter(Contrast %in% contrasts_keep, Scope %in% scopes_keep) %>%
    mutate(col = paste(Contrast, Model, sep = " || "),
           val = ifelse(is.na(r), "-", sprintf("%.2f", r))) %>%
    select(Scope, col, val) %>%
    pivot_wider(names_from = col, values_from = val) %>%
    arrange(Scope) %>%
    mutate(Scope = lab_scope(Scope)) %>%
    rename(` ` = Scope)
}

ct = levels(dat$Contrast)
block1 = wide_for(ct[1:2])                       # contrasts 1 + 2, all scopes
sc     = levels(dat$Scope); sc = sc[sc %in% dat$Scope]
half   = ceiling(length(sc) / 2)
b2a    = wide_for(ct[3], sc[1:half])
b2b    = wide_for(ct[3], sc[(half + 1):length(sc)])
# pad the shorter half so the two can sit side by side
if (nrow(b2a) > nrow(b2b)) b2b[(nrow(b2b) + 1):nrow(b2a), ] = ""
block2 = bind_cols(b2a, tibble(`  ` = rep("", nrow(b2a))),
                   set_names(b2b, paste0(names(b2b), " ")))

### save CSVs ----------------------------------------------------------------
write_csv(dat %>% transmute(Scope = lab_scope(Scope), Contrast, Model, r, n_genes = n),
          file.path(out_dir, paste0(script_ind, "covariate_robustness_long.csv")))
write_csv(block1, file.path(out_dir, paste0(script_ind, "TableBlock1.csv")))
write_csv(block2, file.path(out_dir, paste0(script_ind, "TableBlock2.csv")))

### render -------------------------------------------------------------------
cap = paste0("Supplementary Table 8. Pearson correlation (r) between the log2 fold ",
             "changes of two contrasts, under cumulative covariate adjustment, for the ",
             "pooled astrocyte model and each subcluster. Covariates are added to the ",
             "five-covariate base model (cohort, APOE group, CD33 group, brain region, ",
             "sex). r is computed over ", cfg$note, ". Gene counts per cell are given in ",
             script_ind, "covariate_robustness_long.csv. A dash indicates too few shared ",
             "genes to correlate.")

if (!requireNamespace("flextable", quietly = TRUE)) {
  message("flextable not installed - CSVs written only.")
} else {
  library(flextable)
  # two-level header: contrast spanning its four covariate columns
  mk = function(df, cap = NULL) {
    nm  = names(df)
    top = ifelse(grepl(" \\|\\| ", nm), sub(" \\|\\| .*$", "", nm), "")
    bot = ifelse(grepl(" \\|\\| ", nm), sub("^.* \\|\\| ", "", nm), nm)
    ft = flextable(df) %>%
      set_header_df(mapping = data.frame(keys = nm, top = top, bot = bot),
                    key = "keys") %>%
      merge_h(part = "header") %>%
      bold(part = "header") %>% italic(i = 1, part = "header") %>%
      align(align = "center", part = "header") %>%
      align(j = 1, align = "left", part = "all") %>%
      border_remove() %>%
      hline_top(part = "header", border = fp_border_default(width = 1.5)) %>%
      hline(i = 1, part = "header", border = fp_border_default(width = 0.75)) %>%
      hline_bottom(part = "header", border = fp_border_default(width = 1.5)) %>%
      hline_bottom(part = "body", border = fp_border_default(width = 1.5)) %>%
      fontsize(size = 8, part = "all") %>% autofit()
    if (!is.null(cap)) ft = set_caption(ft, cap)
    ft
  }
  ft1 = mk(block1, cap)
  ft2 = mk(block2)
  ok = tryCatch({ save_as_docx(ft1, ft2,
        path = file.path(out_dir, paste0(script_ind, "covariate_robustness_table.docx"))); TRUE },
        error = function(e) { message("  FAILED docx: ", conditionMessage(e)); FALSE })
  if (ok) message("  wrote covariate_robustness_table.docx")
  try(save_as_image(ft1, path = file.path(out_dir, paste0(script_ind, "Block1.png")), res = 300), silent = TRUE)
  try(save_as_image(ft2, path = file.path(out_dir, paste0(script_ind, "Block2.png")), res = 300), silent = TRUE)
}

message("Done. ", length(sc), " scopes x ", length(ct), " contrasts x ",
        length(cfg$levels), " covariate models. Outputs in: ", out_dir)
