# LD_X01: Cohort summary tables for thesis supplementary.
#   Table S1 - categorical variables, folded into two panels shown side by side.
#   Table S2 - continuous variables (single table).
# Numbers are recomputed from source metadata so the tables always match the data.
# Outputs: CSVs (always) + a Word doc and PNG/HTML (if flextable is installed).

library(tidyverse)

### paths (all on RDS) ------------------------------------------------------
base        = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2"
in_dir      = file.path(base, "data_TREM2_michael/A_input")
samp_path   = file.path(in_dir, "TREM2_Samplesheet_snRNAseq_GliaEnriched_with_cohort_info.tsv")
plaque_path = file.path(in_dir, "TREM2_plaque_data_Sam.csv")
out_dir     = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind  = "X01_"
message("Writing outputs to: ", out_dir)

### load --------------------------------------------------------------------
samp = read_tsv(samp_path, show_col_types = FALSE)
donor_col = "BrainBankNetworkIDFormatted"

# donor-level table: one row per individual
donors = samp %>% distinct(.data[[donor_col]], .keep_all = TRUE)
n_don  = nrow(donors)            # 70
n_samp = nrow(samp)              # 117

### helper: counts for one categorical variable -----------------------------
# returns tibble(Variable, Category, donor_str, sample_str) in given level order
cat_block = function(col, var_label, levels_order, region_only = FALSE) {
  recode_na = function(v) ifelse(is.na(v) | v == "NA", "Missing", as.character(v))
  d = recode_na(donors[[col]]); s = recode_na(samp[[col]])
  tibble(Category = levels_order) %>%
    mutate(
      d_n = map_int(Category, ~ sum(d == .x)),
      s_n = map_int(Category, ~ sum(s == .x)),
      donor_str  = if (region_only) "—"
                   else sprintf("%d (%.1f)", d_n, 100 * d_n / n_don),
      sample_str = as.character(s_n),
      Variable   = c(var_label, rep("", length(Category) - 1))
    ) %>%
    select(Variable, Category, donor_str, sample_str)
}

### build all categorical blocks --------------------------------------------
blocks = list(
  cat_block("NeuropathologicalDiagnosis", "Diagnosis", c("AD", "Control")),
  cat_block("TREM2Variant", "TREM2 variant", c("CV", "R62H", "R47H")),
  cat_block("APOEgroup", "APOE ε4 status", c("APOE4-neg", "APOE4-pos")),
  cat_block("APOE", "APOE genotype",
            c("E2/E2", "E2/E3", "E2/E4", "E3/E3", "E3/E4", "E4/E4")),
  cat_block("CD33Group", "CD33 group", c("CV", "CD33var", "Missing")),
  cat_block("Sex", "Sex", c("F", "M")),
  cat_block("cohort", "Sequencing cohort", c("BiogenInitial", "BiogenExtra")),
  cat_block("Braak", "Braak stage",
            c("0", "I", "I,II", "II", "III", "IV", "V", "V,VI", "VI", "Missing")),
  cat_block("BrainRegion", "Brain region", c("MTG", "SSC"), region_only = TRUE)
)

# tidy display labels
relabel = c("APOE4-neg" = "ε4-negative", "APOE4-pos" = "ε4-positive",
            "CD33var" = "CD33-variant", "F" = "Female", "M" = "Male",
            "I,II" = "I–II", "V,VI" = "V–VI")
tidy_panel = function(df) df %>% mutate(Category = recode(Category, !!!relabel))

# split blocks into two balanced panels (16 rows each)
left  = bind_rows(blocks[1:5]) %>% tidy_panel()   # Diagnosis..CD33  (16 rows)
right = bind_rows(blocks[6:9]) %>% tidy_panel()   # Sex..Region      (16 rows)

# pad to equal length, then place side by side with a spacer column
pad = function(df, n) {
  if (nrow(df) < n) df = bind_rows(df, df[rep(NA, n - nrow(df)), ])
  mutate(df, across(everything(), ~ replace_na(.x, "")))
}
n_rows = max(nrow(left), nrow(right))
left   = pad(left, n_rows); right = pad(right, n_rows)

S1 = bind_cols(
  set_names(left,  c("L1", "L2", "L3", "L4")),
  tibble(SP = rep("", n_rows)),
  set_names(right, c("R1", "R2", "R3", "R4"))
)

### Table S2: continuous variables ------------------------------------------
num_summary = function(x, level, label) {
  x = suppressWarnings(as.numeric(x)); na = sum(is.na(x)); x = x[!is.na(x)]
  tibble(Variable = label, Level = level, n = length(x),
         `Mean ± SD` = sprintf("%.1f ± %.1f", mean(x), sd(x)),
         Median = sprintf("%.2f", median(x)),
         Range = sprintf("%.2f–%.2f", min(x), max(x)),
         Missing = na)
}

# region-matched amyloid plaque density (donor + region join)
plaque = read_csv(plaque_path, show_col_types = FALSE)
samp_pd = samp %>%
  left_join(plaque %>% select(all_of(donor_col), BrainRegion, TotalDensity),
            by = c(donor_col, "BrainRegion"))

S2 = bind_rows(
  num_summary(donors$Age,                "donor",  "Age (years)"),
  num_summary(donors$PostMortemInterval, "donor",  "Post-mortem interval (h)"),
  num_summary(samp_pd$TotalDensity,      "sample", "Amyloid plaque density"),
  num_summary(samp$pctAT8PositiveArea,   "sample", "AT8 % area (phospho-tau)"),
  num_summary(samp$pctPHF1PositiveArea,  "sample", "PHF1 % area (phospho-tau)"),
  num_summary(samp$pct4G8PositiveArea,   "sample", "4G8 % area (amyloid-β)")
)

### save CSVs ---------------------------------------------------------------
write_csv(S1, file.path(out_dir, paste0(script_ind, "TableS1_categorical_sidebyside.csv")))
write_csv(S2, file.path(out_dir, paste0(script_ind, "TableS2_continuous.csv")))

### render (flextable -> Word/PNG/HTML) -------------------------------------
if (requireNamespace("flextable", quietly = TRUE)) {
  library(flextable)
  hdr = c("Variable", "Category", "Donors n (%)", "Samples n")

  ft1 = flextable(S1) %>%
    set_header_labels(L1 = hdr[1], L2 = hdr[2], L3 = hdr[3], L4 = hdr[4],
                      SP = "", R1 = hdr[1], R2 = hdr[2], R3 = hdr[3], R4 = hdr[4]) %>%
    bold(part = "header") %>%
    bold(j = c("L1", "R1"), part = "body") %>%
    align(j = c("L3", "L4", "R3", "R4"), align = "center", part = "all") %>%
    width(j = "SP", width = 0.2) %>%
    border_remove() %>% hline_top(part = "header") %>%
    hline_bottom(part = "header") %>% hline_bottom(part = "body") %>%
    fontsize(size = 9, part = "all") %>% autofit() %>%
    set_caption(sprintf("Table S1. Cohort composition (categorical). %d donors, %d samples.",
                        n_don, n_samp))

  ft2 = flextable(S2) %>%
    bold(part = "header") %>% bold(j = "Variable", part = "body") %>%
    align(j = c("n", "Mean ± SD", "Median", "Range", "Missing"),
          align = "center", part = "all") %>%
    fontsize(size = 9, part = "all") %>% autofit() %>%
    set_caption("Table S2. Cohort composition (continuous).")

  save_as_docx(ft1, ft2, path = file.path(out_dir, paste0(script_ind, "cohort_tables.docx")))
  try(save_as_html(ft1, ft2, path = file.path(out_dir, paste0(script_ind, "cohort_tables.html"))), silent = TRUE)
  try(save_as_image(ft1, path = file.path(out_dir, paste0(script_ind, "TableS1.png")), res = 200), silent = TRUE)
  try(save_as_image(ft2, path = file.path(out_dir, paste0(script_ind, "TableS2.png")), res = 200), silent = TRUE)
  message("Rendered Word/HTML/PNG via flextable.")
} else {
  message("flextable not installed - CSVs written. install.packages('flextable') to render Word/PNG.")
  print(knitr::kable(S1)); print(knitr::kable(S2))
}

message("Done. Outputs in: ", out_dir)
