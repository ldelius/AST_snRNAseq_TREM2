# LD_X01: Cohort summary tables for thesis supplementary.
#   Table S1 - demographics, donor level, TREM2 genotype (columns) x diagnosis (sub-rows).
#   Table S2 - neuropathology, sample level, same layout.
# Layout follows the published cohort-table style: variables as labelled rows, each split
# into Control / AD sub-rows, one column per TREM2 genotype plus an "All" column.
# Cells are raw counts (n or n/N) and mean +/- SD; no percentages.
# Numbers are recomputed from source metadata so the tables always match the data.
# Outputs: CSVs (always) + a Word doc and PNG/HTML (if flextable is installed).

library(tidyverse)

### paths -------------------------------------------------------------------
# Canonical location is RDS on the HPC. This script only reads two small metadata
# files, so it can also be run locally against the mounted RDS share (same files,
# same output folder) - the first existing path wins.
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")
message("Using base: ", base)
in_dir      = file.path(base, "data_TREM2_michael/A_input")
samp_path   = file.path(in_dir, "TREM2_Samplesheet_snRNAseq_GliaEnriched_with_cohort_info.tsv")
plaque_path = file.path(in_dir, "TREM2_plaque_data_Sam.csv")
out_dir     = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind  = "LD_X01_"
message("Writing outputs to: ", out_dir)

### load --------------------------------------------------------------------
samp = read_tsv(samp_path, show_col_types = FALSE)
donor_col = "BrainBankNetworkIDFormatted"

# Braak stage as a number (comma-coded ranges take the midpoint, e.g. "V,VI" -> 5.5)
braak_num = c("0" = 0, "I" = 1, "I,II" = 1.5, "II" = 2, "III" = 3,
              "IV" = 4, "V" = 5, "V,VI" = 5.5, "VI" = 6)
samp = samp %>% mutate(BraakNum = unname(braak_num[as.character(Braak)]))

# region-matched amyloid plaque density (donor + region join), sample level
plaque = read_csv(plaque_path, show_col_types = FALSE)
samp_pd = samp %>%
  left_join(plaque %>% select(all_of(donor_col), BrainRegion, TotalDensity),
            by = c(donor_col, "BrainRegion"))

# donor-level table: one row per individual
donors = samp %>% distinct(.data[[donor_col]], .keep_all = TRUE)
n_don  = nrow(donors)            # 70
n_samp = nrow(samp)              # 117

### cell formatters ---------------------------------------------------------
# each takes the subset of rows for one genotype x diagnosis cell, returns a string
f_n       = function(d) as.character(nrow(d))
f_mean_sd = function(col) function(d) {
  x = suppressWarnings(as.numeric(d[[col]])); x = x[!is.na(x)]
  if (length(x) == 0) "—" else if (length(x) == 1) sprintf("%.1f", x)
  else sprintf("%.1f ± %.1f", mean(x), sd(x))
}
f_frac    = function(col, val) function(d)                       # carriers out of group
  sprintf("%d/%d", sum(d[[col]] %in% val), nrow(d))
f_split   = function(col, a, b) function(d)                      # two counts, e.g. F / M
  sprintf("%d / %d", sum(d[[col]] == a, na.rm = TRUE), sum(d[[col]] == b, na.rm = TRUE))
f_region  = function(d) sprintf("%d (%d/%d)", nrow(d),
                                sum(d$BrainRegion == "MTG"), sum(d$BrainRegion == "SSC"))

### block builder -----------------------------------------------------------
# one label row (variable name, empty cells) + one row per diagnosis group
variants = c("CV", "R47H", "R62H")
diagnoses = c("Control", "AD")

make_block = function(data, label, f) {
  rows = map_dfr(diagnoses, function(g) {
    cells = map_chr(variants, function(v)
      f(filter(data, TREM2Variant == v, NeuropathologicalDiagnosis == g)))
    tibble(Row = g, CV = cells[1], R47H = cells[2], R62H = cells[3],
           All = f(filter(data, NeuropathologicalDiagnosis == g)))
  })
  bind_rows(tibble(Row = label, CV = "", R47H = "", R62H = "", All = ""), rows)
}

### Table S1: demographics (donor level) ------------------------------------
S1 = bind_rows(
  make_block(donors,  "Number of donors",              f_n),
  make_block(samp,    "Number of samples (MTG/SSC)",   f_region),
  make_block(donors,  "Age at death (years)",          f_mean_sd("Age")),
  make_block(donors,  "Sex (F / M)",                   f_split("Sex", "F", "M")),
  make_block(donors,  "APOE ε4-positive",         f_frac("APOEgroup", "APOE4-pos")),
  make_block(donors,  "APOE ε4/ε4 homozygous", f_frac("APOE", "E4/E4")),
  make_block(donors,  "CD33 variant carriers",         f_frac("CD33Group", "CD33var")),
  make_block(donors,  "Sequencing cohort (initial / extra)",
             f_split("cohort", "BiogenInitial", "BiogenExtra")),
  make_block(donors,  "Post-mortem interval (h)",      f_mean_sd("PostMortemInterval")),
  make_block(donors,  "Braak stage",                   f_mean_sd("BraakNum"))
)
S1_labels = which(S1$CV == "" & S1$All == "")   # label rows, for bold/indent styling

### Table S2: neuropathology (sample level) ---------------------------------
S2 = bind_rows(
  make_block(samp_pd, "Amyloid plaque density",        f_mean_sd("TotalDensity")),
  make_block(samp,    "4G8 % area (amyloid-β)",   f_mean_sd("pct4G8PositiveArea")),
  make_block(samp,    "AT8 % area (phospho-tau)",      f_mean_sd("pctAT8PositiveArea")),
  make_block(samp,    "PHF1 % area (phospho-tau)",     f_mean_sd("pctPHF1PositiveArea"))
)
S2_labels = which(S2$CV == "" & S2$All == "")

### missingness (reported as footnotes, not as table rows) ------------------
miss_note = function(d, cols) {
  n = map_int(cols, ~ sum(is.na(d[[.x]]) | d[[.x]] == "NA"))
  paste(sprintf("%s: %d missing", cols, n)[n > 0], collapse = "; ")
}
note1 = miss_note(donors, c("Braak", "PostMortemInterval", "CD33Group", "APOE"))
note2 = miss_note(samp_pd, c("TotalDensity", "pct4G8PositiveArea",
                             "pctAT8PositiveArea", "pctPHF1PositiveArea"))

### Table S3: Braak stage distribution (donor level) ------------------------
# The mean in Table S1 collapses an ordinal scale; this shows which stages are
# actually represented in each genotype x diagnosis group. Counts of donors.
braak_levels  = c("0", "I", "I,II", "II", "III", "IV", "V", "V,VI", "VI")
braak_display = c("0", "I", "I–II", "II", "III", "IV", "V", "V–VI", "VI")
f_stage   = function(lv) function(d) as.character(sum(d$Braak == lv, na.rm = TRUE))
f_missing = function(d) as.character(sum(is.na(d$Braak) | d$Braak == "NA"))

# 10 blocks of 3 rows would run to 30 rows, so fold into two panels of 5 blocks
# shown side by side (15 rows). Both halves have label rows at the same positions.
blocks3 = c(map2(braak_levels, braak_display,
                 ~ make_block(donors, paste("Braak stage", .y), f_stage(.x))),
            list(make_block(donors, "Not recorded", f_missing)))

left3  = bind_rows(blocks3[1:5])    # stages 0 .. III
right3 = bind_rows(blocks3[6:10])   # stages IV .. VI, plus Not recorded
stopifnot(nrow(left3) == nrow(right3))

S3 = bind_cols(
  set_names(left3, paste0("L_", names(left3))),
  tibble(SP = rep("", nrow(left3))),
  set_names(right3, paste0("R_", names(right3)))
)
S3_labels = which(S3$L_CV == "" & S3$L_All == "")

### captions ----------------------------------------------------------------
cap1 = sprintf(paste("Supplementary Table 1. Demographic characteristics across TREM2",
                     "genotypes and disease groups (%d donors, %d samples). Values are",
                     "counts (n, or n/N of the group) or mean ± SD. Braak stage was",
                     "converted to a numeric scale, with comma-coded stages taking the",
                     "midpoint (I,II = 1.5; V,VI = 5.5)."), n_don, n_samp)
cap2 = sprintf(paste("Supplementary Table 2. Neuropathological findings across TREM2",
                     "genotypes and disease groups, quantified per tissue section",
                     "(n = %d samples). Values are mean ± SD."), n_samp)
cap3 = sprintf(paste("Supplementary Table 3. Distribution of Braak stages across TREM2",
                     "genotypes and disease groups (%d donors). Values are donor counts.",
                     "Stages recorded as a range in the source records are kept as such",
                     "(I–II, V–VI) rather than assigned to a single stage."), n_don)

### save CSVs ---------------------------------------------------------------
write_csv(S1, file.path(out_dir, paste0(script_ind, "TableS1_demographics.csv")))
write_csv(S2, file.path(out_dir, paste0(script_ind, "TableS2_neuropathology.csv")))
write_csv(S3, file.path(out_dir, paste0(script_ind, "TableS3_braak_distribution.csv")))

### render: styled HTML (no extra packages, opens in Word) ------------------
# Base-R writer so a thesis-ready table always exists, even without flextable.
# Open the .html in Word (or copy-paste into the document) - it pastes as a real table.
esc = function(x) { x = gsub("&", "&amp;", x); x = gsub("<", "&lt;", x); gsub(">", "&gt;", x) }

# row_cols = column indices holding row labels (2 of them for a two-panel table);
# spacer = index of the blank separator column, if any.
html_table = function(df, label_rows, caption,
                      hdr = c("TREM2 genotype", "CV", "R47H", "R62H", "All"),
                      row_cols = 1, spacer = integer(0)) {
  hcls = ifelse(seq_along(hdr) %in% row_cols, "v", ifelse(seq_along(hdr) %in% spacer, "sp", "c"))
  th   = paste0("<th class='", hcls, "'>", esc(hdr), "</th>", collapse = "")
  body = map_chr(seq_len(nrow(df)), function(i) {
    lab_cls = if (i %in% label_rows) "lab" else "sub"
    cells = map_chr(seq_along(df), function(j) {
      cls = if (j %in% row_cols) lab_cls else if (j %in% spacer) "sp" else "c"
      paste0("<td class='", cls, "'>", esc(as.character(df[[j]][i])), "</td>")
    })
    paste0("<tr>", paste(cells, collapse = ""), "</tr>")
  })
  paste0("<p class='cap'>", esc(caption), "</p>\n<table>\n<thead><tr>", th,
         "</tr></thead>\n<tbody>\n", paste(body, collapse = "\n"), "\n</tbody>\n</table>")
}

css = paste(
  "body{font-family:'Times New Roman',Georgia,serif;font-size:10pt;margin:24px;}",
  "table{border-collapse:collapse;margin:0 0 28px 0;}",
  "th,td{padding:3px 14px 3px 4px;text-align:left;vertical-align:top;}",
  "thead th{border-top:1.5px solid #000;border-bottom:1.5px solid #000;font-style:italic;font-weight:bold;}",
  "tbody tr:last-child td{border-bottom:1.5px solid #000;}",
  "td.lab{font-weight:bold;}", "td.sub{padding-left:22px;}",
  "th.sp,td.sp{padding:0 10px;border:none;}",
  "p.cap{font-size:9pt;margin:0 0 6px 0;max-width:760px;}",
  "p.note{font-size:8pt;font-style:italic;margin:-22px 0 28px 0;}", sep = "\n")

html_doc = paste0(
  "<!DOCTYPE html>\n<html><head><meta charset='utf-8'>\n<style>\n", css, "\n</style>\n",
  "</head><body>\n",
  html_table(S1, S1_labels, cap1), if (nzchar(note1)) paste0("<p class='note'>", esc(note1), "</p>") else "",
  "\n", html_table(S2, S2_labels, cap2), if (nzchar(note2)) paste0("<p class='note'>", esc(note2), "</p>") else "",
  "\n", html_table(S3, S3_labels, cap3,
                   hdr = c("TREM2 genotype", "CV", "R47H", "R62H", "All", "",
                           "TREM2 genotype", "CV", "R47H", "R62H", "All"),
                   row_cols = c(1, 7), spacer = 6),
  "\n</body></html>\n")

html_path = file.path(out_dir, paste0(script_ind, "cohort_tables.html"))
con = file(html_path, open = "w", encoding = "UTF-8")
writeLines(html_doc, con, useBytes = FALSE); close(con)
message("Wrote styled HTML (open in Word): ", html_path)

### render (flextable -> Word/PNG) ------------------------------------------
# flextable is the primary renderer for the thesis tables. Failures are reported
# with their actual error message rather than swallowed, so a broken render is
# obvious in the console instead of looking like "only CSVs were written".
if (!requireNamespace("flextable", quietly = TRUE)) {
  stop("flextable is not available to this R session. Check .libPaths() / the ",
       "loaded R module, or install it. The styled HTML above is already written ",
       "and can be used in the meantime.")
} else {
  library(flextable)

  say = function(label, expr) {
    ok = tryCatch({ force(expr); TRUE },
                  error = function(e) { message("  FAILED ", label, ": ", conditionMessage(e)); FALSE })
    if (ok) message("  wrote ", label)
    invisible(ok)
  }

  style_tab = function(df, label_rows, caption, note) {
    ft = flextable(df) %>%
      set_header_labels(Row = "TREM2 genotype", CV = "CV", R47H = "R47H",
                        R62H = "R62H", All = "All") %>%
      bold(part = "header") %>%
      italic(j = "Row", part = "header") %>%
      italic(j = c("CV", "R47H", "R62H", "All"), part = "header") %>%
      bold(i = label_rows, j = "Row", part = "body") %>%
      padding(i = setdiff(seq_len(nrow(df)), label_rows), j = "Row",
              padding.left = 18, part = "body") %>%
      align(j = c("CV", "R47H", "R62H", "All"), align = "left", part = "all") %>%
      border_remove() %>%
      hline_top(part = "header", border = fp_border_default(width = 1.5)) %>%
      hline_bottom(part = "header", border = fp_border_default(width = 1.5)) %>%
      hline_bottom(part = "body", border = fp_border_default(width = 1.5)) %>%
      fontsize(size = 9, part = "all") %>%
      set_caption(caption) %>%
      autofit()
    if (nzchar(note)) ft = add_footer_lines(ft, note) %>% fontsize(size = 8, part = "footer")
    ft
  }

  # two-panel variant for S3: same look, headers repeated either side of a spacer
  style_tab_wide = function(df, label_rows, caption) {
    lab_l = "L_Row"; lab_r = "R_Row"
    num_cols = setdiff(names(df), c(lab_l, lab_r, "SP"))
    hdr = set_names(as.list(c("TREM2 genotype", "CV", "R47H", "R62H", "All", "",
                              "TREM2 genotype", "CV", "R47H", "R62H", "All")), names(df))
    flextable(df) %>%
      set_header_labels(values = hdr) %>%
      bold(part = "header") %>% italic(part = "header") %>%
      bold(i = label_rows, j = c(lab_l, lab_r), part = "body") %>%
      padding(i = setdiff(seq_len(nrow(df)), label_rows), j = c(lab_l, lab_r),
              padding.left = 18, part = "body") %>%
      align(j = num_cols, align = "left", part = "all") %>%
      border_remove() %>%
      hline_top(part = "header", border = fp_border_default(width = 1.5)) %>%
      hline_bottom(part = "header", border = fp_border_default(width = 1.5)) %>%
      hline_bottom(part = "body", border = fp_border_default(width = 1.5)) %>%
      # the spacer column carries no rules, so the two panels read as separate
      border(j = "SP", border = fp_border_default(width = 0), part = "all") %>%
      width(j = "SP", width = 0.25) %>%
      fontsize(size = 9, part = "all") %>%
      set_caption(caption) %>%
      autofit()
  }

  ft1 = style_tab(S1, S1_labels, cap1, note1)
  ft2 = style_tab(S2, S2_labels, cap2, note2)
  ft3 = style_tab_wide(S3, S3_labels, cap3)

  message("flextable ", as.character(packageVersion("flextable")), " - rendering:")
  say("cohort_tables.docx",
      save_as_docx(ft1, ft2, ft3, path = file.path(out_dir, paste0(script_ind, "cohort_tables.docx"))))
  # PNGs need webshot2/chromote (headless Chrome); the docx is the thesis deliverable
  say("TableS1.png",
      save_as_image(ft1, path = file.path(out_dir, paste0(script_ind, "TableS1.png")), res = 300))
  say("TableS2.png",
      save_as_image(ft2, path = file.path(out_dir, paste0(script_ind, "TableS2.png")), res = 300))
  say("TableS3.png",
      save_as_image(ft3, path = file.path(out_dir, paste0(script_ind, "TableS3.png")), res = 300))
}

message("Done. Outputs in: ", out_dir)
