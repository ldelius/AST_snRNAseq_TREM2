# LD_X21: supplementary table - do the WGCNA modules survive additional
# correction for pathology?
#
# The main network (LD_F03c_v02, DEG-seeded, 7-covariate corrected) was rebuilt
# from scratch in LD_F03d_v01 after further group-protected correction for:
#   Model A  Age + Braak
#   Model B  Age + Braak + plaque density + AT8 + PHF1 + 4G8
# Each model re-ran the whole pipeline independently, so its module numbering is
# its own and carries no relation to the main one. LD_F03d never compared the
# resulting module sets, which is what this script does.
#
# MATCHING RULE: each main module is matched to the rebuilt module SHARING THE
# LARGEST NUMBER OF ITS GENES (largest intersection), as stated in the Methods -
# not by highest Jaccard. The two rules pick the same partner for 28 of the 30
# module-model pairs; they differ only for main M3 under Model A and main M10
# under Model B, both weak matches either way (J <= 0.34).
#
# Reported per module:
#   Retained  % of ALL genes in the original module found in the matched module.
#             STRICT denominator: a gene that the rebuilt network pushes into the
#             unassigned set (M0) counts as lost. NB this is a harsher measure
#             than "% of the genes still assigned to any module", which is what
#             the Results text quotes - the two differ most for modules that shed
#             many genes to M0 (e.g. M12: 71% strict vs 86% among-assigned).
#   J         Jaccard index of that match, kept as a supporting statistic since
#             it also penalises the matched module being much larger.
#
# M0 is excluded: it is WGCNA's unassigned set, not a module, so "preservation"
# of it is not meaningful.
#
# Cheap: reads three gene-membership CSVs. Nothing recomputed, no .qs loaded.

library(tidyverse)

### paths -------------------------------------------------------------------
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")

w_dir   = file.path(base, "LD_F_DESeq_pseudobulk_WGCNA")
main_csv = file.path(w_dir, "LD_F03c_v02/LD_F03c_v02_Module_genes.csv")
a_csv    = file.path(w_dir, "LD_F03d_v01/ModelA/LD_F03d_v01_ModelA_Module_genes.csv")
b_csv    = file.path(w_dir, "LD_F03d_v01/ModelB/LD_F03d_v01_ModelB_Module_genes.csv")
out_dir  = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X21_"
for (p in c(main_csv, a_csv, b_csv)) if (!file.exists(p)) stop("Missing input: ", p)
message("Using base: ", base)

### read module memberships --------------------------------------------------
# each CSV is one column per module, padded with NA; the *_TF columns are the
# transcription-factor subsets of the same modules and are not separate modules
read_modules = function(p) {
  x = read_csv(p, show_col_types = FALSE)
  x = x[, !grepl("_TF$", names(x)), drop = FALSE]
  m = map(x, ~ .x[!is.na(.x) & .x != ""])
  m[order(as.integer(sub("^M", "", names(m))))]
}
main = read_modules(main_csv)
mods = list(A = read_modules(a_csv), B = read_modules(b_csv))
message("  modules: main ", length(main), " | Model A ", length(mods$A),
        " | Model B ", length(mods$B), " (each including M0)")

### best match per main module ----------------------------------------------
jacc = function(a, b) length(intersect(a, b)) / length(union(a, b))

best_match = function(genes, other) {
  other = other[names(other) != "M0"]                 # never match onto unassigned
  ov = map_int(other, ~ length(intersect(genes, .x))) # Methods rule: largest overlap
  i  = which.max(ov)
  j  = map_dbl(other, ~ jacc(genes, .x))
  tibble(match = names(other)[i],
         J = j[i],                                    # Jaccard OF THE MATCH
         retained = ov[i] / length(genes),            # strict denominator
         jacc_match = names(other)[which.max(j)])     # what Jaccard would have picked
}

main_mods = main[names(main) != "M0"]
res = imap_dfr(main_mods, function(g, m) {
  a = best_match(g, mods$A); b = best_match(g, mods$B)
  tibble(Module = m, N = length(g),
         A_match = a$match, A_retained = a$retained, A_J = a$J, A_jacc_match = a$jacc_match,
         B_match = b$match, B_retained = b$retained, B_J = b$J, B_jacc_match = b$jacc_match)
}) %>% arrange(as.integer(sub("^M", "", Module)))

# report where the two matching rules disagree
diff_a = res %>% filter(A_match != A_jacc_match)
diff_b = res %>% filter(B_match != B_jacc_match)
if (nrow(diff_a) || nrow(diff_b)) {
  message("Matching rules disagree (overlap vs Jaccard):")
  if (nrow(diff_a)) for (i in seq_len(nrow(diff_a)))
    message("  Model A, main ", diff_a$Module[i], ": overlap -> ", diff_a$A_match[i],
            " | Jaccard -> ", diff_a$A_jacc_match[i])
  if (nrow(diff_b)) for (i in seq_len(nrow(diff_b)))
    message("  Model B, main ", diff_b$Module[i], ": overlap -> ", diff_b$B_match[i],
            " | Jaccard -> ", diff_b$B_jacc_match[i])
} else message("Overlap and Jaccard matching agree for every module.")

fmt_pct = function(x) sprintf("%.0f%%", 100 * x)
fmt_j   = function(x) sprintf("%.2f", x)

tab = res %>%
  transmute(Module, `Genes` = N,
            `Best match` = A_match, `Retained` = fmt_pct(A_retained), `J` = fmt_j(A_J),
            `Best match ` = B_match, `Retained ` = fmt_pct(B_retained), `J ` = fmt_j(B_J))

write_csv(res, file.path(out_dir, paste0(script_ind, "pathology_robustness_long.csv")))
write_csv(tab, file.path(out_dir, paste0(script_ind, "pathology_robustness_table.csv")))

### render -------------------------------------------------------------------
cap = paste0("Supplementary Table 10. Stability of the co-expression modules under ",
             "additional correction for pathology. The network was rebuilt from scratch ",
             "after further correction for Age and Braak stage (Model A) and for Age, ",
             "Braak stage, amyloid plaque density and AT8, PHF1 and 4G8 percentage area ",
             "(Model B). Module numbering is independent in each run, so each module of ",
             "the main network is matched to the rebuilt module sharing the largest number ",
             "of its genes. Retained is the percentage of all genes in the original module ",
             "found in that matched module; genes reassigned to the unassigned set count as ",
             "lost. J is the Jaccard index of the match. The unassigned set (M0) is not a ",
             "module and is excluded.")

if (!requireNamespace("flextable", quietly = TRUE)) {
  message("flextable not installed - CSVs written only.")
} else {
  library(flextable)
  nm = names(tab)
  top = c("", "", rep("Model A: + Age + Braak", 3),
          rep("Model B: + Age + Braak + plaque + tau/amyloid area", 3))
  ft = flextable(tab) %>%
    set_header_df(mapping = data.frame(keys = nm, top = top, bot = nm), key = "keys") %>%
    merge_h(part = "header") %>%
    bold(part = "header") %>% italic(i = 1, part = "header") %>%
    bold(j = "Module", part = "body") %>%
    align(j = nm[-1], align = "center", part = "all") %>%
    border_remove() %>%
    hline_top(part = "header", border = fp_border_default(width = 1.5)) %>%
    hline(i = 1, part = "header", border = fp_border_default(width = 0.75)) %>%
    hline_bottom(part = "header", border = fp_border_default(width = 1.5)) %>%
    hline_bottom(part = "body", border = fp_border_default(width = 1.5)) %>%
    fontsize(size = 9, part = "all") %>%
    set_caption(cap) %>% autofit()

  ok = tryCatch({ save_as_docx(ft, path = file.path(out_dir, paste0(script_ind, "pathology_robustness_table.docx"))); TRUE },
                error = function(e) { message("  FAILED docx: ", conditionMessage(e)); FALSE })
  if (ok) message("  wrote pathology_robustness_table.docx")
  try(save_as_image(ft, path = file.path(out_dir, paste0(script_ind, "pathology_robustness_table.png")), res = 300),
      silent = TRUE)
}

well_preserved = sum(res$A_J >= 0.5 & res$B_J >= 0.5)
message("Done. ", nrow(res), " modules; ", well_preserved,
        " with J >= 0.5 in both models. Outputs in: ", out_dir)
