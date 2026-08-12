# LD_X13: Supplementary table of software packages and versions.
#
# TWO STAGES, because the analysis env and the table-rendering env differ:
#   1. On the HPC (conda R43_240426, the environment the thesis analyses ran in):
#      queries the installed version of every package and writes/updates a CSV.
#   2. On the laptop: fills in anything still missing (flextable, which is only
#      installed locally) and renders the Word/PNG table.
# Versions are only ever FILLED IN, never overwritten, so the HPC values win for
# every analysis package. Run on the HPC first, then locally.
#
# Rationale: package versions must be read from the live library, not from old job
# logs - the conda env has been updated over time, so different logs record
# different versions of the same package.

library(tidyverse)

### paths -------------------------------------------------------------------
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")
out_dir = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X13_"
cache = file.path(out_dir, paste0(script_ind, "package_versions.csv"))
message("Using base: ", base, " | ", R.version.string)

### the packages to report --------------------------------------------------
# Order is alphabetical, matching how the table is read rather than load order.
pkgs = tribble(
  ~Package,             ~Purpose,
  "AnnotationDbi",      "Annotation backend",
  "BiocParallel",       "Parallelisation",
  "circlize",           "Colour scales for heatmaps",
  "clusterProfiler",    "GO over-representation",
  "colorRamps",         "Colour palettes",
  "ComplexHeatmap",     "Heatmaps",
  "DESeq2",             "Pseudobulk differential expression, variance stabilisation",
  "DOSE",               "clusterProfiler dependency",
  "enrichplot",         "Enrichment visualisation",
  "fgsea",              "Pre-ranked gene set enrichment",
  "flextable",          "Rendering of supplementary tables",
  "ggrepel",            "Plot labelling",
  "harmony",            "Batch integration across samples",
  "limma",              "Covariate removal from expression matrices",
  "lme4",               "Mixed-model backend",
  "lmerTest",           "Linear mixed models",
  "matrixStats",        "Matrix utilities",
  "msigdbr",            "MSigDB Hallmark gene sets",
  "org.Hs.eg.db",       "GO annotations",
  "patchwork",          "Figure composition",
  "pheatmap",           "Heatmaps",
  "presto",             "Fast Wilcoxon implementation used by Seurat",
  "qs",                 "Object serialisation",
  "sccomp",             "Differential abundance",
  "Seurat",             "Single-nucleus handling, normalisation, clustering, marker testing, module scores",
  "tidyverse",          "Data handling",
  "variancePartition",  "Covariate variance decomposition",
  "viridis",            "Colour palettes",
  "WGCNA",              "Co-expression network construction"
)

# Packages that ARE used somewhere in the scripts but are not in the list above.
# Set include_extras = TRUE to add them to the table; leave FALSE to keep the
# table restricted to what the thesis text actually describes.
include_extras = FALSE
extras = tribble(
  ~Package,        ~Purpose,
  "CellChat",      "Cell-cell communication analysis",
  "future",        "Parallelisation backend for Seurat",
  "ggrastr",       "Rasterisation of large scatter layers",
  "scales",        "Palette and axis scale helpers",
  "ggsignif",      "Significance brackets on plots",
  "SeuratObject",  "Data structure underlying Seurat"
)
if (include_extras) pkgs = bind_rows(pkgs, extras) %>% arrange(tolower(Package))

### stage 1: fill in versions from the CURRENT library -----------------------
ver_of = function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)

prev = if (file.exists(cache)) read_csv(cache, show_col_types = FALSE) else
  tibble(Package = character(), Version = character(), Source = character())

# The HPC env is the one the thesis analyses ran in, so it is authoritative: a run
# there OVERWRITES whatever is on record. A laptop run only fills gaps (in practice
# just flextable), so it can never silently replace an analysis package's version
# with the laptop's copy. Run order then does not matter.
on_hpc = startsWith(base, "/rds")
message(if (on_hpc) "HPC run: versions found here overwrite the cache."
        else        "Local run: only filling gaps, HPC values kept.")

tab = pkgs %>%
  left_join(prev %>% select(Package, Version, Source), by = "Package") %>%
  mutate(found   = map_chr(Package, ver_of),
         use_new = !is.na(found) & (on_hpc | is.na(Version)),
         Source  = if_else(use_new, R.version.string, Source),
         Version = if_else(use_new, found, Version)) %>%
  select(Package, Version, Purpose, Source)

write_csv(tab, cache)
message("Cache updated: ", cache)

missing = tab$Package[is.na(tab$Version)]
if (length(missing))
  message("Still missing (not installed in this R): ", paste(missing, collapse = ", "),
          "\n  -> run this script in the other environment to fill them in.")

# MSigDB data release behind msigdbr, if the installed version exposes it
msig_release = tryCatch(as.character(msigdbr::msigdbr_version()), error = function(e) NA_character_)

### stage 2: render (needs flextable, i.e. the laptop) ----------------------
# The analysis environment is the one to cite: report the R version that most
# packages came from (the HPC env), not every environment that touched the table.
# flextable is recorded from the laptop but is a rendering tool, not an analysis one.
r_analysis = names(sort(table(na.omit(tab$Source)), decreasing = TRUE))[1]
cap = paste0("Supplementary Table 5. R packages used in this thesis, with the ",
             "version of each. Analyses were run under ", r_analysis, ".",
             if (!is.na(msig_release)) paste0(" MSigDB data release ", msig_release, ".") else "")

if (!requireNamespace("flextable", quietly = TRUE)) {
  message("flextable not available here - CSV written; render on the laptop.")
} else {
  library(flextable)
  ft = tab %>%
    select(Package, Version, Purpose) %>%
    mutate(Version = replace_na(Version, "not recorded")) %>%
    flextable() %>%
    bold(part = "header") %>% italic(part = "header") %>%
    italic(j = "Package", part = "body") %>%   # package names as code/proper nouns
    border_remove() %>%
    hline_top(part = "header", border = fp_border_default(width = 1.5)) %>%
    hline_bottom(part = "header", border = fp_border_default(width = 1.5)) %>%
    hline_bottom(part = "body", border = fp_border_default(width = 1.5)) %>%
    fontsize(size = 9, part = "all") %>%
    set_caption(cap) %>%
    width(j = c("Package", "Version", "Purpose"), width = c(1.6, 1.0, 4.2))

  out_docx = file.path(out_dir, paste0(script_ind, "package_versions.docx"))
  ok = tryCatch({ save_as_docx(ft, path = out_docx); TRUE },
                error = function(e) { message("  FAILED docx: ", conditionMessage(e)); FALSE })
  if (ok) message("  wrote ", basename(out_docx))
  try(save_as_image(ft, path = file.path(out_dir, paste0(script_ind, "package_versions.png")),
                    res = 300), silent = TRUE)
}

message("Done. ", sum(!is.na(tab$Version)), "/", nrow(tab), " versions recorded.")
