# LD_X17b: Recompute LD_E04c concordance under shared, nominal-significance and
# adjusted-significance gene-selection rules without refitting models.

library(tidyverse)

### paths -------------------------------------------------------------------
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")
e_out    = file.path(base, "LD_E_DESeq_pseudobulk")
ckpt_dir = file.path(e_out, "LD_E04c_ckpt")
deg_csv  = file.path(e_out, "LD_E02c/LD_E02c_v01_DEGs_by_cluster_genes.csv")
out_csv  = file.path(e_out, "LD_E04c_effect_robustness_summary_ALLMEASURES.csv")
if (!dir.exists(ckpt_dir)) stop("Missing checkpoints: ", ckpt_dir)
if (!file.exists(deg_csv)) stop("Missing input: ", deg_csv)

### settings -----------------------------------------------------------------
PADJ_CUT = 0.10    # for r_sig_both (matches the LD_X07 figure)
PNOM_CUT = 0.05    # for r_nom_both (nominal, unadjusted)
MIN_N    = 3       # fewer shared genes than this -> r is not computed

# the three pairs of the concordance figure (y, x)
pairs = list(
  P_R62H_vs_AD   = c("R62H_vs_CV", "CV_AD_vs_Control"),
  P_R47H_vs_AD   = c("R47H_vs_CV", "CV_AD_vs_Control"),
  P_R62H_vs_R47H = c("R62H_vs_CV", "R47H_vs_CV"))
levels_keep = c("M0_base", "M1_Age", "M2_Age_PMI", "M3_Age_PMI_Braak")

### gene universe (same as LD_E04c) -----------------------------------------
deg_universe = read_csv(deg_csv, show_col_types = FALSE) %>% unlist(use.names = FALSE) %>% unique()
deg_universe = deg_universe[!is.na(deg_universe) & deg_universe != ""]
message("DEG universe: ", length(deg_universe), " genes")

### read checkpoints ---------------------------------------------------------
# filenames are <scope>__<level>__<contrast>.rds, written by LD_E04c's ckpt_file()
files = list.files(ckpt_dir, pattern = "rds$", full.names = TRUE)
message("Checkpoints: ", length(files))
key_of = function(f) sub("\\.rds$", "", basename(f))

res = new.env(parent = emptyenv())
for (f in files) {
  o = readRDS(f)
  if (is.null(o$res)) next
  r = as.data.frame(o$res)
  g = if ("gene" %in% names(r)) as.character(r$gene) else rownames(r)
  keep = g %in% deg_universe
  assign(paste(o$scope, o$level, o$contrast, sep = "|"),
         tibble(gene = g[keep], log2FC = r$log2FoldChange[keep],
                padj = r$padj[keep], pval = r$pvalue[keep]),
         envir = res)
}
keys = ls(res)
scopes = unique(sub("\\|.*$", "", keys))
message("Scopes: ", length(scopes))

### compute the three measures ----------------------------------------------
r_of = function(x, y) if (length(x) >= MIN_N) suppressWarnings(cor(x, y)) else NA_real_

out = NULL
for (sc in scopes) for (pn in names(pairs)) for (lv in levels_keep) {
  ry = res[[paste(sc, lv, pairs[[pn]][1], sep = "|")]]
  rx = res[[paste(sc, lv, pairs[[pn]][2], sep = "|")]]
  if (is.null(ry) || is.null(rx)) next
  d = inner_join(ry, rx, by = "gene", suffix = c("_y", "_x")) %>%
    filter(!is.na(log2FC_y), !is.na(log2FC_x))
  d_sig = d %>% filter(!is.na(padj_y), !is.na(padj_x), padj_y < PADJ_CUT, padj_x < PADJ_CUT)
  d_nom = d %>% filter(!is.na(pval_y), !is.na(pval_x), pval_y < PNOM_CUT, pval_x < PNOM_CUT)
  out = bind_rows(out, tibble(
    scope = sc, pair = pn, level = lv,
    n_shared   = nrow(d),     r_shared   = r_of(d$log2FC_x,     d$log2FC_y),
    n_sig_both = nrow(d_sig), r_sig_both = r_of(d_sig$log2FC_x, d_sig$log2FC_y),
    n_nom_both = nrow(d_nom), r_nom_both = r_of(d_nom$log2FC_x, d_nom$log2FC_y)))
}

out = out %>% mutate(level = factor(level, levels = levels_keep)) %>% arrange(scope, pair, level)
write_csv(out, out_csv)

cov = out %>% summarise(cells = n(),
                        shared = sum(!is.na(r_shared)),
                        nominal = sum(!is.na(r_nom_both)),
                        sig_both = sum(!is.na(r_sig_both)))
message("Coverage of ", cov$cells, " cells:  r_shared ", cov$shared,
        " | r_nom_both ", cov$nominal, " | r_sig_both ", cov$sig_both)
message("Median genes per cell:  shared ", median(out$n_shared),
        " | nominal ", median(out$n_nom_both), " | sig_both ", median(out$n_sig_both))
message("Wrote: ", out_csv)
