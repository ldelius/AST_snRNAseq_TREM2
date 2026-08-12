message("\n\n##########################################################################\n",
        "# Start LD_F03e v01: WGCNA hub genes and donor-aware LMMs ", Sys.time(),
        "\n##########################################################################\n\n")

# Donor and sample random intercepts account for paired brain regions within
# donors. TREM2 associations remain descriptive because modules are DEG-seeded.
library(qs)
library(tidyverse)
library(WGCNA)
library(lme4)
library(lmerTest)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

in_dir  = paste0(main_dir, "LD_F_DESeq_pseudobulk_WGCNA/LD_F03c_v02/")
out_dir = paste0(main_dir, "LD_F_DESeq_pseudobulk_WGCNA/LD_F03e_v01/")
dir.create(out_dir, showWarnings = FALSE)

script_ind = "LD_F03e_v01_"


### load F03c v02 checkpoint

bulk_data = qread(file = paste0(in_dir, "LD_F03c_v02_bulk_data.qs"))

mod_gene_tab = bulk_data$wgcna$mod_gene_tab
mods         = unique(mod_gene_tab$module)
me_mat       = bulk_data$wgcna$mod_eigengene_mat     # modules x cluster_samples
meta         = bulk_data$meta
meta_lmm     = meta[match(colnames(me_mat), meta$cluster_sample), ]


### kME and hub genes --------------------------------------------------------

message("\n\n   *Compute kME and hub-gene ranking \n")

input_mat = bulk_data$wgcna$input_mat        # cluster_samples x genes
MEs       = bulk_data$wgcna$network$MEs      # cluster_samples x modules (ME0, ME1, ...)

kME_df = as.data.frame(WGCNA::signedKME(input_mat, MEs, outputColumnName = "kME"))
colnames(kME_df) = gsub("^kME", "kME_M", colnames(kME_df))
kME_df$gene = rownames(kME_df)

kME_mat_full = as.matrix(kME_df[, setdiff(colnames(kME_df), "gene")])
rownames(kME_mat_full) = kME_df$gene

write_csv(cbind(gene = rownames(kME_mat_full), as.data.frame(kME_mat_full)),
         file = paste0(out_dir, script_ind, "Module_membership_kME_full_matrix.csv"))

# each gene's kME to its OWN assigned module
own_col = paste0("kME_", mod_gene_tab$module)
mod_gene_tab$kME_own = kME_mat_full[cbind(match(mod_gene_tab$gene, rownames(kME_mat_full)),
                                          match(own_col, colnames(kME_mat_full)))]

write_csv(mod_gene_tab,
         file = paste0(out_dir, script_ind, "Module_membership_kME.csv"))

# top-10 hub genes per module by own-module kME (M0 = unassigned/grey, excluded)
hub_genes = mod_gene_tab %>%
  filter(module != "M0") %>%
  group_by(module) %>%
  arrange(desc(kME_own), .by_group = TRUE) %>%
  mutate(hub_rank = row_number()) %>%
  filter(hub_rank <= 10) %>%
  ungroup() %>%
  select(module, hub_rank, gene, kME_own, colors)

write_csv(hub_genes,
         file = paste0(out_dir, script_ind, "Module_hub_genes_top10.csv"))


### donor-aware pathology LMMs -----------------------------------------------
# Sample identifies donor-region combinations, so donor and sample intercepts
# account for both cross-region donor correlation and sample-level variation.

message("\n\n   *Donor-aware LMM per module x pathology trait \n")

path_traits = c("plaque_dens", "pct4G8PositiveArea", "pctAT8PositiveArea",
                "pctPHF1PositiveArea", "Braak_numeric", "Age", "PostMortemInterval")

lmm_path_list = list()

for (mod1 in rownames(me_mat)){
  for (tr1 in path_traits){

    df1 = data.frame(
      eigengene    = as.numeric(me_mat[mod1, ]),
      trait        = as.numeric(meta_lmm[[tr1]]),
      cluster_name = factor(meta_lmm$cluster_name),
      donor        = factor(meta_lmm$BrainBankNetworkIDFormatted),
      sample       = factor(meta_lmm$sample),
      stringsAsFactors = FALSE
    )
    df1 = df1[complete.cases(df1), ]
    df1$cluster_name = droplevels(df1$cluster_name)
    df1$donor        = droplevels(df1$donor)
    df1$sample       = droplevels(df1$sample)

    if (nrow(df1) < 10 || nlevels(df1$sample) < 3 || nlevels(df1$donor) < 3 || sd(df1$trait) == 0) next

    fit = tryCatch(
      lmerTest::lmer(eigengene ~ trait + (1 | cluster_name) + (1 | donor) + (1 | sample),
                     data = df1, REML = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    sm = tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(sm) || !"trait" %in% rownames(sm)) next

    lmm_path_list[[paste(mod1, tr1)]] = data.frame(
      module     = mod1,
      trait      = tr1,
      n_obs      = nrow(df1),
      n_clusters = nlevels(df1$cluster_name),
      n_donors   = nlevels(df1$donor),
      n_samples  = nlevels(df1$sample),
      estimate   = sm["trait", "Estimate"],
      SE         = sm["trait", "Std. Error"],
      df         = sm["trait", "df"],
      t_value    = sm["trait", "t value"],
      pvalue     = sm["trait", "Pr(>|t|)"],
      singular   = isSingular(fit),
      stringsAsFactors = FALSE
    )
  }
}

lmm_path = do.call(rbind, lmm_path_list)
# BH-FDR across all module x pathology-trait pairs, same convention as the
# F03c v02 Pearson screen (BH across all module x trait pairs combined)
lmm_path$padj_BH = p.adjust(lmm_path$pvalue, method = "BH")

write_csv(lmm_path,
         file = paste0(out_dir, script_ind, "Module_eigengene_pathology_LMM.csv"))


### donor-aware TREM2-variant LMMs -------------------------------------------

message("\n\n   *TREM2Variant LMM, corrected random-effect structure \n")

lmm_trem2_omni_list  = list()
lmm_trem2_pairs_list = list()

for (mod1 in rownames(me_mat)){

  df1 = data.frame(
    eigengene    = as.numeric(me_mat[mod1, ]),
    TREM2Variant = factor(meta_lmm$TREM2Variant),
    cluster_name = factor(meta_lmm$cluster_name),
    donor        = factor(meta_lmm$BrainBankNetworkIDFormatted),
    sample       = factor(meta_lmm$sample),
    stringsAsFactors = FALSE
  )
  df1 = df1[complete.cases(df1), ]
  df1$TREM2Variant = droplevels(df1$TREM2Variant)
  df1$cluster_name = droplevels(df1$cluster_name)
  df1$donor        = droplevels(df1$donor)
  df1$sample       = droplevels(df1$sample)

  if (nlevels(df1$TREM2Variant) < 2 || nlevels(df1$donor) < 3) next

  fit = tryCatch(
    lmerTest::lmer(eigengene ~ TREM2Variant + (1 | cluster_name) + (1 | donor) + (1 | sample),
                   data = df1, REML = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) next

  aov1 = tryCatch(anova(fit, type = 2), error = function(e) NULL)
  get_aov = function(col){
    if (!is.null(aov1) && "TREM2Variant" %in% rownames(aov1) && col %in% colnames(aov1)){
      aov1["TREM2Variant", col]
    } else NA_real_
  }

  lmm_trem2_omni_list[[mod1]] = data.frame(
    module     = mod1,
    n_obs      = nrow(df1),
    n_clusters = nlevels(df1$cluster_name),
    n_donors   = nlevels(df1$donor),
    n_samples  = nlevels(df1$sample),
    F_value    = get_aov("F value"),
    num_df     = get_aov("NumDF"),
    den_df     = get_aov("DenDF"),
    pvalue     = get_aov("Pr(>F)"),
    singular   = isSingular(fit),
    stringsAsFactors = FALSE
  )

  lv       = levels(df1$TREM2Variant)
  ref      = lv[1]
  fe_names = names(fixef(fit))
  n_fe     = length(fe_names)
  all_pairs = combn(lv, 2, simplify = FALSE)

  pairs_rows = list()
  for (pp in all_pairs){
    g1 = pp[1]; g2 = pp[2]
    L = rep(0, n_fe)
    if (g1 != ref) L[which(fe_names == paste0("TREM2Variant", g1))] = -1
    if (g2 != ref) L[which(fe_names == paste0("TREM2Variant", g2))] =  1

    ct = tryCatch(
      lmerTest::contest(fit, L = matrix(L, nrow = 1), joint = FALSE),
      error = function(e) NULL
    )
    if (is.null(ct) || nrow(ct) == 0) next

    pairs_rows[[length(pairs_rows) + 1]] = data.frame(
      module   = mod1,
      contrast = paste0(g2, " - ", g1),
      estimate = ct[1, "Estimate"],
      SE       = ct[1, "Std. Error"],
      df       = ct[1, "df"],
      t.ratio  = ct[1, "t value"],
      p.value  = ct[1, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }
  if (length(pairs_rows) > 0){
    lmm_trem2_pairs_list[[mod1]] = do.call(rbind, pairs_rows)
  }
}

lmm_trem2_omni  = do.call(rbind, lmm_trem2_omni_list)
lmm_trem2_pairs = do.call(rbind, lmm_trem2_pairs_list)

if (!is.null(lmm_trem2_omni) && nrow(lmm_trem2_omni) > 0){
  lmm_trem2_omni$padj_BH = p.adjust(lmm_trem2_omni$pvalue, method = "BH")
  write_csv(lmm_trem2_omni,
           file = paste0(out_dir, script_ind, "Module_eigengene_LMM_TREM2_corrected_omnibus.csv"))
}

if (!is.null(lmm_trem2_pairs) && nrow(lmm_trem2_pairs) > 0){
  lmm_trem2_pairs$padj_BH = p.adjust(lmm_trem2_pairs$p.value, method = "BH")
  write_csv(lmm_trem2_pairs,
           file = paste0(out_dir, script_ind, "Module_eigengene_LMM_TREM2_corrected_pairwise.csv"))
}


### save replot bundle -------------------------------------------------------

message("\n\n   *Save small replot bundle \n")

replot_bundle = list(
  me_mat        = me_mat,                                            # modules x cluster_samples
  meta          = meta[, c("cluster_sample", "cluster_name", "sample",
                           "BrainBankNetworkIDFormatted", "TREM2Variant", path_traits)],
  mod_gene_tab  = mod_gene_tab[, c("gene", "module", "module_number", "colors")],
  mods          = mods,
  GO_results_tab = bulk_data$GO_results$by_comp_GO_res
)

qsave(replot_bundle, file = paste0(out_dir, script_ind, "replot_bundle.qs"))

message("\n\nDone. Written to ", out_dir, "\n")
