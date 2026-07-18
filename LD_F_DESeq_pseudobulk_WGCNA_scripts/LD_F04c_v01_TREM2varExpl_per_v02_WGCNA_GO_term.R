message("\n\n##########################################################################\n",
        "# Start LD_F04c v01: TREM2 variance-explained per WGCNA GO term: ", Sys.time(),
        "\n##########################################################################\n",
        "\n   For each GO-BP term enriched in the v02 WGCNA astrocyte modules, computes",
        "\n   per-gene TREM2 variance (from H01 variancePartition) under two setups:",
        "\n",
        "\n   Setup A  -- genome-level (FAIR, addresses circularity):",
        "\n     Foreground: ALL genes annotated to the GO term (org.Hs.eg.db GOALL),",
        "\n     Background: all expressed genes (H01 genome-wide).",
        "\n     Tests whether the biological process is intrinsically TREM2-driven",
        "\n     independently of the DESeq filter. Each GO term shown once (deduplicated",
        "\n     by ID); n_modules column records how many included modules contained it.",
        "\n",
        "\n   Setup B  -- within-filter (fair within-set comparison):",
        "\n     Foreground: module genes annotated to the GO term (DESeq-filtered set),",
        "\n     Background: all genes across included modules (DESeq-filtered set).",
        "\n     Tests which processes stand out within the TREM2-associated gene set.",
        "\n     Module x term combinations kept (module genes differ per module).",
        "\n",
        "\n   Included modules: M1 M5 M8 M11 M12 M14 M15 (astrocyte-relevant).",
        "\n   Excluded: M0 (grey); M2 M3 M6 M7 (contamination);",
        "\n             M4 M9 M10 M13 (no interpretable GO signature).",
        "\n##########################################################################\n\n")

library(qs)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)


### directories / index --------------------------------------------------------

main_dir   = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)
out_dir    = paste0(main_dir, "LD_F_DESeq_pseudobulk_WGCNA/LD_F04c_v01/")
dir.create(out_dir, showWarnings = FALSE)
script_ind = "LD_F04c_v01_"
prim_var   = "TREM2Variant"

mods_include = c("M1","M5","M8","M11","M12","M14","M15")


### 1. per-gene TREM2 variance from H01 variancePartition ---------------------

bulk_h01 = qread(paste0(main_dir, "LD_H_VarPartition_output/LD_H01_v02_bulk_data.qs"))
varPart   = as.data.frame(bulk_h01$varPart_analysis$varPart)

if (!(prim_var %in% colnames(varPart)))
  stop("'", prim_var, "' not a column of the H01 varPart table.")

gene_var         = setNames(varPart[[prim_var]], rownames(varPart))
gene_var         = gene_var[!is.na(gene_var)]
bg_genome_median = median(gene_var)

message(sprintf("\n   H01 genome-wide TREM2 variance: median = %.4f  (n = %d genes)\n",
                bg_genome_median, length(gene_var)))


### 2. v02 GO results filtered to included modules ----------------------------

go_all  = read_csv(paste0(main_dir,
                          "LD_F_DESeq_pseudobulk_WGCNA/LD_F03c_v02/LD_F03c_v02_GO_results_by_comp.csv"),
                   show_col_types = FALSE)
go_filt = go_all[go_all$module %in% mods_include, ]
message(sprintf("   GO terms in included modules: %d  (modules: %s)\n",
                nrow(go_filt), paste(sort(unique(go_filt$module)), collapse = ", ")))


### 3. v02 bulk_data: Setup B background = all included-module genes ----------

bd = qread(paste0(main_dir,
                  "LD_F_DESeq_pseudobulk_WGCNA/LD_F03c_v02/LD_F03c_v02_bulk_data.qs"))

bg_B_genes  = intersect(unique(unlist(bd$wgcna$mod_genes[mods_include])), names(gene_var))
bg_B_var    = gene_var[bg_B_genes]
bg_B_median = median(bg_B_var)

message(sprintf("   Setup B background: %d module genes with H01 variance; median = %.4f\n",
                length(bg_B_genes), bg_B_median))


### 4. Build GO-term -> all annotated SYMBOL genes for Setup A ----------------

message("   Building GO-to-SYMBOL mapping (GOALL, hierarchy-aware)...")

go_ids_unique = unique(go_filt$ID)

go_anno = suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db,
                        keys    = go_ids_unique,
                        columns = "SYMBOL",
                        keytype = "GOALL")
)
go_anno = go_anno[!is.na(go_anno$SYMBOL), ]

# first column = the key (named "GOALL"); split to list
go_to_genes_A = split(go_anno$SYMBOL, go_anno[[1]])

message(sprintf("   GO-to-gene mapping built for %d term IDs\n", length(go_to_genes_A)))


###############################################################################
# 5. SETUP A: all annotated genes per GO term vs genome-wide background
#    One test per unique GO term (deduplicate across modules).
###############################################################################

message("\n   Running Setup A (all annotated genes vs genome-wide)...\n")

# unique GO terms with module membership recorded
go_unique = go_filt %>%
  group_by(ID, Description) %>%
  summarise(modules_containing = paste(sort(unique(module)), collapse = ", "),
            n_modules          = n_distinct(module),
            .groups            = "drop")

res_A       = list()
term_vals_A = list()

for (i in seq_len(nrow(go_unique))){

  go_id   = go_unique$ID[i]
  genes_i = go_to_genes_A[[go_id]]
  if (is.null(genes_i)) next
  genes_i = intersect(genes_i, names(gene_var))
  if (length(genes_i) < 3) next

  v_i  = gene_var[genes_i]
  rest = gene_var[setdiff(names(gene_var), genes_i)]
  wp   = suppressWarnings(wilcox.test(v_i, rest, alternative = "greater")$p.value)

  term_vals_A[[go_id]] = as.numeric(v_i)
  res_A[[go_id]] = tibble(
    ID                 = go_id,
    Description        = go_unique$Description[i],
    modules_containing = go_unique$modules_containing[i],
    n_modules          = go_unique$n_modules[i],
    n_all_annotated    = length(genes_i),
    median_TREM2_var   = median(v_i),
    mean_TREM2_var     = mean(v_i),
    max_TREM2_var      = max(v_i),
    frac_above_0.10    = mean(v_i > 0.10),
    median_minus_bg    = median(v_i) - bg_genome_median,
    wilcox_p           = wp
  )
}

res_A = bind_rows(res_A)
res_A$wilcox_padj = p.adjust(res_A$wilcox_p, method = "BH")
res_A = res_A %>% arrange(desc(median_TREM2_var))

write_csv(res_A,
          paste0(out_dir, script_ind, "SetupA_all_annotated_vs_genome.csv"))
message(sprintf("   Setup A: %d unique GO terms tested\n", nrow(res_A)))


###############################################################################
# 6. SETUP B: module genes per GO term vs all included-module genes
#    Module x term combinations kept (different module gene sets).
###############################################################################

message("   Running Setup B (module genes vs module-gene background)...\n")

res_B       = list()
term_vals_B = list()

for (i in seq_len(nrow(go_filt))){

  genes_i = unlist(strsplit(go_filt$geneID[i], "/"))
  genes_i = intersect(genes_i, names(gene_var))
  if (length(genes_i) < 3) next

  v_i  = gene_var[genes_i]
  rest = bg_B_var[setdiff(names(bg_B_var), genes_i)]
  wp   = suppressWarnings(wilcox.test(v_i, rest, alternative = "greater")$p.value)

  key = paste(go_filt$module[i], go_filt$ID[i], sep = "|")
  term_vals_B[[key]] = as.numeric(v_i)
  res_B[[key]] = tibble(
    key              = key,
    module           = go_filt$module[i],
    ID               = go_filt$ID[i],
    Description      = go_filt$Description[i],
    n_module_genes   = length(genes_i),
    median_TREM2_var = median(v_i),
    mean_TREM2_var   = mean(v_i),
    max_TREM2_var    = max(v_i),
    frac_above_0.10  = mean(v_i > 0.10),
    median_minus_bg  = median(v_i) - bg_B_median,
    wilcox_p         = wp
  )
}

res_B = bind_rows(res_B)
res_B$wilcox_padj = p.adjust(res_B$wilcox_p, method = "BH")
res_B = res_B %>% arrange(desc(median_TREM2_var))

write_csv(res_B %>% dplyr::select(-key),
          paste0(out_dir, script_ind, "SetupB_module_genes_vs_module_bg.csv"))
message(sprintf("   Setup B: %d module x term combinations tested\n", nrow(res_B)))


###############################################################################
# 7. PLOTS
###############################################################################

top_n      = 30   # ranked dotplot
top_n_dist = 10   # violin/distribution plot

# ── Setup A: ranked dotplot ──────────────────────────────────────────────────

pl_A = res_A %>%
  slice_head(n = top_n) %>%
  mutate(
    Description   = factor(Description, levels = rev(unique(Description))),
    n_modules_fct = factor(ifelse(n_modules == 1, "1 module", "2+ modules"))
  )

p_A_dot = ggplot(pl_A, aes(x = median_TREM2_var, y = Description)) +
  geom_point(aes(size = n_all_annotated, colour = n_modules_fct)) +
  geom_vline(xintercept = bg_genome_median, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("1 module" = "dodgerblue", "2+ modules" = "orange"),
                      name   = "module presence") +
  scale_size_continuous(name = "n annotated genes", range = c(2, 7)) +
  labs(
    title    = "Setup A: TREM2 variance per GO term (all annotated genes)",
    subtitle = paste0("Background: genome-wide median = ", round(bg_genome_median, 4),
                      " | included modules: ", paste(mods_include, collapse = ", ")),
    x = "median per-gene TREM2 variance fraction", y = NULL
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))

pdf(paste0(out_dir, script_ind, "SetupA_ranked_dotplot.pdf"), width = 11, height = 9)
print(p_A_dot)
dev.off()


# ── Setup A: Jaccard-deduplicated dotplot (gene-overlap-based, representative) ─
#
# Cluster GO terms by pairwise Jaccard on their genome-wide annotated gene sets.
# Within each cluster keep the term with the highest median TREM2 variance.
# Threshold = 0.5: terms sharing >50% of annotated genes -> same cluster.

jaccard = function(a, b) length(intersect(a, b)) / length(union(a, b))

ids_ordered = res_A$ID   # already sorted by desc(median_TREM2_var)

assigned    = rep(NA_integer_, length(ids_ordered))
cluster_rep = character(0)
cluster_idx = 0L

for (i in seq_along(ids_ordered)){
  if (!is.na(assigned[i])) next          # already absorbed into a cluster
  cluster_idx          = cluster_idx + 1L
  cluster_rep[cluster_idx] = ids_ordered[i]
  assigned[i]          = cluster_idx
  genes_i = go_to_genes_A[[ids_ordered[i]]]
  for (j in seq_along(ids_ordered)){
    if (!is.na(assigned[j])) next
    genes_j = go_to_genes_A[[ids_ordered[j]]]
    if (jaccard(genes_i, genes_j) >= 0.5) assigned[j] = cluster_idx
  }
}

rep_ids_jaccard = ids_ordered[!duplicated(assigned)]   # one rep per cluster

# build grouping table: one row per term showing its representative
rep_lookup = tibble(
  ID           = ids_ordered,
  cluster_idx  = assigned,
  rep_ID       = cluster_rep[assigned]
) %>%
  left_join(res_A %>% dplyr::select(ID, Description, median_TREM2_var, n_modules,
                                     modules_containing, n_all_annotated),
            by = "ID") %>%
  left_join(res_A %>% dplyr::select(ID, rep_Description = Description),
            by = c("rep_ID" = "ID")) %>%
  mutate(is_representative = ID == rep_ID) %>%
  arrange(cluster_idx, desc(median_TREM2_var)) %>%
  dplyr::select(cluster_idx, rep_ID, rep_Description, is_representative,
                ID, Description, median_TREM2_var, n_all_annotated,
                n_modules, modules_containing)

write_csv(rep_lookup,
          paste0(out_dir, script_ind, "SetupA_jaccard0.5_grouping.csv"))

res_A_jac = res_A %>%
  filter(ID %in% rep_ids_jaccard, wilcox_padj < 0.05) %>%
  arrange(desc(median_TREM2_var)) %>%
  slice_head(n = top_n) %>%
  mutate(
    Description   = factor(Description, levels = rev(unique(Description))),
    n_modules_fct = factor(ifelse(n_modules == 1, "1 module", "2+ modules"))
  )

message(sprintf("   Jaccard deduplication: %d terms -> %d representatives\n",
                length(ids_ordered), length(rep_ids_jaccard)))

p_A_jac = ggplot(res_A_jac, aes(x = median_TREM2_var, y = Description)) +
  geom_point(aes(size = n_all_annotated, colour = n_modules_fct)) +
  geom_vline(xintercept = bg_genome_median, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("1 module" = "dodgerblue", "2+ modules" = "orange"),
                      name   = "module presence") +
  scale_size_continuous(name = "n annotated genes", range = c(2, 7)) +
  labs(
    title    = "Setup A: TREM2 variance per GO term (Jaccard-deduplicated, Jaccard >= 0.5)",
    subtitle = paste0("Representative = highest median TREM2 variance per gene-overlap cluster; ",
                      "background median = ", round(bg_genome_median, 4)),
    x = "median per-gene TREM2 variance fraction", y = NULL
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))

pdf(paste0(out_dir, script_ind, "SetupA_jaccard_dedup_dotplot.pdf"), width = 11, height = 9)
print(p_A_jac)
dev.off()


# ── Setup B: ranked dotplot ──────────────────────────────────────────────────

pl_B = res_B %>%
  slice_head(n = top_n) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

p_B_dot = ggplot(pl_B, aes(x = median_TREM2_var, y = Description)) +
  geom_point(aes(size = n_module_genes, colour = module)) +
  geom_vline(xintercept = bg_B_median, linetype = "dashed", colour = "grey40") +
  scale_size_continuous(name = "n module genes", range = c(2, 7)) +
  labs(
    title    = "Setup B: TREM2 variance per GO term (module genes vs module background)",
    subtitle = paste0("Background: all included-module genes, median = ", round(bg_B_median, 4),
                      " | included modules: ", paste(mods_include, collapse = ", ")),
    x = "median per-gene TREM2 variance fraction", y = NULL
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))

pdf(paste0(out_dir, script_ind, "SetupB_ranked_dotplot.pdf"), width = 11, height = 9)
print(p_B_dot)
dev.off()


# ── Setup A: distribution (violin + box) top 10 terms vs genome-wide ────────

top_ids_A   = head(res_A$ID, top_n_dist)
top_desc_A  = res_A$Description[match(top_ids_A, res_A$ID)]

dist_A = bind_rows(
  lapply(seq_along(top_ids_A), function(j){
    tibble(set = top_desc_A[j], TREM2_var = term_vals_A[[top_ids_A[j]]])
  })
)
dist_A = bind_rows(dist_A,
                   tibble(set = "(genome-wide background)", TREM2_var = as.numeric(gene_var)))

set_order_A = c(top_desc_A, "(genome-wide background)")
dist_A$set  = factor(dist_A$set, levels = rev(unique(set_order_A)))

p_A_dist = ggplot(dist_A, aes(x = TREM2_var, y = set,
                               fill = set == "(genome-wide background)")) +
  geom_violin(scale = "width", colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.6) +
  geom_vline(xintercept = bg_genome_median, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("TRUE" = "grey70", "FALSE" = "dodgerblue"), guide = "none") +
  labs(
    title    = "Setup A: per-gene TREM2 variance — top 10 GO terms vs genome-wide",
    subtitle = "All genes annotated to each term; dashed = genome-wide median",
    x = "per-gene TREM2 variance fraction", y = NULL
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))

pdf(paste0(out_dir, script_ind, "SetupA_top10_distribution.pdf"), width = 9, height = 7)
print(p_A_dist)
dev.off()


# ── Setup A: distribution top 10 Jaccard-deduplicated (Jaccard >= 0.5) ───────

top_ids_A_jac  = head(res_A$ID[res_A$ID %in% rep_ids_jaccard], top_n_dist)
top_desc_A_jac = res_A$Description[match(top_ids_A_jac, res_A$ID)]

dist_A_jac = bind_rows(
  lapply(seq_along(top_ids_A_jac), function(j){
    tibble(set = top_desc_A_jac[j], TREM2_var = term_vals_A[[top_ids_A_jac[j]]])
  })
)
dist_A_jac = bind_rows(dist_A_jac,
                        tibble(set = "(genome-wide background)", TREM2_var = as.numeric(gene_var)))
dist_A_jac$set = factor(dist_A_jac$set,
                         levels = rev(c(top_desc_A_jac, "(genome-wide background)")))

p_A_jac_dist = ggplot(dist_A_jac, aes(x = TREM2_var, y = set,
                                        fill = set == "(genome-wide background)")) +
  geom_violin(scale = "width", colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.6) +
  geom_vline(xintercept = bg_genome_median, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("TRUE" = "grey70", "FALSE" = "dodgerblue"), guide = "none") +
  labs(
    title    = "Setup A: per-gene TREM2 variance — top 10 GO terms (Jaccard-deduplicated)",
    subtitle = "Jaccard >= 0.5; all genes annotated to each term; dashed = genome-wide median",
    x = "per-gene TREM2 variance fraction", y = NULL
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))

pdf(paste0(out_dir, script_ind, "SetupA_top10_jaccard_distribution.pdf"), width = 9, height = 7)
print(p_A_jac_dist)
dev.off()


# ── Setup B: distribution (violin + box) top 10 terms vs module background ──

top_keys_B  = head(res_B$key, top_n_dist)
top_desc_B  = res_B$Description[match(top_keys_B, res_B$key)]

dist_B = bind_rows(
  lapply(seq_along(top_keys_B), function(j){
    tibble(set = paste0(top_desc_B[j], "\n(", res_B$module[res_B$key == top_keys_B[j]][1], ")"),
           TREM2_var = term_vals_B[[top_keys_B[j]]])
  })
)
dist_B = bind_rows(dist_B,
                   tibble(set = "(module background)", TREM2_var = as.numeric(bg_B_var)))

set_order_B = c(paste0(top_desc_B, "\n(",
                        res_B$module[match(top_keys_B, res_B$key)], ")"),
                "(module background)")
dist_B$set = factor(dist_B$set, levels = rev(unique(set_order_B)))

p_B_dist = ggplot(dist_B, aes(x = TREM2_var, y = set,
                               fill = set == "(module background)")) +
  geom_violin(scale = "width", colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.6) +
  geom_vline(xintercept = bg_B_median, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("TRUE" = "grey70", "FALSE" = "steelblue"), guide = "none") +
  labs(
    title    = "Setup B: per-gene TREM2 variance — top 10 GO terms vs module background",
    subtitle = "Module genes per term; dashed = all-module-gene median; module in label",
    x = "per-gene TREM2 variance fraction", y = NULL
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7))

pdf(paste0(out_dir, script_ind, "SetupB_top10_distribution.pdf"), width = 9, height = 7)
print(p_B_dist)
dev.off()


###############################################################################
# 8. M1 sub-theme clustering: Jaccard on M1 module genes per term
#    One dot per GO term; y-axis = data-driven cluster (named by representative).
#    Representative = term with highest median TREM2 variance in each cluster.
###############################################################################

message("\n   Clustering M1 GO terms by Jaccard on module gene sets...\n")

m1_go = go_filt[go_filt$module == "M1", ]

# gene sets: M1 module genes annotated to each term, intersected with H01 gene_var
m1_gene_sets = lapply(seq_len(nrow(m1_go)), function(i){
  intersect(unlist(strsplit(m1_go$geneID[i], "/")), names(gene_var))
})
names(m1_gene_sets) = paste0("M1|", m1_go$ID)

# pull per-term stats from res_B, ordered desc median TREM2 var; cap to top 50
m1_keys  = paste0("M1|", m1_go$ID)
m1_stats = res_B[res_B$key %in% m1_keys, ]
m1_stats = head(m1_stats[order(-m1_stats$median_TREM2_var), ], 50)

ids_m1      = m1_stats$key
assigned_m1 = rep(NA_integer_, length(ids_m1))
k = 0L

for (i in seq_along(ids_m1)){
  if (!is.na(assigned_m1[i])) next
  k = k + 1L
  assigned_m1[i] = k
  gi = m1_gene_sets[[ids_m1[i]]]
  for (j in seq_along(ids_m1)){
    if (!is.na(assigned_m1[j])) next
    gj = m1_gene_sets[[ids_m1[j]]]
    if (length(gi) > 0 && length(gj) > 0 && jaccard(gi, gj) >= 0.2)
      assigned_m1[j] = k
  }
}

rep_keys_m1 = ids_m1[!duplicated(assigned_m1)]
rep_desc_m1 = m1_stats$Description[match(rep_keys_m1, m1_stats$key)]

m1_cluster_tab = tibble(
  key             = ids_m1,
  cluster_idx     = assigned_m1,
  rep_key         = rep_keys_m1[assigned_m1],
  rep_Description = rep_desc_m1[assigned_m1]
) %>%
  left_join(m1_stats %>% dplyr::select(key, Description, median_TREM2_var, n_module_genes),
            by = "key")

message(sprintf("   M1: top %d terms -> %d clusters (Jaccard >= 0.2 on module genes)\n",
                nrow(m1_stats), length(rep_keys_m1)))

write_csv(m1_cluster_tab %>%
            dplyr::select(cluster_idx, rep_key, rep_Description,
                          key, Description, median_TREM2_var, n_module_genes),
          paste0(out_dir, script_ind, "M1_jaccard0.4_subtheme_grouping.csv"))

# order clusters by their representative's median TREM2 var (already top of cluster)
rep_order_m1 = m1_cluster_tab %>%
  filter(key == rep_key) %>%
  arrange(desc(median_TREM2_var)) %>%
  pull(rep_Description)

# keep only clusters with >= 2 terms OR the cluster is the top representative
# (single-term clusters still show, they just have one dot)
m1_cluster_tab$rep_Description = factor(m1_cluster_tab$rep_Description,
                                         levels = rev(rep_order_m1))

# add n_terms per cluster for annotation
m1_cluster_tab = m1_cluster_tab %>%
  group_by(cluster_idx) %>%
  mutate(n_terms_cluster = n()) %>%
  ungroup()

# crossbar data: median per cluster
m1_crossbar = m1_cluster_tab %>%
  group_by(rep_Description) %>%
  summarise(med = median(median_TREM2_var),
            n   = n(),
            .groups = "drop")

p_m1 = ggplot(m1_cluster_tab,
               aes(x = median_TREM2_var, y = rep_Description)) +
  geom_jitter(aes(size = n_module_genes), height = 0.2, alpha = 0.6,
              colour = "steelblue", shape = 16) +
  geom_crossbar(data = m1_crossbar,
                aes(x = med, xmin = med, xmax = med, y = rep_Description),
                width = 0.5, colour = "black", linewidth = 0.6) +
  geom_text(data = m1_crossbar,
            aes(x = Inf, y = rep_Description,
                label = sprintf("n=%d  %.3fx", n, med / bg_B_median)),
            hjust = 1.05, size = 2.8, colour = "grey20") +
  geom_vline(xintercept = bg_B_median, linetype = "dashed", colour = "grey40") +
  scale_size_continuous(name = "n module genes", range = c(1.5, 6)) +
  labs(
    title    = "M1 GO sub-themes: Jaccard clustering on module gene sets",
    subtitle = "Top 50 M1 terms; Jaccard >= 0.2; representative = highest median TREM2 var per cluster; dashed = all-module median",
    x = "median per-gene TREM2 variance fraction (Setup B)", y = NULL
  ) +
  theme_bw() +
  theme(axis.text.y  = element_text(size = 7),
        plot.margin  = margin(5, 80, 5, 5))

pdf(paste0(out_dir, script_ind, "M1_jaccard_subtheme_dotplot.pdf"), width = 11, height = 9)
print(p_m1)
dev.off()


sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed LD_F04c v01 ", Sys.time(),
        "\n##########################################################################\n\n")
