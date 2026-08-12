# LD_X20: Supplementary table of the top 10 WGCNA hub genes per module.

library(tidyverse)

### paths -------------------------------------------------------------------
base_candidates = c("/rds/general/user/lvd25/home/AST_scRNAseq_TREM2",   # HPC
                    "/Volumes/lvd25/home/AST_scRNAseq_TREM2")            # RDS mounted locally
base = base_candidates[dir.exists(base_candidates)][1]
if (is.na(base)) stop("Neither RDS path is reachable - is the share mounted?")

hub_csv    = file.path(base, "LD_F_DESeq_pseudobulk_WGCNA/LD_F03e_v01",
                       "LD_F03e_v01_Module_hub_genes_top10.csv")
out_dir    = file.path(base, "LD_X_Thesis_Presentation_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_ind = "LD_X20_"
if (!file.exists(hub_csv)) stop("Missing input: ", hub_csv)
message("Using base: ", base)

### data ---------------------------------------------------------------------
hub = read_csv(hub_csv, show_col_types = FALSE)

tab = hub %>%
  arrange(as.integer(sub("^M", "", module)), hub_rank) %>%   # M1, M2, ... M15, not M1, M10, M11
  group_by(module) %>%
  summarise(`Top 10 hub genes (by kME)` = paste(gene[order(hub_rank)], collapse = ", "),
            .groups = "drop") %>%
  mutate(mod_n = as.integer(sub("^M", "", module))) %>%
  arrange(mod_n) %>%
  select(Module = module, `Top 10 hub genes (by kME)`)

write_csv(tab, file.path(out_dir, paste0(script_ind, "hub_genes_top10.csv")))
# long form too, with the kME values, for anyone who wants the numbers
write_csv(hub %>% arrange(as.integer(sub("^M", "", module)), hub_rank) %>%
            select(module, hub_rank, gene, kME_own),
          file.path(out_dir, paste0(script_ind, "hub_genes_top10_long.csv")))

### render -------------------------------------------------------------------
cap = paste0("Supplementary Table 9. The ten hub genes of each co-expression module, ",
             "ranked by module membership (kME) within their own module and listed in ",
             "rank order. Modules were defined by WGCNA on the DEG-seeded gene set. ",
             "The unassigned set (M0) is not a module and is not shown. kME values are ",
             "given in ", script_ind, "hub_genes_top10_long.csv.")

if (!requireNamespace("flextable", quietly = TRUE)) {
  message("flextable not installed - CSVs written only.")
} else {
  library(flextable)
  ft = flextable(tab) %>%
    bold(part = "header") %>% italic(part = "header") %>%
    bold(j = "Module", part = "body") %>%
    italic(j = "Top 10 hub genes (by kME)", part = "body") %>%   # gene symbols
    valign(valign = "top", part = "body") %>%
    border_remove() %>%
    hline_top(part = "header", border = fp_border_default(width = 1.5)) %>%
    hline_bottom(part = "header", border = fp_border_default(width = 1.5)) %>%
    hline_bottom(part = "body", border = fp_border_default(width = 1.5)) %>%
    fontsize(size = 9, part = "all") %>%
    set_caption(cap) %>%
    width(j = c("Module", "Top 10 hub genes (by kME)"), width = c(0.8, 5.9))

  ok = tryCatch({ save_as_docx(ft, path = file.path(out_dir, paste0(script_ind, "hub_genes_table.docx"))); TRUE },
                error = function(e) { message("  FAILED docx: ", conditionMessage(e)); FALSE })
  if (ok) message("  wrote hub_genes_table.docx")
  try(save_as_image(ft, path = file.path(out_dir, paste0(script_ind, "hub_genes_table.png")), res = 300),
      silent = TRUE)
}

message("Done. ", nrow(tab), " modules x 10 hub genes. Outputs in: ", out_dir)
