message("\n\n##########################################################################\n",
        "# Start G03a: CellChat per-run visualisation - ungrouped ", Sys.time(),
        "\n##########################################################################\n",
        "\n   Loads the ungrouped CellChat object from G02a and produces per-run\n",
        "   visualisations:\n",
        "     - aggregated network circle plots (count + weight)\n",
        "     - interaction heatmaps (count + weight)\n",
        "     - signalling-role heatmaps (full clusters + MIC->AST pathways)\n",
        "     - signalling-role scatter (2D senders vs receivers)\n",
        "     - per-pathway plots looped over MIC->AST relevant pathways\n",
        "         (chord, role network heatmap, gene expression violin)\n",
        "     - MIC->AST L-R interaction table\n",
        "     - MIC->AST L-R bubble plot (CellChat-native + custom ggplot)\n",
        "     - pathway-gene Z-score heatmaps (uses G01a merged pseudobulk)\n",
        "\n   Single CellChat object - no across-run comparisons possible.\n",
        "   computeCommunProbPathway() and aggregateNet() were already run in G02a;\n",
        "   netAnalysis_computeCentrality() runs here (G02a skipped it).\n",
        "\n##########################################################################\n\n")


#set environment/load packages
library(qs)
library(tidyverse)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(pheatmap)
library(viridis)
library(colorRamps)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "LD_G03a_v001_"

#specify input + output directories
in_dir  = paste0(main_dir, "LD_G_MIC_AST_communication_analysis_output/")
out_dir = paste0(in_dir, "LD_G03a/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}


#input files
in_cellchat = paste0(in_dir, "LD_G02a_v001_cellchat_ungrouped.qs")
in_bulkdata = paste0(in_dir, "LD_G01a_bulk_data_AST_MIC.qs")



###########################################################
# functions
###########################################################

#custom colour palette for variable values defined in vector v
pal = function(v){
  v2 = length(unique(v))
  if (v2 == 2){
    p2 = c("grey20", "dodgerblue")
  } else if (v2 == 3){
    p2 = c("dodgerblue", "grey20", "orange")
  } else if (v2 == 4){
    p2 = c("dodgerblue", "green4", "grey20", "orange")
  } else if (v2 < 6){
    p2 = matlab.like(6)[1:v2]
  } else {
    p2 = matlab.like(v2)
  }
  return(p2)
}



###########################################################
# 1. load CellChat object and merged pseudobulk
###########################################################

message("\n\n          *** Load CellChat object and merged pseudobulk... ", Sys.time(), "\n\n")


cellchat = qread(file = in_cellchat)
bulk_data = qread(file = in_bulkdata)

cat("CellChat object - cells: ", length(cellchat@idents),
    "  clusters: ", length(levels(cellchat@idents)),
    "  pathways: ", length(cellchat@netP$pathways), "\n", sep = "")



###########################################################
# 2. define cluster ordering (biology-driven)
###########################################################

# AST first (SLC1A2, then GFAP, then CHI3L1), then MIC (HOM, DAM, IRM, HLA, CRM).
# Within each subtype, alphanumeric order. Applied to plots so figures from
# different runs/comparisons share the same column/row order.

cluster_levels_present = levels(cellchat@idents)

cat("\nCluster levels present in CellChat object:\n")
cat("  ", paste(cluster_levels_present, collapse = "\n  "), "\n", sep = "")

ast_order = c(
  sort(grep("^AST_SLC1A2",  cluster_levels_present, value = TRUE)),
  sort(grep("^AST_GFAP",    cluster_levels_present, value = TRUE)),
  sort(grep("^AST_CHI3L1",  cluster_levels_present, value = TRUE))
)
mic_order = c(
  sort(grep("^HOM_",        cluster_levels_present, value = TRUE)),
  sort(grep("^DAM_",        cluster_levels_present, value = TRUE)),
  sort(grep("^IRM_",        cluster_levels_present, value = TRUE)),
  sort(grep("^HLA_",        cluster_levels_present, value = TRUE)),
  sort(grep("^CRM_",        cluster_levels_present, value = TRUE))
)
cluster_order = c(ast_order, mic_order)


### sanity check: cluster_order must cover all clusters in cellchat
unmatched = setdiff(cluster_levels_present, cluster_order)
if (length(unmatched) > 0){
  warning("Clusters not matched by the AST/MIC prefix patterns - appended at end: ",
          paste(unmatched, collapse = ", "))
  cluster_order = c(cluster_order, sort(unmatched))
}


### reorder CellChat idents to apply this order to all plots
# updateClusterLabels() only re-orders @idents factor levels, not @net or
# centrality matrices - so plot functions that read precomputed structures
# can show the old order. Reorder the factor directly, then re-run
# aggregateNet() so @net is in the new order. Centrality runs in section 3, so it picks up the new order automatically.
cellchat@idents = factor(cellchat@idents, levels = cluster_order)
cellchat = aggregateNet(cellchat)


### MIC and AST cluster lists for downstream filtering
clusters_AST = ast_order
clusters_MIC = mic_order

cat("\nCluster ordering (n = ", length(cluster_order), "):\n", sep = "")
cat("  AST (n = ", length(clusters_AST), "): ", paste(clusters_AST, collapse = ", "), "\n", sep = "")
cat("  MIC (n = ", length(clusters_MIC), "): ", paste(clusters_MIC, collapse = ", "), "\n", sep = "")



###########################################################
# 3. compute centrality (signalling-role analysis)
###########################################################

#  Computes sender / receiver / mediator / influencer scores per pathway per cluster.
# centrality is needed for role heatmaps (sections 6 + 8),
# role scatter (section 9), and the signallingRole_network sub-plot in section 10.

message("\n\n          *** Compute centrality... ", Sys.time(), "\n\n")


cellchat = netAnalysis_computeCentrality(cellchat, slot.name = "netP")



###########################################################
# 4. aggregated network circle plots (count + weight)
###########################################################

# Two circle plots: edge count (number of significant interactions) and edge
# weight (sum of communication probability).
# Node size refelects number of cells; the amount of clusters is to high to really see sth in the circle plot,
# but it's a useful overview of the overall communication landscape and a starting point for the more detailed plots that follow.

message("\n\n          *** Aggregated network circle plots... ", Sys.time(), "\n\n")


pdf(file = paste0(out_dir, script_ind, "aggregated_network_circle.pdf"),
    width = 14, height = 7)
{
  groupSize = as.numeric(table(cellchat@idents))
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = "Interaction weights / strength")
}
dev.off()



###########################################################
# 5. interaction heatmaps (count + weight)
###########################################################

# CellChat-native heatmap of the same data as the circle plots. Easier to read
# at 38 clusters than the circle plot.

message("\n\n          *** Interaction heatmaps... ", Sys.time(), "\n\n")


pdf(file = paste0(out_dir, script_ind, "aggregated_network_heatmap.pdf"),
    width = 16, height = 8)
{
  ht1 = netVisual_heatmap(cellchat, measure = "count",  color.heatmap = "Reds")
  ht2 = netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds")
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
}
dev.off()



###########################################################
# 6. signalling-role heatmaps - full cluster set
###########################################################

# Outgoing and incoming signalling strength per pathway per cluster, across
# all pathways.

message("\n\n          *** Signalling-role heatmaps (full)... ", Sys.time(), "\n\n")


n_pathways_full = length(cellchat@netP$pathways)
heatmap_height_full = max(15, 0.25 * n_pathways_full + 4)

pdf(file = paste0(out_dir, script_ind, "signallingRole_heatmap_full.pdf"),
    width = 14, height = heatmap_height_full)
{
  ht1 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",
                                          width = 7, height = heatmap_height_full - 2,
                                          color.heatmap = "BuGn")
  ht2 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                          width = 7, height = heatmap_height_full - 2,
                                          color.heatmap = "GnBu")
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
}
dev.off()



###########################################################
# 7. identify MIC->AST relevant pathways
###########################################################

# A pathway is "MIC->AST relevant" if at least one significant L-R interaction
# has source in MIC clusters and target in AST clusters. Used to scope per-
# pathway plots (sections 9, 10, 13) and the pathway-gene heatmap (section 14)
# to the question of interest, instead of looping over all ~100+ pathways.

message("\n\n          *** Identify MIC->AST relevant pathways... ", Sys.time(), "\n\n")


net_tab = subsetCommunication(cellchat)

net_tab_MIC_AST = net_tab[net_tab$source %in% clusters_MIC &
                          net_tab$target %in% clusters_AST, ]

pathways_MIC_AST = sort(unique(net_tab_MIC_AST$pathway_name))

cat("Total pathways in CellChat:           ", n_pathways_full, "\n",
    "MIC->AST relevant pathways:           ", length(pathways_MIC_AST), "\n",
    "Total L-R interactions:               ", nrow(net_tab), "\n",
    "MIC->AST L-R interactions:            ", nrow(net_tab_MIC_AST), "\n",
    sep = "")

write_csv(net_tab,
          file = paste0(out_dir, script_ind, "interactions_all.csv"))
write_csv(net_tab_MIC_AST,
          file = paste0(out_dir, script_ind, "interactions_MIC_to_AST.csv"))



###########################################################
# 8. signalling-role heatmaps - MIC->AST pathways only
###########################################################

# Same as section 6 but restricted to MIC->AST relevant pathways.
# should highlight which MIC clusters dominate the outgoing signal and which AST clusters receive it.

message("\n\n          *** Signalling-role heatmaps (MIC->AST pathways)... ", Sys.time(), "\n\n")


if (length(pathways_MIC_AST) > 0){

  heatmap_height_sub = max(8, 0.3 * length(pathways_MIC_AST) + 2)

  pdf(file = paste0(out_dir, script_ind, "signallingRole_heatmap_MIC_AST_pathways.pdf"),
      width = 14, height = heatmap_height_sub)
  {
    ht1 = netAnalysis_signalingRole_heatmap(cellchat,
                                            signaling = pathways_MIC_AST,
                                            pattern = "outgoing",
                                            width = 7, height = heatmap_height_sub - 2,
                                            color.heatmap = "BuGn")
    ht2 = netAnalysis_signalingRole_heatmap(cellchat,
                                            signaling = pathways_MIC_AST,
                                            pattern = "incoming",
                                            width = 7, height = heatmap_height_sub - 2,
                                            color.heatmap = "GnBu")
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  }
  dev.off()

} else {
  message("    No MIC->AST pathways - skipping section 8")
}



###########################################################
# 9. signalling-role scatter (2D senders vs receivers)
###########################################################

# Each cluster as a dot in (outgoing strength, incoming strength) space.
# Aggregated across all pathways.

message("\n\n          *** Signalling-role scatter... ", Sys.time(), "\n\n")


p_scatter = netAnalysis_signalingRole_scatter(cellchat)

pdf(file = paste0(out_dir, script_ind, "signallingRole_scatter.pdf"),
    width = 8, height = 6)
print(p_scatter)
dev.off()



###########################################################
# 10. per-pathway plots: chord + role network (MIC->AST pathways)
###########################################################

# Looped over MIC->AST relevant pathways. Two plot types per pathway, each
# saved to its own multi-page PDF (one page per pathway).
#  - netVisual_aggregate as chord: shows L-R flow through the pathway across
#    all clusters
#  - netAnalysis_signalingRole_network: heatmap of sender / receiver / mediator
#    / influencer scores per cluster for that pathway


message("\n\n          *** Per-pathway chord + role network (MIC->AST)... ", Sys.time(), "\n\n")


### chord plots
pdf(file = paste0(out_dir, script_ind, "per_pathway_chord_MIC_AST.pdf"),
    width = 9, height = 9)
{
  for (pw in pathways_MIC_AST){
    tryCatch({
      netVisual_aggregate(cellchat, signaling = pw, layout = "chord")
      title(main = pw, line = -1)
    }, error = function(e){
      message("    chord skipped for ", pw, ": ", conditionMessage(e))
    })
  }
}
dev.off()


### role network heatmaps
pdf(file = paste0(out_dir, script_ind, "per_pathway_signallingRole_network_MIC_AST.pdf"),
    width = 12, height = 4)
{
  for (pw in pathways_MIC_AST){
    tryCatch({
      netAnalysis_signalingRole_network(cellchat, signaling = pw,
                                        width = 10, height = 2.5, font.size = 9)
    }, error = function(e){
      message("    role network skipped for ", pw, ": ", conditionMessage(e))
    })
  }
}
dev.off()



###########################################################
# 11. MIC->AST L-R interactions: ranked table + summary
###########################################################

# Already exported the full MIC->AST table in section 7. Here add a pathway-
# level summary (sum prob per source-target-pathway combination, ranked).

message("\n\n          *** MIC->AST L-R rankings... ", Sys.time(), "\n\n")


if (length(pathways_MIC_AST) > 0){

  net_tab_MIC_AST_ranked = net_tab_MIC_AST[order(-net_tab_MIC_AST$prob), ]
  write_csv(net_tab_MIC_AST_ranked,
            file = paste0(out_dir, script_ind, "interactions_MIC_to_AST_ranked.csv"))


  path_tab_sel = net_tab_MIC_AST %>%
    group_by(source, target, pathway_name) %>%
    summarise(sum_prob = sum(prob),
              n_LR     = n(),
              .groups  = "drop") %>%
    arrange(desc(sum_prob))

  write_csv(path_tab_sel,
            file = paste0(out_dir, script_ind, "pathway_by_cluster_pair_MIC_to_AST.csv"))

  cat("Top 10 MIC->AST source-target-pathway combinations by summed probability:\n")
  print(head(path_tab_sel, 10))

} else {
  message("    No MIC->AST pathways - skipping section 11")
  net_tab_MIC_AST_ranked = net_tab_MIC_AST  # empty - keep variable defined
  path_tab_sel = tibble(source = character(), target = character(),
                        pathway_name = character(), sum_prob = numeric(),
                        n_LR = integer())
}



###########################################################
# 12. MIC->AST L-R bubble plot (CellChat-native + custom)
###########################################################

# Two bubble plots:
#  - netVisual_bubble: standard CellChat L-R bubble across MIC sources / AST
#    targets, restricted to MIC->AST pathways
#  - custom ggplot bubble: same data but with pathway-coloured dots, matches
#    Michael's F02

message("\n\n          *** MIC->AST bubble plots... ", Sys.time(), "\n\n")


### CellChat-native bubble
n_LR_MIC_AST = length(unique(net_tab_MIC_AST$interaction_name))
bubble_height = max(8, 0.25 * n_LR_MIC_AST + 2)

if (length(pathways_MIC_AST) > 0){
  pdf(file = paste0(out_dir, script_ind, "bubble_MIC_to_AST.pdf"),
      width = 12, height = bubble_height)
  {
    p = netVisual_bubble(cellchat,
                         sources.use   = clusters_MIC,
                         targets.use   = clusters_AST,
                         signaling     = pathways_MIC_AST,
                         remove.isolate = TRUE,
                         angle.x        = 45)
    print(p)
  }
  dev.off()
}


### custom ggplot bubble (Michael's F02 style)
if (length(pathways_MIC_AST) > 0 && nrow(net_tab_MIC_AST_ranked) > 0){

  t1 = net_tab_MIC_AST_ranked
  t1$source_target = paste0(t1$source, " => ", t1$target)
  t1 = t1[order(match(t1$pathway_name, unique(path_tab_sel$pathway_name)), -t1$prob), ]

  p_custom = ggplot(t1, aes(x = source_target, y = interaction_name_2,
                             size = prob, color = pathway_name)) +
    geom_point() +
    scale_color_manual(limits = unique(t1$pathway_name),
                       values = pal(unique(t1$pathway_name))) +
    scale_x_discrete(limits = unique(t1$source_target)) +
    scale_y_discrete(limits = unique(t1$interaction_name_2)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6))

  pdf(file = paste0(out_dir, script_ind, "bubble_MIC_to_AST_custom.pdf"),
      width = 14, height = bubble_height)
  print(p_custom)
  dev.off()

} else {
  message("    No MIC->AST pathways - skipping custom bubble plot")
}



###########################################################
# 13. per-pathway gene expression (CellChat-native violin)
###########################################################

# CellChat-native plotGeneExpression: violin plots of ligand/receptor genes
# split by cluster identity. Lightweight screen of which genes drive each
# pathway. Multi-page PDF, one page per MIC->AST pathway.

message("\n\n          *** Per-pathway gene expression violins... ", Sys.time(), "\n\n")


if (length(pathways_MIC_AST) > 0){

  pdf(file = paste0(out_dir, script_ind, "per_pathway_gene_expression_MIC_AST.pdf"),
      width = 10, height = 6)
  {
    for (pw in pathways_MIC_AST){
      tryCatch({
        p = plotGeneExpression(cellchat, signaling = pw)
        print(p + plot_annotation(title = pw))
      }, error = function(e){
        message("    gene expression skipped for ", pw, ": ", conditionMessage(e))
      })
    }
  }
  dev.off()

} else {
  message("    No MIC->AST pathways - skipping section 13")
}



###########################################################
# 14. pathway-gene Z-score heatmaps (uses G01a merged pseudobulk)
###########################################################

# For each MIC->AST pathway, plot expression Z-score (corrected VST) of the
# pathway's ligand / receptor / cofactor genes across all cluster_samples.
# Genes pulled from cellchat@DB; expression matrix from G01a merged pseudobulk.
# Annotated by lineage / TREM2Variant / NeuropathologicalDiagnosis / cluster.

message("\n\n          *** Pathway-gene Z-score heatmaps (MIC->AST)... ", Sys.time(), "\n\n")


if (length(pathways_MIC_AST) > 0){

  ### extract pathway-gene table for MIC->AST pathways
  # Loop pulls ligand subunits, receptor subunits, agonists, antagonists, co-A
  # and co-I receptors for each L-R pair in each pathway. Matches Michael's F02
  # extraction.
  # Note: cellchat@DB is the DEG-filtered DB from G02a, so pathway_genes_tab
  # only contains interactions that survived DEG filtering. This is consistent
  # with what was actually tested by computeCommunProb.

  pathway_genes_tab = tibble(pathway = character(),
                             interaction_name = character(),
                             gene_type = character(),
                             gene = character())

  for (pw in pathways_MIC_AST){

    path_ints = unique(net_tab_MIC_AST$interaction_name[net_tab_MIC_AST$pathway_name == pw])

    db = cellchat@DB
    t2 = db$interaction
    t3 = t2[t2$interaction_name %in% path_ints, ]

    for (j in seq_len(nrow(t3))){

      ### cofactors
      for (cf in c("agonist", "antagonist", "co_A_receptor", "co_I_receptor")){
        cf_name = t3[[cf]][j]
        if (cf_name %in% rownames(db$cofactor)){
          t4 = db$cofactor[rownames(db$cofactor) == cf_name, ]
          t5 = tibble(pathway          = pw,
                      interaction_name = t3$interaction_name[j],
                      gene_type        = cf,
                      gene             = unlist(t4))
          pathway_genes_tab = rbind(pathway_genes_tab, t5)
        }
      }

      ### ligand / ligand subunits
      if (t3$ligand[j] %in% rownames(db$complex)){
        t4 = db$complex[rownames(db$complex) == t3$ligand[j], ]
        t5 = tibble(pathway          = pw,
                    interaction_name = t3$interaction_name[j],
                    gene_type        = "ligand_subunit",
                    gene             = unlist(t4))
      } else {
        t5 = tibble(pathway          = pw,
                    interaction_name = t3$interaction_name[j],
                    gene_type        = "ligand",
                    gene             = t3$ligand[j])
      }
      pathway_genes_tab = rbind(pathway_genes_tab, t5)

      ### receptor / receptor subunits
      if (t3$receptor[j] %in% rownames(db$complex)){
        t4 = db$complex[rownames(db$complex) == t3$receptor[j], ]
        t5 = tibble(pathway          = pw,
                    interaction_name = t3$interaction_name[j],
                    gene_type        = "receptor_subunit",
                    gene             = unlist(t4))
      } else {
        t5 = tibble(pathway          = pw,
                    interaction_name = t3$interaction_name[j],
                    gene_type        = "receptor",
                    gene             = t3$receptor[j])
      }
      pathway_genes_tab = rbind(pathway_genes_tab, t5)
    }
  }

  pathway_genes_tab = pathway_genes_tab[pathway_genes_tab$gene != "", ]
  pathway_genes_tab = pathway_genes_tab[!is.na(pathway_genes_tab$gene), ]

  write_csv(pathway_genes_tab,
            file = paste0(out_dir, script_ind, "pathway_genes_table_MIC_AST.csv"))

  cat("Pathway-gene table: ",
      nrow(pathway_genes_tab), " rows across ",
      length(unique(pathway_genes_tab$pathway)), " pathways\n", sep = "")


  ### plot pathway-gene heatmaps
  # One page per pathway, columns = cluster_samples ordered by lineage then
  # cluster_name, rows = pathway genes labelled by gene_type. Z-score from the
  # merged G01a pseudobulk corrected VST.

  z_mat = bulk_data$gene_Z_scores$clusters_combined
  meta  = bulk_data$meta


  ### order columns: lineage > cluster_name > sample
  meta_ord = meta[order(match(meta$lineage, c("MIC", "AST")),
                         match(meta$cluster_name, cluster_order),
                         meta$sample), ]
  z_mat = z_mat[, meta_ord$cluster_sample, drop = FALSE]


  ### column annotation
  annot_col = data.frame(
    lineage                   = meta_ord$lineage,
    cluster_name              = meta_ord$cluster_name,
    TREM2Variant              = meta_ord$TREM2Variant,
    NeuropathologicalDiagnosis = meta_ord$NeuropathologicalDiagnosis,
    row.names                 = meta_ord$cluster_sample
  )

  annot_colors = list(
    lineage                    = c(AST = "dodgerblue", MIC = "orange"),
    TREM2Variant               = setNames(pal(levels(meta_ord$TREM2Variant)),
                                          levels(meta_ord$TREM2Variant)),
    NeuropathologicalDiagnosis = setNames(pal(levels(meta_ord$NeuropathologicalDiagnosis)),
                                          levels(meta_ord$NeuropathologicalDiagnosis))
  )


  pdf(file = paste0(out_dir, script_ind, "pathway_gene_heatmaps_MIC_AST.pdf"),
      width = 18, height = 8)
  {
    for (pw in pathways_MIC_AST){

      t1 = pathway_genes_tab[pathway_genes_tab$pathway == pw, ]
      t1 = t1[!duplicated(paste0(t1$gene, t1$gene_type)), ]
      t1 = t1[t1$gene %in% rownames(z_mat), ]

      if (nrow(t1) < 2){
        message("    skip ", pw, " - <2 genes available in merged pseudobulk")
        next
      }

      ### pull the Z-score submatrix for this pathway's genes
      pl_mat = z_mat[t1$gene, , drop = FALSE]
      rownames(pl_mat) = paste0(t1$gene, " (", t1$gene_type, ")")

      ### guard against all-NaN matrix (from constant genes in G01a)
      # all-NaN -> max(abs(...), na.rm=TRUE) returns -Inf -> lims invalid ->
      # pheatmap errors. tryCatch would silently drop the pathway, this is
      # explicit.
      if (all(is.nan(pl_mat))){
        message("    skip ", pw, " - all genes constant (NaN) in merged pseudobulk")
        next
      }

      lims = 0.7 * c(-max(abs(pl_mat), na.rm = TRUE),
                      max(abs(pl_mat), na.rm = TRUE))

      tryCatch({
        pheatmap::pheatmap(
          pl_mat,
          cluster_rows       = FALSE,
          cluster_cols       = FALSE,
          show_rownames      = TRUE,
          show_colnames      = FALSE,
          annotation_col     = annot_col,
          annotation_colors  = annot_colors,
          color              = colorRampPalette(c("blue", "white", "red"))(250),
          breaks             = seq(lims[1], lims[2], length.out = 251),
          border_color       = NA,
          fontsize           = 8,
          fontsize_row       = 6,
          main               = paste0(pw, " - ", nrow(pl_mat), " genes")
        )
      }, error = function(e){
        message("    heatmap skipped for ", pw, ": ", conditionMessage(e))
      })
    }
  }
  dev.off()

} else {
  message("    No MIC->AST pathways - skipping section 14")
}



###########################################################
# 15. save updated CellChat object (with centrality)
###########################################################

message("\n\n          *** Save updated CellChat object... ", Sys.time(), "\n\n")


qsave(cellchat,
      file = paste0(out_dir, script_ind, "cellchat_ungrouped_with_centrality.qs"))


message("\n\n##########################################################################\n",
        "# Finished G03a ", Sys.time(),
        "\n##########################################################################\n\n")