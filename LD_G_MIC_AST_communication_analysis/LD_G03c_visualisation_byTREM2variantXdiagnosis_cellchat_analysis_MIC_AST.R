message("\n\n##########################################################################\n",
        "# Start LD_G03c: CellChat visualisation by TREM2 variant and diagnosis ", Sys.time(),
        "\n##########################################################################\n\n")

library(qs)
library(tidyverse)
library(CellChat)
library(Seurat)
library(patchwork)
library(ComplexHeatmap)
library(pheatmap)
library(viridis)
library(colorRamps)
library(circlize)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

script_ind = "LD_G03c_v001_"

in_dir  = paste0(main_dir, "LD_G_MIC_AST_communication_analysis_output/")
out_dir = paste0(in_dir, "LD_G03c/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}


#input CellChat objects (G02c output)
# Control variant groups have few donors; interpret thresholded results cautiously.
group_levels = c("CV_Control", "R47H_Control", "R62H_Control",
                 "CV_AD",      "R47H_AD",      "R62H_AD")
in_cellchat = setNames(
  paste0(in_dir, "LD_G02c_v001_cellchat_TREM2Variant_x_Diagnosis_", group_levels, ".qs"),
  group_levels
)
in_bulkdata = paste0(in_dir, "LD_G01a_bulk_data_AST_MIC.qs")


#significance cutoff for custom MIC->AST plots and ranked CSVs (BH-FDR padj
#within MIC->AST hypothesis space, applied at L-R level)
padj_cutoff = 0.1



###########################################################
# functions
###########################################################

#custom colour palette
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
# 1. load CellChat objects + merged pseudobulk
###########################################################

message("\n\n          *** Load CellChat objects and merged pseudobulk... ", Sys.time(), "\n\n")


cc_list = lapply(in_cellchat, qread)
bulk_data = qread(file = in_bulkdata)

for (g in group_levels){
  cat(g, ": cells=", length(cc_list[[g]]@idents),
      "  clusters=", length(levels(cc_list[[g]]@idents)),
      "  pathways=", length(cc_list[[g]]@netP$pathways), "\n", sep = "")
}



###########################################################
# 2. define cluster ordering, apply per group
###########################################################

# Build cluster order from union across groups. Each group CellChat object may
# be missing some clusters (small subgroups + filterCommunication can drop them).

message("\n\n          *** Apply cluster ordering... ", Sys.time(), "\n\n")


all_clusters = sort(unique(unlist(lapply(cc_list, function(x) levels(x@idents)))))

ast_order = c(
  sort(grep("^AST_SLC1A2",  all_clusters, value = TRUE)),
  sort(grep("^AST_GFAP",    all_clusters, value = TRUE)),
  sort(grep("^AST_CHI3L1",  all_clusters, value = TRUE))
)
mic_order = c(
  sort(grep("^HOM_",        all_clusters, value = TRUE)),
  sort(grep("^DAM_",        all_clusters, value = TRUE)),
  sort(grep("^IRM_",        all_clusters, value = TRUE)),
  sort(grep("^HLA_",        all_clusters, value = TRUE)),
  sort(grep("^CRM_",        all_clusters, value = TRUE))
)
cluster_order = c(ast_order, mic_order)

unmatched = setdiff(all_clusters, cluster_order)
if (length(unmatched) > 0){
  warning("Clusters not matched by AST/MIC prefix - appended at end: ",
          paste(unmatched, collapse = ", "))
  cluster_order = c(cluster_order, sort(unmatched))
}

clusters_AST = ast_order
clusters_MIC = mic_order


### apply order to each CellChat object's idents + re-aggregate net
for (g in group_levels){
  this_order = intersect(cluster_order, levels(cc_list[[g]]@idents))
  cc_list[[g]]@idents = factor(cc_list[[g]]@idents, levels = this_order)
  cc_list[[g]] = aggregateNet(cc_list[[g]])
}


cat("\nCluster ordering (n = ", length(cluster_order), "):\n", sep = "")
cat("  AST (n = ", length(clusters_AST), "): ", paste(clusters_AST, collapse = ", "), "\n", sep = "")
cat("  MIC (n = ", length(clusters_MIC), "): ", paste(clusters_MIC, collapse = ", "), "\n", sep = "")


###########################################################
# 3. compute centrality per group
###########################################################

message("\n\n          *** Compute centrality per group... ", Sys.time(), "\n\n")


for (g in group_levels){
  cc_list[[g]] = netAnalysis_computeCentrality(cc_list[[g]], slot.name = "netP")
}



###########################################################
# 3b. compute BH-FDR padj at L-R level within MIC->AST hypothesis space
###########################################################

# CellChat stores L-R-level permutation p-values in @net$pval (3D array
# [src, tgt, LR]). @netP$pval does NOT exist - that slot only holds prob and
# pathways post-aggregation. We compute padj at L-R level, masked to MIC->AST.

message("\n\n          *** Compute BH-FDR padj for MIC->AST L-R edges... ", Sys.time(), "\n\n")


for (g in group_levels){
  pval_arr = cc_list[[g]]@net$pval
  if (is.null(pval_arr)){
    stop("@net$pval is NULL for ", g, " - cannot compute padj. ",
         "Check that computeCommunProb() ran with nboot > 0.")
  }

  src_idx = which(dimnames(pval_arr)[[1]] %in% clusters_MIC)
  tgt_idx = which(dimnames(pval_arr)[[2]] %in% clusters_AST)

  if (length(src_idx) == 0 || length(tgt_idx) == 0){
    cc_list[[g]]@net$padj = array(NA, dim = dim(pval_arr), dimnames = dimnames(pval_arr))
    cat(g, ": no MIC sources or AST targets - padj all NA\n", sep = "")
    next
  }

  mic_ast_pvals = as.vector(pval_arr[src_idx, tgt_idx, , drop = FALSE])
  mic_ast_padj  = p.adjust(mic_ast_pvals, method = "BH")

  padj_arr = array(NA, dim = dim(pval_arr), dimnames = dimnames(pval_arr))
  padj_arr[src_idx, tgt_idx, ] = array(mic_ast_padj,
                                        dim = c(length(src_idx), length(tgt_idx),
                                                dim(pval_arr)[3]))

  cc_list[[g]]@net$padj = padj_arr

  cat(g, ": MIC->AST L-R edges tested = ", length(mic_ast_pvals),
      "  | padj < ", padj_cutoff, " = ", sum(mic_ast_padj < padj_cutoff, na.rm = TRUE),
      "  | padj < 0.05 = ", sum(mic_ast_padj < 0.05, na.rm = TRUE), "\n", sep = "")
}



###########################################################
# 4. aggregated network circle plots
###########################################################

# 6 panels per page (2 rows x 3 cols), separate pages for count and weight.

message("\n\n          *** Aggregated network circle plots... ", Sys.time(), "\n\n")


pdf(file = paste0(out_dir, script_ind, "aggregated_network_circle_pval.pdf"),
    width = 18, height = 14)
{
  # page 1: count
  par(mfrow = c(2, 3), xpd = TRUE)
  for (g in group_levels){
    cc = cc_list[[g]]
    netVisual_circle(cc@net$count, vertex.weight = as.numeric(table(cc@idents)),
                     weight.scale = TRUE, label.edge = FALSE,
                     title.name = paste0(g, " - count"))
  }

  # page 2: weight
  par(mfrow = c(2, 3), xpd = TRUE)
  for (g in group_levels){
    cc = cc_list[[g]]
    netVisual_circle(cc@net$weight, vertex.weight = as.numeric(table(cc@idents)),
                     weight.scale = TRUE, label.edge = FALSE,
                     title.name = paste0(g, " - weight"))
  }
}
dev.off()



###########################################################
# 5. interaction heatmaps (count + weight)
###########################################################

# 6 panels in one row gets unreadable - split: row 1 = CV_Ctrl/CV_AD/R47H_Ctrl,
# row 2 = R47H_AD/R62H_Ctrl/R62H_AD. Two pages per measure (count, weight).

message("\n\n          *** Interaction heatmaps... ", Sys.time(), "\n\n")


# One group per page (12 pages total, 6 groups x 2 measures). Avoids mismatched-
# nrow errors when small groups (R47H_Control especially) have dropped clusters.
pdf(file = paste0(out_dir, script_ind, "aggregated_network_heatmap_pval.pdf"),
    width = 10, height = 9)
{
  for (g in group_levels){
    tryCatch({
      ht_count = netVisual_heatmap(cc_list[[g]], measure = "count", color.heatmap = "Reds",
                                    title.name = paste0(g, " - count"))
      draw(ht_count)
    }, error = function(e){ message("    section 5 count skipped for ", g, ": ", conditionMessage(e)) })

    tryCatch({
      ht_weight = netVisual_heatmap(cc_list[[g]], measure = "weight", color.heatmap = "Reds",
                                     title.name = paste0(g, " - weight"))
      draw(ht_weight)
    }, error = function(e){ message("    section 5 weight skipped for ", g, ": ", conditionMessage(e)) })
  }
}
dev.off()



###########################################################
# 6. signalling-role heatmaps - full pathway set
###########################################################

# Outgoing pages first, then incoming. Two pages per pattern (3 groups each)
# to keep panels readable.

message("\n\n          *** Signalling-role heatmaps (full)... ", Sys.time(), "\n\n")


# Pathway sets differ across groups; pathway dims can't be padded easily.
# One group per page; outgoing then incoming (12 pages total).
n_pw_max = max(sapply(cc_list, function(x) length(x@netP$pathways)))
heatmap_height_full = max(15, 0.25 * n_pw_max + 4)


pdf(file = paste0(out_dir, script_ind, "signallingRole_heatmap_full_pval.pdf"),
    width = 10, height = heatmap_height_full)
{
  for (g in group_levels){
    tryCatch({
      ht_out = netAnalysis_signalingRole_heatmap(cc_list[[g]], pattern = "outgoing",
                                                 width = 8, height = heatmap_height_full - 2,
                                                 color.heatmap = "BuGn",
                                                 title = paste0(g, " - outgoing"))
      draw(ht_out)
    }, error = function(e){ message("    section 6 outgoing skipped for ", g, ": ", conditionMessage(e)) })

    tryCatch({
      ht_in = netAnalysis_signalingRole_heatmap(cc_list[[g]], pattern = "incoming",
                                                width = 8, height = heatmap_height_full - 2,
                                                color.heatmap = "GnBu",
                                                title = paste0(g, " - incoming"))
      draw(ht_in)
    }, error = function(e){ message("    section 6 incoming skipped for ", g, ": ", conditionMessage(e)) })
  }
}
dev.off()



###########################################################
# 7. identify MIC->AST relevant pathways per group + union
###########################################################

message("\n\n          *** Identify MIC->AST relevant pathways... ", Sys.time(), "\n\n")


net_tabs           = list()
net_tabs_MIC_AST   = list()
pathways_MIC_AST   = list()


### helper: pull padj for each row of net_tab from cc@net$padj at L-R level
add_padj_column = function(tab, cc){
  padj_arr = cc@net$padj
  if (nrow(tab) == 0 || is.null(padj_arr)){
    tab$padj = if (nrow(tab) == 0) numeric(0) else NA_real_
    return(tab)
  }
  tab$padj = mapply(function(s, t, lr){
    if (s %in% dimnames(padj_arr)[[1]] &&
        t %in% dimnames(padj_arr)[[2]] &&
        lr %in% dimnames(padj_arr)[[3]]){
      padj_arr[s, t, lr]
    } else NA_real_
  }, tab$source, tab$target, tab$interaction_name)
  tab
}


for (g in group_levels){
  net_tabs[[g]]         = subsetCommunication(cc_list[[g]])
  net_tabs[[g]]         = add_padj_column(net_tabs[[g]], cc_list[[g]])

  net_tabs_MIC_AST[[g]] = net_tabs[[g]][net_tabs[[g]]$source %in% clusters_MIC &
                                         net_tabs[[g]]$target %in% clusters_AST, ]
  pathways_MIC_AST[[g]] = sort(unique(net_tabs_MIC_AST[[g]]$pathway_name))

  cat(g, ": total pathways=", length(cc_list[[g]]@netP$pathways),
      "  MIC->AST pathways=", length(pathways_MIC_AST[[g]]),
      "  MIC->AST L-Rs=", nrow(net_tabs_MIC_AST[[g]]), "\n", sep = "")
}

pathways_MIC_AST_union = sort(unique(unlist(pathways_MIC_AST)))
cat("\nUnion MIC->AST pathways across groups: ", length(pathways_MIC_AST_union), "\n")


### save tables
for (g in group_levels){
  write_csv(net_tabs[[g]],
            file = paste0(out_dir, script_ind, "interactions_all_", g, ".csv"))
  write_csv(net_tabs_MIC_AST[[g]],
            file = paste0(out_dir, script_ind, "interactions_MIC_to_AST_", g, ".csv"))
}



###########################################################
# 8. signalling-role heatmaps - MIC->AST pathways only
###########################################################

# Same as section 6 but restricted to union MIC->AST pathways.
# Shows total outgoing/incoming per cluster - not direction-restricted.
# For MIC->AST flow specifically, see section 8b.

message("\n\n          *** Signalling-role heatmaps (MIC->AST pathways)... ", Sys.time(), "\n\n")


if (length(pathways_MIC_AST_union) > 0){

  heatmap_height_sub = max(8, 0.3 * length(pathways_MIC_AST_union) + 2)

  # one group per page; pathway dimensions differ across groups
  pdf(file = paste0(out_dir, script_ind, "signallingRole_heatmap_MIC_AST_pathways_pval.pdf"),
      width = 10, height = heatmap_height_sub)
  {
    for (g in group_levels){
      pw_g = intersect(pathways_MIC_AST_union, cc_list[[g]]@netP$pathways)
      if (length(pw_g) < 2){
        message("    section 8 skipping ", g, ": <2 MIC->AST pathways present")
        next
      }
      tryCatch({
        ht_out = netAnalysis_signalingRole_heatmap(cc_list[[g]],
                                                   signaling = pw_g,
                                                   pattern = "outgoing",
                                                   width = 8, height = heatmap_height_sub - 2,
                                                   color.heatmap = "BuGn",
                                                   title = paste0(g, " - outgoing"))
        draw(ht_out)
      }, error = function(e){ message("    section 8 outgoing skipped for ", g, ": ", conditionMessage(e)) })

      tryCatch({
        ht_in = netAnalysis_signalingRole_heatmap(cc_list[[g]],
                                                  signaling = pw_g,
                                                  pattern = "incoming",
                                                  width = 8, height = heatmap_height_sub - 2,
                                                  color.heatmap = "GnBu",
                                                  title = paste0(g, " - incoming"))
        draw(ht_in)
      }, error = function(e){ message("    section 8 incoming skipped for ", g, ": ", conditionMessage(e)) })
    }
  }
  dev.off()
}



###########################################################
# 8b. NEW: MIC->AST direction-restricted signalling heatmap (custom)
###########################################################

# Custom heatmap built from cellchat@net$prob masked to (MIC sources x AST
# targets x significant L-Rs), then aggregated by pathway via cc@LR$LRsig.
# Rows = pathway; cols = MIC clusters (outgoing) or AST clusters (incoming);
# value = sum prob across significant MIC->AST L-Rs in that pathway.

message("\n\n          *** MIC->AST direction-restricted signalling heatmap... ", Sys.time(), "\n\n")


### sanity check: cc@LR$LRsig must exist
for (g in group_levels){
  if (is.null(cc_list[[g]]@LR$LRsig)){
    stop("@LR$LRsig is NULL for ", g, " - cannot map L-Rs to pathways. ",
         "Check G02c workflow.")
  }
}


build_MIC_AST_role_mats = function(cc, pathways_to_use, clusters_MIC, clusters_AST,
                                   padj_cutoff_use){
  # cc@net$prob and cc@net$padj are 3D arrays: [source, target, L-R]
  prob = cc@net$prob
  padj = cc@net$padj
  if (is.null(prob) || is.null(padj)) return(list(out = NULL, inc = NULL))

  lr_sig = cc@LR$LRsig
  lr_to_pw = setNames(lr_sig$pathway_name, lr_sig$interaction_name)

  src = intersect(clusters_MIC, dimnames(prob)[[1]])
  tgt = intersect(clusters_AST, dimnames(prob)[[2]])
  if (length(src) == 0 || length(tgt) == 0) return(list(out = NULL, inc = NULL))

  # mask non-significant L-Rs to zero
  mask = is.na(padj) | padj >= padj_cutoff_use
  prob[mask] = 0

  lr_names = dimnames(prob)[[3]]
  lr_pw    = lr_to_pw[lr_names]
  pw_keep  = intersect(pathways_to_use, unique(lr_pw))
  if (length(pw_keep) == 0) return(list(out = NULL, inc = NULL))

  prob_sub = prob[src, tgt, , drop = FALSE]

  prob_pw = array(0, dim = c(length(src), length(tgt), length(pw_keep)),
                  dimnames = list(src, tgt, pw_keep))
  for (pw in pw_keep){
    lrs_in_pw = lr_names[lr_pw == pw & !is.na(lr_pw)]
    if (length(lrs_in_pw) == 0) next
    prob_pw[, , pw] = apply(prob_sub[, , lrs_in_pw, drop = FALSE], c(1, 2), sum)
  }

  out_mat = apply(prob_pw, c(1, 3), sum); out_mat = t(out_mat)
  inc_mat = apply(prob_pw, c(2, 3), sum); inc_mat = t(inc_mat)

  return(list(out = out_mat, inc = inc_mat))
}


if (length(pathways_MIC_AST_union) > 0){

  role_mats = lapply(group_levels, function(g){
    build_MIC_AST_role_mats(cc_list[[g]], pathways_MIC_AST_union,
                            clusters_MIC, clusters_AST,
                            padj_cutoff_use = padj_cutoff)
  })
  names(role_mats) = group_levels

  ### DIAGNOSTIC: report total mass in each group's out/inc matrix.
  # Zero sums => padj masking filtered everything out for that group.
  # Non-zero sums but apparently empty plots => colour scale washed out by
  # very few hot cells (consider per-panel scale or quantile capping).
  cat("\n=== DIAGNOSTIC: role_mats sums per group ===\n")
  for (g in group_levels){
    out_sum = if (!is.null(role_mats[[g]]$out)) sum(role_mats[[g]]$out, na.rm = TRUE) else NA
    inc_sum = if (!is.null(role_mats[[g]]$inc)) sum(role_mats[[g]]$inc, na.rm = TRUE) else NA
    cat(g, ": out_sum =", out_sum, " | inc_sum =", inc_sum, "\n")
  }

  # split 6 groups into two rows of 3 for readable side-by-side panels
  groups_row1 = group_levels[1:3]   # CV_Control, CV_AD, R47H_Control
  groups_row2 = group_levels[4:6]   # R47H_AD, R62H_Control, R62H_AD


  heatmap_h = max(8, 0.3 * length(pathways_MIC_AST_union) + 2)

  pdf(file = paste0(out_dir, script_ind, "signallingRole_heatmap_MIC_to_AST_only_padj.pdf"),
      width = 22, height = heatmap_h)
  {
    # outgoing - row of 6 panels (3 + 3 across two pages for readability)
    for (groups_row in list(groups_row1, groups_row2)){
      hts_out = lapply(groups_row, function(g){
        m = role_mats[[g]]$out
        if (is.null(m) || nrow(m) < 2) return(NULL)
        Heatmap(m,
                name             = paste0("MIC->AST out\n", g),
                col              = colorRamp2(c(0, max(m, na.rm = TRUE)),
                                              c("white", "darkgreen")),
                cluster_rows     = FALSE,
                cluster_columns  = FALSE,
                row_names_side   = "left",
                column_names_rot = 45,
                column_title     = paste0(g, " - MIC->AST outgoing\n(per MIC cluster, summed across AST targets)\nSignificance: BH-FDR padj < ", padj_cutoff, " (MIC->AST L-R scope)"),
                row_names_gp     = gpar(fontsize = 8),
                column_names_gp  = gpar(fontsize = 8))
      })
      hts_out = hts_out[!sapply(hts_out, is.null)]
      if (length(hts_out) >= 1) draw(Reduce(`+`, hts_out), ht_gap = unit(0.5, "cm"))
    }

    # incoming
    for (groups_row in list(groups_row1, groups_row2)){
      hts_in = lapply(groups_row, function(g){
        m = role_mats[[g]]$inc
        if (is.null(m) || nrow(m) < 2) return(NULL)
        Heatmap(m,
                name             = paste0("MIC->AST inc\n", g),
                col              = colorRamp2(c(0, max(m, na.rm = TRUE)),
                                              c("white", "darkblue")),
                cluster_rows     = FALSE,
                cluster_columns  = FALSE,
                row_names_side   = "left",
                column_names_rot = 45,
                column_title     = paste0(g, " - MIC->AST incoming\n(per AST cluster, summed across MIC sources)\nSignificance: BH-FDR padj < ", padj_cutoff, " (MIC->AST L-R scope)"),
                row_names_gp     = gpar(fontsize = 8),
                column_names_gp  = gpar(fontsize = 8))
      })
      hts_in = hts_in[!sapply(hts_in, is.null)]
      if (length(hts_in) >= 1) draw(Reduce(`+`, hts_in), ht_gap = unit(0.5, "cm"))
    }
  }
  dev.off()
}



###########################################################
# 9. signalling-role scatter
###########################################################

# 6 panels per page (2 rows x 3 cols).

message("\n\n          *** Signalling-role scatter... ", Sys.time(), "\n\n")


pdf(file = paste0(out_dir, script_ind, "signallingRole_scatter_pval.pdf"),
    width = 18, height = 12)
{
  ps = lapply(group_levels, function(g){
    netAnalysis_signalingRole_scatter(cc_list[[g]]) +
      ggtitle(g)
  })
  print(wrap_plots(ps, nrow = 2, ncol = 3))
}
dev.off()



###########################################################
# 10. per-pathway plots: chord + role network (MIC->AST pathways)
###########################################################

# Multi-page PDF, one chord per page per group. Role network panel layout in
# 2x3 grid via base mfrow (base graphics, not circlize).

message("\n\n          *** Per-pathway chord + role network (MIC->AST pathways)... ", Sys.time(), "\n\n")


### chord plots - one per page per group (circlize doesn't honor par(mfrow))
pdf(file = paste0(out_dir, script_ind, "per_pathway_chord_MIC_AST_pathways_pval.pdf"),
    width = 9, height = 9)
{
  for (pw in pathways_MIC_AST_union){
    for (g in group_levels){
      cc = cc_list[[g]]
      if (pw %in% cc@netP$pathways){
        tryCatch({
          netVisual_aggregate(cc, signaling = pw, layout = "chord")
          title(main = paste0(g, " - ", pw), line = -1)
        }, error = function(e){
          plot.new(); title(main = paste0(g, " - ", pw, "\nchord error: ", conditionMessage(e)))
          message("    chord skipped for ", g, " ", pw, ": ", conditionMessage(e))
        })
      } else {
        plot.new(); title(main = paste0(g, " - ", pw, "\n(absent in this group)"))
      }
    }
  }
}
dev.off()


### role network heatmaps - 6 groups per pathway, 2x3 base graphics layout
pdf(file = paste0(out_dir, script_ind, "per_pathway_signallingRole_network_MIC_AST_pathways_pval.pdf"),
    width = 24, height = 8)
{
  for (pw in pathways_MIC_AST_union){
    par(mfrow = c(2, 3), xpd = TRUE)
    for (g in group_levels){
      cc = cc_list[[g]]
      if (pw %in% cc@netP$pathways){
        tryCatch({
          netAnalysis_signalingRole_network(cc, signaling = pw,
                                            width = 10, height = 2.5, font.size = 9)
          mtext(paste0(g, " - ", pw), side = 3, line = 0.5, cex = 0.9)
        }, error = function(e){
          plot.new(); title(main = paste0(g, " - ", pw, "\nrole network error"))
          message("    role network skipped for ", g, " ", pw, ": ", conditionMessage(e))
        })
      } else {
        plot.new(); title(main = paste0(g, " - ", pw, "\n(absent in this group)"))
      }
    }
  }
}
dev.off()



###########################################################
# 10b. NEW: MIC->AST direction-restricted per-pathway chord (custom)
###########################################################

# Custom chord: MIC senders, AST receivers, no autocrine.
# One page per group, two side-by-side panels per page:
#   left  = filtered by raw pval < pval_cutoff_use
#   right = filtered by BH-FDR padj < padj_cutoff_use
# Lets you see what's lost to FDR correction.

message("\n\n          *** MIC->AST-only per-pathway chord (custom)... ", Sys.time(), "\n\n")


pval_cutoff = 0.05


draw_MIC_AST_chord = function(cc, pw, clusters_MIC, clusters_AST,
                              mode = c("padj", "pval"),
                              ribbon_threshold = 0,
                              padj_cutoff_use = 0.1, pval_cutoff_use = 0.05){
  # Title and caption are added by the caller, not here.
  mode = match.arg(mode)

  prob = cc@net$prob
  if (is.null(prob)){
    plot.new()
    return(invisible(NULL))
  }

  # choose significance array based on mode
  if (mode == "padj"){
    sig_arr = cc@net$padj
    cutoff  = padj_cutoff_use
  } else {
    sig_arr = cc@net$pval
    cutoff  = pval_cutoff_use
  }
  if (is.null(sig_arr)){
    plot.new()
    return(invisible(NULL))
  }

  lr_sig = cc@LR$LRsig
  lr_in_pw = lr_sig$interaction_name[lr_sig$pathway_name == pw]
  lr_in_pw = intersect(lr_in_pw, dimnames(prob)[[3]])
  if (length(lr_in_pw) == 0){
    plot.new()
    mtext("(pathway absent)", side = 3, line = -2, cex = 0.7)
    return(invisible(NULL))
  }

  src = intersect(clusters_MIC, dimnames(prob)[[1]])
  tgt = intersect(clusters_AST, dimnames(prob)[[2]])
  if (length(src) == 0 || length(tgt) == 0){
    plot.new()
    mtext("(no MIC src / AST tgt)", side = 3, line = -2, cex = 0.7)
    return(invisible(NULL))
  }

  prob_sub = prob[src, tgt, lr_in_pw, drop = FALSE]
  sig_sub  = sig_arr[src, tgt, lr_in_pw, drop = FALSE]
  mask = is.na(sig_sub) | sig_sub >= cutoff
  prob_sub[mask] = 0

  m = apply(prob_sub, c(1, 2), sum)
  m[m <= ribbon_threshold] = 0
  if (sum(m) == 0){
    plot.new()
    msg = if (mode == "padj") paste0("(no L-R survives padj < ", padj_cutoff_use, ")")
          else                paste0("(no L-R survives pval < ", pval_cutoff_use, ")")
    mtext(msg, side = 3, line = -2, cex = 0.7)
    return(invisible(NULL))
  }

  df = expand.grid(from = rownames(m), to = colnames(m), stringsAsFactors = FALSE)
  df$value = as.vector(m)
  df = df[df$value > 0, ]
  if (nrow(df) == 0){
    plot.new()
    mtext("(no MIC->AST ribbons)", side = 3, line = -2, cex = 0.7)
    return(invisible(NULL))
  }

  # per-cluster colours: gradient within each lineage so individual sub-clusters
  # are distinguishable while still keeping the AST/MIC visual distinction
  ast_colors = colorRampPalette(c("dodgerblue4", "lightskyblue1"))(length(clusters_AST))
  mic_colors = colorRampPalette(c("darkorange3", "lightyellow"))(length(clusters_MIC))
  names(ast_colors) = clusters_AST
  names(mic_colors) = clusters_MIC

  grid.col = c(mic_colors[src], ast_colors[tgt])

  circos.clear()
  circos.par(gap.after = c(rep(2, length(src) - 1), 10,
                           rep(2, length(tgt) - 1), 10))
  chordDiagram(df,
               grid.col       = grid.col,
               directional    = 1,
               direction.type = c("diffHeight", "arrows"),
               link.arr.type  = "big.arrow",
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.05))
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y){
    sector.name = get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.5, sector.name,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                cex = 0.7)
  }, bg.border = NA)
  circos.clear()
  invisible(NULL)
}


if (length(pathways_MIC_AST_union) > 0){
  pdf(file = paste0(out_dir, script_ind, "per_pathway_chord_MIC_to_AST_only_pval_vs_padj.pdf"),
      width = 16, height = 9)
  {
    for (pw in pathways_MIC_AST_union){
      for (g in group_levels){
        # 1x2 layout: pval (left) vs padj (right)
        par(mfrow = c(1, 2), oma = c(3, 0, 3, 0))

        # left: raw pval
        tryCatch({
          draw_MIC_AST_chord(cc_list[[g]], pw, clusters_MIC, clusters_AST,
                             mode = "pval", pval_cutoff_use = pval_cutoff)
        }, error = function(e){
          plot.new()
          message("    pval chord error for ", g, " ", pw, ": ", conditionMessage(e))
        })
        title(main = paste0("raw p < ", pval_cutoff), line = 1, cex.main = 1)

        # right: padj
        tryCatch({
          draw_MIC_AST_chord(cc_list[[g]], pw, clusters_MIC, clusters_AST,
                             mode = "padj", padj_cutoff_use = padj_cutoff)
        }, error = function(e){
          plot.new()
          message("    padj chord error for ", g, " ", pw, ": ", conditionMessage(e))
        })
        title(main = paste0("BH-FDR padj < ", padj_cutoff), line = 1, cex.main = 1)

        # page-level title at top, caption at bottom
        mtext(paste0(g, " - ", pw), side = 3, line = 0.5, outer = TRUE, cex = 1.1, font = 2)
        mtext("MIC->AST L-R scope. Left: raw permutation p-value. Right: BH-FDR within MIC->AST hypothesis space.",
              side = 1, line = 1, outer = TRUE, cex = 0.7)
      }
    }
  }
  dev.off()
}



###########################################################
# 11. MIC->AST L-R rankings tables
###########################################################

# Per-group ranked tables + pathway-level summary. Ranked by padj ascending,
# then prob descending.

message("\n\n          *** MIC->AST L-R rankings... ", Sys.time(), "\n\n")


for (g in group_levels){

  net_g = net_tabs_MIC_AST[[g]]

  if (nrow(net_g) > 0){
    net_g_ranked = net_g[order(net_g$padj, -net_g$prob, na.last = TRUE), ]
    write_csv(net_g_ranked,
              file = paste0(out_dir, script_ind, "interactions_MIC_to_AST_ranked_", g, ".csv"))

    path_tab = net_g %>%
      group_by(source, target, pathway_name) %>%
      summarise(sum_prob = sum(prob), n_LR = n(), .groups = "drop") %>%
      arrange(desc(sum_prob))

    write_csv(path_tab,
              file = paste0(out_dir, script_ind, "pathway_by_cluster_pair_MIC_to_AST_", g, ".csv"))

    cat(g, " - top 5 MIC->AST source-target-pathway combinations:\n", sep = "")
    print(head(path_tab, 5))
  } else {
    message("    No MIC->AST L-Rs for ", g)
  }
}



###########################################################
# 12. MIC->AST L-R bubble plots
###########################################################

# Wider PDF, smaller text. One page per group.

message("\n\n          *** MIC->AST bubble plots... ", Sys.time(), "\n\n")


pdf(file = paste0(out_dir, script_ind, "bubble_MIC_to_AST_pval.pdf"),
    width = 24, height = 18)
{
  for (g in group_levels){
    if (length(pathways_MIC_AST[[g]]) == 0){
      message("    skipping bubble for ", g, " (no MIC->AST pathways)")
      next
    }
    tryCatch({
      p = netVisual_bubble(cc_list[[g]],
                           sources.use    = intersect(clusters_MIC, levels(cc_list[[g]]@idents)),
                           targets.use    = intersect(clusters_AST, levels(cc_list[[g]]@idents)),
                           signaling      = pathways_MIC_AST[[g]],
                           remove.isolate = TRUE,
                           angle.x        = 90,
                           font.size      = 6,
                           font.size.title = 9) +
        ggtitle(paste0("MIC -> AST - ", g))
      print(p)
    }, error = function(e){
      message("    bubble skipped for ", g, ": ", conditionMessage(e))
    })
  }
}
dev.off()


### custom ggplot bubble per group
pdf(file = paste0(out_dir, script_ind, "bubble_MIC_to_AST_custom_pval.pdf"),
    width = 26, height = 18)
{
  for (g in group_levels){
    net_g = net_tabs_MIC_AST[[g]]
    if (nrow(net_g) == 0) next

    t1 = net_g[order(-net_g$prob), ]
    t1$source_target = paste0(t1$source, " => ", t1$target)
    t1 = t1[order(t1$pathway_name, -t1$prob), ]

    p_custom = ggplot(t1, aes(x = source_target, y = interaction_name_2,
                              size = prob, color = pathway_name)) +
      geom_point() +
      scale_color_manual(limits = unique(t1$pathway_name),
                         values = pal(unique(t1$pathway_name))) +
      scale_x_discrete(limits = unique(t1$source_target)) +
      scale_y_discrete(limits = unique(t1$interaction_name_2)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
            axis.text.y = element_text(size = 5),
            legend.position = "right") +
      ggtitle(paste0("MIC -> AST - ", g))
    print(p_custom)
  }
}
dev.off()



###########################################################
# 13. per-pathway gene expression violins (custom, from G01b Seurat)
###########################################################

# Custom violin plots, independent of plotGeneExpression() which fails with
# 'S4SXP': should not happen errors under ggplot2 4.0+ via Seurat::VlnPlot.
#
# Design (decided 2026-05-03):
#   - one PDF page per (pathway, gene) combination
#   - x-axis = clusters relevant to the gene type:
#       * ligand / ligand_subunit / agonist / antagonist -> MIC clusters only
#       * receptor / receptor_subunit / co_A_receptor / co_I_receptor -> AST clusters only
#   - violins split/coloured by group (CV_Control, CV_AD, R47H_Control,
#     R47H_AD, R62H_Control, R62H_AD)
#   - genes per pathway = UNION across groups of L-Rs surviving padj < cutoff
#   - expression pulled from G01b merged Seurat (RNA log-normalised layer)

message("\n\n          *** Per-pathway gene expression violins (custom)... ", Sys.time(), "\n\n")


### load merged Seurat from G01b for gene expression
in_seur_merged = paste0(in_dir, "LD_G01b_seur_merged.qs")
if (!file.exists(in_seur_merged)){
  stop("Cannot find G01b merged Seurat for violins: ", in_seur_merged)
}
seur_for_vln = qread(in_seur_merged)


### fail fast if expected meta columns are missing
required_meta = c("TREM2Variant", "NeuropathologicalDiagnosis", "cluster_name")
missing_meta  = setdiff(required_meta, colnames(seur_for_vln@meta.data))
if (length(missing_meta) > 0){
  stop("seur_for_vln missing meta columns: ", paste(missing_meta, collapse = ", "))
}


### add a 'group' column matching G02c's split (TREM2Variant_x_Diagnosis)
seur_for_vln$group = paste0(as.character(seur_for_vln$TREM2Variant), "_",
                            as.character(seur_for_vln$NeuropathologicalDiagnosis))
seur_for_vln$group = factor(seur_for_vln$group, levels = group_levels)


### palette for groups (consistent across all violin pages)
group_palette = setNames(
  c("#9ecae1", "#2171b5",   # CV: light blue, dark blue
    "#fdae6b", "#d94801",   # R47H: light orange, dark orange
    "#a1d99b", "#238b45"),  # R62H: light green, dark green
  group_levels
)


### helper: gene-type -> which cluster set (MIC for ligand-side, AST for receptor-side)
gene_type_to_clusters = function(gene_type){
  if (gene_type %in% c("ligand", "ligand_subunit", "agonist", "antagonist")){
    return(clusters_MIC)
  } else if (gene_type %in% c("receptor", "receptor_subunit", "co_A_receptor", "co_I_receptor")){
    return(clusters_AST)
  } else {
    return(c(clusters_MIC, clusters_AST))   # fallback
  }
}


### build pathway -> union of significant L-Rs across groups
build_sig_lrs_per_pathway = function(net_tabs_MIC_AST, padj_cutoff){
  # returns named list: pathway_name -> character vector of interaction_name
  out = list()
  for (g in names(net_tabs_MIC_AST)){
    nt = net_tabs_MIC_AST[[g]]
    nt = nt[!is.na(nt$padj) & nt$padj < padj_cutoff, ]
    for (pw in unique(nt$pathway_name)){
      out[[pw]] = unique(c(out[[pw]],
                            nt$interaction_name[nt$pathway_name == pw]))
    }
  }
  out
}

sig_lrs_per_pw = build_sig_lrs_per_pathway(net_tabs_MIC_AST, padj_cutoff)


### build pathway -> gene list (using same gene-extraction logic as section 14)
### but restricted to the SIGNIFICANT L-Rs only (per pathway)
build_pathway_genes_from_lrs = function(lrs, db){
  # lrs: character vector of interaction_name within one pathway
  # db: cellchat@DB
  out = tibble(interaction_name = character(),
               gene_type        = character(),
               gene             = character())
  if (length(lrs) == 0) return(out)

  t3 = db$interaction[db$interaction$interaction_name %in% lrs, ]
  for (j in seq_len(nrow(t3))){
    ### cofactors
    for (cf in c("agonist", "antagonist", "co_A_receptor", "co_I_receptor")){
      cf_name = t3[[cf]][j]
      if (cf_name %in% rownames(db$cofactor)){
        t4 = db$cofactor[rownames(db$cofactor) == cf_name, ]
        out = rbind(out, tibble(interaction_name = t3$interaction_name[j],
                                gene_type = cf, gene = unlist(t4)))
      }
    }
    ### ligand
    if (t3$ligand[j] %in% rownames(db$complex)){
      t4 = db$complex[rownames(db$complex) == t3$ligand[j], ]
      out = rbind(out, tibble(interaction_name = t3$interaction_name[j],
                              gene_type = "ligand_subunit", gene = unlist(t4)))
    } else {
      out = rbind(out, tibble(interaction_name = t3$interaction_name[j],
                              gene_type = "ligand", gene = t3$ligand[j]))
    }
    ### receptor
    if (t3$receptor[j] %in% rownames(db$complex)){
      t4 = db$complex[rownames(db$complex) == t3$receptor[j], ]
      out = rbind(out, tibble(interaction_name = t3$interaction_name[j],
                              gene_type = "receptor_subunit", gene = unlist(t4)))
    } else {
      out = rbind(out, tibble(interaction_name = t3$interaction_name[j],
                              gene_type = "receptor", gene = t3$receptor[j]))
    }
  }
  out = out[out$gene != "" & !is.na(out$gene), ]
  out = out[!duplicated(paste0(out$gene, "|", out$gene_type)), ]
  out
}


### pull a gene's per-cell expression from RNA log-norm data layer
get_gene_expr = function(seur, gene){
  if (!(gene %in% rownames(seur))) return(NULL)
  as.numeric(GetAssayData(seur, assay = "RNA", layer = "data")[gene, ])
}


### plot one violin per (pathway, gene)
### use any group's @DB as reference (DBs are identical across groups)
db_ref = cc_list[[1]]@DB

cairo_pdf(filename = paste0(out_dir, script_ind, "per_pathway_gene_expression_MIC_AST_custom.pdf"),
          width = 14, height = 5, onefile = TRUE)
{
  for (pw in pathways_MIC_AST_union){
    sig_lrs = sig_lrs_per_pw[[pw]]
    if (is.null(sig_lrs) || length(sig_lrs) == 0){
      grid::grid.newpage()
      grid::grid.text(paste0(pw, "\n(no L-R survives padj < ", padj_cutoff, " in any group)"),
                      gp = grid::gpar(fontsize = 14))
      next
    }

    pg = build_pathway_genes_from_lrs(sig_lrs, db_ref)
    if (nrow(pg) == 0){
      grid::grid.newpage()
      grid::grid.text(paste0(pw, "\n(no genes resolved from significant L-Rs)"),
                      gp = grid::gpar(fontsize = 14))
      next
    }

    for (j in seq_len(nrow(pg))){
      gene      = pg$gene[j]
      gene_type = pg$gene_type[j]

      expr = get_gene_expr(seur_for_vln, gene)
      if (is.null(expr)){
        grid::grid.newpage()
        grid::grid.text(paste0(pw, " - ", gene, " (", gene_type, ")\n(gene not in Seurat)"),
                        gp = grid::gpar(fontsize = 12))
        next
      }

      clust_subset = gene_type_to_clusters(gene_type)
      meta = seur_for_vln@meta.data
      keep_cells = which(as.character(meta$cluster_name) %in% clust_subset)
      if (length(keep_cells) == 0){
        grid::grid.newpage()
        grid::grid.text(paste0(pw, " - ", gene, " (", gene_type, ")\n(no cells in relevant cluster set)"),
                        gp = grid::gpar(fontsize = 12))
        next
      }

      df = data.frame(
        expression   = expr[keep_cells],
        cluster_name = factor(as.character(meta$cluster_name[keep_cells]),
                              levels = intersect(cluster_order, clust_subset)),
        group        = meta$group[keep_cells]
      )
      df = df[!is.na(df$group), ]

      if (nrow(df) == 0 || all(df$expression == 0)){
        grid::grid.newpage()
        grid::grid.text(paste0(pw, " - ", gene, " (", gene_type, ")\n(no expression in relevant clusters)"),
                        gp = grid::gpar(fontsize = 12))
        next
      }

      p = ggplot(df, aes(x = cluster_name, y = expression, fill = group)) +
        geom_violin(scale = "width", trim = TRUE, position = position_dodge(width = 0.8),
                    alpha = 0.8, linewidth = 0.2) +
        scale_fill_manual(values = group_palette, drop = FALSE) +
        labs(x = NULL, y = "log-norm expression",
             title = paste0(pw, " - ", gene, " (", gene_type, ")"),
             fill  = "group") +
        theme_classic(base_size = 10) +
        theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
              axis.text.y     = element_text(size = 8),
              legend.position = "right",
              plot.title      = element_text(size = 11, face = "bold"))

      tryCatch(print(p),
               error = function(e){
                 grid::grid.newpage()
                 grid::grid.text(paste0(pw, " - ", gene, " (", gene_type, ")\nrender error:\n",
                                        conditionMessage(e)),
                                 gp = grid::gpar(fontsize = 11))
                 message("    violin print failed for ", pw, " ", gene, ": ",
                         conditionMessage(e))
               })
    }
  }
}
dev.off()


### free Seurat memory before section 14
rm(seur_for_vln); gc()



###########################################################
# 14. pathway-gene Z-score heatmaps (uses G01a merged pseudobulk)
###########################################################

# Z-score matrix is joint across all cluster_samples (from G01a) - not group-
# specific - so this section produces ONE PDF per pathway, looped over union.

message("\n\n          *** Pathway-gene Z-score heatmaps... ", Sys.time(), "\n\n")


if (length(pathways_MIC_AST_union) > 0){

  pathway_genes_tab = tibble(pathway = character(),
                             interaction_name = character(),
                             gene_type = character(),
                             gene = character())

  for (pw in pathways_MIC_AST_union){
    # find first group containing this pathway - DBs identical across groups
    g_use = group_levels[sapply(group_levels, function(g) pw %in% cc_list[[g]]@netP$pathways)][1]
    if (is.na(g_use)) next

    cc = cc_list[[g_use]]
    net_g = net_tabs_MIC_AST[[g_use]]
    path_ints = unique(net_g$interaction_name[net_g$pathway_name == pw])

    db = cc@DB
    t2 = db$interaction
    t3 = t2[t2$interaction_name %in% path_ints, ]

    for (j in seq_len(nrow(t3))){
      ### cofactors
      for (cf in c("agonist", "antagonist", "co_A_receptor", "co_I_receptor")){
        cf_name = t3[[cf]][j]
        if (cf_name %in% rownames(db$cofactor)){
          t4 = db$cofactor[rownames(db$cofactor) == cf_name, ]
          pathway_genes_tab = rbind(pathway_genes_tab,
                                    tibble(pathway = pw,
                                           interaction_name = t3$interaction_name[j],
                                           gene_type = cf, gene = unlist(t4)))
        }
      }
      ### ligand
      if (t3$ligand[j] %in% rownames(db$complex)){
        t4 = db$complex[rownames(db$complex) == t3$ligand[j], ]
        pathway_genes_tab = rbind(pathway_genes_tab,
                                  tibble(pathway = pw,
                                         interaction_name = t3$interaction_name[j],
                                         gene_type = "ligand_subunit", gene = unlist(t4)))
      } else {
        pathway_genes_tab = rbind(pathway_genes_tab,
                                  tibble(pathway = pw,
                                         interaction_name = t3$interaction_name[j],
                                         gene_type = "ligand", gene = t3$ligand[j]))
      }
      ### receptor
      if (t3$receptor[j] %in% rownames(db$complex)){
        t4 = db$complex[rownames(db$complex) == t3$receptor[j], ]
        pathway_genes_tab = rbind(pathway_genes_tab,
                                  tibble(pathway = pw,
                                         interaction_name = t3$interaction_name[j],
                                         gene_type = "receptor_subunit", gene = unlist(t4)))
      } else {
        pathway_genes_tab = rbind(pathway_genes_tab,
                                  tibble(pathway = pw,
                                         interaction_name = t3$interaction_name[j],
                                         gene_type = "receptor", gene = t3$receptor[j]))
      }
    }
  }

  pathway_genes_tab = pathway_genes_tab[pathway_genes_tab$gene != "" & !is.na(pathway_genes_tab$gene), ]

  write_csv(pathway_genes_tab,
            file = paste0(out_dir, script_ind, "pathway_genes_table_MIC_AST.csv"))

  cat("Pathway-gene table: ", nrow(pathway_genes_tab), " rows across ",
      length(unique(pathway_genes_tab$pathway)), " pathways\n", sep = "")


  ### plot heatmaps
  z_mat = bulk_data$gene_Z_scores$clusters_combined
  meta  = bulk_data$meta

  meta_ord = meta[order(match(meta$lineage, c("MIC", "AST")),
                        match(meta$cluster_name, cluster_order),
                        meta$sample), ]
  z_mat = z_mat[, meta_ord$cluster_sample, drop = FALSE]

  annot_col = data.frame(
    lineage                    = meta_ord$lineage,
    cluster_name               = meta_ord$cluster_name,
    TREM2Variant               = meta_ord$TREM2Variant,
    NeuropathologicalDiagnosis = meta_ord$NeuropathologicalDiagnosis,
    row.names                  = meta_ord$cluster_sample
  )

  annot_colors = list(
    lineage                    = c(AST = "dodgerblue", MIC = "orange"),
    TREM2Variant               = setNames(pal(levels(meta_ord$TREM2Variant)),
                                          levels(meta_ord$TREM2Variant)),
    NeuropathologicalDiagnosis = setNames(pal(levels(meta_ord$NeuropathologicalDiagnosis)),
                                          levels(meta_ord$NeuropathologicalDiagnosis))
  )


  pdf(file = paste0(out_dir, script_ind, "pathway_gene_heatmaps_MIC_AST.pdf"),
      width = 22, height = 8)
  {
    for (pw in pathways_MIC_AST_union){
      t1 = pathway_genes_tab[pathway_genes_tab$pathway == pw, ]
      t1 = t1[!duplicated(paste0(t1$gene, t1$gene_type)), ]
      t1 = t1[t1$gene %in% rownames(z_mat), ]

      if (nrow(t1) < 2){
        message("    skip ", pw, " - <2 genes available")
        next
      }

      pl_mat = z_mat[t1$gene, , drop = FALSE]
      rownames(pl_mat) = paste0(t1$gene, " (", t1$gene_type, ")")

      if (all(is.nan(pl_mat))){
        message("    skip ", pw, " - all genes constant")
        next
      }

      lims = 0.7 * c(-max(abs(pl_mat), na.rm = TRUE),
                      max(abs(pl_mat), na.rm = TRUE))

      tryCatch({
        pheatmap::pheatmap(
          pl_mat,
          cluster_rows = FALSE, cluster_cols = FALSE,
          show_rownames = TRUE, show_colnames = FALSE,
          annotation_col = annot_col, annotation_colors = annot_colors,
          color = colorRampPalette(c("blue", "white", "red"))(250),
          breaks = seq(lims[1], lims[2], length.out = 251),
          border_color = NA,
          fontsize = 8, fontsize_row = 6,
          main = paste0(pw, " - ", nrow(pl_mat), " genes")
        )
      }, error = function(e){
        message("    heatmap skipped for ", pw, ": ", conditionMessage(e))
      })
    }
  }
  dev.off()
}



###########################################################
# 15. save updated CellChat objects (with centrality + padj)
###########################################################

message("\n\n          *** Save updated CellChat objects... ", Sys.time(), "\n\n")


for (g in group_levels){
  qsave(cc_list[[g]],
        file = paste0(out_dir, script_ind, "cellchat_TREM2Variant_x_Diagnosis_", g, "_with_centrality.qs"))
}


message("\n\n##########################################################################\n",
        "# Finished G03c ", Sys.time(),
        "\n##########################################################################\n\n")
