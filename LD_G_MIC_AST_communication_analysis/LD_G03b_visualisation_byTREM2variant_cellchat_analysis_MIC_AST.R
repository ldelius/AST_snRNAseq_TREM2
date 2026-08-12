message("\n\n##########################################################################\n",
        "# Start LD_G03b: CellChat visualisation by TREM2 variant ", Sys.time(),
        "\n##########################################################################\n\n")

library(qs)
library(tidyverse)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(pheatmap)
library(viridis)
library(colorRamps)
library(circlize)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "LD_G03b_v001_"

#specify input + output directories
in_dir  = paste0(main_dir, "LD_G_MIC_AST_communication_analysis_output/")
out_dir = paste0(in_dir, "LD_G03b/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}


#input CellChat objects (G02b output)
variant_levels = c("CV", "R47H", "R62H")
in_cellchat = setNames(
  paste0(in_dir, "LD_G02b_v001_cellchat_TREM2Variant_", variant_levels, ".qs"),
  variant_levels
)
in_bulkdata = paste0(in_dir, "LD_G01a_bulk_data_AST_MIC.qs")


#significance cutoff for custom MIC->AST plots and ranked CSVs (BH-FDR padj
#within MIC->AST hypothesis space, applied at L-R level)
padj_cutoff = 0.1



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
# 1. load CellChat objects + merged pseudobulk
###########################################################

message("\n\n          *** Load CellChat objects and merged pseudobulk... ", Sys.time(), "\n\n")


cc_list = lapply(in_cellchat, qread)
bulk_data = qread(file = in_bulkdata)

for (v in variant_levels){
  cat(v, ": cells=", length(cc_list[[v]]@idents),
      "  clusters=", length(levels(cc_list[[v]]@idents)),
      "  pathways=", length(cc_list[[v]]@netP$pathways), "\n", sep = "")
}



###########################################################
# 2. define cluster ordering, apply per variant
###########################################################

# Build cluster order from the union of clusters across all variants. Each
# variant CellChat object may be missing some clusters (small subgroups +
# filterCommunication can drop them) - apply cluster_order intersected with
# each variant's idents.

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
for (v in variant_levels){
  this_order = intersect(cluster_order, levels(cc_list[[v]]@idents))
  cc_list[[v]]@idents = factor(cc_list[[v]]@idents, levels = this_order)
  cc_list[[v]] = aggregateNet(cc_list[[v]])
}


cat("\nCluster ordering (n = ", length(cluster_order), "):\n", sep = "")
cat("  AST (n = ", length(clusters_AST), "): ", paste(clusters_AST, collapse = ", "), "\n", sep = "")
cat("  MIC (n = ", length(clusters_MIC), "): ", paste(clusters_MIC, collapse = ", "), "\n", sep = "")


###########################################################
# 3. compute centrality per variant
###########################################################

message("\n\n          *** Compute centrality per variant... ", Sys.time(), "\n\n")


for (v in variant_levels){
  cc_list[[v]] = netAnalysis_computeCentrality(cc_list[[v]], slot.name = "netP")
}



###########################################################
# 3b. compute BH-FDR padj at L-R level within MIC->AST hypothesis space
###########################################################

# CellChat stores L-R-level permutation p-values in @net$pval (4D array:
# [source, target, L-R, ...] - usually [src, tgt, LR]). @netP$pval does NOT
# exist - that slot only holds prob and pathways post-aggregation.
#
# We compute BH-FDR padj at the L-R level: extract pvals for (MIC source,
# AST target, all L-Rs), adjust within that subset, write back to a new
# @net$padj array. Used by sections 8b, 10b, and the CSVs only.

message("\n\n          *** Compute BH-FDR padj for MIC->AST L-R edges... ", Sys.time(), "\n\n")


for (v in variant_levels){
  pval_arr = cc_list[[v]]@net$pval
  if (is.null(pval_arr)){
    stop("@net$pval is NULL for ", v, " - cannot compute padj. ",
         "Check that computeCommunProb() ran with nboot > 0.")
  }

  src_idx = which(dimnames(pval_arr)[[1]] %in% clusters_MIC)
  tgt_idx = which(dimnames(pval_arr)[[2]] %in% clusters_AST)

  if (length(src_idx) == 0 || length(tgt_idx) == 0){
    cc_list[[v]]@net$padj = array(NA, dim = dim(pval_arr), dimnames = dimnames(pval_arr))
    cat(v, ": no MIC sources or AST targets - padj all NA\n", sep = "")
    next
  }

  mic_ast_pvals = as.vector(pval_arr[src_idx, tgt_idx, , drop = FALSE])
  mic_ast_padj  = p.adjust(mic_ast_pvals, method = "BH")

  padj_arr = array(NA, dim = dim(pval_arr), dimnames = dimnames(pval_arr))
  padj_arr[src_idx, tgt_idx, ] = array(mic_ast_padj,
                                        dim = c(length(src_idx), length(tgt_idx),
                                                dim(pval_arr)[3]))

  cc_list[[v]]@net$padj = padj_arr

  cat(v, ": MIC->AST L-R edges tested = ", length(mic_ast_pvals),
      "  | padj < ", padj_cutoff, " = ", sum(mic_ast_padj < padj_cutoff, na.rm = TRUE),
      "  | padj < 0.05 = ", sum(mic_ast_padj < 0.05, na.rm = TRUE), "\n", sep = "")
}



###########################################################
# 4. aggregated network circle plots
###########################################################

# Three panels per page (CV / R47H / R62H), separate pages for count vs weight.

message("\n\n          *** Aggregated network circle plots... ", Sys.time(), "\n\n")


pdf(file = paste0(out_dir, script_ind, "aggregated_network_circle_pval.pdf"),
    width = 18, height = 7)
{
  # page 1: count
  par(mfrow = c(1, 3), xpd = TRUE)
  for (v in variant_levels){
    cc = cc_list[[v]]
    netVisual_circle(cc@net$count, vertex.weight = as.numeric(table(cc@idents)),
                     weight.scale = TRUE, label.edge = FALSE,
                     title.name = paste0(v, " - Number of interactions"))
  }

  # page 2: weight
  par(mfrow = c(1, 3), xpd = TRUE)
  for (v in variant_levels){
    cc = cc_list[[v]]
    netVisual_circle(cc@net$weight, vertex.weight = as.numeric(table(cc@idents)),
                     weight.scale = TRUE, label.edge = FALSE,
                     title.name = paste0(v, " - Interaction weights / strength"))
  }
}
dev.off()



###########################################################
# 5. interaction heatmaps (count + weight)
###########################################################

# CellChat-native heatmap. Three panels per row (CV/R47H/R62H), one row for
# count, one row for weight.

message("\n\n          *** Interaction heatmaps... ", Sys.time(), "\n\n")


# cluster sets differ across variants (filterCommunication can drop clusters
# with too few cells), so heatmap concatenation by `+` errors out on row mismatch.
# Plot one variant per page (count first, then weight).
pdf(file = paste0(out_dir, script_ind, "aggregated_network_heatmap_pval.pdf"),
    width = 10, height = 10)
{
  for (v in variant_levels){
    tryCatch({
      ht = netVisual_heatmap(cc_list[[v]], measure = "count",  color.heatmap = "Reds",
                             title.name = paste0(v, " - count"))
      draw(ht)
    }, error = function(e){ message("    section 5 count skipped for ", v, ": ", conditionMessage(e)) })
  }

  for (v in variant_levels){
    tryCatch({
      ht = netVisual_heatmap(cc_list[[v]], measure = "weight", color.heatmap = "Reds",
                             title.name = paste0(v, " - weight"))
      draw(ht)
    }, error = function(e){ message("    section 5 weight skipped for ", v, ": ", conditionMessage(e)) })
  }
}
dev.off()



###########################################################
# 6. signalling-role heatmaps - full pathway set
###########################################################

# Outgoing + incoming per variant. Outgoing pages first, then incoming pages.
# Variants side-by-side as separate panels per page.

message("\n\n          *** Signalling-role heatmaps (full)... ", Sys.time(), "\n\n")


# pathway counts vary by variant - so do cluster sets after filterCommunication
# Concatenating heatmaps requires matching dims; pathway dimensions can't be
# easily padded across variants, so plot one variant per page (separate pages
# for outgoing and incoming, two pages per variant).
n_pw_max = max(sapply(cc_list, function(x) length(x@netP$pathways)))
heatmap_height_full = max(15, 0.25 * n_pw_max + 4)


pdf(file = paste0(out_dir, script_ind, "signallingRole_heatmap_full_pval.pdf"),
    width = 10, height = heatmap_height_full)
{
  for (v in variant_levels){
    tryCatch({
      ht_out = netAnalysis_signalingRole_heatmap(cc_list[[v]], pattern = "outgoing",
                                                 width = 8, height = heatmap_height_full - 2,
                                                 color.heatmap = "BuGn",
                                                 title = paste0(v, " - outgoing"))
      draw(ht_out)
    }, error = function(e){ message("    section 6 outgoing skipped for ", v, ": ", conditionMessage(e)) })

    tryCatch({
      ht_in = netAnalysis_signalingRole_heatmap(cc_list[[v]], pattern = "incoming",
                                                width = 8, height = heatmap_height_full - 2,
                                                color.heatmap = "GnBu",
                                                title = paste0(v, " - incoming"))
      draw(ht_in)
    }, error = function(e){ message("    section 6 incoming skipped for ", v, ": ", conditionMessage(e)) })
  }
}
dev.off()



###########################################################
# 7. identify MIC->AST relevant pathways per variant + union
###########################################################

# Per-variant pathway set + union across variants. The union is used for
# section 8/10/12 panels so the same row order is comparable across variants.

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


for (v in variant_levels){
  net_tabs[[v]]         = subsetCommunication(cc_list[[v]])
  net_tabs[[v]]         = add_padj_column(net_tabs[[v]], cc_list[[v]])

  net_tabs_MIC_AST[[v]] = net_tabs[[v]][net_tabs[[v]]$source %in% clusters_MIC &
                                         net_tabs[[v]]$target %in% clusters_AST, ]
  pathways_MIC_AST[[v]] = sort(unique(net_tabs_MIC_AST[[v]]$pathway_name))

  cat(v, ": total pathways=", length(cc_list[[v]]@netP$pathways),
      "  MIC->AST pathways=", length(pathways_MIC_AST[[v]]),
      "  MIC->AST L-Rs=", nrow(net_tabs_MIC_AST[[v]]), "\n", sep = "")
}

pathways_MIC_AST_union = sort(unique(unlist(pathways_MIC_AST)))
cat("\nUnion MIC->AST pathways across variants: ", length(pathways_MIC_AST_union), "\n")


### save tables
for (v in variant_levels){
  write_csv(net_tabs[[v]],
            file = paste0(out_dir, script_ind, "interactions_all_", v, ".csv"))
  write_csv(net_tabs_MIC_AST[[v]],
            file = paste0(out_dir, script_ind, "interactions_MIC_to_AST_", v, ".csv"))
}



###########################################################
# 8. signalling-role heatmaps - MIC->AST pathways only
###########################################################

# Restricted to the union of MIC-to-AST pathways. Values are total outgoing and
# incoming signalling per cluster; direction-specific flow is shown in section 8b.

message("\n\n          *** Signalling-role heatmaps (MIC->AST pathways)... ", Sys.time(), "\n\n")


if (length(pathways_MIC_AST_union) > 0){

  heatmap_height_sub = max(8, 0.3 * length(pathways_MIC_AST_union) + 2)

  # one variant per page; pathway dimensions differ across variants
  pdf(file = paste0(out_dir, script_ind, "signallingRole_heatmap_MIC_AST_pathways_pval.pdf"),
      width = 10, height = heatmap_height_sub)
  {
    for (v in variant_levels){
      pw_v = intersect(pathways_MIC_AST_union, cc_list[[v]]@netP$pathways)
      if (length(pw_v) < 2){
        message("    section 8 skipping ", v, ": <2 MIC->AST pathways present")
        next
      }
      tryCatch({
        ht_out = netAnalysis_signalingRole_heatmap(cc_list[[v]],
                                                   signaling = pw_v,
                                                   pattern = "outgoing",
                                                   width = 8, height = heatmap_height_sub - 2,
                                                   color.heatmap = "BuGn",
                                                   title = paste0(v, " - outgoing"))
        draw(ht_out)
      }, error = function(e){ message("    section 8 outgoing skipped for ", v, ": ", conditionMessage(e)) })

      tryCatch({
        ht_in = netAnalysis_signalingRole_heatmap(cc_list[[v]],
                                                  signaling = pw_v,
                                                  pattern = "incoming",
                                                  width = 8, height = heatmap_height_sub - 2,
                                                  color.heatmap = "GnBu",
                                                  title = paste0(v, " - incoming"))
        draw(ht_in)
      }, error = function(e){ message("    section 8 incoming skipped for ", v, ": ", conditionMessage(e)) })
    }
  }
  dev.off()
}



###########################################################
# 8b. NEW: MIC->AST direction-restricted signalling heatmap (custom)
###########################################################

# Custom heatmap built from cellchat@netP$prob masked to MIC sources x AST
# targets, summed across (source, target) pairs to give per-pathway-per-cluster
# signal. Two panels per variant: outgoing (rows = pathways, cols = MIC clusters,
# value = sum prob to any AST target) and incoming (rows = pathways, cols = AST
# clusters, value = sum prob from any MIC source).

message("\n\n          *** MIC->AST direction-restricted signalling heatmap... ", Sys.time(), "\n\n")


### sanity check: cc@LR$LRsig must exist (interaction_name -> pathway_name lookup)
### populated by computeCommunProb() in standard CellChat workflows
for (v in variant_levels){
  if (is.null(cc_list[[v]]@LR$LRsig)){
    stop("@LR$LRsig is NULL for ", v, " - cannot map L-Rs to pathways. ",
         "Check G02b workflow.")
  }
}


build_MIC_AST_role_mats = function(cc, pathways_to_use, clusters_MIC, clusters_AST,
                                   padj_cutoff_use){
  # cc@net$prob and cc@net$padj are 3D arrays at L-R level: [source, target, L-R]
  # We aggregate L-Rs to pathway via cc@LR$LRsig (interaction_name -> pathway_name)
  prob = cc@net$prob
  padj = cc@net$padj
  if (is.null(prob) || is.null(padj)) return(list(out = NULL, inc = NULL))

  # interaction_name -> pathway_name lookup
  lr_sig = cc@LR$LRsig
  lr_to_pw = setNames(lr_sig$pathway_name, lr_sig$interaction_name)

  # restrict source/target dims to MIC sources / AST targets present
  src = intersect(clusters_MIC, dimnames(prob)[[1]])
  tgt = intersect(clusters_AST, dimnames(prob)[[2]])
  if (length(src) == 0 || length(tgt) == 0) return(list(out = NULL, inc = NULL))

  # mask non-significant edges to zero
  mask = is.na(padj) | padj >= padj_cutoff_use
  prob[mask] = 0

  # per-pathway: sum L-R probs that map to that pathway
  lr_names = dimnames(prob)[[3]]
  lr_pw    = lr_to_pw[lr_names]
  pw_keep  = intersect(pathways_to_use, unique(lr_pw))
  if (length(pw_keep) == 0) return(list(out = NULL, inc = NULL))

  prob_sub = prob[src, tgt, , drop = FALSE]

  # aggregate over L-Rs per pathway to get [src, tgt, pathway]
  prob_pw = array(0, dim = c(length(src), length(tgt), length(pw_keep)),
                  dimnames = list(src, tgt, pw_keep))
  for (pw in pw_keep){
    lrs_in_pw = lr_names[lr_pw == pw & !is.na(lr_pw)]
    if (length(lrs_in_pw) == 0) next
    prob_pw[, , pw] = apply(prob_sub[, , lrs_in_pw, drop = FALSE], c(1, 2), sum)
  }

  # outgoing: per MIC source per pathway, summed across AST targets
  out_mat = apply(prob_pw, c(1, 3), sum)   # [src, pathway]
  out_mat = t(out_mat)                       # rows = pathway, cols = MIC source

  # incoming: per AST target per pathway, summed across MIC sources
  inc_mat = apply(prob_pw, c(2, 3), sum)   # [tgt, pathway]
  inc_mat = t(inc_mat)                       # rows = pathway, cols = AST target

  return(list(out = out_mat, inc = inc_mat))
}


if (length(pathways_MIC_AST_union) > 0){

  # build per-variant matrices
  role_mats = lapply(variant_levels, function(v){
    build_MIC_AST_role_mats(cc_list[[v]], pathways_MIC_AST_union,
                            clusters_MIC, clusters_AST,
                            padj_cutoff_use = padj_cutoff)
  })
  names(role_mats) = variant_levels

  # plot - one row of panels for outgoing, one for incoming
  heatmap_h = max(8, 0.3 * length(pathways_MIC_AST_union) + 2)

  pdf(file = paste0(out_dir, script_ind, "signallingRole_heatmap_MIC_to_AST_only_padj.pdf"),
      width = 20, height = heatmap_h)
  {
    # outgoing: one heatmap per variant
    hts_out = lapply(variant_levels, function(v){
      m = role_mats[[v]]$out
      if (is.null(m) || nrow(m) < 2) return(NULL)
      Heatmap(m,
              name             = paste0("MIC->AST out\n", v),
              col              = colorRamp2(c(0, max(m, na.rm = TRUE)),
                                            c("white", "darkgreen")),
              cluster_rows     = FALSE,
              cluster_columns  = FALSE,
              row_names_side   = "left",
              column_names_rot = 45,
              column_title     = paste0(v, " - MIC->AST outgoing\n(per MIC cluster, summed across AST targets)\nSignificance: BH-FDR padj < ", padj_cutoff, " (MIC->AST L-R scope)"),
              row_names_gp     = gpar(fontsize = 8),
              column_names_gp  = gpar(fontsize = 8))
    })
    hts_out = hts_out[!sapply(hts_out, is.null)]
    if (length(hts_out) >= 1) draw(Reduce(`+`, hts_out), ht_gap = unit(0.5, "cm"))

    # incoming: one heatmap per variant
    hts_in = lapply(variant_levels, function(v){
      m = role_mats[[v]]$inc
      if (is.null(m) || nrow(m) < 2) return(NULL)
      Heatmap(m,
              name             = paste0("MIC->AST inc\n", v),
              col              = colorRamp2(c(0, max(m, na.rm = TRUE)),
                                            c("white", "darkblue")),
              cluster_rows     = FALSE,
              cluster_columns  = FALSE,
              row_names_side   = "left",
              column_names_rot = 45,
              column_title     = paste0(v, " - MIC->AST incoming\n(per AST cluster, summed across MIC sources)\nSignificance: BH-FDR padj < ", padj_cutoff, " (MIC->AST L-R scope)"),
              row_names_gp     = gpar(fontsize = 8),
              column_names_gp  = gpar(fontsize = 8))
    })
    hts_in = hts_in[!sapply(hts_in, is.null)]
    if (length(hts_in) >= 1) draw(Reduce(`+`, hts_in), ht_gap = unit(0.5, "cm"))
  }
  dev.off()
}



###########################################################
# 9. signalling-role scatter
###########################################################

# 2D senders vs receivers per variant. Side-by-side panels.

message("\n\n          *** Signalling-role scatter... ", Sys.time(), "\n\n")


pdf(file = paste0(out_dir, script_ind, "signallingRole_scatter_pval.pdf"),
    width = 18, height = 6)
{
  ps = lapply(variant_levels, function(v){
    netAnalysis_signalingRole_scatter(cc_list[[v]]) +
      ggtitle(v)
  })
  print(wrap_plots(ps, nrow = 1))
}
dev.off()



###########################################################
# 10. per-pathway plots: chord + role network (MIC->AST pathways)
###########################################################

# Multi-page PDF: for each pathway in the union, three side-by-side panels
# (one per variant). Pathways absent in a variant get a placeholder.
# These are the standard CellChat per-pathway plots and show all directions
# (not direction-restricted - see section 10b for MIC->AST-only chords).

message("\n\n          *** Per-pathway chord + role network (MIC->AST pathways)... ", Sys.time(), "\n\n")


### chord plots - one per page per variant (circlize doesn't honor par(mfrow))
pdf(file = paste0(out_dir, script_ind, "per_pathway_chord_MIC_AST_pathways_pval.pdf"),
    width = 9, height = 9)
{
  for (pw in pathways_MIC_AST_union){
    for (v in variant_levels){
      cc = cc_list[[v]]
      if (pw %in% cc@netP$pathways){
        tryCatch({
          netVisual_aggregate(cc, signaling = pw, layout = "chord")
          title(main = paste0(v, " - ", pw), line = -1)
        }, error = function(e){
          plot.new(); title(main = paste0(v, " - ", pw, "\nchord error: ", conditionMessage(e)))
          message("    chord skipped for ", v, " ", pw, ": ", conditionMessage(e))
        })
      } else {
        plot.new(); title(main = paste0(v, " - ", pw, "\n(absent in this variant)"))
      }
    }
  }
}
dev.off()


### role network heatmaps
pdf(file = paste0(out_dir, script_ind, "per_pathway_signallingRole_network_MIC_AST_pathways_pval.pdf"),
    width = 36, height = 4)
{
  for (pw in pathways_MIC_AST_union){
    par(mfrow = c(1, 3), xpd = TRUE)
    for (v in variant_levels){
      cc = cc_list[[v]]
      if (pw %in% cc@netP$pathways){
        tryCatch({
          netAnalysis_signalingRole_network(cc, signaling = pw,
                                            width = 10, height = 2.5, font.size = 9)
          mtext(paste0(v, " - ", pw), side = 3, line = 0.5, cex = 0.9)
        }, error = function(e){
          plot.new(); title(main = paste0(v, " - ", pw, "\nrole network error"))
          message("    role network skipped for ", v, " ", pw, ": ", conditionMessage(e))
        })
      } else {
        plot.new(); title(main = paste0(v, " - ", pw, "\n(absent in this variant)"))
      }
    }
  }
}
dev.off()



###########################################################
# 10b. NEW: MIC->AST direction-restricted per-pathway chord (custom)
###########################################################

# Custom chord using circlize: for each pathway, draw a chord with MIC clusters
# on one side and AST clusters on the other, ribbons = communication probability
# from MIC source to AST target for that pathway. No autocrine, no AST->AST,
# no MIC->MIC. One page per pathway, three panels per page (CV / R47H / R62H).

message("\n\n          *** MIC->AST-only per-pathway chord (custom)... ", Sys.time(), "\n\n")


draw_MIC_AST_chord = function(cc, pw, clusters_MIC, clusters_AST,
                              ribbon_threshold = 0, padj_cutoff_use){
  # Aggregate L-R level prob to pathway level for this single pathway, masking
  # non-significant L-R edges via @net$padj (BH-FDR within MIC->AST).
  # Title is added by the caller, not here, to avoid double-title overwrites.
  prob = cc@net$prob
  padj = cc@net$padj
  if (is.null(prob) || is.null(padj)){
    plot.new()
    return(invisible(NULL))
  }

  # interaction_name -> pathway_name
  lr_sig = cc@LR$LRsig
  lr_in_pw = lr_sig$interaction_name[lr_sig$pathway_name == pw]
  lr_in_pw = intersect(lr_in_pw, dimnames(prob)[[3]])
  if (length(lr_in_pw) == 0){
    plot.new()
    return(invisible(NULL))
  }

  src = intersect(clusters_MIC, dimnames(prob)[[1]])
  tgt = intersect(clusters_AST, dimnames(prob)[[2]])
  if (length(src) == 0 || length(tgt) == 0){
    plot.new()
    return(invisible(NULL))
  }

  # mask non-significant L-Rs to zero, then sum across L-Rs in this pathway
  prob_sub = prob[src, tgt, lr_in_pw, drop = FALSE]
  padj_sub = padj[src, tgt, lr_in_pw, drop = FALSE]
  mask = is.na(padj_sub) | padj_sub >= padj_cutoff_use
  prob_sub[mask] = 0

  m = apply(prob_sub, c(1, 2), sum)   # [src, tgt]
  m[m <= ribbon_threshold] = 0
  if (sum(m) == 0){
    plot.new()
    return(invisible(NULL))
  }

  # convert to long form for circlize::chordDiagram
  df = expand.grid(from = rownames(m), to = colnames(m), stringsAsFactors = FALSE)
  df$value = as.vector(m)
  df = df[df$value > 0, ]
  if (nrow(df) == 0){
    plot.new()
    mtext("(no MIC->AST signal)", side = 3, line = -2, cex = 0.8)
    return(invisible(NULL))
  }

  # colour scheme: MIC sectors orange-ish, AST sectors blue-ish
  grid.col = setNames(
    c(rep("orange", length(src)), rep("dodgerblue", length(tgt))),
    c(src, tgt)
  )

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
  pdf(file = paste0(out_dir, script_ind, "per_pathway_chord_MIC_to_AST_only_padj.pdf"),
      width = 9, height = 9)
  {
    for (pw in pathways_MIC_AST_union){
      for (v in variant_levels){
        cc = cc_list[[v]]
        tryCatch({
          draw_MIC_AST_chord(cc, pw, clusters_MIC, clusters_AST,
                             padj_cutoff_use = padj_cutoff)
        }, error = function(e){
          plot.new()
          message("    chord error for ", v, " ", pw, ": ", conditionMessage(e))
        })
        title(main = paste0(v, " - ", pw), line = -1)
        mtext(paste0("Significance: BH-FDR padj < ", padj_cutoff, " (MIC->AST L-R scope)"),
              side = 1, line = 3, cex = 0.7)
      }
    }
  }
  dev.off()
}



###########################################################
# 11. MIC->AST L-R rankings tables
###########################################################

# Per-variant ranked tables + pathway-level summary.

message("\n\n          *** MIC->AST L-R rankings... ", Sys.time(), "\n\n")


for (v in variant_levels){

  net_v = net_tabs_MIC_AST[[v]]

  if (nrow(net_v) > 0){
    # ranked CSV: sort by padj ascending (NAs last), then prob descending
    net_v_ranked = net_v[order(net_v$padj, -net_v$prob, na.last = TRUE), ]
    write_csv(net_v_ranked,
              file = paste0(out_dir, script_ind, "interactions_MIC_to_AST_ranked_", v, ".csv"))

    path_tab = net_v %>%
      group_by(source, target, pathway_name) %>%
      summarise(sum_prob = sum(prob), n_LR = n(), .groups = "drop") %>%
      arrange(desc(sum_prob))

    write_csv(path_tab,
              file = paste0(out_dir, script_ind, "pathway_by_cluster_pair_MIC_to_AST_", v, ".csv"))

    cat(v, " - top 5 MIC->AST source-target-pathway combinations:\n", sep = "")
    print(head(path_tab, 5))
  } else {
    message("    No MIC->AST L-Rs for ", v)
  }
}



###########################################################
# 12. MIC->AST L-R bubble plots
###########################################################

# Wider PDF, smaller text. Per variant separately because trying to combine
# three netVisual_bubble plots side-by-side at this density is unreadable.

message("\n\n          *** MIC->AST bubble plots... ", Sys.time(), "\n\n")


pdf(file = paste0(out_dir, script_ind, "bubble_MIC_to_AST_pval.pdf"),
    width = 24, height = 18)
{
  for (v in variant_levels){
    if (length(pathways_MIC_AST[[v]]) == 0){
      message("    skipping bubble for ", v, " (no MIC->AST pathways)")
      next
    }
    tryCatch({
      p = netVisual_bubble(cc_list[[v]],
                           sources.use    = intersect(clusters_MIC, levels(cc_list[[v]]@idents)),
                           targets.use    = intersect(clusters_AST, levels(cc_list[[v]]@idents)),
                           signaling      = pathways_MIC_AST[[v]],
                           remove.isolate = TRUE,
                           angle.x        = 90,
                           font.size      = 6,
                           font.size.title = 9) +
        ggtitle(paste0("MIC -> AST - ", v))
      print(p)
    }, error = function(e){
      message("    bubble skipped for ", v, ": ", conditionMessage(e))
    })
  }
}
dev.off()


### custom ggplot bubble per variant
pdf(file = paste0(out_dir, script_ind, "bubble_MIC_to_AST_custom_pval.pdf"),
    width = 26, height = 18)
{
  for (v in variant_levels){
    net_v = net_tabs_MIC_AST[[v]]
    if (nrow(net_v) == 0) next

    t1 = net_v[order(-net_v$prob), ]
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
      ggtitle(paste0("MIC -> AST - ", v))
    print(p_custom)
  }
}
dev.off()



###########################################################
# 13. per-pathway gene expression violins (FIXED for cairo_pdf)
###########################################################

# G03a's section 13 PDF was corrupted - switching to cairo_pdf and isolating
# each page in its own tryCatch so device state stays clean.

message("\n\n          *** Per-pathway gene expression violins... ", Sys.time(), "\n\n")


cairo_pdf(filename = paste0(out_dir, script_ind, "per_pathway_gene_expression_MIC_AST.pdf"),
          width = 14, height = 6, onefile = TRUE)
{
  for (pw in pathways_MIC_AST_union){
    for (v in variant_levels){
      cc = cc_list[[v]]
      if (!(pw %in% cc@netP$pathways)){
        grid::grid.newpage()
        grid::grid.text(paste0(v, " - ", pw, "\n(absent in this variant)"),
                        gp = grid::gpar(fontsize = 14))
        next
      }
      tryCatch({
        p = plotGeneExpression(cc, signaling = pw)
        print(p + plot_annotation(title = paste0(v, " - ", pw)))
      }, error = function(e){
        grid::grid.newpage()
        grid::grid.text(paste0(v, " - ", pw, "\ngene expression error:\n", conditionMessage(e)),
                        gp = grid::gpar(fontsize = 12))
        message("    gene expression skipped for ", v, " ", pw, ": ", conditionMessage(e))
      })
    }
  }
}
dev.off()



###########################################################
# 14. pathway-gene Z-score heatmaps (uses G01a merged pseudobulk)
###########################################################

# As G03a section 14: rows = pathway genes, cols = cluster_samples ordered by
# lineage > cluster > sample. Annotated by lineage, TREM2Variant, diagnosis,
# cluster. The Z-score matrix is NOT variant-specific - it's the joint Z-score
# across all cluster_samples - so this section produces ONE PDF, not one per
# variant. Pathways looped over the union set.

message("\n\n          *** Pathway-gene Z-score heatmaps... ", Sys.time(), "\n\n")


if (length(pathways_MIC_AST_union) > 0){

  # build pathway->genes table from union of variants' DBs
  # (DBs are identical across variants since DEG-filter step was the same)
  pathway_genes_tab = tibble(pathway = character(),
                             interaction_name = character(),
                             gene_type = character(),
                             gene = character())

  for (pw in pathways_MIC_AST_union){
    # use first variant where this pathway exists - all variants share the DB
    v_use = variant_levels[sapply(variant_levels, function(v) pw %in% cc_list[[v]]@netP$pathways)][1]
    if (is.na(v_use)) next

    cc = cc_list[[v_use]]
    net_v = net_tabs_MIC_AST[[v_use]]
    path_ints = unique(net_v$interaction_name[net_v$pathway_name == pw])

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
# 15. save updated CellChat objects (with centrality)
###########################################################

# For re-runs of this script's plotting sections without recomputing centrality.
# G03c does NOT consume these - it loads from G02c independently.

message("\n\n          *** Save updated CellChat objects... ", Sys.time(), "\n\n")


for (v in variant_levels){
  qsave(cc_list[[v]],
        file = paste0(out_dir, script_ind, "cellchat_TREM2Variant_", v, "_with_centrality.qs"))
}


message("\n\n##########################################################################\n",
        "# Finished G03b ", Sys.time(),
        "\n##########################################################################\n\n")
