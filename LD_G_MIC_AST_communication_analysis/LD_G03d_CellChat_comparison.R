message("\n\n##########################################################################\n",
        "# Start LD_G03d: Pairwise CellChat comparisons in AD ", Sys.time(),
        "\n##########################################################################\n\n")

library(qs)
library(tidyverse)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)


### define directories and script index
main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

script_ind = "LD_G03d_v001_"

in_dir  = paste0(main_dir, "LD_G_MIC_AST_communication_analysis_output/")
out_dir = paste0(in_dir, "LD_G03d/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}


#input CellChat objects (G02c output - AD groups only)
ad_groups = c("CV_AD", "R47H_AD", "R62H_AD")
in_cellchat = setNames(
  paste0(in_dir, "LD_G02c_v001_cellchat_TREM2Variant_x_Diagnosis_", ad_groups, ".qs"),
  ad_groups
)


### pairwise comparisons (always (variant of interest) vs CV as reference,
### R47H vs R62H added as exploratory)
# R47H_AD has 8 donors, so comparisons involving it are descriptive.
pairs = list(
  CV_vs_R47H   = c("CV_AD",   "R47H_AD"),
  CV_vs_R62H   = c("CV_AD",   "R62H_AD"),
  R47H_vs_R62H = c("R47H_AD", "R62H_AD")
)



###########################################################
# 1. load CellChat objects and run centrality (if not already)
###########################################################

message("\n\n          *** Load CellChat objects... ", Sys.time(), "\n\n")


cc_list = lapply(in_cellchat, qread)


### ensure centrality is computed in each object
for (g in ad_groups){
  if (length(cc_list[[g]]@netP$centr) == 0){
    cc_list[[g]] = netAnalysis_computeCentrality(cc_list[[g]], slot.name = "netP")
  }
  cat(g, ": cells=", length(cc_list[[g]]@idents),
      "  clusters=", length(levels(cc_list[[g]]@idents)),
      "  pathways=", length(cc_list[[g]]@netP$pathways), "\n", sep = "")
}



###########################################################
# 2. cluster ordering across the 3 AD groups (union)
###########################################################

# CellChat::mergeCellChat requires the same idents levels across objects
# being merged. If small variants dropped clusters, we have to lift the
# smaller objects to the union cluster set via liftCellChat() before merging.

message("\n\n          *** Define union cluster set + lift if needed... ", Sys.time(), "\n\n")


all_clusters_union = sort(unique(unlist(lapply(cc_list, function(x) levels(x@idents)))))

ast_order = c(
  sort(grep("^AST_SLC1A2",  all_clusters_union, value = TRUE)),
  sort(grep("^AST_GFAP",    all_clusters_union, value = TRUE)),
  sort(grep("^AST_CHI3L1",  all_clusters_union, value = TRUE))
)
mic_order = c(
  sort(grep("^HOM_",        all_clusters_union, value = TRUE)),
  sort(grep("^DAM_",        all_clusters_union, value = TRUE)),
  sort(grep("^IRM_",        all_clusters_union, value = TRUE)),
  sort(grep("^HLA_",        all_clusters_union, value = TRUE)),
  sort(grep("^CRM_",        all_clusters_union, value = TRUE))
)
cluster_order = c(ast_order, mic_order)

unmatched = setdiff(all_clusters_union, cluster_order)
if (length(unmatched) > 0){
  warning("Clusters not matched by AST/MIC prefix - appended at end: ",
          paste(unmatched, collapse = ", "))
  cluster_order = c(cluster_order, sort(unmatched))
}

clusters_AST = ast_order
clusters_MIC = mic_order


### lift each object to the union set (no-op if already complete)
### liftCellChat expands @idents and pads internal slots with zeros for
### missing clusters; required for mergeCellChat to align.
for (g in ad_groups){
  present = levels(cc_list[[g]]@idents)
  missing = setdiff(cluster_order, present)
  if (length(missing) > 0){
    cat("  ", g, ": lifting to add ", length(missing), " missing clusters: ",
        paste(missing, collapse = ", "), "\n", sep = "")
    cc_list[[g]] = tryCatch(
      liftCellChat(cc_list[[g]], cluster_order),
      error = function(e){
        message("    liftCellChat failed for ", g, ": ", conditionMessage(e),
                " - leaving object un-lifted; pairs requiring this group ",
                "may fail at mergeCellChat and will be skipped.")
        cc_list[[g]]
      })
  } else {
    cat("  ", g, ": already has full cluster set\n", sep = "")
  }
  ### apply consistent ordering on whatever levels are present
  this_order = intersect(cluster_order, levels(cc_list[[g]]@idents))
  cc_list[[g]]@idents = factor(cc_list[[g]]@idents, levels = this_order)
}



###########################################################
# 3. helpers: rankNet wrapper + per-pair analysis
###########################################################

### plot_rankNet_block: stacked + side-by-side(raw) PDF + raw contribution CSV.
###   sources.use / targets.use = NULL  -> unfiltered (all clusters, default rankNet).
###   sources.use / targets.use = char  -> directional slice (e.g. MIC -> AST).
### File names: <script_ind><pair_label>_rankNet_<label>_pathway.pdf / _data.csv
plot_rankNet_block = function(merged_cc, comparison, label, pair_label,
                              out_dir, script_ind,
                              sources.use = NULL, targets.use = NULL,
                              pdf_width = 9, pdf_height = 14){

  cat("    -> ", pair_label, " | rankNet ", label,
      " (sources=", ifelse(is.null(sources.use), "all", length(sources.use)),
      ", targets=", ifelse(is.null(targets.use), "all", length(targets.use)),
      ")\n", sep = "")

  pdf_path = paste0(out_dir, script_ind, pair_label, "_rankNet_", label, "_pathway.pdf")
  csv_path = paste0(out_dir, script_ind, pair_label, "_rankNet_", label, "_data.csv")

  pdf(file = pdf_path, width = pdf_width, height = pdf_height)
  {
    tryCatch({
      gg1 = rankNet(merged_cc, mode = "comparison", measure = "weight",
                    comparison = comparison,
                    sources.use = sources.use, targets.use = targets.use,
                    stacked = TRUE, do.stat = TRUE) +
            ggtitle(paste0(pair_label, " - ", label, " info flow (stacked)"))
      print(gg1)
    }, error = function(e){ message("       stacked failed: ", conditionMessage(e)) })

    tryCatch({
      gg2 = rankNet(merged_cc, mode = "comparison", measure = "weight",
                    comparison = comparison,
                    sources.use = sources.use, targets.use = targets.use,
                    stacked = FALSE, do.stat = TRUE,
                    show.raw = TRUE) +
            ggtitle(paste0(pair_label, " - ", label, " info flow (side-by-side, raw)"))
      print(gg2)
    }, error = function(e){ message("       side-by-side failed: ", conditionMessage(e)) })
  }
  dev.off()

  tryCatch({
    rn_data = rankNet(merged_cc, mode = "comparison", measure = "weight",
                      comparison = comparison,
                      sources.use = sources.use, targets.use = targets.use,
                      stacked = FALSE, do.stat = TRUE,
                      show.raw = TRUE, return.data = TRUE)
    rn_df = if (is.data.frame(rn_data)) rn_data else rn_data$signaling.contribution
    write.csv(rn_df, file = csv_path, row.names = FALSE)
  }, error = function(e){ message("       data CSV failed: ", conditionMessage(e)) })

  invisible(NULL)
}


### plot_directional_compareInteractions: barplot of total count + total weight
###   per dataset, restricted to sources.use -> targets.use slice. CellChat's
###   compareInteractions() has no source/target filter, so we compute manually
###   from @net[[i]]$count and @net[[i]]$weight ([source, target] matrices).
plot_directional_compareInteractions = function(merged_cc, sources.use, targets.use,
                                                group_names, label, pair_label,
                                                out_dir, script_ind){

  cat("    -> ", pair_label, " | compareInteractions ", label, "\n", sep = "")

  pdf_path = paste0(out_dir, script_ind, pair_label, "_compareInteractions_", label, ".pdf")
  csv_path = paste0(out_dir, script_ind, pair_label, "_compareInteractions_", label, ".csv")

  df = tryCatch({
    src_use = intersect(sources.use, rownames(merged_cc@net[[1]]$count))
    tgt_use = intersect(targets.use, colnames(merged_cc@net[[1]]$count))

    if (length(src_use) == 0 | length(tgt_use) == 0){
      message("       no matching sources/targets - skipping ", label)
      return(invisible(NULL))
    }

    do.call(rbind, lapply(seq_along(merged_cc@net), function(i){
      data.frame(
        dataset = group_names[i],
        count   = sum(merged_cc@net[[i]]$count [src_use, tgt_use]),
        weight  = sum(merged_cc@net[[i]]$weight[src_use, tgt_use])
      )
    }))
  }, error = function(e){
    message("       directional compareInteractions data failed: ", conditionMessage(e)); NULL
  })

  if (is.null(df)) return(invisible(NULL))
  df$dataset = factor(df$dataset, levels = group_names)

  pdf(file = pdf_path, width = 8, height = 5)
  tryCatch({
    g1 = ggplot(df, aes(x = dataset, y = count, fill = dataset)) +
           geom_col() + theme_classic() +
           theme(legend.position = "none",
                 axis.text.x = element_text(angle = 45, hjust = 1)) +
           labs(x = NULL, y = "# interactions",
                title = paste0(pair_label, " - ", label, " # interactions"))
    g2 = ggplot(df, aes(x = dataset, y = weight, fill = dataset)) +
           geom_col() + theme_classic() +
           theme(legend.position = "none",
                 axis.text.x = element_text(angle = 45, hjust = 1)) +
           labs(x = NULL, y = "interaction weight",
                title = paste0(pair_label, " - ", label, " interaction weight"))
    print(g1 + g2)
  }, error = function(e){ message("       directional compareInteractions plot failed: ", conditionMessage(e)) })
  dev.off()

  tryCatch(write.csv(df, file = csv_path, row.names = FALSE),
           error = function(e){ message("       CSV failed: ", conditionMessage(e)) })

  invisible(NULL)
}


### plot_3group_rankNet_stats:
###   For >2 datasets CellChat's do.stat does not annotate significance in the
###   plot. This helper recreates that idea: per-pathway Kruskal-Wallis test
###   on per-cluster-pair contributions (sum over L-Rs in pathway, restricted
###   to sources.use x targets.use slice), BH-corrected across pathways, plus
###   pairwise post-hoc Wilcoxon (BH per pair). Outputs:
###     <pair_label>_rankNet_<label>_3group_stats.csv
###     <pair_label>_rankNet_<label>_3group_annotated.pdf  - bar plot with
###       y-axis labels coloured by dominant group when KW q < q_thresh.
plot_3group_rankNet_stats = function(merged_cc, group_names, label, pair_label,
                                     out_dir, script_ind,
                                     sources.use = NULL, targets.use = NULL,
                                     q_thresh = 0.05){

  cat("    -> ", pair_label, " | 3-group stats ", label, "\n", sep = "")

  csv_path = paste0(out_dir, script_ind, pair_label, "_rankNet_", label, "_3group_stats.csv")
  pdf_path = paste0(out_dir, script_ind, pair_label, "_rankNet_", label, "_3group_annotated.pdf")

  ### Use pathway-level prob (@netP$prob) - same source as CellChat's rankNet,
  ### so this helper's pathway set matches the rankNet PDF exactly.
  ### @netP[[i]]$prob is [source, target, pathway], already filtered by
  ### pathway-level permutation testing in computeCommunProbPathway().
  pathways = sort(unique(unlist(lapply(merged_cc@netP, function(x){
    if (!is.null(dimnames(x$prob))) dimnames(x$prob)[[3]] else x$pathways
  }))))
  if (length(pathways) == 0){
    message("       no pathways in @netP - skipping 3-group stats")
    return(invisible(NULL))
  }

  per_pathway_vectors = function(P){
    lapply(seq_along(merged_cc@netP), function(i){
      prob_arr = merged_cc@netP[[i]]$prob   # [source, target, pathway]
      pathways_avail = dimnames(prob_arr)[[3]]
      if (!(P %in% pathways_avail)) return(numeric(0))
      src_use = if (is.null(sources.use)) rownames(prob_arr) else intersect(sources.use, rownames(prob_arr))
      tgt_use = if (is.null(targets.use)) colnames(prob_arr) else intersect(targets.use, colnames(prob_arr))
      if (length(src_use) == 0 || length(tgt_use) == 0) return(numeric(0))
      mat = prob_arr[src_use, tgt_use, P, drop = FALSE]
      as.vector(mat)
    })
  }

  rows = lapply(pathways, function(P){
    vecs = per_pathway_vectors(P)
    names(vecs) = group_names

    per_group_sum    = sapply(vecs, function(v) if (length(v) == 0) 0 else sum(v))
    per_group_median = sapply(vecs, function(v) if (length(v) == 0) NA_real_ else median(v))
    n_pairs          = sapply(vecs, length)

    long_df = do.call(rbind, lapply(seq_along(vecs), function(i){
      if (length(vecs[[i]]) == 0) return(NULL)
      data.frame(value = vecs[[i]], group = names(vecs)[i])
    }))
    kw_p = NA_real_
    if (!is.null(long_df) && length(unique(long_df$group)) >= 2 && var(long_df$value) > 0){
      kw_p = tryCatch(kruskal.test(value ~ group, data = long_df)$p.value,
                      error = function(e) NA_real_)
    }

    grp_pairs = combn(group_names, 2, simplify = FALSE)
    pw = sapply(grp_pairs, function(gp){
      a = vecs[[gp[1]]]; b = vecs[[gp[2]]]
      if (length(a) == 0 || length(b) == 0) return(NA_real_)
      if (all(c(a, b) == 0)) return(NA_real_)
      tryCatch(wilcox.test(a, b)$p.value, error = function(e) NA_real_)
    })
    names(pw) = sapply(grp_pairs, function(gp) paste0("p_", gp[1], "_vs_", gp[2]))

    out = data.frame(pathway = P,
                     n_pairs_per_group = paste(n_pairs, collapse = ";"),
                     check.names = FALSE)
    for (g in group_names){ out[[paste0("sum_",    g)]] = per_group_sum[g] }
    for (g in group_names){ out[[paste0("median_", g)]] = per_group_median[g] }
    out$kw_p = kw_p
    for (nm in names(pw)) out[[nm]] = unname(pw[nm])
    out
  })

  stats_df = do.call(rbind, rows)
  stats_df$kw_q = p.adjust(stats_df$kw_p, method = "BH")
  pw_cols = grep("^p_", colnames(stats_df), value = TRUE)
  for (pc in pw_cols){
    stats_df[[sub("^p_", "q_", pc)]] = p.adjust(stats_df[[pc]], method = "BH")
  }

  ### dominant_group = group with the highest summed contribution (matches bar height)
  sum_cols = paste0("sum_", group_names)
  stats_df$dominant_group = apply(stats_df[, sum_cols, drop = FALSE], 1, function(r){
    if (all(is.na(r)) || all(r == 0, na.rm = TRUE)) return(NA_character_)
    group_names[which.max(replace(r, is.na(r), -Inf))]
  })

  star = function(q){
    if (is.na(q)) return("")
    if (q < 0.001) return("***")
    if (q < 0.01)  return("**")
    if (q < 0.05)  return("*")
    ""
  }
  stats_df$kw_sig_label = sapply(stats_df$kw_q, star)
  stats_df$total_sum = rowSums(stats_df[, paste0("sum_", group_names), drop = FALSE], na.rm = TRUE)
  stats_df = stats_df[order(stats_df$kw_q, -stats_df$total_sum), ]

  write.csv(stats_df, file = csv_path, row.names = FALSE)

  plot_df = stats_df[stats_df$total_sum > 0, ]
  if (nrow(plot_df) == 0){
    message("       no non-zero pathways - skipping plot")
    return(invisible(stats_df))
  }

  long = do.call(rbind, lapply(group_names, function(g){
    data.frame(pathway = plot_df$pathway,
               group   = g,
               value   = plot_df[[paste0("sum_", g)]])
  }))
  long$pathway = factor(long$pathway, levels = rev(plot_df$pathway))
  long$group   = factor(long$group, levels = group_names)

  group_cols = scales::hue_pal()(length(group_names))
  names(group_cols) = group_names

  label_cols = sapply(levels(long$pathway), function(P){
    row = plot_df[plot_df$pathway == P, ]
    if (nrow(row) == 0) return("black")
    if (is.na(row$kw_q) || row$kw_q >= q_thresh) return("black")
    if (is.na(row$dominant_group)) return("black")
    group_cols[[row$dominant_group]]
  })
  pathway_labels = sapply(levels(long$pathway), function(P){
    row = plot_df[plot_df$pathway == P, ]
    if (nrow(row) == 0) return(P)
    s = row$kw_sig_label
    if (s == "") P else paste0(P, " ", s)
  })

  p = ggplot(long, aes(x = value, y = pathway, fill = group)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7) +
        scale_fill_manual(values = group_cols) +
        scale_y_discrete(labels = pathway_labels) +
        theme_classic() +
        theme(axis.text.y = element_text(colour = label_cols, size = 7),
              legend.position = "right") +
        labs(title = paste0(pair_label, " - ", label, " info flow (3-group, KW BH-q)"),
             subtitle = paste0("y-label colour = dominant group when KW q < ", q_thresh,
                               "; * q<0.05, ** q<0.01, *** q<0.001"),
             x = "info flow (sum of probabilities)", y = NULL,
             fill = "dataset")

  pdf_height = max(6, min(40, 0.18 * nrow(plot_df) + 2))
  pdf(file = pdf_path, width = 9, height = pdf_height)
  tryCatch(print(p), error = function(e){ message("       3-group plot failed: ", conditionMessage(e)) })
  dev.off()

  invisible(stats_df)
}


### plot_by_MIC_bubble_3group:
###   Per-MIC-subtype bubble plot across all 3 AD groups in one PDF.
###   One page per MIC source subtype. Each page = one panel with three
###   datasets shown side-by-side per (source, target) cluster pair.
###   No max.dataset filter - shows every L-R with non-zero signal in any group.
plot_by_MIC_bubble_3group = function(merged_cc, group_names,
                                     clusters_MIC, clusters_AST,
                                     out_dir, script_ind, pair_label = "all_AD"){

  cat("    -> ", pair_label, " by-MIC bubble (3 groups, ", length(group_names), " datasets)\n", sep = "")
  pdf_path = paste0(out_dir, script_ind, pair_label, "_bubble_MIC_to_AST_by_MIC.pdf")

  ### Restrict to MIC->AST pathways present in any dataset (avoids dumping
  ### the entire L-R database into the bubble).
  pw_signaling = tryCatch({
    pw_lists = lapply(seq_along(merged_cc@net), function(i){
      ### derive pathways from each netP entry
      merged_cc@netP[[i]]$pathways
    })
    sort(unique(unlist(pw_lists)))
  }, error = function(e){ NULL })

  src_use = intersect(clusters_MIC, rownames(merged_cc@net[[1]]$prob))
  tgt_use = intersect(clusters_AST, colnames(merged_cc@net[[1]]$prob))

  if (length(src_use) == 0 || length(tgt_use) == 0){
    message("       no MIC sources or AST targets found - skipping")
    return(invisible(NULL))
  }
  if (is.null(pw_signaling) || length(pw_signaling) == 0){
    message("       no signalling pathways found in @netP - skipping (would dump full database)")
    return(invisible(NULL))
  }

  pdf(pdf_path, width = 4 + 1.2 * length(tgt_use) * length(group_names) / 3,
                height = 14)

  for (src in src_use){
    tryCatch({
      p = netVisual_bubble(merged_cc,
                           sources.use     = src,
                           targets.use     = tgt_use,
                           signaling       = pw_signaling,
                           comparison      = seq_along(group_names),
                           title.name      = paste0(pair_label, " | source=", src,
                                                    ": MIC->AST L-Rs across groups"),
                           angle.x         = 90,
                           remove.isolate  = TRUE,
                           font.size       = 7,
                           font.size.title = 9)
      print(p)
    }, error = function(e){ message("       bubble (", src, ") failed: ", conditionMessage(e)) })
  }

  dev.off()
  invisible(NULL)
}


### plot_aggregated_bubble_MIC_to_AST:
###   Aggregated MIC->AST L-R bubble plot across N groups (one column per group).
###   Each row = one L-R; bubble value = max(prob) over all (MIC subtype x AST
###   subtype) edges in that group. Same logic as the per-pair aggregated bubble,
###   generalised to >2 groups.
plot_aggregated_bubble_MIC_to_AST = function(cc_list, group_names,
                                             clusters_MIC, clusters_AST,
                                             out_dir, script_ind,
                                             pair_label = "all_AD"){

  cat("    -> ", pair_label, " aggregated bubble MIC->AST (", length(group_names), " groups)\n", sep = "")
  pdf_path = paste0(out_dir, script_ind, pair_label, "_bubble_MIC_to_AST_aggr.pdf")
  csv_path = paste0(out_dir, script_ind, pair_label, "_bubble_MIC_to_AST_aggr.csv")

  agg_long = tryCatch({
    do.call(rbind, lapply(group_names, function(g){
      cc = cc_list[[g]]
      nt = subsetCommunication(cc)
      nt = nt[nt$source %in% clusters_MIC & nt$target %in% clusters_AST, ]
      if (nrow(nt) == 0) return(NULL)
      lr_label = if ("interaction_name_2" %in% colnames(nt)) "interaction_name_2" else "interaction_name"
      a = aggregate(nt[, "prob"],
                    by = list(LR = nt[[lr_label]], pathway_name = nt$pathway_name),
                    FUN = function(x) max(x, na.rm = TRUE))
      colnames(a)[3] = "prob"
      a$dataset = g
      a
    }))
  }, error = function(e){ message("       aggregated bubble data failed: ", conditionMessage(e)); NULL })

  if (is.null(agg_long) || nrow(agg_long) == 0){
    message("       no MIC->AST L-Rs across groups - skipping")
    return(invisible(NULL))
  }
  agg_long$dataset = factor(agg_long$dataset, levels = group_names)

  ### order rows by max prob across groups, descending (top of plot = strongest L-R anywhere)
  lr_order = aggregate(agg_long$prob, by = list(LR = agg_long$LR), FUN = max)
  lr_order = lr_order[order(-lr_order$x), "LR"]
  agg_long$LR = factor(agg_long$LR, levels = rev(lr_order))

  write.csv(agg_long, file = csv_path, row.names = FALSE)

  ### dynamic height so long lists stay readable
  pdf_height = max(8, min(40, 0.18 * length(unique(agg_long$LR)) + 2))
  pdf(pdf_path, width = 2 + 1.5 * length(group_names), height = pdf_height)
  tryCatch({
    p = ggplot(agg_long, aes(x = dataset, y = LR, size = prob, colour = prob)) +
          geom_point() +
          scale_colour_viridis_c(option = "magma", end = 0.85) +
          scale_size_continuous(range = c(1, 5)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                axis.text.y = element_text(size = 6)) +
          labs(title = paste0(pair_label, ": MIC->AST L-Rs (aggregated, max prob)"),
               x = NULL, y = NULL,
               size = "max prob", colour = "max prob")
    print(p)
  }, error = function(e){ message("       aggregated bubble plot failed: ", conditionMessage(e)) })
  dev.off()

  invisible(agg_long)
}


### plot_per_group_MIC_AST_circles:
###   Per-group circle plots, all 3 AD groups side-by-side, restricted to
###   MIC<->AST edges. One PDF per measure (count or weight), two pages:
###     page 1: MIC -> AST direction (3 panels)
###     page 2: AST -> MIC direction (3 panels)
###   Shared edge.weight.max across panels within a page so chord widths
###   are visually comparable between groups.
plot_per_group_MIC_AST_circles = function(cc_list, group_names,
                                          clusters_MIC, clusters_AST,
                                          out_dir, script_ind, measure = "count"){

  cat("    -> all_AD per-group MIC<->AST circle (", measure, ")\n", sep = "")
  pdf_path = paste0(out_dir, script_ind, "all_AD_circle_MIC_AST_", measure, ".pdf")

  directional_mat = function(cc, src_set, tgt_set){
    M = cc@net[[measure]]
    out = matrix(0, nrow = nrow(M), ncol = ncol(M), dimnames = dimnames(M))
    keep_src = intersect(src_set, rownames(M))
    keep_tgt = intersect(tgt_set, colnames(M))
    if (length(keep_src) > 0 && length(keep_tgt) > 0){
      out[keep_src, keep_tgt] = M[keep_src, keep_tgt]
    }
    out
  }

  mats_MtoA = lapply(cc_list, directional_mat, src_set = clusters_MIC, tgt_set = clusters_AST)
  mats_AtoM = lapply(cc_list, directional_mat, src_set = clusters_AST, tgt_set = clusters_MIC)

  pdf(pdf_path, width = 6 * length(group_names), height = 7)

  ### Per-panel auto-scaling so within-panel variation is visible.
  ### edge.width.max bumped to 15 (default 8) for a wider visible width range.
  ### Cross-panel chord widths are NOT directly comparable - read the
  ### compareInteractions bar plot for absolute magnitudes between groups.
  par(mfrow = c(1, length(group_names)), xpd = TRUE, mar = c(2, 2, 4, 2))
  for (g in group_names){
    tryCatch(
      netVisual_circle(mats_MtoA[[g]],
                       weight.scale   = TRUE,
                       edge.width.max = 15,
                       title.name     = paste0(g, " - MIC->AST ", measure,
                                               " (auto-scaled per panel)")),
      error = function(e) message("       ", g, " MIC->AST circle failed: ", conditionMessage(e)))
  }

  par(mfrow = c(1, length(group_names)), xpd = TRUE, mar = c(2, 2, 4, 2))
  for (g in group_names){
    tryCatch(
      netVisual_circle(mats_AtoM[[g]],
                       weight.scale   = TRUE,
                       edge.width.max = 15,
                       title.name     = paste0(g, " - AST->MIC ", measure,
                                               " (auto-scaled per panel)")),
      error = function(e) message("       ", g, " AST->MIC circle failed: ", conditionMessage(e)))
  }

  dev.off()
  invisible(NULL)
}


run_pair_comparison = function(pair_name, group_a, group_b, cc_a, cc_b,
                               clusters_MIC, clusters_AST, out_dir, script_ind){
  cat("\n\n--- Pair: ", pair_name, " (", group_a, " vs ", group_b, ") ---\n", sep = "")

  ### 3a. merge - 'add.names' becomes the dataset names.
  ### If lift failed earlier and idents levels still differ, mergeCellChat
  ### will error - skip this pair rather than abort the whole script.
  obj_list = list(cc_a, cc_b)
  names(obj_list) = c(group_a, group_b)
  merged_cc = tryCatch(
    mergeCellChat(obj_list, add.names = names(obj_list)),
    error = function(e){
      message("    mergeCellChat failed for ", pair_name, ": ", conditionMessage(e),
              " - skipping this pair.")
      NULL
    })
  if (is.null(merged_cc)) return(invisible(NULL))


  ### 3b. (compareInteractions moved to all-AD section as a single 3-condition plot)


  ### 3c. netVisual_diffInteraction: circle plot of differential network
  ###     red edge = up in group_b vs group_a, blue = down.
  pdf(file = paste0(out_dir, script_ind, pair_name, "_diffInteraction_circle.pdf"),
      width = 16, height = 8)
  {
    par(mfrow = c(1, 2), xpd = TRUE)
    netVisual_diffInteraction(merged_cc, weight.scale = TRUE, measure = "count",
                              title.name = paste0(pair_name, " - diff count"))
    netVisual_diffInteraction(merged_cc, weight.scale = TRUE, measure = "weight",
                              title.name = paste0(pair_name, " - diff weight"))
  }
  dev.off()


  ### 3d. netVisual_heatmap: same diff visualisation as heatmap
  pdf(file = paste0(out_dir, script_ind, pair_name, "_diffInteraction_heatmap.pdf"),
      width = 14, height = 10)
  {
    tryCatch({
      ht1 = netVisual_heatmap(merged_cc, measure = "count",
                              title.name = paste0(pair_name, " - diff count"))
      draw(ht1)
    }, error = function(e){ message("    diff heatmap count failed: ", conditionMessage(e)) })

    tryCatch({
      ht2 = netVisual_heatmap(merged_cc, measure = "weight",
                              title.name = paste0(pair_name, " - diff weight"))
      draw(ht2)
    }, error = function(e){ message("    diff heatmap weight failed: ", conditionMessage(e)) })
  }
  dev.off()


  ### 3e. rankNet: pathway-level differential information flow.
  ###     do.stat = TRUE applies a Wilcoxon rank-sum test per pathway.
  ###     Three views per pair:
  ###       (i)   all_pathways  - whole network (MIC<->MIC, AST<->AST, both directions)
  ###       (ii)  MIC_to_AST    - directional slice: MIC sources -> AST targets
  ###       (iii) AST_to_MIC    - reverse direction (lower priority, completeness)
  plot_rankNet_block(merged_cc, comparison = c(1, 2),
                     label = "all_pathways", pair_label = pair_name,
                     out_dir = out_dir, script_ind = script_ind)

  plot_rankNet_block(merged_cc, comparison = c(1, 2),
                     label = "MIC_to_AST", pair_label = pair_name,
                     out_dir = out_dir, script_ind = script_ind,
                     sources.use = clusters_MIC, targets.use = clusters_AST)

  plot_rankNet_block(merged_cc, comparison = c(1, 2),
                     label = "AST_to_MIC", pair_label = pair_name,
                     out_dir = out_dir, script_ind = script_ind,
                     sources.use = clusters_AST, targets.use = clusters_MIC)


  ### 3f. bubble plots of differentially up-regulated MIC->AST L-Rs
  ###     Two output PDFs:
  ###       (i)  *_bubble_MIC_to_AST_aggr.pdf  - one column per dataset,
  ###            L-Rs aggregated across MIC/AST subtype combos via max(prob).
  ###       (ii) *_bubble_MIC_to_AST_by_MIC.pdf - one page per MIC subtype,
  ###            columns = AST subtype x dataset.
  src_use = intersect(clusters_MIC, levels(cc_a@idents))
  tgt_use = intersect(clusters_AST, levels(cc_a@idents))

  nt_a_pw = subsetCommunication(cc_a)
  nt_b_pw = subsetCommunication(cc_b)
  nt_a_pw = nt_a_pw[nt_a_pw$source %in% clusters_MIC & nt_a_pw$target %in% clusters_AST, ]
  nt_b_pw = nt_b_pw[nt_b_pw$source %in% clusters_MIC & nt_b_pw$target %in% clusters_AST, ]
  pw_signaling = sort(unique(c(nt_a_pw$pathway_name, nt_b_pw$pathway_name)))

  if (length(pw_signaling) == 0){
    message("    no MIC->AST pathways in either group - skipping diff bubbles")
  } else {

    ### (i) Aggregated bubble: one column per dataset, max(prob) across subtype combos
    pdf(file = paste0(out_dir, script_ind, pair_name, "_bubble_MIC_to_AST_aggr.pdf"),
        width = 6, height = 14)
    tryCatch({
      lr_label = "interaction_name_2"
      if (!(lr_label %in% colnames(nt_a_pw))) lr_label = "interaction_name"

      agg_a = aggregate(nt_a_pw[, "prob"],
                        by = list(LR = nt_a_pw[[lr_label]],
                                  pathway_name = nt_a_pw$pathway_name),
                        FUN = function(x) max(x, na.rm = TRUE))
      agg_b = aggregate(nt_b_pw[, "prob"],
                        by = list(LR = nt_b_pw[[lr_label]],
                                  pathway_name = nt_b_pw$pathway_name),
                        FUN = function(x) max(x, na.rm = TRUE))
      colnames(agg_a)[3] = "prob"; colnames(agg_b)[3] = "prob"
      agg_a$dataset = group_a
      agg_b$dataset = group_b

      agg_long = rbind(agg_a, agg_b)
      agg_long$dataset = factor(agg_long$dataset, levels = c(group_a, group_b))

      ### order rows by max prob across datasets, descending
      lr_order = aggregate(agg_long$prob, by = list(LR = agg_long$LR), FUN = max)
      lr_order = lr_order[order(-lr_order$x), "LR"]
      agg_long$LR = factor(agg_long$LR, levels = rev(lr_order))

      p_aggr = ggplot(agg_long, aes(x = dataset, y = LR,
                                    size = prob, colour = prob)) +
        geom_point() +
        scale_colour_viridis_c(option = "magma", end = 0.85) +
        scale_size_continuous(range = c(1, 5)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_text(size = 6)) +
        labs(title = paste0(pair_name, ": MIC->AST L-Rs (aggregated, max prob)"),
             x = NULL, y = NULL,
             size = "max prob", colour = "max prob")
      print(p_aggr)
    }, error = function(e){ message("    bubble aggregated failed: ", conditionMessage(e)) })
    dev.off()


    ### (ii) Per-MIC-subtype bubble: one page per source cluster
    pdf(file = paste0(out_dir, script_ind, pair_name, "_bubble_MIC_to_AST_by_MIC.pdf"),
        width = 12, height = 14)
    for (src in src_use){
      tryCatch({
        p_up = netVisual_bubble(merged_cc,
                                sources.use     = src,
                                targets.use     = tgt_use,
                                signaling       = pw_signaling,
                                comparison      = c(1, 2),
                                max.dataset     = 2,
                                title.name      = paste0(pair_name, " | source=", src,
                                                         ": L-Rs UP in ", group_b),
                                angle.x         = 90,
                                remove.isolate  = TRUE,
                                font.size       = 7,
                                font.size.title = 9)
        print(p_up)
      }, error = function(e){ message("    bubble UP (", src, ") failed: ", conditionMessage(e)) })

      tryCatch({
        p_down = netVisual_bubble(merged_cc,
                                  sources.use     = src,
                                  targets.use     = tgt_use,
                                  signaling       = pw_signaling,
                                  comparison      = c(1, 2),
                                  max.dataset     = 1,
                                  title.name      = paste0(pair_name, " | source=", src,
                                                           ": L-Rs UP in ", group_a),
                                  angle.x         = 90,
                                  remove.isolate  = TRUE,
                                  font.size       = 7,
                                  font.size.title = 9)
        print(p_down)
      }, error = function(e){ message("    bubble DOWN (", src, ") failed: ", conditionMessage(e)) })
    }
    dev.off()
  }


  ### 3g. extract differential L-R tables for both directions
  tryCatch({
    nt_a = subsetCommunication(cc_a)
    nt_b = subsetCommunication(cc_b)

    nt_a = nt_a[nt_a$source %in% clusters_MIC & nt_a$target %in% clusters_AST, ]
    nt_b = nt_b[nt_b$source %in% clusters_MIC & nt_b$target %in% clusters_AST, ]

    nt_a$edge = paste0(nt_a$source, "|", nt_a$target, "|", nt_a$interaction_name)
    nt_b$edge = paste0(nt_b$source, "|", nt_b$target, "|", nt_b$interaction_name)

    keep_cols = c("source", "target", "ligand", "receptor", "pathway_name",
                  "interaction_name", "interaction_name_2")
    info = rbind(nt_a[, c("edge", keep_cols)], nt_b[, c("edge", keep_cols)])
    info = info[!duplicated(info$edge), ]

    a_lookup = setNames(nt_a$prob, nt_a$edge)
    b_lookup = setNames(nt_b$prob, nt_b$edge)
    pa_lookup = setNames(nt_a$pval, nt_a$edge)
    pb_lookup = setNames(nt_b$pval, nt_b$edge)

    diff_tab = info
    diff_tab$prob_a   = ifelse(diff_tab$edge %in% names(a_lookup),
                                a_lookup[diff_tab$edge], 0)
    diff_tab$prob_b   = ifelse(diff_tab$edge %in% names(b_lookup),
                                b_lookup[diff_tab$edge], 0)
    diff_tab$pval_a   = ifelse(diff_tab$edge %in% names(pa_lookup),
                                pa_lookup[diff_tab$edge], 1)
    diff_tab$pval_b   = ifelse(diff_tab$edge %in% names(pb_lookup),
                                pb_lookup[diff_tab$edge], 1)
    diff_tab$prob_diff = diff_tab$prob_b - diff_tab$prob_a
    diff_tab$abs_diff  = abs(diff_tab$prob_diff)
    diff_tab = diff_tab[order(-diff_tab$abs_diff), ]

    names(diff_tab)[names(diff_tab) == "prob_a"] = paste0("prob_", group_a)
    names(diff_tab)[names(diff_tab) == "prob_b"] = paste0("prob_", group_b)
    names(diff_tab)[names(diff_tab) == "pval_a"] = paste0("pval_", group_a)
    names(diff_tab)[names(diff_tab) == "pval_b"] = paste0("pval_", group_b)

    diff_tab = diff_tab[, !colnames(diff_tab) %in% "edge"]

    write_csv(diff_tab,
              file = paste0(out_dir, script_ind, pair_name,
                            "_diffMIC_to_AST_LRs_ranked_by_abs_diff.csv"))

    cat("  -> ", pair_name, ": diff CSV saved (", nrow(diff_tab), " edges)\n", sep = "")
    cat("     top 5 increased in ", group_b, ":\n", sep = "")
    print(head(diff_tab[diff_tab$prob_diff > 0, c("source", "target", "ligand", "receptor",
                                                   paste0("prob_", group_a), paste0("prob_", group_b),
                                                   "prob_diff")], 5))
    cat("     top 5 increased in ", group_a, ":\n", sep = "")
    decreased = diff_tab[diff_tab$prob_diff < 0, ]
    decreased = decreased[order(decreased$prob_diff), ]
    print(head(decreased[, c("source", "target", "ligand", "receptor",
                              paste0("prob_", group_a), paste0("prob_", group_b),
                              "prob_diff")], 5))
  }, error = function(e){ message("    diff CSV failed: ", conditionMessage(e)) })


  ### 3h. save merged object for re-use
  qsave(merged_cc,
        file = paste0(out_dir, script_ind, pair_name, "_merged.qs"))

  invisible(merged_cc)
}



###########################################################
# 4. run all three pairs
###########################################################

message("\n\n          *** Running pairwise comparisons... ", Sys.time(), "\n\n")


for (pair_name in names(pairs)){
  pair_groups = pairs[[pair_name]]
  group_a = pair_groups[1]
  group_b = pair_groups[2]

  cc_a = cc_list[[group_a]]
  cc_b = cc_list[[group_b]]

  run_pair_comparison(pair_name, group_a, group_b, cc_a, cc_b,
                      clusters_MIC, clusters_AST, out_dir, script_ind)
}



###########################################################
# 5. all-three merged object (for global rankNet across all 3 AD groups)
###########################################################

# Useful for one-shot pathway-level overview comparing all 3 variants in AD.
# rankNet handles >2 datasets (do.stat applied per pair via Kruskal-Wallis or
# pairwise wilcoxon depending on CellChat version).

message("\n\n          *** All-three merged + global rankNet... ", Sys.time(), "\n\n")


all_merged = tryCatch(
  mergeCellChat(cc_list, add.names = ad_groups),
  error = function(e){
    message("    all-AD mergeCellChat failed: ", conditionMessage(e),
            " - skipping all-three section.")
    NULL
  })

if (!is.null(all_merged)){
  ### Single 3-condition compareInteractions barplot (replaces per-pair version)
  pdf(file = paste0(out_dir, script_ind, "all_AD_compareInteractions.pdf"),
      width = 8, height = 5)
  {
    tryCatch({
      g1 = compareInteractions(all_merged, show.legend = FALSE,
                               group = seq_along(ad_groups),
                               measure = "count")  +
           ggtitle("All AD groups - total # interactions")
      g2 = compareInteractions(all_merged, show.legend = FALSE,
                               group = seq_along(ad_groups),
                               measure = "weight") +
           ggtitle("All AD groups - total interaction weight")
      print(g1 + g2)
    }, error = function(e){ message("    all-AD compareInteractions failed: ", conditionMessage(e)) })
  }
  dev.off()

  ### Directional compareInteractions: MIC -> AST and AST -> MIC slices
  plot_directional_compareInteractions(all_merged,
                                       sources.use = clusters_MIC,
                                       targets.use = clusters_AST,
                                       group_names = ad_groups,
                                       label = "MIC_to_AST",
                                       pair_label = "all_AD",
                                       out_dir = out_dir, script_ind = script_ind)

  plot_directional_compareInteractions(all_merged,
                                       sources.use = clusters_AST,
                                       targets.use = clusters_MIC,
                                       group_names = ad_groups,
                                       label = "AST_to_MIC",
                                       pair_label = "all_AD",
                                       out_dir = out_dir, script_ind = script_ind)

  ### Three rankNet views across all 3 AD groups
  plot_rankNet_block(all_merged, comparison = seq_along(ad_groups),
                     label = "all_pathways", pair_label = "all_AD",
                     out_dir = out_dir, script_ind = script_ind)

  plot_rankNet_block(all_merged, comparison = seq_along(ad_groups),
                     label = "MIC_to_AST", pair_label = "all_AD",
                     out_dir = out_dir, script_ind = script_ind,
                     sources.use = clusters_MIC, targets.use = clusters_AST)

  plot_rankNet_block(all_merged, comparison = seq_along(ad_groups),
                     label = "AST_to_MIC", pair_label = "all_AD",
                     out_dir = out_dir, script_ind = script_ind,
                     sources.use = clusters_AST, targets.use = clusters_MIC)

  ### 3-group statistics + annotated plot for each slice (KW + pairwise Wilcoxon, BH-corrected)
  plot_3group_rankNet_stats(all_merged, group_names = ad_groups,
                            label = "all_pathways", pair_label = "all_AD",
                            out_dir = out_dir, script_ind = script_ind)

  plot_3group_rankNet_stats(all_merged, group_names = ad_groups,
                            label = "MIC_to_AST", pair_label = "all_AD",
                            out_dir = out_dir, script_ind = script_ind,
                            sources.use = clusters_MIC, targets.use = clusters_AST)

  plot_3group_rankNet_stats(all_merged, group_names = ad_groups,
                            label = "AST_to_MIC", pair_label = "all_AD",
                            out_dir = out_dir, script_ind = script_ind,
                            sources.use = clusters_AST, targets.use = clusters_MIC)

  ### Per-group MIC<->AST circle plots: 3 panels per page, MIC->AST and AST->MIC pages
  plot_per_group_MIC_AST_circles(cc_list, group_names = ad_groups,
                                 clusters_MIC = clusters_MIC,
                                 clusters_AST = clusters_AST,
                                 out_dir = out_dir, script_ind = script_ind,
                                 measure = "count")

  plot_per_group_MIC_AST_circles(cc_list, group_names = ad_groups,
                                 clusters_MIC = clusters_MIC,
                                 clusters_AST = clusters_AST,
                                 out_dir = out_dir, script_ind = script_ind,
                                 measure = "weight")

  ### Aggregated MIC->AST bubble across all 3 AD groups (one column per group)
  plot_aggregated_bubble_MIC_to_AST(cc_list, group_names = ad_groups,
                                    clusters_MIC = clusters_MIC,
                                    clusters_AST = clusters_AST,
                                    out_dir = out_dir, script_ind = script_ind,
                                    pair_label = "all_AD")

  ### Per-MIC-subtype bubble across all 3 AD groups (one page per MIC subtype)
  plot_by_MIC_bubble_3group(all_merged, group_names = ad_groups,
                            clusters_MIC = clusters_MIC,
                            clusters_AST = clusters_AST,
                            out_dir = out_dir, script_ind = script_ind,
                            pair_label = "all_AD")

  qsave(all_merged,
        file = paste0(out_dir, script_ind, "all_AD_merged.qs"))
}


message("\n\n##########################################################################\n",
        "# Finished G03d ", Sys.time(),
        "\n##########################################################################\n\n")
