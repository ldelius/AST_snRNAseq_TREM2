message("\n\n##########################################################################\n",
        "# Start LD_E03c: Five-covariate pseudobulk GSEA ", Sys.time(),
        "\n##########################################################################\n\n")

# Uses FDR-significant five-covariate E02c results.
library(qs)
library(tidyverse)
library(DESeq2)
library(colorRamps)
library(viridis)
library(pheatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(ggrepel)
library(msigdbr)
library(fgsea)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "LD_E03c_"

#specify output directory
out_dir = paste0(main_dir,"LD_E_DESeq_pseudobulk/")


# load the 5-covariate astrocyte DESeq2 results (E02c)
bulk_data = qread(file = paste0(out_dir, "LD_E02c/LD_E02c_v01_bulk_data.qs"))



### get user-defined gene sets for GSEA

gsea_sets = list()

# Green24 astrocyte subpopulation markers (astrocytes only)
t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/Green24_S2_subpopulation_markers.csv"))
t1 = t1[t1$cell.type == "astrocytes",]

for (cl in unique(t1$state)){
  t2 = t1[t1$state == cl & t1$avg_log2FC>log2(1.2) & t1$p_val_adj<0.05,]
  gsea_sets[[paste0("Green24_",cl)]] = t2$gene
}

# User-defined sets are limited to Green et al. astrocyte-state signatures.

lengths(gsea_sets)



###########################################################
# functions
###########################################################

#custom colour palette for variable values defined in vector v
pal = function(v){
  v2 = length(unique(v))
  if (v2 == 2){
    p2 = c("grey20", "dodgerblue")
  } else if (v2 ==3){
    p2 = c("dodgerblue", "grey20", "orange")
  } else if (v2 ==4){
    p2 = c("dodgerblue", "green4","grey20", "orange")
  } else if (v2<6){
    p2 = matlab.like(6)[1:v2]
  } else {
    p2 = matlab.like(v2)
  }
  return(p2)
}


### function create_path_comp_mat_list: create list of matrices for plotting pathways by subcategory vs sel comparisons

create_path_comp_mat_list = function(gsea_res_tab, comps_sel = unique(gsea_res_tab$comp), 
                                     comp_ref = comps_sel){
  
  path_comp_mat_list = list()
  
  for (cat1 in unique(gsea_res_tab$sub_cat)){
    
    t2 = gsea_res_tab[gsea_res_tab$sub_cat == cat1 & gsea_res_tab$comp %in% comps_sel & 
                        gsea_res_tab$pval<0.05,]
    
    pathways_sel = unique(t2$pathway)
    
    path_comp_mat = matrix(nrow = length(pathways_sel), ncol = length(comps_sel),
                           dimnames = list(pathways_sel, comps_sel))
    
    for (comp1 in comps_sel){
      for (path1 in pathways_sel){
        t3 = t2[t2$comp == comp1 & t2$pathway == path1,]
        if (nrow(t3)>0){
          path_comp_mat[path1, comp1] = t2[["signed_neg_log10_pval"]][t2$comp == comp1 & 
                                                                        t2$pathway == path1]
        }
      }
    }
    
    path_comp_mat[is.na(path_comp_mat)] = 0
    
    #keep only pathways reg in ref_comps
    
    m1 = path_comp_mat != 0
    keep_pathways = apply(as.matrix(m1[,colnames(m1) %in% comp_ref]) ,1, any) 
    names(keep_pathways)[keep_pathways ] 
    m2 = path_comp_mat[keep_pathways,]
    
    if (is.matrix(m2)){
      if (nrow(m2)>1){path_comp_mat_list[[cat1]] = m2}
    }
  }
  
  return(path_comp_mat_list)
}




#### create list of matrices with top N terms by comparison

cat1 = "HALLMARKS"

create_top_path_comp_mat_list = function(gsea_res_tab, comps_sel = unique(gsea_res_tab$comp), 
                                         comp_ref = comps_sel, N_top = 10){
  
  path_comp_mat_list = list()
  
  for (cat1 in unique(gsea_res_tab$sub_cat)){
    
    t2 = gsea_res_tab[gsea_res_tab$sub_cat == cat1 & gsea_res_tab$comp %in% comps_sel & 
                        gsea_res_tab$pval<0.05,]
    
    #extract top pathways by comp
    
    pathways_sel = NULL
    
    for (comp1 in comps_sel){
      
      t3 = t2[t2$comp == comp1,]
      
      if (nrow(t3)>0){
        if (nrow(t3)>N_top){
          t3 = t3[order(t3$pval),]
          t3 = t3[1:N_top,]
        }
        pathways_sel = c(pathways_sel, t3$pathway)
      }
    }
    
    path_comp_mat = matrix(nrow = length(pathways_sel), ncol = length(comps_sel),
                           dimnames = list(pathways_sel, comps_sel))
    
    for (comp1 in comps_sel){
      for (path1 in pathways_sel){
        t3 = t2[t2$comp == comp1 & t2$pathway == path1,]
        if (nrow(t3)>0){
          path_comp_mat[path1, comp1] = t2[["signed_neg_log10_pval"]][t2$comp == comp1 & 
                                                                        t2$pathway == path1]
        }
      }
    }
    
    path_comp_mat[is.na(path_comp_mat)] = 0
    
    #keep only pathways reg in ref_comps
    
    m1 = path_comp_mat != 0
    keep_pathways = apply(as.matrix(m1[,colnames(m1) %in% comp_ref]) ,1, any) 
    names(keep_pathways)[keep_pathways ] 
    m2 = path_comp_mat[keep_pathways,]
    
    if (is.matrix(m2)){
      if (nrow(m2)>1){path_comp_mat_list[[cat1]] = m2}
    }
  }
  
  return(path_comp_mat_list)
}




### plot_path_comp_mat_list: plot list of matrices for plotting pathways by subcategory vs sel comparisons

plot_path_comp_mat_list = function(path_comp_mat_list){
  
  for (cat1 in names(path_comp_mat_list)){
    
    path_comp_mat = path_comp_mat_list[[cat1]]
    
    path_comp_mat[path_comp_mat>5] = 5
    path_comp_mat[path_comp_mat< -5] = -5
    
    lim_c = max(abs(path_comp_mat)) 
    
    p1 = pheatmap::pheatmap(path_comp_mat, 
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            show_rownames = TRUE, show_colnames = TRUE,
                            color = colorRampPalette(c("blue", "white", "red"))(250),
                            breaks = seq(-lim_c, lim_c, length.out = 251),
                            border_color = NA, cellwidth = 10, cellheight = 10, 
                            fontsize = 10, 
                            main = paste0(cat1))
  }
  return(p1)
}




#################################################
# GSEA of DEGs by comp
#################################################

deseq_res_list = bulk_data$deseq_results
#deseq_res_list = bulk_data$deseq_results[1:3]

category_codes = c("H")  # Hallmark only (Green added separately as user-defined sets)

gsea_res_tab = NULL

for (cat1 in category_codes){
  
  message("\n   *run GSEA for category ", cat1, " - ", Sys.time(),"\n")
  
  pathways = msigdbr(species = "Homo sapiens", collection = cat1) %>%
    dplyr::select(gs_name, gene_symbol) %>%
    group_by(gs_name) %>%
    summarize(gene_symbol = list(gene_symbol)) %>%
    deframe()
  
  for (comp1 in names(deseq_res_list)){
    
    message("   *run GSEA for comp ", comp1)
    
    t1 = deseq_res_list[[comp1]]
    
    t1$signed_log10p = -sign(t1$log2FoldChange)*log10(t1$pvalue)
    t2 = t1[!is.na(t1$signed_log10p) & t1$signed_log10p != 0,]
    t2 = t2[order(t2$signed_log10p),]
    ranks = t2$signed_log10p
    names(ranks) = t2$gene
    
    set.seed(42)
    
    t3 <- fgsea(pathways = pathways, 
                stats    = ranks,
                minSize  = 15,
                maxSize  = 500)
    
    t4 = cbind(category_code = cat1, comp = comp1, t3)
    
    gsea_res_tab = rbind(gsea_res_tab, t4)
    
    gc()
    
  }
}

### add GSEA for user-defined gene sets 

message("\n   *run GSEA for user defined sets - ", Sys.time(), "\n")

for (comp1 in names(deseq_res_list)){
  
  message("   *run GSEA for comp ", comp1)
  
  t1 = deseq_res_list[[comp1]]
  
  t1$signed_log10p = -sign(t1$log2FoldChange)*log10(t1$pvalue)
  t2 = t1[!is.na(t1$signed_log10p),]
  t2 = t2[order(t2$signed_log10p),]
  ranks = t2$signed_log10p
  names(ranks) = t2$gene
  
  set.seed(42)
  
  t3 <- fgsea(pathways = gsea_sets,
              stats    = ranks,
              minSize  = 5,
              maxSize  = 500)
  
  t4 = cbind(category_code = "user_def_sets", comp = comp1, t3)
  
  gsea_res_tab = rbind(gsea_res_tab, t4)
  
  gc()
  
}


###define subcategories for plotting 

t1 = gsea_res_tab

t1$signed_neg_log10_padj = -sign(t1$NES)*log10(t1$padj)
t1$signed_neg_log10_pval = -sign(t1$NES)*log10(t1$pval)

t1$sub_cat = t1$category_code
t1$sub_cat[grepl("HALLMARK_", t1$pathway)] = "HALLMARKS"
t1$sub_cat[grepl("REACTOME_", t1$pathway)] = "REACTOME"
t1$sub_cat[grepl("KEGG_", t1$pathway)] = "KEGG"
t1$sub_cat[grepl("WP_", t1$pathway)] = "WP"
t1$sub_cat[grepl("GOBP_", t1$pathway)] = "GOBP"
t1$sub_cat[grepl("GOMF_", t1$pathway)] = "GOMF"
t1$sub_cat[grepl("GOCC_", t1$pathway)] = "GOCC"
t1$sub_cat[grepl("HP_", t1$pathway)] = "HP"

gsea_res_tab = t1

### extract leading edge genes

message("\n   *Extract leading edge genes - ", Sys.time(),"\n")

t1 = gsea_res_tab
t2 = t1$leadingEdge
t3 = lapply(t2, paste0, collapse = "/")
gsea_res_tab$leadingEdge_collapsed = unlist(t3)


### save results

bulk_data$gsea_res_tab = gsea_res_tab

write_csv(gsea_res_tab, file = paste0(out_dir,script_ind, "GSEA_results.csv") )

qsave(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.qs"))



#################################################
# plot heatmap GSEA pathways vs comps by category
#################################################

# bulk_data = qread(file = paste0(out_dir,script_ind, "bulk_data.qs"))

create_path_comp_mat_list = function(gsea_res_tab, comps_sel = unique(gsea_res_tab$comp),
                                     comp_ref = comps_sel, sig_cut = 0.05){

  path_comp_mat_list = list()

  for (cat1 in unique(gsea_res_tab$sub_cat)){

    t2 = gsea_res_tab[gsea_res_tab$sub_cat == cat1 & gsea_res_tab$comp %in% comps_sel &
                          !is.na(gsea_res_tab$padj) & gsea_res_tab$padj < sig_cut,]

    pathways_sel = unique(t2$pathway)
    
    path_comp_mat = matrix(nrow = length(pathways_sel), ncol = length(comps_sel),
                           dimnames = list(pathways_sel, comps_sel))
    
    for (comp1 in comps_sel){
      for (path1 in pathways_sel){
        t3 = t2[t2$comp == comp1 & t2$pathway == path1,]
        if (nrow(t3)>0){
          path_comp_mat[path1, comp1] = t2[["NES"]][t2$comp == comp1 & t2$pathway == path1]
        }
      }
    }
    
    path_comp_mat[is.na(path_comp_mat)] = 0
    
    #keep only pathways reg in ref_comps
    
    m1 = path_comp_mat != 0
    keep_pathways = apply(as.matrix(m1[,colnames(m1) %in% comp_ref]) ,1, any) 
    names(keep_pathways)[keep_pathways ] 
    m2 = path_comp_mat[keep_pathways,]
    
    if (is.matrix(m2)){
      if (nrow(m2)>1){path_comp_mat_list[[cat1]] = m2}
    }
  }
  
  return(path_comp_mat_list)
}




#### create list of matrices with top N terms by comparison

cat1 = "HALLMARKS"

create_top_path_comp_mat_list = function(gsea_res_tab, comps_sel = unique(gsea_res_tab$comp),
                                     comp_ref = comps_sel, N_top = 10, sig_cut = 0.05){

  path_comp_mat_list = list()

  for (cat1 in unique(gsea_res_tab$sub_cat)){

    t2 = gsea_res_tab[gsea_res_tab$sub_cat == cat1 & gsea_res_tab$comp %in% comps_sel &
                        !is.na(gsea_res_tab$padj) & gsea_res_tab$padj < sig_cut,]

    #extract top pathways by comp
    
    pathways_sel = NULL
    
    for (comp1 in comps_sel){
      
      t3 = t2[t2$comp == comp1,]
      
      if (nrow(t3)>0){
        if (nrow(t3)>N_top){
          t3 = t3[order(t3$pval),]
          t3 = t3[1:N_top,]
        }
        pathways_sel = c(pathways_sel, t3$pathway)
      }
    }
    
    path_comp_mat = matrix(nrow = length(pathways_sel), ncol = length(comps_sel),
                           dimnames = list(pathways_sel, comps_sel))
    
    for (comp1 in comps_sel){
      for (path1 in pathways_sel){
        t3 = t2[t2$comp == comp1 & t2$pathway == path1,]
        if (nrow(t3)>0){
          path_comp_mat[path1, comp1] = t2[["NES"]][t2$comp == comp1 & t2$pathway == path1]
        }
      }
    }
    
    path_comp_mat[is.na(path_comp_mat)] = 0
    
    #keep only pathways reg in ref_comps
    
    m1 = path_comp_mat != 0
    keep_pathways = apply(as.matrix(m1[,colnames(m1) %in% comp_ref]) ,1, any) 
    names(keep_pathways)[keep_pathways ] 
    m2 = path_comp_mat[keep_pathways,]
    
    if (is.matrix(m2)){
      if (nrow(m2)>1){path_comp_mat_list[[cat1]] = m2}
    }
  }
  
  return(path_comp_mat_list)
}




### plot_path_comp_mat_list: plot list of matrices for plotting pathways by subcategory vs sel comparisons

plot_path_comp_mat_list = function(path_comp_mat_list, prefix = ""){

  for (cat1 in names(path_comp_mat_list)){

    path_comp_mat = path_comp_mat_list[[cat1]]

    lim_c = max(abs(path_comp_mat))

    p1 = pheatmap::pheatmap(path_comp_mat,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            show_rownames = TRUE, show_colnames = TRUE,
                            color = colorRampPalette(c("blue", "white", "red"))(250),
                            breaks = seq(-lim_c, lim_c, length.out = 251),
                            border_color = NA, cellwidth = 10, cellheight = 10,
                            fontsize = 10,
                            main = paste0(prefix, cat1))
  }
  return(p1)
}



###define subcategories and comparisons for plotting 

gsea_res_tab = bulk_data$gsea_res_tab

# comps are named "<cluster>_<comp_tag>". The "largest7" variant keeps only comps
# whose subcluster is among the 7 largest by nuclei: SLC1A2 s0/s3/s4, GFAP s1/s2,
# CHI3L1 s6/s9 (s0 naturally drops from the cluster-vs-ref column, being the ref).
largest7 = c("AST_SLC1A2_s0", "AST_SLC1A2_s3", "AST_SLC1A2_s4",
             "AST_GFAP_s1", "AST_GFAP_s2", "AST_CHI3L1_s6", "AST_CHI3L1_s9")
all_comps     = unique(gsea_res_tab$comp)
comps_largest = all_comps[vapply(all_comps,
                  function(c1) any(startsWith(c1, paste0(largest7, "_"))), logical(1))]

subcl_variants = list(all_subclusters = all_comps, largest7 = comps_largest)
cutoffs        = c(FDR05 = 0.05, FDR10 = 0.10)   # plotted at BOTH; fgsea computed once


### heatmaps (DECIDED settings): all-significant pathways at FDR < 0.10, for the
### two subcluster variants x two collections = 4 heatmaps, each as PDF and PNG.
### Hallmark rows are ordered by functional family (not clustered); Green clustered.

SIG_CUT = 0.10   # FDR < 0.1, consistent with the DEG threshold used elsewhere

# MSigDB Hallmark process categories (Liberzon et al. 2015) used to order rows.
hallmark_family = c(
  HALLMARK_ALLOGRAFT_REJECTION="immune", HALLMARK_COAGULATION="immune",
  HALLMARK_COMPLEMENT="immune", HALLMARK_INTERFERON_ALPHA_RESPONSE="immune",
  HALLMARK_INTERFERON_GAMMA_RESPONSE="immune", HALLMARK_IL6_JAK_STAT3_SIGNALING="immune",
  HALLMARK_INFLAMMATORY_RESPONSE="immune",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB="signaling", HALLMARK_ANDROGEN_RESPONSE="signaling",
  HALLMARK_ESTROGEN_RESPONSE_EARLY="signaling", HALLMARK_ESTROGEN_RESPONSE_LATE="signaling",
  HALLMARK_IL2_STAT5_SIGNALING="signaling", HALLMARK_KRAS_SIGNALING_UP="signaling",
  HALLMARK_KRAS_SIGNALING_DN="signaling", HALLMARK_MTORC1_SIGNALING="signaling",
  HALLMARK_NOTCH_SIGNALING="signaling", HALLMARK_PI3K_AKT_MTOR_SIGNALING="signaling",
  HALLMARK_HEDGEHOG_SIGNALING="signaling", HALLMARK_TGF_BETA_SIGNALING="signaling",
  HALLMARK_WNT_BETA_CATENIN_SIGNALING="signaling",
  HALLMARK_APOPTOSIS="pathway", HALLMARK_HYPOXIA="pathway",
  HALLMARK_PROTEIN_SECRETION="pathway", HALLMARK_UNFOLDED_PROTEIN_RESPONSE="pathway",
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY="pathway",
  HALLMARK_BILE_ACID_METABOLISM="metabolic", HALLMARK_CHOLESTEROL_HOMEOSTASIS="metabolic",
  HALLMARK_FATTY_ACID_METABOLISM="metabolic", HALLMARK_GLYCOLYSIS="metabolic",
  HALLMARK_HEME_METABOLISM="metabolic", HALLMARK_OXIDATIVE_PHOSPHORYLATION="metabolic",
  HALLMARK_XENOBIOTIC_METABOLISM="metabolic",
  HALLMARK_E2F_TARGETS="proliferation", HALLMARK_G2M_CHECKPOINT="proliferation",
  HALLMARK_MYC_TARGETS_V1="proliferation", HALLMARK_MYC_TARGETS_V2="proliferation",
  HALLMARK_P53_PATHWAY="proliferation", HALLMARK_MITOTIC_SPINDLE="proliferation",
  HALLMARK_DNA_REPAIR="DNA_damage", HALLMARK_UV_RESPONSE_DN="DNA_damage",
  HALLMARK_UV_RESPONSE_UP="DNA_damage",
  HALLMARK_ADIPOGENESIS="development", HALLMARK_ANGIOGENESIS="development",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION="development", HALLMARK_MYOGENESIS="development",
  HALLMARK_SPERMATOGENESIS="development", HALLMARK_PANCREAS_BETA_CELLS="development",
  HALLMARK_APICAL_JUNCTION="cellular_component", HALLMARK_APICAL_SURFACE="cellular_component",
  HALLMARK_PEROXISOME="cellular_component")
fam_levels = c("immune", "signaling", "pathway", "metabolic", "proliferation",
               "DNA_damage", "development", "cellular_component", "other")

# draw one NES heatmap (a single matrix) to PDF and PNG
draw_heatmap = function(mat, collection, fbase){
  if (collection == "Hallmark"){
    fam = hallmark_family[rownames(mat)]; fam[is.na(fam)] = "other"
    ord = order(factor(fam, levels = fam_levels), rownames(mat))   # family-ordered
    mat = mat[ord, , drop = FALSE]; fam = fam[ord]
    gaps   = head(cumsum(rle(fam)$lengths), -1)                    # gap between families
    annot  = data.frame(family = factor(fam, levels = fam_levels)); rownames(annot) = rownames(mat)
    crows  = FALSE
  } else {                                                         # cluster Green terms
    gaps = NULL; annot = NA; crows = TRUE
  }
  lim = max(abs(mat))
  W = ncol(mat) * 0.18 + 7
  H = nrow(mat) * 0.16 + 3
  args = list(mat, cluster_rows = crows, cluster_cols = FALSE,
              color = colorRampPalette(c("blue", "white", "red"))(250),
              breaks = seq(-lim, lim, length.out = 251),
              border_color = NA, cellwidth = 10, cellheight = 10, fontsize = 8,
              gaps_row = gaps, main = collection)
  if (!identical(annot, NA)) args$annotation_row = annot
  pdf(paste0(fbase, ".pdf"), width = W, height = H); do.call(pheatmap::pheatmap, args); dev.off()
  png(paste0(fbase, ".png"), width = W, height = H, units = "in", res = 300); do.call(pheatmap::pheatmap, args); dev.off()
}

for (vn in names(subcl_variants)){
  ml = create_path_comp_mat_list(gsea_res_tab, comps_sel = subcl_variants[[vn]], sig_cut = SIG_CUT)
  for (sc in names(ml)){
    coll  = if (sc == "HALLMARKS") "Hallmark" else if (sc == "user_def_sets") "Green" else sc
    fbase = paste0(out_dir, script_ind, "GSEA_heatmap_", vn, "_", coll, "_FDR10_allsig")
    draw_heatmap(ml[[sc]], coll, fbase)
    message("   wrote: ", basename(fbase), " (.pdf + .png)")
  }
}


### union row-counts: how many DISTINCT pathways (heatmap rows) per cutoff/variant

row_counts = NULL
for (vn in names(subcl_variants)) for (cn in names(cutoffs)) for (sc in unique(gsea_res_tab$sub_cat)){
  t2 = gsea_res_tab[gsea_res_tab$sub_cat == sc & gsea_res_tab$comp %in% subcl_variants[[vn]] &
                      !is.na(gsea_res_tab$padj) & gsea_res_tab$padj < cutoffs[[cn]], ]
  row_counts = rbind(row_counts, data.frame(subcluster_set = vn, cutoff = cn, sub_cat = sc,
                                            n_distinct_pathways = length(unique(t2$pathway))))
}
write_csv(row_counts, file = paste0(out_dir, script_ind, "GSEA_heatmap_row_counts.csv"))
message("\n##### distinct pathways (heatmap rows) per cutoff #####\n")
print(row_counts, row.names = FALSE)


### dropped report: pathways nominal-significant (p<0.05) but NOT FDR-significant

dropped = NULL
for (cn in names(cutoffs)){
  d = gsea_res_tab[!is.na(gsea_res_tab$pval) & gsea_res_tab$pval < 0.05 &
                     (is.na(gsea_res_tab$padj) | gsea_res_tab$padj >= cutoffs[[cn]]),
                   c("sub_cat", "comp", "pathway", "NES", "pval", "padj")]
  if (nrow(d) > 0){ d$cutoff = cn; dropped = rbind(dropped, d) }
}
write_csv(dropped, file = paste0(out_dir, script_ind, "GSEA_dropped_nominal_not_FDR.csv"))
message("   dropped (nominal-sig but not FDR-sig) rows written: ",
        ifelse(is.null(dropped), 0, nrow(dropped)))


#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed LD_E03c ", Sys.time(),
        "\n##########################################################################\n\n\n")
