message("\n\n##########################################################################\n",
        "# Start LD_E03b2: DEG characterisation - GSEA of astrocyte pseudobulk results", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n")

# what this script does:
## - runs GSEA on each DESeq2 comparison made in the previous script (Are whole groups of biologically related genes shifting together?)
## - saves and plots GSEA results

# Open packages necessary for analysis.
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
library(msigdbr) # Michael used msigdf; replaced with msigdbr (CRAN, actively maintained)
library(fgsea)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "LD_E03b2_v02_"

#specify output directory
out_dir = paste0(main_dir,"LD_E_DESeq_pseudobulk/")


# load astrocyte DEseq2 dataset
bulk_data = qread(file = paste0(out_dir, "LD_E02a2_v02_bulk_data.qs"))



### get user-defined gene sets for GSEA
# GSEA needs predefined gene sets. we use the standard MSigDB library plus additional custom signatures

gsea_sets = list()

# Green24 astrocyte subpopulation markers (astrocytes only)
t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/Green24_S2_subpopulation_markers.csv"))
t1 = t1[t1$cell.type == "astrocytes",]

for (cl in unique(t1$state)){
  t2 = t1[t1$state == cl & t1$avg_log2FC>log2(1.2) & t1$p_val_adj<0.05,]
  gsea_sets[[paste0("Green24_",cl)]] = t2$gene
}

t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/GOI_sets_251020.csv"))

for (set1 in unique(t1$gene_set)){
  t2 = t1[t1$gene_set == set1,]
  gsea_sets[[set1]] = t2$gene
}


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

category_codes = c("H", "C2", "C5", "C8") # msigdbr v10+ uses uppercase collection names
#category_codes = c("H", "C2")

gsea_res_tab = NULL

for (cat1 in category_codes){
  
  message("\n   *run GSEA for category ", cat1, " - ", Sys.time(),"\n")
  
  pathways = msigdbr(species = "Homo sapiens", collection = cat1) %>% # Michael used msigdf.human; replaced with msigdbr(); msigdbr v10+ uses collection= instead of category=
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

plot_path_comp_mat_list = function(path_comp_mat_list){
  
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
                            main = paste0(cat1))
  }
  return(p1)
}



###define subcategories and comparisons for plotting 

gsea_res_tab = bulk_data$gsea_res_tab


### extract NES matrix and plot (all pathways, all comps)

path_comp_mat_list = create_path_comp_mat_list(gsea_res_tab = gsea_res_tab)

pdf(file = paste0(out_dir,script_ind, 
                  "GSEA_heatmaps_by_sub_cat.pdf"), 
    width = 40, height = 70)
{
  plot_path_comp_mat_list(path_comp_mat_list)
}
dev.off()


### extract NES matrix and plot (all pathways, all comps)

path_comp_mat_list = create_top_path_comp_mat_list(gsea_res_tab = gsea_res_tab)

pdf(file = paste0(out_dir,script_ind, 
                  "GSEA_heatmaps_by_sub_cat_top10_pathways.pdf"), 
    width = 40, height = 50)
{
  plot_path_comp_mat_list(path_comp_mat_list)
}
dev.off()





#get info on version of R, used packages etc
sessionInfo()


message("\n\n##########################################################################\n",
        "# Completed LD_E03b2 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


