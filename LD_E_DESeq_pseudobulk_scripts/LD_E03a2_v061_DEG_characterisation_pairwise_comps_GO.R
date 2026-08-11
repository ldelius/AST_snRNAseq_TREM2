message("\n\n##########################################################################\n",
        "# Start LD_E03a2: Pairwise log2FC comparisons and GO ", Sys.time(),
        "\n##########################################################################\n\n")

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
library(fgsea)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "LD_E03a2_v02_"

#specify output directory
out_dir = paste0(main_dir, "LD_E_DESeq_pseudobulk/")


# load astrocyte DESeq2 dataset
bulk_data = qread(file = paste0(out_dir, "LD_E02a2_v02_bulk_data.qs"))


### define clusters and reference cluster

comp_clusters = unique(bulk_data$meta$cluster_name)

ref_cluster = "AST_SLC1A2_s0"


### define comparisons for log2FC comparison

deseq_res_list = bulk_data$deseq_results

names(deseq_res_list)

log2FC_comps = list()

l1 = list()

# Option 1: TREM2 variant effect vs general AD effect (per cluster)

# R62H vs CV (AD only) vs CV AD vs Control
for (cl in comp_clusters){
  l1[[paste0(cl, "_AD_R62H_vs_CV_vs_CV_AD_vs_Control")]] =
    c(paste0(cl, "_AD_TREM2_R62H_vs_CV"), paste0(cl, "_TREM2_CV_AD_vs_Control"))
}

# R47H vs CV (AD only) vs CV AD vs Control
for (cl in comp_clusters){
  l1[[paste0(cl, "_AD_R47H_vs_CV_vs_CV_AD_vs_Control")]] =
    c(paste0(cl, "_AD_TREM2_R47H_vs_CV"), paste0(cl, "_TREM2_CV_AD_vs_Control"))
}

# R47H vs R62H (AD only) vs CV AD vs Control
for (cl in comp_clusters){
  l1[[paste0(cl, "_AD_R47H_vs_R62H_vs_CV_AD_vs_Control")]] =
    c(paste0(cl, "_AD_TREM2_R47H_vs_R62H"), paste0(cl, "_TREM2_CV_AD_vs_Control"))
}


# Option 2: R62H vs R47H variant effects compared against each other (per cluster)

# R62H vs CV vs R47H vs CV
for (cl in comp_clusters){
  l1[[paste0(cl, "_AD_R62H_vs_CV_vs_R47H_vs_CV")]] =
    c(paste0(cl, "_AD_TREM2_R62H_vs_CV"), paste0(cl, "_AD_TREM2_R47H_vs_CV"))
}


# Option 3: astrocyte subtype identity vs disease/variant effects (per non-reference cluster)

# cluster vs ref vs CV AD vs Control
for (cl in comp_clusters[comp_clusters != ref_cluster]){
  l1[[paste0(cl, "_CV_vs_", ref_cluster, "_vs_CV_AD_vs_Control")]] =
    c(paste0(cl, "_CV_vs_", ref_cluster), paste0(cl, "_TREM2_CV_AD_vs_Control"))
}

# cluster vs ref vs R62H vs CV
for (cl in comp_clusters[comp_clusters != ref_cluster]){
  l1[[paste0(cl, "_CV_vs_", ref_cluster, "_vs_AD_R62H_vs_CV")]] =
    c(paste0(cl, "_CV_vs_", ref_cluster), paste0(cl, "_AD_TREM2_R62H_vs_CV"))
}

# cluster vs ref vs R47H vs CV
for (cl in comp_clusters[comp_clusters != ref_cluster]){
  l1[[paste0(cl, "_CV_vs_", ref_cluster, "_vs_AD_R47H_vs_CV")]] =
    c(paste0(cl, "_CV_vs_", ref_cluster), paste0(cl, "_AD_TREM2_R47H_vs_CV"))
}


# Option 4: R47H vs R62H vs each individual variant vs CV (per cluster)

# R47H vs R62H vs R62H vs CV
for (cl in comp_clusters){
  l1[[paste0(cl, "_AD_R47H_vs_R62H_vs_R62H_vs_CV")]] =
    c(paste0(cl, "_AD_TREM2_R47H_vs_R62H"), paste0(cl, "_AD_TREM2_R62H_vs_CV"))
}

# R47H vs R62H vs R47H vs CV
for (cl in comp_clusters){
  l1[[paste0(cl, "_AD_R47H_vs_R62H_vs_R47H_vs_CV")]] =
    c(paste0(cl, "_AD_TREM2_R47H_vs_R62H"), paste0(cl, "_AD_TREM2_R47H_vs_CV"))
}


# check and remove invalid comparisons (e.g. clusters with too few samples for DESeq)
for (comp1 in names(l1)){
  v1 = l1[[comp1]]
  if (is.null(deseq_res_list[[v1[1]]]) | is.null(deseq_res_list[[v1[2]]])){
    warning("Incomplete comparison ", comp1, " will be removed")
    l1 = l1[names(l1) != comp1]
  }
}

log2FC_comps = l1

names(log2FC_comps)

### get gene sets for alternative labels of log2FC plots

goi_sets = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/GOI_sets_251020.csv"))

t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/cell_type_markers_241219_w_astr_subtype_markers.csv"))
subtype_markers = t1[t1$lineage %in% c("Astrocyte", "AST", "AST_RG"),] # Michael used "Microglia"; replaced with astrocyte lineages


###get TF gene panels
GOI = list()
t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol



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


###function: plot bulk gene expression heatmap with annotations

bulkdata_heatmap = function(pl_mat, pl_meta, x_col, 
                            p_mat = FALSE,
                            pl_genes = NULL, 
                            meta_annot_cols = NULL,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            show_rownames = TRUE, show_colnames = TRUE,
                            color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                            lims = NULL,  cellwidth = 15, cellheight = 10, 
                            fontsize = 10, number_color = "grey70",
                            title = "Z-score vst-norm gene expression"){
  
  if (is.null(pl_genes)){pl_genes = rownames(pl_mat)}
  
  pl_mat = pl_mat[match(pl_genes, rownames(pl_mat), nomatch = 0),]
  pl_meta = pl_meta[match(pl_meta[[x_col]], colnames(pl_mat), nomatch = 0),]
  
  #define annotation bars
  
  if (!is.null(meta_annot_cols)){
    
    annot_col = data.frame(row.names = pl_meta[[x_col]]) #if not defined cbind converts factor values to factor levels
    
    for (col1 in meta_annot_cols){
      v1 = pl_meta[match(colnames(pl_mat), pl_meta[[x_col]]),][[col1]]
      if (is.numeric(v1)){v1 = factor(v1, levels = unique(pl_meta[[col1]][order(pl_meta[[col1]])]))} else {
        v1 = factor(v1, levels = unique(pl_meta[[col1]]))
      }
      annot_col = as.data.frame(cbind(annot_col, v1))
    }
    colnames(annot_col) = meta_annot_cols
    rownames(annot_col) = colnames(pl_mat)
    
    annot_colors = lapply(meta_annot_cols, function(x){
      v1 = pal(unique(annot_col[[x]]) )
      names(v1) = levels(annot_col[[x]])
      return(v1)
    })
    names(annot_colors) = meta_annot_cols
    
  }else{
    annot_col = NULL
    annot_colors = NULL
  }
  
  # define plot limits
  
  if (is.null(lims)){lims = c(-0.7*max(abs(na.omit(pl_mat))), 0.7*max(abs(na.omit(pl_mat))))}
  
  #create plot
  
  p1 = pheatmap::pheatmap(pl_mat, 
                          display_numbers = p_mat,
                          cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                          show_rownames = show_rownames, show_colnames = show_colnames,
                          color = color,
                          breaks = seq(lims[1], lims[2], length.out = length(color)+1),
                          annotation_col = annot_col, annotation_colors = annot_colors,
                          border_color = NA, cellwidth = cellwidth, cellheight = cellheight, 
                          fontsize = fontsize, number_color = number_color, fontsize_number = fontsize, 
                          main = title
  )
  
  return(p1)
  
}





###########################################################
# extract log2FC corr for pairwise comparison of log2FC_comps pairs 
###########################################################

###generate table with log2FC and p values, extract genes extract genes diff reg in comps

keep_genes = unique(unlist(bulk_data$DEGs))

log2FC_tab = NULL

comps_all = names(deseq_res_list)

for (comp in comps_all){
  
  t2 = deseq_res_list[[comp]]
  
  t3 = tibble(comp = comp, gene = keep_genes, 
              log2FC = t2$log2FoldChange[match(keep_genes, rownames(t2))],
              pval = t2$pvalue[match(keep_genes, rownames(t2))],
              padj = t2$padj[match(keep_genes, rownames(t2))])
  
  log2FC_tab = rbind(log2FC_tab, t3)
  
}


###extract genes diff expressed in reference comparison, check diff expression in comparisons

reg_groups = c("down_down", "up_up","down_up", "up_down", "down_nreg",  "up_nreg", 
         "nreg_down", "nreg_up", "nreg_nreg")

log2FC_comps_DEGs = list()

log2FC_comps_tab = NULL

for (comp in names(log2FC_comps)){
  
  comp1 = log2FC_comps[[comp]][1]
  comp_ref = log2FC_comps[[comp]][2]
  
  t1 = log2FC_tab[log2FC_tab$comp == comp1, ]
  t2 = log2FC_tab[log2FC_tab$comp == comp_ref, ]
  
  t3 = t1
  t3$reg = "nreg"
  t3$reg[(t3$pval<0.05 &!is.na(t3$pval) & t3$log2FC>0)] = "up"
  t3$reg[(t3$pval<0.05 &!is.na(t3$pval) & t3$log2FC<0)] = "down"
  
  t3$comp_ref = comp_ref
  t3$log2FC_ref = t2$log2FC
  t3$pval_ref = t2$pval
  t3$padj_ref = t2$padj
  t3$reg_ref = "nreg"
  t3$reg_ref[t3$pval_ref<0.05 &!is.na(t3$pval_ref) & t3$log2FC_ref>0] = "up"
  t3$reg_ref[t3$pval_ref<0.05 &!is.na(t3$pval_ref) & t3$log2FC_ref<0] = "down"
  
  #keep genes regulated in one or the other condition, group by reg in comp vs comp_ref
  
  t3 = t3[(t3$reg != "nreg" | t3$reg_ref!="nreg"), ]
  
  t3$reg_comp_ref[t3$reg != "nreg"] = "reg_comp"
  t3$reg_comp_ref[t3$reg_ref != "nreg"] = "reg_ref"
  t3$reg_comp_ref[t3$reg != "nreg" & t3$reg_ref != "nreg"] = "reg_comp_ref"
  
  t3$reg_group = paste0(t3$reg, "_", t3$reg_ref)
  
  #remove genes without determined pval in either comp or ref comp, collect all others in combined table
  
  t3 = t3[!is.na(t3$pval) & !is.na(t3$pval_ref),]
  
  log2FC_comps_tab = rbind(log2FC_comps_tab, t3)
  
  for (reg_group in reg_groups){
    log2FC_comps_DEGs[[paste0(comp, "_", reg_group)]] = t3$gene[t3$reg_group == reg_group]
  }
}  



### save table with genes with different response comparison vs reference comparison including TFs

l1 = log2FC_comps_DEGs
l2 = lapply(l1, function(x){x = x[x%in% GOI$TF]})
names(l2) = paste0(names(l1), "_TF")

l3 = c(l1, l2)

m1 = matrix(nrow = max(lengths(l3)), ncol = length(l3))
colnames(m1) = names(l3)

for (i in names(l3)){
  v1 = l3[[i]]
  if(length(v1)>0){
    m1[1:length(v1),i] = v1
  }
}
m1[is.na(m1)] = ""

write_csv(as_tibble(m1), file = paste0(out_dir,script_ind, "Pairwise_comps_DEGs.csv"))

t1 = tibble(gene_set = names(l3), N_genes = lengths(l3))

write_csv(t1, file = paste0(out_dir,script_ind, "Pairwise_comps_DEGs_N.csv"))


bulk_data$deseq_log2FC_comps$comps = log2FC_comps
bulk_data$deseq_log2FC_comps$log2FC_comps_tab = log2FC_comps_tab
bulk_data$deseq_log2FC_comps$log2FC_comps_DEGs = log2FC_comps_DEGs

qsave(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.qs"))



###########################################################
# plot log2FC corr for different comparisons vs reference comparisons 
#   colour by regulation in reference vs comparison, label genes with max reg (distance from origin) 
###########################################################

pl = list()

for (comp in names(log2FC_comps)){
  
  comp1 = log2FC_comps[[comp]][1]
  comp_ref = log2FC_comps[[comp]][2]
  
  t1 = log2FC_comps_tab
  
  t3 = t1[t1$comp == comp1 & t1$comp_ref == comp_ref,]
  
  #for each group, label 10 most regulated genes (furthest from origin)
  
  t3$label = ""
  
  for (reg_group in reg_groups){
    
    t4 = t3[t3$reg_group == reg_group, ]
    t4$dist = sqrt((t4$log2FC)^2 + (t4$log2FC_ref)^2)
    t5 = t4[order(-t4$dist),]
    
    if (nrow(t5)>10){t5 = t5[1:10,]}
    
    t3$label[t3$gene %in% t5$gene] = t3$gene[t3$gene %in% t5$gene]
    
  }
  
  #order to put ref reg to background, then comp reg, then on top both reg
  t3 = t3[order(match(t3$reg_group, reg_groups)),]
  
  xmax = max(abs(t3$log2FC_ref))
  ymax = max(abs(t3$log2FC))
  
  pl[[paste0(comp1, " log2FC vs ", comp_ref)]] = 
    ggplot(t3, aes(x = log2FC_ref, y = log2FC, color = reg_group))+
    geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_hline(yintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_vline(xintercept = c(0), linewidth = 0.3, color = "grey30", linetype = 1)+
    geom_hline(yintercept = c(0), linewidth = 0.3, color = "grey30", linetype = 1)+
    geom_smooth(color = "grey30", method = "lm", formula=y~x)+
    geom_point(alpha = 0.8, size = 2)+
    geom_label_repel(aes(label = label), seed = 42, min.segment.length = 0.2, max.overlaps = Inf,
                     max.time = 5)+
    scale_color_manual(limits = c("down_down", "up_up","down_up", "up_down", 
                                  "down_nreg",  "up_nreg", "nreg_down", "nreg_up", "nreg_nreg"), 
                       values = c("blue", "red", "magenta3", "magenta3", 
                                  "grey40", "grey40","grey40","grey40", "grey"))+
    coord_cartesian(xlim = c(-xmax,xmax), ylim = c(-ymax, ymax))+
    theme_minimal()+
    labs(title = paste0(comp1, " log2FC vs ", comp_ref),
         x = paste0("log2FC_", comp_ref),
         y = paste0("log2FC_", comp1),
         colour = paste0("comp pval<0.05"))
  
}
  

pdf(file = paste0(out_dir,script_ind, "Log2FC_corr_sel_comps.pdf"), 
    width = 9, height = 8)
{
  lapply(pl, function(x){x})
}
dev.off()



###########################################################
# plot log2FC corr for different comparisons vs reference comparisons (only DEGs in both comparisons)
#   colour by regulation in reference vs comparison, label genes with max reg (distance from origin) 
###########################################################

pl = list()

for (comp in names(log2FC_comps)){
  
  comp1 = log2FC_comps[[comp]][1]
  comp_ref = log2FC_comps[[comp]][2]
  
  t1 = log2FC_comps_tab
  
  t3 = t1[t1$comp == comp1 & t1$comp_ref == comp_ref,]
  
  t3 = t3[t3$reg != "nreg" & t3$reg_ref != "nreg",]
  
  #for each group, label 10 most regulated genes (furthest from origin)
  
  t3$label = ""
  
  for (reg_group in reg_groups){
    
    t4 = t3[t3$reg_group == reg_group, ]
    t4$dist = sqrt((t4$log2FC)^2 + (t4$log2FC_ref)^2)
    t5 = t4[order(-t4$dist),]
    
    if (nrow(t5)>10){t5 = t5[1:10,]}
    
    t3$label[t3$gene %in% t5$gene] = t3$gene[t3$gene %in% t5$gene]
    
  }
  
  #order to put concordantly regulated genes in foreground
  t3 = t3[order(match(t3$reg_group, reg_groups)),]
  
  xmax = max(abs(t3$log2FC_ref))
  ymax = max(abs(t3$log2FC))
  
  pl[[paste0(comp1, " log2FC vs ", comp_ref)]] = 
    ggplot(t3, aes(x = log2FC_ref, y = log2FC, color = reg_group))+
    geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_hline(yintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_vline(xintercept = c(0), linewidth = 0.3, color = "grey30", linetype = 1)+
    geom_hline(yintercept = c(0), linewidth = 0.3, color = "grey30", linetype = 1)+
    geom_smooth(color = "grey30", method = "lm", formula=y~x)+
    geom_point(alpha = 0.8, size = 2)+
    geom_label_repel(aes(label = label), seed = 42, min.segment.length = 0.2, max.overlaps = Inf,
                     max.time = 5)+
    scale_color_manual(limits = c("down_down", "up_up","down_up", "up_down", 
                                  "down_nreg",  "up_nreg", "nreg_down", "nreg_up", "nreg_nreg"), 
                       values = c("blue", "red", "magenta3", "magenta3", 
                                  "grey40", "grey40","grey40","grey40", "grey"))+
    coord_cartesian(xlim = c(-xmax,xmax), ylim = c(-ymax, ymax))+
    theme_minimal()+
    labs(title = paste0(comp1, " log2FC vs ", comp_ref),
         x = paste0("log2FC_", comp_ref),
         y = paste0("log2FC_", comp1),
         colour = paste0("comp pval<0.05"))
  
}


pdf(file = paste0(out_dir,script_ind, "Log2FC_corr_sel_comps_reg_both.pdf"), 
    width = 9, height = 8)
{
  lapply(pl, function(x){x})
}
dev.off()



###########################################################
# plot log2FC corr for different comparisons vs reference comparisons 
#   colour and label by subtype markers 
###########################################################

pl = list()

for (comp in names(log2FC_comps)){
  
  comp1 = log2FC_comps[[comp]][1]
  comp_ref = log2FC_comps[[comp]][2]
  
  t1 = log2FC_comps_tab
  
  t3 = t1[t1$comp == comp1 & t1$comp_ref == comp_ref,]
  
  t3 = t3[t3$reg != "nreg" & t3$reg_ref != "nreg",]
  
  t3$marker = subtype_markers$subtype[match(t3$gene, subtype_markers$gene)]
  
  #label subtype markers
  
  t3$label = ""
  t3$label[!is.na(t3$marker)] = t3$gene[!is.na(t3$marker)]
  
  #order to put labelled to foreground
  t3 = t3[order(t3$label, na.last = FALSE),]
  
  xmax = max(abs(t3$log2FC_ref))
  ymax = max(abs(t3$log2FC))
  
  pl[[paste0(comp1, " log2FC vs ", comp_ref)]] = 
    ggplot(t3, aes(x = log2FC_ref, y = log2FC, color = marker))+
    geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_hline(yintercept = c(-log2(1.2), log2(1.2)), linewidth = 0.3, color = "grey30", linetype = 2)+
    geom_vline(xintercept = c(0), linewidth = 0.3, color = "grey30", linetype = 1)+
    geom_hline(yintercept = c(0), linewidth = 0.3, color = "grey30", linetype = 1)+
    geom_smooth(color = "grey30", method = "lm", formula=y~x)+
    geom_point(alpha = 0.8, size = 2)+
    geom_label_repel(aes(label = label), seed = 42, min.segment.length = 0, max.overlaps = Inf,
                     max.time = 5)+
    scale_color_manual(limits = unique(subtype_markers$subtype), 
                       values = pal(unique(subtype_markers$subtype)) )+
    coord_cartesian(xlim = c(-xmax,xmax), ylim = c(-ymax, ymax))+
    theme_minimal()+
    labs(title = paste0(comp1, " log2FC vs ", comp_ref),
         x = paste0("log2FC_", comp_ref),
         y = paste0("log2FC_", comp1),
         colour = paste0("subtype markers"))
  
}


pdf(file = paste0(out_dir,script_ind, "Log2FC_corr_sel_comps_subtype_markers_highlighted.pdf"), 
    width = 9, height = 8)
{
  lapply(pl, function(x){x})
}
dev.off()




#################################################
### GO-BP over-representation analysis of genes by comparison and reg group (only reg in both comparisons)
#################################################

regs = c("up_up", "down_down", "up_down", "down_up")

GO_list = list()
GO_results_tab = NULL

l1 = log2FC_comps_DEGs
l2 = l1[!grepl("_nreg", names(l1))]

N_comps = length(l2)

for (i in 1:N_comps){
  
  comp = names(l2)[i]
  
  message("\n          *** GO analysis DEGs by comparison and cluster ", comp, 
          " (", i, " of ",N_comps,  ") - ", Sys.time(),"\n")
  
  ego = NULL
  
  go_genes =  l2[[comp]]
  
  if (length(go_genes)>2){
    
    ego = enrichGO(gene         = go_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
    
    GO_list[[comp]] = ego
  }
  
  if (!is.null(ego)){
    
    t1 = ego@result[ego@result$p.adjust<=0.05,]
    
    if (nrow(t1)>0){
      t2 = cbind(comp = comp, t1)
      GO_results_tab = rbind(GO_results_tab, t2)
    }
  }
}

### extract/add comparison, cluster, direction of change

t1 = GO_results_tab
t1 = t1[t1$Count>1,]

for (reg in regs){
  t1$up_down[grepl(reg, t1$comp)] = reg
}

for (i in 1:nrow(t1)){

  v1 = str_remove(t1$comp[i], paste0("_", t1$up_down[i]))
  t1$comp_group[i] = v1

}

GO_results_tab = t1


###save GO results

write_csv(GO_results_tab, file = paste0(out_dir,script_ind, "GO_results_pairwise_comps.csv"))

bulk_data$deseq_log2FC_comps$GO_results_tab = GO_results_tab

qsave(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.qs"))



########################################################################
# visualise GO analysis (dotplot by cluster for top 10 terms by comp)
########################################################################

t1 = GO_results_tab

t2 = NULL

for (cl in unique(t1$comp)){
  t3 = t1[t1$comp == cl,]
  if (nrow(t3)>10){t3 = t3[1:10,]}
  t2 = rbind(t2, t3)
}

pl = list()

for (comp_group in unique(t2$comp_group)){
  
  t3 = t2[t2$comp_group == comp_group,]
  
  pl[[comp_group]] = ggplot(t3, aes(x = paste0(comp_group, "_", up_down), 
                                    y = Description, size = Count, colour = up_down))+
    geom_point()+
    scale_color_manual (limits = regs, values =pal(regs))+
    scale_size_continuous(limits = c(0, max(t2$Count)))+
    scale_x_discrete(limits = unique(paste0(t3$comp_group, "_", t3$up_down)))+
    scale_y_discrete(limits = unique(t3$Description))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    labs(title = comp_group)
  
}

pdf(file = paste0(out_dir,script_ind, "GO_results_pairwise_comps_dotplot_top_terms.pdf"), 
    width = 10, height = 10)
lapply(pl, function(x){x})
dev.off()



### visualise GO analysis (dotplot by cluster for top 5 terms by comp combined)

t1 = GO_results_tab

t2 = NULL

for (cl in unique(t1$comp)){
  t3 = t1[t1$comp == cl,]
  if (nrow(t3)>10){t3 = t3[1:5,]}
  t2 = rbind(t2, t3)
}

terms_comb = unique(t2$Description)

t3 = t1[t1$Description %in% terms_comb, ]

p1 = ggplot(t3, aes(x = comp, y = Description, size = Count, colour = up_down))+
  geom_point()+
  scale_color_manual (limits = regs, values =pal(regs))+
  scale_size_continuous(limits = c(0, max(t3$Count)))+
  scale_x_discrete(limits = unique(t3$comp))+
  scale_y_discrete(limits = unique(t3$Description))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(title = "Top10 enriched GO terms for each comparison")


pdf(file = paste0(out_dir,script_ind, "GO_results_pairwise_comps_dotplot_top_terms_comb.pdf"), 
    width = 12, height = 10)
plot(p1)
dev.off()




#get info on version of R, used packages etc
sessionInfo()


message("\n\n##########################################################################\n",
        "# Completed LD_E03a2 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")
