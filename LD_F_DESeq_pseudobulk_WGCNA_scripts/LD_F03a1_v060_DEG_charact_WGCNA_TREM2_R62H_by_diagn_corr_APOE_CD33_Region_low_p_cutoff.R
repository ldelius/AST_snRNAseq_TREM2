message("\n\n##########################################################################\n",
        "# Start LD_F03a1: DEG characterisation WGCNA ", Sys.time(),
        "\n##########################################################################\n\n")

library(qs)
library(tidyverse)
library(DESeq2)
library(colorRamps)
library(viridis)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(ggrepel)
library(WGCNA)
library(enrichplot)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify output directory
in_dir  = paste0(main_dir,"LD_F_DESeq_pseudobulk_WGCNA/")
out_dir = paste0(main_dir,"LD_F_DESeq_pseudobulk_WGCNA/LD_F03a1_v03/")
dir.create(out_dir, showWarnings = FALSE)

#specify script/output index as prefix for file names
script_ind = "LD_F03a1_v04_"


### load DEseq2 dataset
bulk_data = qread(file = paste0(in_dir,"LD_F02a1_v03_bulk_data.rda"))



### get gene of interest gene sets

GOI = list()
t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol

t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/GOI_sets_251020.csv"))

for (goi_set in unique(t1$gene_set)){
  GOI[[goi_set]] = t1$gene[t1$gene_set == goi_set]
}

### get marker gene sets

t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/cell_type_markers_241219_w_astr_subtype_markers.csv"))
t2 = t1[t1$level == "Astrocyte_subtypes",]

subtype_marker_list = list()

for (subtype in unique(t2$subtype)){
  subtype_marker_list[[subtype]] = t2$gene[t2$subtype == subtype]
}


### define covariates for module - covariate correlation analysis

t1 = bulk_data$meta
names(t1)

t1$RNA_counts = apply(bulk_data$counts, 2, sum)

#define numeric measure for Braak staging

unique(t1$Braak)

Braak_conv = tibble(Braak_stage = c("0", "I", "II", "III", "IV", "V", "V,VI", "VI"),
                    Braak_numeric = c(0, 1, 2,3,4,5, 5.5, 6))

t1$Braak_numeric = Braak_conv$Braak_numeric[match(t1$Braak, Braak_conv$Braak_stage)]

bulk_data$meta = t1

names(t1)
covars1 = c("group", "Age", "Sex", "APOEgroup","CD33Group", 
            "Braak_numeric","BrainRegion",
            "plaque_dens", "pct4G8PositiveArea", "pctAT8PositiveArea", "pctPHF1PositiveArea",
            "PostMortemInterval", "N_cells", "RNA_counts", "cohort")
covars2 = c("all","NeuropathologicalDiagnosis", "TREM2Variant")



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
    p2 = c("dodgerblue", "grey20", "orange", "green4")
  } else if (v2<6){
    p2 = matlab.like(6)[1:v2]
  } else {
    p2 = matlab.like(v2)
  }
  return(p2)
}



###function: plot bulk gene expression heatmap with annotations (incl gene annotations)

bulkdata_heatmap = function(pl_mat, pl_meta, x_col, pl_genes = NULL, 
                            meta_annot_cols = NULL, gene_annot = NULL,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            show_rownames = T, show_colnames = T,
                            color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                            lims = NULL,  cellwidth = 15, cellheight = 10, 
                            fontsize = 10, title = "Z-score vst-norm gene expression"){
  
  if (is.null(pl_genes)){pl_genes = rownames(pl_mat)}
  
  pl_mat = pl_mat[match(pl_genes, rownames(pl_mat), nomatch = 0),]
  pl_meta = pl_meta[match(pl_meta[[x_col]], colnames(pl_mat), nomatch = 0),]
  
  #define annotation bars
  
  if (!is.null(meta_annot_cols)){
    
    annot_row = NULL
    annot_col = data.frame(row.names = pl_meta[[x_col]]) #if not defined cbind converts factor values to factor levels
    
    for (col1 in meta_annot_cols){
      v1 = pl_meta[match(colnames(pl_mat), pl_meta[[x_col]]),][[col1]]
      v1 = factor(v1, levels = unique(pl_meta[[col1]]))
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
    annot_row = NULL
    annot_colors = NULL
  }
  
  if (!is.null(gene_annot)){
    
    annot_row = data.frame(row.names = row.names(gene_annot))
    
    for (col1 in colnames(gene_annot)){
      v1 = gene_annot[[col1]]
      v1 = factor(v1, levels = unique(v1))
      annot_row = as.data.frame(cbind(annot_row, v1))
    }
    colnames(annot_row) = colnames(gene_annot)
    rownames(annot_row) = rownames(gene_annot)
    
    annot_colors_row = lapply(colnames(gene_annot), function(x){
      v1 = pal(unique(annot_row[[x]]) )
      names(v1) = levels(annot_row[[x]])
      return(v1)
    })
    names(annot_colors_row) = colnames(gene_annot)
    annot_colors = c(annot_colors, annot_colors_row)
  }
  
  # define plot limits
  
  if (is.null(lims)){lims = c(-0.7*max(abs(na.omit(pl_mat))), 0.7*max(abs(na.omit(pl_mat))))}
  
  #create plot
  
  p1 = pheatmap::pheatmap(pl_mat, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                          color = color,
                          breaks = seq(lims[1], lims[2], length.out = length(color)+1),
                          show_rownames = show_rownames, show_colnames = show_colnames,
                          annotation_col = annot_col, annotation_row = annot_row,
                          annotation_colors = annot_colors,
                          border_color = NA, cellwidth = cellwidth, cellheight = cellheight, 
                          fontsize = fontsize, main = title
  )
  
  return(p1)
  
}


### create heatmap annotations for Complex Heatmap

#annot_df = row_anno_df

create_heatmap_annot = function(annot_df, annot_dim = c("column", "row"), 
                                annotation_name_side = c("left", "top")){
  
  #define colors
  
  annot_colors = list()
  
  for (col1 in names(annot_df)){
    
    v1 = annot_df[[col1]]
    
    if (is.numeric(v1)){
      v2 = colorRamp2(breaks = seq(from = min(na.omit(v1)), to = max(na.omit(v1)), length.out = 100),
                      colors = viridis(100))
      annot_colors = c(annot_colors, list(v2))
      
    } else {
      v2 = pal(unique(v1))
      names(v2) = unique(v1)
      annot_colors = c(annot_colors, list(v2))
    }
  }
  names(annot_colors) = names(annot_df)
  
  # Create the Column Annotation Object
  annot <- HeatmapAnnotation(
    df = annot_df,
    col = annot_colors,
    which = annot_dim,
    annotation_name_side = annotation_name_side[1]
  )
  
  return(annot)
}



#######################################
# calculate WGCNA soft-thresholding power
#######################################

message("\n\n          *** Calculate soft-thresholding power... ", Sys.time(),"\n\n")

input_mat = t(bulk_data$vst_mat[unique(unlist(bulk_data$DEGs)),])
bulk_data$wgcna$input_mat = input_mat


allowWGCNAThreads()          # allow multi-threading (optional)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

### plot thresholding power graphs

t1 = as_tibble(sft$fitIndices)
t1$signed_R2 = -sign(t1$slope)*t1$SFT.R.sq

wgcna_power = t1$Power[t1$signed_R2>0.8][1]
if (wgcna_power<6){wgcna_power=6}

bulk_data$wgcna$power_table = t1
bulk_data$wgcna$power = wgcna_power

p1 = ggplot(t1, aes(x = Power, y = signed_R2))+geom_point()+
  geom_line(color = "grey")+geom_label(aes(label = Power))+
  geom_hline(yintercept = 0.8, color = "red")+
  geom_vline(xintercept = wgcna_power, color = "red")+
  labs(title = "Scale independence", x = "Soft Threshold (power)", 
       y = "Scale Free Topology Model Fit, signed R^2")

p2 = ggplot(t1, aes(x = Power, y = mean.k.))+geom_point()+
  geom_line()+geom_label(aes(label = Power))+
  labs(title = "Mean connectivity", x = "Soft Threshold (power)", 
       y = "Mean connectivity")


pdf(file = paste0(out_dir,script_ind, "Power_thresholding_tests.pdf"), 
    width = 9, height = 4)
{
  plot(p1+p2)
}
dev.off()



qsave(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.qs")) 



#######################################
# calculate WGCNA with selected soft-thresholding power
#######################################

temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)

netwk <- blockwiseModules(input_mat,                
                          power = wgcna_power,                
                          networkType = "signed",
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 20000,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          saveTOMs = FALSE,
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor

bulk_data$wgcna$network = netwk

qsave(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.qs")) 



#######################################
# Basic network characterisation and extraction of modules
#######################################

netwk = bulk_data$wgcna$network

### extract gene - module assignment

t1 <- data.frame(
  gene = names(netwk$colors),
  module = paste0("M", netwk$colors),
  module_number = netwk$colors,
  colors = labels2colors(netwk$colors)
)

mod_gene_tab = t1[order(t1$module_number),] 

mods = unique(mod_gene_tab$module)

bulk_data$wgcna$mod_gene_tab = mod_gene_tab

for (mod1 in mods){
  bulk_data$wgcna$mod_genes[[mod1]] = mod_gene_tab$gene[mod_gene_tab$module == mod1]
}

qsave(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.qs")) 



### save table with module genes and TFs

l1 = bulk_data$wgcna$mod_genes
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

write_csv(as_tibble(m1), file = paste0(out_dir,script_ind, "Module_genes.csv"))

t1 = tibble(gene_set = names(l3), N_genes = lengths(l3))

write_csv(t1, file = paste0(out_dir,script_ind, "Module_genes_N.csv"))



###plot network dendrogramme

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)

# Plot the dendrogram and the module colors underneath

pdf(file = paste0(out_dir,script_ind, "WGCNA_network_dendrogram.pdf"), 
    width = 10, height = 3)
{
  plotDendroAndColors(
    netwk$dendrograms[[1]],
    mergedColors,
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05 )
}
dev.off()


###extract module eigengene matrix

m1 = t(netwk$MEs)
rownames(m1) = str_replace_all(rownames(m1), "ME", "M")
m2 = m1[mods,]

bulk_data$wgcna$mod_eigengene_mat = m2



#######################################################
###plot module eigengene heatmap by sample and module
#######################################################

meta = bulk_data$meta

pl_mat_X = bulk_data$wgcna$mod_eigengene_mat
lims_X = 0.2*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


#Create Column Annotation Data (Top Bars)

col_anno_df <- as.data.frame(meta[,c("cluster_name", "APOEgroup", "CD33Group", "BrainRegion",
                                     "plaque_dens",
                                     "NeuropathologicalDiagnosis", "TREM2Variant")])
for (col1 in names(col_anno_df)){
  if (!is.numeric(col_anno_df[[col1]])){
    col_anno_df[[col1]] = factor(col_anno_df[[col1]], levels = unique(col_anno_df[[col1]]))}
  }
rownames(col_anno_df) = meta$cluster_sample

col_annot = create_heatmap_annot(annot_df = col_anno_df)


#plot heatmap

ht_list <- Heatmap(
  pl_mat_X,
  name = "Z-score",
  col = colorRamp2(breaks = seq(from = lims_X[1], to = lims_X[2], length.out = 100),
                   colors = viridis(100)),
  cluster_columns = FALSE,
  column_title = "Cluster_samples",
  top_annotation = col_annot,
  cluster_rows = TRUE,
  row_title = "Modules",
  left_annotation = NULL, 
  show_row_names = TRUE, 
  show_column_names = FALSE,
  width = 15, height = 5
)


pdf(file = paste0(out_dir,script_ind, "Module_eigengene_heatmap_by_sample.pdf"), 
    width = 10, height = 10)
{
  draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
}
dev.off()


### Exploratory module-trait correlation screen (WGCNA approach)
# Cluster-sample pseudobulks pool clusters, so trait correlations may reflect
# cluster identity as well as the biology of interest; treat them as exploratory.

message("\n\n   *Compute exploratory module-trait correlation screen \n")

# --- Build trait table aligned to eigengene columns ---

me_mat   = bulk_data$wgcna$mod_eigengene_mat   # modules x samples
meta_mt  = bulk_data$meta
meta_mt$RNA_counts = apply(bulk_data$counts[, meta_mt$cluster_sample], 2, sum)

# ensure row order matches eigengene column order
meta_mt = meta_mt[match(colnames(me_mat), meta_mt$cluster_sample), ]

# Numeric traits (direct)
numeric_traits = c("plaque_dens", "pct4G8PositiveArea", "pctAT8PositiveArea",
                   "pctPHF1PositiveArea", "Braak_numeric", "Age",
                   "PostMortemInterval", "N_cells", "RNA_counts")
trait_df = meta_mt[, numeric_traits, drop = FALSE]

# Categorical traits as binary dummies
# TREM2Variant: CV is reference; one dummy per non-reference level.
# Each risk-variant dummy compares that variant vs CV only — samples
# carrying *other* non-CV variants are set to NA so they are dropped
# from the per-trait correlation (pairwise complete cases).
for (lv in setdiff(levels(meta_mt$TREM2Variant), "CV")) {
  dummy = rep(NA_integer_, nrow(meta_mt))
  dummy[meta_mt$TREM2Variant == "CV"] = 0L
  dummy[meta_mt$TREM2Variant == lv]   = 1L
  trait_df[[paste0("TREM2_", lv)]] = dummy
}

# APOEgroup: APOE4-pos vs APOE4-neg
trait_df[["APOE4_pos"]]  = as.integer(meta_mt$APOEgroup  == "APOE4-pos")

# CD33Group: CD33var vs CV
trait_df[["CD33var"]]    = as.integer(meta_mt$CD33Group   == "CD33var")

# BrainRegion: SSC vs MTG
trait_df[["Region_SSC"]] = as.integer(meta_mt$BrainRegion == "SSC")

# cohort: dummies for each non-reference level
cohort_levels = levels(meta_mt$cohort)
for (lv in cohort_levels[-1]) {
  trait_df[[paste0("Cohort_", lv)]] = as.integer(meta_mt$cohort == lv)
}

# Drop constant columns (zero variance after removing NAs)
trait_df = trait_df[, sapply(trait_df, function(x) {
  x2 = x[!is.na(x)]
  length(x2) > 1 && var(x2) > 0
}), drop = FALSE]

trait_mat = as.matrix(trait_df)

# --- Pearson correlation per module x trait pair ---

n_mods      = nrow(me_mat)
n_traits    = ncol(trait_mat)
mod_names   = rownames(me_mat)
trait_names = colnames(trait_mat)

cor_mat  = matrix(NA, nrow = n_mods, ncol = n_traits,
                  dimnames = list(mod_names, trait_names))
pval_mat = matrix(NA, nrow = n_mods, ncol = n_traits,
                  dimnames = list(mod_names, trait_names))
n_mat    = matrix(NA, nrow = n_mods, ncol = n_traits,
                  dimnames = list(mod_names, trait_names))

for (mod1 in mod_names) {
  eg = me_mat[mod1, ]
  for (tr1 in trait_names) {
    tr   = trait_mat[, tr1]
    ok   = complete.cases(eg, tr)
    n_ok = sum(ok)
    if (n_ok >= 5) {
      r = cor(eg[ok], tr[ok], method = "pearson")
      cor_mat[mod1,  tr1] = r
      n_mat[mod1,    tr1] = n_ok
      pval_mat[mod1, tr1] = corPvalueStudent(r, n_ok)
    }
  }
}

# BH FDR across all module x trait pairs
padj_mat = matrix(p.adjust(as.vector(pval_mat), method = "BH"),
                  nrow = n_mods, ncol = n_traits,
                  dimnames = list(mod_names, trait_names))

# --- Save long-format results as CSV ---

res_long = expand.grid(module = mod_names, trait = trait_names,
                       stringsAsFactors = FALSE) %>%
  mutate(
    r      = as.vector(cor_mat),
    pvalue = as.vector(pval_mat),
    padj   = as.vector(padj_mat),
    n      = as.vector(n_mat)
  )

write_csv(res_long,
          file = paste0(out_dir, script_ind, "Module_trait_correlation_results.csv"))

# --- Heatmap helpers ---

make_cor_col = function(mat) {
  mx = max(abs(mat), na.rm = TRUE)
  colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))
}

# significance overlay (BH-FDR only): * < 0.05, ** < 0.01, *** < 0.001
sig_label = function(p_raw, p_adj) {
  ifelse(is.na(p_adj), "",
    ifelse(p_adj < 0.001, "***",
      ifelse(p_adj < 0.01, "**",
        ifelse(p_adj < 0.05, "*", ""))))
}

cell_labels = matrix(sig_label(as.vector(pval_mat), as.vector(padj_mat)),
                     nrow = n_mods, ncol = n_traits,
                     dimnames = list(mod_names, trait_names))

# --- Full heatmap (all traits) ---

ht_cor_full = Heatmap(
  cor_mat,
  name              = "Pearson r",
  col               = make_cor_col(cor_mat),
  cell_fun          = function(j, i, x, y, width, height, fill) {
    grid.text(cell_labels[i, j], x, y, gp = gpar(fontsize = 8))
  },
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  row_title         = "Modules",
  column_title      = "Exploratory module-trait screen  (BH-FDR: * <0.05, ** <0.01, *** <0.001)",
  show_row_names    = TRUE,
  show_column_names = TRUE,
  na_col            = "grey90",
  column_names_rot  = 45,
  column_names_gp   = gpar(fontsize = 11),
  row_names_gp      = gpar(fontsize = 11),
  column_title_gp   = gpar(fontsize = 10, fontface = "plain")
)

pdf(file  = paste0(out_dir, script_ind, "Module_trait_correlation_heatmap.pdf"),
    width  = max(n_traits * 0.5 + 3, 8),
    height = max(n_mods   * 0.4 + 3, 6))
{ draw(ht_cor_full, heatmap_legend_side = "right") }
dev.off()

# --- Filtered heatmap: traits with at least one BH-FDR significant module ---

sig_traits = trait_names[apply(padj_mat, 2, function(p) any(p < 0.05, na.rm = TRUE))]

if (length(sig_traits) >= 2) {

  cor_mat_sig  = cor_mat[,      sig_traits, drop = FALSE]
  cell_lbl_sig = cell_labels[,  sig_traits, drop = FALSE]

  ht_cor_sig = Heatmap(
    cor_mat_sig,
    name              = "Pearson r",
    col               = make_cor_col(cor_mat_sig),
    cell_fun          = function(j, i, x, y, width, height, fill) {
      grid.text(cell_lbl_sig[i, j], x, y, gp = gpar(fontsize = 8))
    },
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    row_title         = "Modules",
    column_title      = "Significant traits  (BH-FDR: * <0.05, ** <0.01, *** <0.001)",
    show_row_names    = TRUE,
    show_column_names = TRUE,
    na_col            = "grey90",
    column_names_rot  = 45,
    column_names_gp   = gpar(fontsize = 11),
    row_names_gp      = gpar(fontsize = 11),
    column_title_gp   = gpar(fontsize = 10, fontface = "plain")
  )

  pdf(file  = paste0(out_dir, script_ind, "Module_trait_correlation_heatmap_sig_traits.pdf"),
      width  = max(length(sig_traits) * 0.5 + 3, 6),
      height = max(n_mods * 0.4 + 3, 6))
  { draw(ht_cor_sig, heatmap_legend_side = "right") }
  dev.off()

} else {
  message("   *Fewer than 2 BH-FDR significant traits found; skipping filtered heatmap")
}


######################################################################
### Per-cluster multivariable LM block — commented out (2026-04-21)
### Supervisor not interested; keeping for future reference.
######################################################################
# message("\n\n   *Fit per-cluster multivariable LM per module (omnibus F-tests) \n")
#
# me_mat_lm = bulk_data$wgcna$mod_eigengene_mat
# meta_lm   = bulk_data$meta
# meta_lm   = meta_lm[match(colnames(me_mat_lm), meta_lm$cluster_sample), ]
#
# lm_formula_rhs = "TREM2Variant + cohort + BrainRegion + APOEgroup + CD33Group"
# lm_terms_keep  = c("TREM2Variant", "cohort", "BrainRegion", "APOEgroup", "CD33Group")
#
# mod_names_lm     = rownames(me_mat_lm)
# cluster_names_lm = as.character(unique(meta_lm$cluster_name))
#
# lm_res = list()
#
# for (cl in cluster_names_lm) {
#   idx = which(meta_lm$cluster_name == cl)
#   if (length(idx) < 10) {
#     message("   *Skipping cluster ", cl, " (n=", length(idx), " < 10)")
#     next
#   }
#
#   meta_cl = droplevels(meta_lm[idx, ])
#
#   for (mod1 in mod_names_lm) {
#     y = me_mat_lm[mod1, idx]
#     df_fit = data.frame(eigengene = y, meta_cl, check.names = FALSE)
#
#     fit = try(lm(as.formula(paste("eigengene ~", lm_formula_rhs)), data = df_fit),
#               silent = TRUE)
#     if (inherits(fit, "try-error")) next
#
#     # omnibus per-term F-tests (Type-II-like, partial F for each predictor)
#     a = try(drop1(fit, test = "F"), silent = TRUE)
#     if (inherits(a, "try-error")) next
#     a = a[rownames(a) != "<none>", , drop = FALSE]
#     if (nrow(a) == 0) next
#
#     rss_full = sum(residuals(fit)^2)
#     ss_eff   = a[, "Sum of Sq"]
#
#     lm_res[[length(lm_res) + 1]] = data.frame(
#       cluster      = cl,
#       module       = mod1,
#       predictor    = rownames(a),
#       df_num       = a[, "Df"],
#       df_den       = fit$df.residual,
#       ss_effect    = ss_eff,
#       F_stat       = a[, "F value"],
#       pvalue       = a[, "Pr(>F)"],
#       partial_eta2 = ss_eff / (ss_eff + rss_full),
#       n            = nrow(df_fit),
#       stringsAsFactors = FALSE
#     )
#   }
# }
#
# lm_long = do.call(rbind, lm_res)
#
# # BH FDR across all (cluster x module x predictor) tests
# lm_long$padj = p.adjust(lm_long$pvalue, method = "BH")
#
# write_csv(lm_long,
#           file = paste0(out_dir, script_ind, "Module_trait_lm_results.csv"))
#
# lm_long$sig_label = ifelse(is.na(lm_long$pvalue), "",
#                     ifelse(lm_long$padj   < 0.05, "**",
#                     ifelse(lm_long$pvalue < 0.05, "*", "")))
#
# predictor_order = lm_terms_keep
#
# ht_list_lm = list()
#
# for (cl in unique(lm_long$cluster)) {
#   d_cl = lm_long[lm_long$cluster == cl, ]
#   eta_mat = matrix(NA_real_, nrow = length(mod_names_lm), ncol = length(predictor_order),
#                    dimnames = list(mod_names_lm, predictor_order))
#   lbl_mat = matrix("",       nrow = length(mod_names_lm), ncol = length(predictor_order),
#                    dimnames = list(mod_names_lm, predictor_order))
#   for (i in seq_len(nrow(d_cl))) {
#     if (d_cl$predictor[i] %in% predictor_order) {
#       eta_mat[d_cl$module[i], d_cl$predictor[i]] = d_cl$partial_eta2[i]
#       lbl_mat[d_cl$module[i], d_cl$predictor[i]] = d_cl$sig_label[i]
#     }
#   }
#
#   mx = max(eta_mat, na.rm = TRUE)
#   if (!is.finite(mx) || mx == 0) mx = 1
#   col_fun = colorRamp2(c(0, mx), c("white", "red"))
#
#   ht_list_lm[[cl]] = Heatmap(
#     eta_mat,
#     name              = "partial eta^2",
#     col               = col_fun,
#     cell_fun          = function(j, i, x, y, width, height, fill) {
#       grid.text(lbl_mat[i, j], x, y, gp = gpar(fontsize = 8))
#     },
#     cluster_rows      = FALSE,
#     cluster_columns   = FALSE,
#     row_title         = "Modules",
#     column_title      = paste0(cl, " - LM omnibus (eigengene ~ ", lm_formula_rhs, ")"),
#     show_row_names    = TRUE,
#     show_column_names = TRUE,
#     na_col            = "grey90",
#     column_names_rot  = 45,
#     column_names_gp   = gpar(fontsize = 9),
#     row_names_gp      = gpar(fontsize = 9)
#   )
# }
#
# largest_lm <- names(ht_list_lm)[which.max(sapply(ht_list_lm, function(h) nrow(h@matrix)))]
# pdf(NULL)
# ht_drawn <- draw(ht_list_lm[[largest_lm]], heatmap_legend_side = "right")
# lm_pdf_height <- convertHeight(ComplexHeatmap:::height(ht_drawn), "inches", valueOnly = TRUE)
# lm_pdf_width  <- convertWidth(ComplexHeatmap:::width(ht_drawn),  "inches", valueOnly = TRUE)
# dev.off()
#
# pdf(file  = paste0(out_dir, script_ind, "Module_trait_lm_heatmap_by_cluster.pdf"),
#     width  = lm_pdf_width,
#     height = lm_pdf_height)
# for (cl in names(ht_list_lm)) {
#   draw(ht_list_lm[[cl]], heatmap_legend_side = "right")
# }
# dev.off()


######################################################################
###plot sample module eigengene vs covars by cluster and module
######################################################################

message("\n\n   *Generate module eigengene vs covariate plots \n")

m1 = bulk_data$wgcna$mod_eigengene_mat

meta = bulk_data$meta

gr = unique(meta$group)

cluster_names = unique(meta$cluster_name)

pl = list()

t1 = NULL

for (mod1 in rownames(m1)){
  
  t2 = cbind(meta, module = mod1, eigengene = m1[mod1, ])
  
  t1 = rbind(t1, t2)
}
  
t1$cluster_name = factor(t1$cluster_name, levels = cluster_names)
t1$all = "all"

meta_me = t1

mod1 = "M0"

# cov2 in covars2
# cov1 in covars1
# mod1 in mods

for (cov2 in covars2){
  
  message("   *Generate covariate plots for module ", mod1)
  
  for (cov1 in covars1){
    
    t1 = meta_me
    
    p1 = ggplot(data = t1, aes(x = .data[[cov1]], y = eigengene))+
      theme_bw()
    
    if (is.numeric(t1[[cov1]]) & !is.numeric(t1[[cov2]])){
      
      p1 = p1+
        geom_smooth(aes(color = .data[[cov2]], group = .data[[cov2]]), 
                    method = "lm")+
        geom_point(aes(fill = .data[[cov2]], group = .data[[cov2]]), 
                   shape = 21, size = 2, stroke = 0.3, color = "grey40")+
        scale_fill_manual(limits = unique(t1[[cov2]]), values = pal(unique(t1[[cov2]])))+
        scale_color_manual(limits = unique(t1[[cov2]]), values = pal(unique(t1[[cov2]])))
    }
    
    if (is.numeric(t1[[cov1]]) & is.numeric(t1[[cov2]])){
      
      lims_cov2 = range(t1[[cov2]])
      
      p1 = p1+
        geom_point(aes(fill = .data[[cov2]], group = .data[[cov2]]), 
                   shape = 21, size = 2, stroke = 0.3, color = "grey40")+
        scale_fill_viridis(limits = lims_cov2)
    }
    
    if (!is.numeric(t1[[cov1]]) & !is.numeric(t1[[cov2]])){
      
      p1 = p1+
        geom_point(aes(fill = .data[[cov2]]), shape = 21, 
                   position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
                   size = 2, stroke = 0.3, color = "grey50")+
        scale_x_discrete(limits = unique(t1[[cov1]]))+
        scale_fill_manual(limits = unique(t1[[cov2]]), values = pal(unique(t1[[cov2]])))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    }
    
    if (!is.numeric(t1[[cov1]]) & is.numeric(t1[[cov2]])){
      
      lims_cov2 = range(t1[[cov2]])
      
      p1 = p1+
        geom_point(aes(fill = .data[[cov2]], group = .data[[cov2]]), 
                   shape = 21, size = 2, stroke = 0.3, color = "grey40", 
                   position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+
        scale_x_discrete(limits = unique(t1[[cov1]]))+
        scale_fill_viridis(limits = lims_cov2)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    }
    
    pl[[paste0(cov1, "_by_", cov2)]] = p1 + 
      facet_grid(module~cluster_name, scales = "free")+
      labs(title = paste0("Module eigengene vs ", cov1, "_by_", cov2))
  }
}


message("\n\n   *Plot module eigengene vs covariate plots (total ", length(pl)," plots) \n")

pdf(file = paste0(out_dir,script_ind, "Module_eigengene_vs_covars_scatterplot_by_cluster.pdf"), 
    width = length(cluster_names)*1.5+1, height = length(mods)*1.5+2)
{
  lapply(pl, function(x){x})
}
dev.off()



######################################################################
###plot module eigengene vs group scatterplot by sample and module
### + Kruskal-Wallis p-value per panel
### + pairwise Wilcoxon results saved to CSV with BH-FDR
######################################################################

m1 = bulk_data$wgcna$mod_eigengene_mat

meta = bulk_data$meta

gr = unique(meta$group)

cluster_names = unique(meta$cluster_name)

pl = list()
pw_res = list()
kw_res = list()

# helper: format p-value short
fmt_p = function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  formatC(p, digits = 3, format = "f")
}

for (mod1 in rownames(m1)){

  for (cl in cluster_names){

    meta_cl = meta[meta$cluster_name == cl,]

    t1 = cbind(meta_cl, module = mod1, eigengene = m1[mod1, meta_cl$cluster_sample])

    grp_present = intersect(gr, unique(t1$group))

    # overall Kruskal-Wallis across groups present
    kw_p = NA_real_
    if (length(grp_present) >= 2) {
      grp_ns = table(factor(t1$group, levels = grp_present))
      if (sum(grp_ns >= 2) >= 2) {
        kw_p = suppressWarnings(
          kruskal.test(eigengene ~ factor(group, levels = grp_present), data = t1)$p.value
        )
      }
    }
    kw_res[[length(kw_res) + 1]] = data.frame(
      cluster = cl, module = mod1,
      n_groups = length(grp_present),
      n_total  = nrow(t1),
      kw_pvalue = kw_p,
      stringsAsFactors = FALSE
    )

    # pairwise Wilcoxon across all group pairs present (saved to CSV only)
    pairs = combn(grp_present, 2, simplify = FALSE)

    pair_df = data.frame(
      cluster = cl, module = mod1,
      g1 = sapply(pairs, `[`, 1), g2 = sapply(pairs, `[`, 2),
      n1 = NA_integer_, n2 = NA_integer_,
      pvalue = NA_real_, stringsAsFactors = FALSE
    )
    for (pi in seq_along(pairs)) {
      x = t1$eigengene[t1$group == pairs[[pi]][1]]
      y = t1$eigengene[t1$group == pairs[[pi]][2]]
      pair_df$n1[pi] = length(x)
      pair_df$n2[pi] = length(y)
      if (length(x) >= 2 && length(y) >= 2) {
        pair_df$pvalue[pi] = suppressWarnings(
          wilcox.test(x, y, exact = FALSE)$p.value
        )
      }
    }
    pw_res[[length(pw_res) + 1]] = pair_df

    pl[[paste0(mod1, "_", cl)]] = ggplot(t1, aes(x = group, y = eigengene))+
      geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.6,
                   color = "grey30", width = 0.6) +
      geom_point(aes(fill = group), shape = 21,
                 position = position_jitter(width = 0.15, height = 0),
                 size = 1, stroke = 0.3, color = "grey30")+
      geom_hline(yintercept = 0)+
      scale_x_discrete(limits = gr)+
      scale_color_manual(limits = gr, values = pal(gr))+
      scale_fill_manual(limits = gr, values = pal(gr))+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            legend.position = "none")+
      labs(title = paste0(mod1, " - ", cl, " - Module eigengene by group"))

  }
}

# combine Kruskal-Wallis results, BH across all (cluster x module)
kw_long = do.call(rbind, kw_res)
kw_long$kw_padj = p.adjust(kw_long$kw_pvalue, method = "BH")

write_csv(kw_long,
          file = paste0(out_dir, script_ind, "Module_eigengene_vs_group_kruskal_wallis.csv"))

# combine pairwise results, BH across all (cluster x module x pair)
pw_long = do.call(rbind, pw_res)
pw_long$padj = p.adjust(pw_long$pvalue, method = "BH")

write_csv(pw_long,
          file = paste0(out_dir, script_ind, "Module_eigengene_vs_group_pairwise_wilcoxon.csv"))


pdf(file = paste0(out_dir,script_ind, "Module_eigengene_vs_group_scatterplot_by_cluster.pdf"),
    width = 6, height = 3)
{
  lapply(pl, function(x){x})
}
dev.off()



######################################################################
### Across-cluster LMM per module (TREM2Variant effect)
### Pools pseudobulks across clusters to gain power; controls for
### cluster-level mean offsets and donor-level pseudoreplication via
### random intercepts. One omnibus test per module + pairwise contrasts.
######################################################################

message("\n\n          *** Across-cluster LMM per module (TREM2Variant) - ", Sys.time(),"\n\n")

library(lme4)
library(lmerTest)

  me_mat   = bulk_data$wgcna$mod_eigengene_mat
  meta_lmm = bulk_data$meta
  meta_lmm = meta_lmm[match(colnames(me_mat), meta_lmm$cluster_sample), ]

  lmm_omni_list  = list()
  lmm_pairs_list = list()

  for (mod1 in rownames(me_mat)){

    df1 = data.frame(
      eigengene    = as.numeric(me_mat[mod1, ]),
      TREM2Variant = factor(meta_lmm$TREM2Variant),
      cluster_name = factor(meta_lmm$cluster_name),
      sample       = factor(meta_lmm$sample),
      stringsAsFactors = FALSE
    )
    df1 = df1[complete.cases(df1), ]
    df1$TREM2Variant = droplevels(df1$TREM2Variant)
    df1$cluster_name = droplevels(df1$cluster_name)
    df1$sample       = droplevels(df1$sample)

    if (nlevels(df1$TREM2Variant) < 2) next

    fit = tryCatch(
      lmerTest::lmer(eigengene ~ TREM2Variant + (1 | cluster_name) + (1 | sample),
                     data = df1, REML = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    # Omnibus Type-II F-test (Satterthwaite df via lmerTest) for TREM2Variant
    aov1 = tryCatch(anova(fit, type = 2), error = function(e) NULL)
    get_aov = function(col){
      if (!is.null(aov1) && "TREM2Variant" %in% rownames(aov1) && col %in% colnames(aov1)){
        aov1["TREM2Variant", col]
      } else NA_real_
    }

    lmm_omni_list[[mod1]] = data.frame(
      module     = mod1,
      n_obs      = nrow(df1),
      n_clusters = nlevels(df1$cluster_name),
      n_samples  = nlevels(df1$sample),
      F_value    = get_aov("F value"),
      num_df     = get_aov("NumDF"),
      den_df     = get_aov("DenDF"),
      pvalue     = get_aov("Pr(>F)"),
      stringsAsFactors = FALSE
    )

    # Pairwise contrasts between TREM2 variants (Satterthwaite df via lmerTest::contest)
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
      lmm_pairs_list[[mod1]] = do.call(rbind, pairs_rows)
    }
  }

  lmm_omni  = do.call(rbind, lmm_omni_list)
  lmm_pairs = do.call(rbind, lmm_pairs_list)

  if (!is.null(lmm_omni) && nrow(lmm_omni) > 0){
    lmm_omni$padj_BH = p.adjust(lmm_omni$pvalue, method = "BH")
    write_csv(lmm_omni,
              file = paste0(out_dir, script_ind, "Module_eigengene_LMM_across_clusters_omnibus.csv"))
  }

  if (!is.null(lmm_pairs) && nrow(lmm_pairs) > 0){
    lmm_pairs$padj_BH = p.adjust(lmm_pairs$p.value, method = "BH")
    write_csv(lmm_pairs,
              file = paste0(out_dir, script_ind, "Module_eigengene_LMM_across_clusters_pairwise.csv"))

    # Forest plot: pairwise contrasts per module, estimate ± 95% CI
    fp = lmm_pairs
    fp$ci_lo     = fp$estimate - 1.96 * fp$SE
    fp$ci_hi     = fp$estimate + 1.96 * fp$SE
    fp$sig_label = ifelse(is.na(fp$padj_BH), "",
                   ifelse(fp$padj_BH < 0.001, "***",
                   ifelse(fp$padj_BH < 0.01,  "**",
                   ifelse(fp$padj_BH < 0.05,  "*", ""))))

    p_fp = ggplot(fp, aes(x = estimate, y = module, color = contrast))+
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60")+
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25,
                     position = position_dodge(width = 0.6))+
      geom_point(position = position_dodge(width = 0.6), size = 2)+
      geom_text(aes(label = sig_label),
                position = position_dodge(width = 0.6),
                hjust = -0.4, size = 4, show.legend = FALSE)+
      scale_y_discrete(limits = rev(rownames(me_mat)))+
      labs(title = "Across-cluster LMM: TREM2 variant contrasts per module",
           subtitle = "Estimate ± 95% CI; BH-FDR: * <0.05, ** <0.01, *** <0.001",
           x = "Eigengene difference (LMM estimate)", y = "Module")+
      theme_classic()

    pdf(file = paste0(out_dir, script_ind, "Module_eigengene_LMM_across_clusters_forestplot.pdf"),
        width = 8,
        height = max(4, length(unique(fp$module)) * 0.4 + 2))
    plot(p_fp)
    dev.off()
  }



#########################################
#plot gene z-scores all DEGs ordered by module
#########################################

pl_genes = mod_gene_tab$gene

meta = bulk_data$meta

pl_mat_X = bulk_data$gene_Z_scores$clusters_combined[pl_genes,]
lims_X = 0.1*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


#Create Column Annotation Data (Top Bars)

names(meta)

col_anno_df <- as.data.frame(meta[,c("cluster_name", "APOEgroup", "CD33Group", "BrainRegion",
                                     "plaque_dens",
                                     "NeuropathologicalDiagnosis", "TREM2Variant")])
for (col1 in names(col_anno_df)){
  if (!is.numeric(col_anno_df[[col1]])){
    col_anno_df[[col1]] = factor(col_anno_df[[col1]], levels = unique(col_anno_df[[col1]]))}
}
rownames(col_anno_df) = meta$cluster_sample

col_annot = create_heatmap_annot(annot_df = col_anno_df)


#Create Gene Module Annotation Data (Side Bars)

t1 = bulk_data$wgcna$mod_gene_tab

row_anno_df <- as_tibble(t1)[,c("module")]
row_anno_df$module = factor(row_anno_df$module, levels = unique(row_anno_df$module))
rownames(row_anno_df) = t1$gene

row_annot = create_heatmap_annot(annot_df = row_anno_df[pl_genes,], annot_dim = "row", 
                                 annotation_name_side = "bottom")


#plot heatmap

ht_list <- Heatmap(
  pl_mat_X,
  name = "Z-score",
  col = colorRamp2(breaks = seq(from = lims_X[1], to = lims_X[2], length.out = 100),
                   colors = viridis(100)),
  cluster_columns = FALSE,
  column_title = "cluster_sample",
  top_annotation = col_annot,
  cluster_rows = FALSE,
  row_title = "gene/module",
  left_annotation = row_annot, 
  show_row_names = FALSE, 
  show_column_names = FALSE,
  width = 10, height = 10
)


pdf(file = paste0(out_dir,script_ind, "Gene_z_score_all_DEGs_by_module.pdf"), 
    width = 10, height = 10)
{
  draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "right")
}
dev.off()



#################################################
### GO-BP over-representation analysis of DEGs by module
#################################################

GO_list = list()
GO_results_tab = NULL

N_comps = length(bulk_data$wgcna$mod_genes)

for (i in 1:N_comps){
  
  comp = names(bulk_data$wgcna$mod_genes)[i]
  
  message("\n          *** GO analysis DEGs by module ", comp, 
          " (", i, " of ",N_comps,  ") - ", Sys.time(),"\n")
  
  ego = NULL
  
  go_genes =  bulk_data$wgcna$mod_genes[[comp]]
  
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
      t2 = cbind(module = comp, t1)
      GO_results_tab = rbind(GO_results_tab, t2)
    }
  }
}

GO_results_tab = GO_results_tab[GO_results_tab$Count>1,]

write_csv(GO_results_tab, file = paste0(out_dir,script_ind, "GO_results_by_comp.csv"))

bulk_data$GO_results$by_comp_GO_list = GO_list
bulk_data$GO_results$by_comp_GO_res = GO_results_tab

qsave(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.qs")) 



############################################################################
### visualise GO analysis (dotplot by cluster for top 10 terms by comp, GO term network plot)
############################################################################

t1 = GO_results_tab

t2 = NULL

for (mod1 in unique(t1$module)){
  t3 = t1[t1$module == mod1,]
  if (nrow(t3)>10){t3 = t3[1:10,]}
  t2 = rbind(t2, t3)
}

mod_colors = unique(mod_gene_tab$colors)

p1 = ggplot(t2, aes(x = module, y = Description, size = Count, colour = module))+
  geom_point()+
  scale_color_manual (limits = mods, values =mod_colors)+
  scale_size_continuous(limits = c(0, max(t2$Count)))+
  scale_x_discrete(limits = mods)+
  scale_y_discrete(limits = unique(t2$Description))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(title = "Top10 GO terms per module")

pdf(file = paste0(out_dir,script_ind, "GO_results_by_module_dotplot_top_terms.pdf"), 
    width = 11, height = 12)
plot(p1)
dev.off()



###plot GO term network plot

pl = list()

N_comps = length(GO_list)

for (i in 1:N_comps){
  
  comp = names(GO_list)[i]
  
  message("\n          *** GO analysis DEGs by module ", comp, 
          " (", i, " of ",N_comps,  ") - ", Sys.time(),"\n")
  
  edo = GO_list[[comp]]
  t1 = edo@result
  t1 = t1[t1$p.adjust<0.05,]

  if (nrow(t1)>2){
    tryCatch({
      edo <- pairwise_termsim(edo)
      p1 = emapplot(edo, showCategory = 100)+labs(title = comp)
      pl[[comp]] <- p1
    }, error = function(e){
      message("    skipping emapplot for ", comp, ": ", conditionMessage(e))
    })
  }
}

pdf(file = paste0(out_dir,script_ind, "GO_results_by_module_network_plot.pdf"), 
    width = 15, height = 12)
lapply(pl, function(x){x})
dev.off()



############################################################################
### Plot heatmaps of top20 GO-term-linked genes
############################################################################

###extract GO to plot and linked genes

t1 = GO_results_tab

t2 = NULL

for (mod1 in unique(t1$module)){
  t3 = t1[t1$module == mod1,]
  if (nrow(t3)>20){t3 = t3[1:20,]}
  t2 = rbind(t2, t3)
}


meta = bulk_data$meta

go_genes_list = str_split(t2$geneID, "/")
names(go_genes_list) = paste0(t2$module, "_", t2$Description)

m1 = bulk_data$gene_Z_scores$clusters_combined
pl_mat_X = m1[rownames(m1) %in% unique(unlist(go_genes_list)),]

lims_X = 0.3*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


### Create Column Annotation Data (Top Bars)

col_anno_df <- as.data.frame(meta[,c("cluster_name", "APOEgroup", "CD33Group", "BrainRegion",
                                     "plaque_dens",
                                     "NeuropathologicalDiagnosis", "TREM2Variant")])

for (col1 in names(col_anno_df)){
  if (!is.numeric(col_anno_df[[col1]])){
    col_anno_df[[col1]] = factor(col_anno_df[[col1]], levels = unique(col_anno_df[[col1]]))}
}
rownames(col_anno_df) = meta$cluster_sample

col_annot = create_heatmap_annot(annot_df = col_anno_df)


### Create Gene Module Annotation Data (Side Bars)

t1 = bulk_data$wgcna$mod_gene_tab

row_anno_df <- as_tibble(t1)[,c("module")]
row_anno_df$module = factor(row_anno_df$module, levels = unique(row_anno_df$module))
rownames(row_anno_df) = t1$gene


###plot genes for all GO terms

ht_list = list()

for (pl_set in names(go_genes_list)){
  
  pl_genes = go_genes_list[[pl_set]]
  
  
  row_annot = create_heatmap_annot(annot_df = row_anno_df[pl_genes,], annot_dim = "row", 
                                   annotation_name_side = "bottom")
  
  ht_list[[pl_set]] <- Heatmap(
    pl_mat_X[pl_genes,],
    name = "Z-score",
    col = colorRamp2(breaks = seq(from = lims_X[1], to = lims_X[2], length.out = 100),
                     colors = viridis(100)),
    cluster_columns = FALSE,
    column_title = paste0(pl_set, " - Expression Z-score"),
    top_annotation = col_annot,
    cluster_rows = FALSE,
    row_title = NULL,
    left_annotation = row_annot, 
    show_row_names = TRUE, 
    show_column_names = FALSE,
    width = 10, height = unit(length(pl_genes)*0.5, "cm")
  )
}


largest_go_set <- names(which.max(sapply(go_genes_list, length)))
pdf(NULL)
ht_drawn <- draw(ht_list[[largest_go_set]], heatmap_legend_side = "right", annotation_legend_side = "right")
go_pdf_height <- convertHeight(ComplexHeatmap:::height(ht_drawn), "inches", valueOnly = TRUE)
go_pdf_width  <- convertWidth(ComplexHeatmap:::width(ht_drawn),  "inches", valueOnly = TRUE)
dev.off()

pdf(file = paste0(out_dir, script_ind, "Gene_expression_heatmap_module_GO_genes_by_sample.pdf"),
    width = go_pdf_width, height = go_pdf_height)
for (pl_set in names(go_genes_list)){
  draw(ht_list[[pl_set]], heatmap_legend_side = "right", annotation_legend_side = "right")
}
dev.off()



############################################################################
### Plot heatmaps of GOI genes
############################################################################

m1 = bulk_data$gene_Z_scores$clusters_combined
pl_mat_X = m1[rownames(m1) %in% unique(unlist(GOI)),]

lims_X = 0.1*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


###plot genes for all GO terms

ht_list = list()

for (pl_set in names(GOI)){
  
  pl_genes = intersect(GOI[[pl_set]], rownames(pl_mat_X))
  
  ht_list[[pl_set]] <- Heatmap(
    pl_mat_X[pl_genes,],
    name = "Z-score",
    col = colorRamp2(breaks = seq(from = lims_X[1], to = lims_X[2], length.out = 100),
                     colors = viridis(100)),
    cluster_columns = FALSE,
    column_title = paste0(pl_set, " - Expression Z-score"),
    top_annotation = col_annot,
    cluster_rows = FALSE,
    row_title = NULL,
    left_annotation = NULL,
    show_row_names = TRUE,
    show_column_names = FALSE,
    width = 10, height = unit(length(pl_genes)*0.15, "cm")
  )
}

largest_goi_set <- names(which.max(sapply(names(GOI), function(x)
  length(intersect(GOI[[x]], rownames(pl_mat_X))))))
pdf(NULL)
ht_drawn <- draw(ht_list[[largest_goi_set]], heatmap_legend_side = "right", annotation_legend_side = "right")
goi_pdf_height <- convertHeight(ComplexHeatmap:::height(ht_drawn), "inches", valueOnly = TRUE)
goi_pdf_width  <- convertWidth(ComplexHeatmap:::width(ht_drawn),  "inches", valueOnly = TRUE)
dev.off()

pdf(file = paste0(out_dir, script_ind, "Gene_expression_heatmap_GOI_by_sample.pdf"),
    width = goi_pdf_width, height = goi_pdf_height)
for (pl_set in names(GOI)){
  draw(ht_list[[pl_set]], heatmap_legend_side = "right", annotation_legend_side = "right")
}
dev.off()



############################################################################
### Plot heatmaps of subtype markers
############################################################################

pl_gene_list = subtype_marker_list

m1 = bulk_data$gene_Z_scores$clusters_combined
pl_mat_X = m1[rownames(m1) %in% unique(unlist(pl_gene_list)),]

lims_X = 0.2*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))


###plot genes for all GO terms

ht_list = list()

for (pl_set in names(pl_gene_list )){
  
  pl_genes = intersect(pl_gene_list[[pl_set]], rownames(pl_mat_X))
  
  if (length(pl_genes)>1){
    
    ht_list[[pl_set]] <- Heatmap(
      pl_mat_X[pl_genes,],
      name = "Z-score",
      col = colorRamp2(breaks = seq(from = lims_X[1], to = lims_X[2], length.out = 100),
                       colors = viridis(100)),
      cluster_columns = FALSE,
      column_title = paste0(pl_set, " - Expression Z-score"),
      top_annotation = col_annot,
      cluster_rows = FALSE,
      row_title = NULL,
      left_annotation = NULL, 
      show_row_names = TRUE, 
      show_column_names = FALSE,
      width = 10, height = unit(length(pl_genes)*0.5, "cm")
    )
  }
}


largest_marker_set <- names(which.max(sapply(names(ht_list), function(x)
  nrow(ht_list[[x]]@matrix))))
pdf(NULL)
ht_drawn <- draw(ht_list[[largest_marker_set]], heatmap_legend_side = "right", annotation_legend_side = "right")
marker_pdf_height <- convertHeight(ComplexHeatmap:::height(ht_drawn), "inches", valueOnly = TRUE)
marker_pdf_width  <- convertWidth(ComplexHeatmap:::width(ht_drawn),  "inches", valueOnly = TRUE)
dev.off()

pdf(file = paste0(out_dir, script_ind, "Gene_expression_heatmap_subtype_marker_list_by_sample.pdf"),
    width = marker_pdf_width, height = marker_pdf_height)
for (pl_set in names(ht_list)){
  draw(ht_list[[pl_set]], heatmap_legend_side = "right", annotation_legend_side = "right")
}
dev.off()




######################################################################
###plot module eigengene vs plaque_dens scatterplot by TREM2Variant and module
######################################################################

m1 = bulk_data$wgcna$mod_eigengene_mat

meta = bulk_data$meta

gr = unique(meta$group)

pl = list()
comb_plot_tab = NULL

for (mod1 in rownames(m1)){
  
  for (cl in cluster_names){
    
    meta_cl = meta[meta$cluster_name == cl,]
    
    t1 = cbind(meta_cl, module = mod1, eigengene = m1[mod1, meta_cl$cluster_sample])
    
    comb_plot_tab = rbind(comb_plot_tab,t1)
    
  }
}


###plot combined by cluster vs module

t1 = comb_plot_tab
t1$cluster_name = factor(t1$cluster_name, levels = cluster_names)
t1$module = factor(t1$module, levels = mods)

p1 = ggplot(t1)+
  geom_smooth(aes(x = plaque_dens, y = eigengene, color = TREM2Variant, group = TREM2Variant), 
              method = "lm")+
  geom_point(aes(x = plaque_dens, y = eigengene, fill = TREM2Variant, group = TREM2Variant), 
             shape = 21, size = 2, stroke = 0.3, color = "grey40")+
  facet_grid(rows = vars(module), cols = vars(cluster_name), scales = "free")+
  scale_fill_manual(limits = unique(t1$TREM2Variant), values = pal(unique(t1$TREM2Variant)))+
  scale_color_manual(limits = unique(t1$TREM2Variant), values = pal(unique(t1$TREM2Variant)))+
  theme_bw()+
  labs(title = paste0("Module eigengene by TREM2 variant by cluster vs module"))

pdf(file = paste0(out_dir,script_ind, "Module_eigengene_vs_plaque_dens_scatterplot_by_TREM2Variant_cluster_combined.pdf"),
    width = length(cluster_names) * 1.5 + 2, height = length(mods) * 1.5 + 2)
{
  plot(p1)
}
dev.off()





#get info on version of R, used packages etc
sessionInfo()


message("\n\n##########################################################################\n",
        "# Completed LD_F03a1 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")
