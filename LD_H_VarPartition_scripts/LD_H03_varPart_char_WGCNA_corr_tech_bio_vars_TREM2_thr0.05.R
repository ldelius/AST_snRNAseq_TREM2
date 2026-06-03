message("\n\n##########################################################################\n",
        "# Start LD_H03: variable-gene WGCNA characterisation", Sys.time(),
        "\n##########################################################################\n",
        "\n   WGCNA on the corrected vst matrix (TREM2-driven variable genes, thr 0.05)",
        "\n   modules -> activity -> covariate/group association -> GO characterisation",
        "\n##########################################################################\n\n")

# what this script does (adapted from Michael's G04a31 microglia WGCNA script):
# - builds on LD_H02 output (LD_H02_v02_bulk_data.qs, hybrid-corrected)
# - v03: WGCNA soft-power step aligned to my F03a1/F03a2 scripts (was v02)
# - runs WGCNA on the TREM2-driven variable genes (thr 0.05) using the
#   covariate-CORRECTED vst matrix from H02 (vst_mat_corr)
# - extracts modules, module activity scores, module-covariate associations,
#   module-vs-group plots, and per-module + combined GO-BP characterisation
#
# DATA-WIRING / NAMING changes only - WGCNA logic + parameters kept as Michael's,
# EXCEPT the soft-thresholding power, which is chosen data-driven from MY fit:
#  - paths/naming -> mine; script_ind "LD_H03_v03_"
#  - input -> LD_H02_v02_bulk_data.qs (not G03a3_bulk_data.qs)
#  - marker subtype filter -> "Astrocyte_subtypes" (Michael: "Microglia_subtypes")
#  - GOI / TF / marker tables still read from data_TREM2_michael/A_input (my convention,
#    same as my LD_F03a2 script)
#  - soft power: chosen EXACTLY as in my F03a1/F03a2 scripts (unsigned
#    pickSoftThreshold; first power with signed R^2 > 0.8, floored at 6) so the
#    WGCNA construction matches F for a clean F-vs-H comparison; see below
#  - added runtime confirmation block (gene-set name/size, vst_mat_corr slot)

# Open packages necessary for analysis.
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
library(variancePartition)
library(BiocParallel)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify output directory
out_dir = paste0(main_dir,"LD_H_VarPartition_output/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

#specify script/output index as prefix for file names
script_ind = "LD_H03_v03_"


### load dataset (LD_H02 output: var_genes + corrected vst matrix + z-scores)
bulk_data = qread(file = paste0(out_dir,"LD_H02_v02_bulk_data.qs"))


#######################################
# confirm required H02 slots + gene set (points 3 & 4 of the brief)
#######################################

message("\n   var_genes sets available in H02 output:")
print(data.frame(gene_set = names(bulk_data$varPart_analysis$var_genes),
                 N_genes  = lengths(bulk_data$varPart_analysis$var_genes),
                 row.names = NULL))

# corrected vst matrix must exist (WGCNA input)
if (is.null(bulk_data$vst_mat_corr)){
  stop("   !ERROR: bulk_data$vst_mat_corr not found - rerun LD_H02 first")
}

# the exact gene set Michael used
gene_set_name = "TREM2Variant_thr_0.05"
if (!gene_set_name %in% names(bulk_data$varPart_analysis$var_genes)){
  stop("   !ERROR: '", gene_set_name, "' not in var_genes; available: ",
       paste(names(bulk_data$varPart_analysis$var_genes), collapse = ", "))
}

n_var_genes = length(bulk_data$varPart_analysis$var_genes[[gene_set_name]])
message("\n   Using gene set '", gene_set_name, "' with ", n_var_genes, " genes\n")

# WGCNA needs enough genes to form modules (minModuleSize = 30 below)
if (n_var_genes < 30){
  stop("   !ERROR: only ", n_var_genes,
       " genes - too few for WGCNA (minModuleSize = 30). Lower the threshold in H02.")
}
if (n_var_genes < 100){
  message("   !! FLAG: only ", n_var_genes,
          " genes (<100) - WGCNA modules may be few/small; interpret with caution\n")
}


### specify WGCNA input gene matrix, re-order to focus on TREM2 variants

keep_genes = bulk_data$varPart_analysis$var_genes[[gene_set_name]]

t1 = bulk_data$meta
t2 = t1[order(t1$TREM2Variant),]
t2$group = t2$TREM2Variant
bulk_data$meta = t2

bulk_data$gene_Z_scores_uncorr = bulk_data$gene_Z_scores_uncorr[,rownames(t2)]
t3 = t2[rownames(t2) %in% colnames(bulk_data$gene_Z_scores_corr),]
bulk_data$gene_Z_scores_corr = bulk_data$gene_Z_scores_corr[,rownames(t3)]

input_mat = t(bulk_data$vst_mat_corr[keep_genes, rownames(t3)])
bulk_data$wgcna$input_mat = input_mat



### get gene of interest gene sets

GOI = list()
t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/Transcription Factors hg19 - Fantom5_21-12-21.csv"))
GOI$TF = t1$Symbol

t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/GOI_sets_251020.csv"))

for (goi_set in unique(t1$gene_set)){
  GOI[[goi_set]] = t1$gene[t1$gene_set == goi_set]
}

### get marker gene sets (astrocyte subtypes; Michael used Microglia_subtypes)

t1 = read_csv(paste0(main_dir,"data_TREM2_michael/A_input/cell_type_markers_241219_w_astr_subtype_markers.csv"))
t2 = t1[t1$level == "Astrocyte_subtypes",]

subtype_marker_list = list()

for (subtype in unique(t2$subtype)){
  subtype_marker_list[[subtype]] = t2$gene[t2$subtype == subtype]
}



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



### create heatmap annotations for Complex Heatmap

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

    } else if (is.factor(v1)) {
      v2 = pal(levels(na.omit(v1)))
      names(v2) = levels(na.omit(v1))
      annot_colors = c(annot_colors, list(v2))
    } else {
      stop("   !ERROR: categorical annotations need to be factors")
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
# calculate WGCNA soft-thresholding power (ALIGNED TO MY F WGCNA SCRIPTS)
#######################################
# Soft-power selection is matched EXACTLY to my F03a1/F03a2 WGCNA, so the
# corrected (H) and less-corrected (F) analyses differ ONLY by the correction
# level - not by the network-construction method:
#  - pickSoftThreshold uses the DEFAULT (unsigned) networkType, as in F
#  - power = first power with signed scale-free R^2 > 0.8, floored at 6 (F's rule)
#  - blockwiseModules below already matches F (networkType = "signed", deepSplit 2,
#    minModuleSize 30, maxBlockSize 20000, mergeCutHeight 0.25, numericLabels)
# (Earlier H03 versions used a signed pickSoftThreshold + a sample-size floor of 12;
#  that was changed to match F for a clean F-vs-H comparison.)

message("\n\n          *** Calculate soft-thresholding power... ", Sys.time(),"\n\n")

allowWGCNAThreads()          # allow multi-threading (optional)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function on MY corrected matrix
sft = pickSoftThreshold(
  input_mat,                 # <= Input data (my vst_mat_corr, TREM2 thr0.05 genes)
  powerVector = powers,      # networkType left at DEFAULT (unsigned) to match my F scripts
  verbose = 5
)

### plot thresholding power graphs

t1 = as_tibble(sft$fitIndices)
t1$signed_R2 = -sign(t1$slope)*t1$SFT.R.sq

### power choice ALIGNED TO MY F WGCNA scripts (F03a1/F03a2) for comparability:
### first power reaching signed scale-free R^2 > 0.8, floored at 6.
n_samples = nrow(input_mat)
wgcna_power = t1$Power[t1$signed_R2 > 0.8][1]
if (wgcna_power < 6){ wgcna_power = 6 }

bulk_data$wgcna$power_table = t1
bulk_data$wgcna$power = wgcna_power

message("\n   >>> Chosen soft-thresholding power = ", wgcna_power,
        " (", n_samples, " samples) <<<\n")
print(t1[, c("Power", "SFT.R.sq", "slope", "signed_R2", "mean.k.")])

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
# blockwiseModules parameters kept EXACTLY as Michael's G04a31 (for comparability)

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



### extract module activity score (mean module gene z-score by cluster_sample)

z_mat = bulk_data$gene_Z_scores_corr

m2 = matrix(nrow = length(mods), ncol = ncol(z_mat),
            dimnames = list(mods, colnames(z_mat)))

for (mod1 in mods){

  mod_genes = mod_gene_tab$gene[mod_gene_tab$module == mod1]

  for (cs in colnames(z_mat)){
    m2[mod1,cs] = mean(z_mat[mod_genes, cs])
  }
}

bulk_data$wgcna$mod_activity_mat = m2

qsave(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.qs"))



#######################################################
###plot module activity heatmap by sample and module
#######################################################

prim_vars = bulk_data$varPart_analysis$primary_var

pl_mat_X = bulk_data$wgcna$mod_activity_mat

lims_X = 0.1*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pl_meta = bulk_data$meta[colnames(pl_mat_X),]

pl_mat_X = pl_mat_X[, rownames(pl_meta)]

#Create Column Annotation Data (Top Bars)
annot_vars = c(bulk_data$varPart_analysis$candidate_covars, prim_vars)

col_anno_df <- as.data.frame(pl_meta[,annot_vars])

for (col1 in names(col_anno_df)){
  if (!(is.numeric(col_anno_df[[col1]]) | is.factor(col_anno_df[[col1]]))){
    col_anno_df[[col1]] = factor(col_anno_df[[col1]], levels = unique(col_anno_df[[col1]]))}
}
rownames(col_anno_df) = pl_meta$cluster_sample

col_annot = create_heatmap_annot(annot_df = col_anno_df)


#plot heatmap

ht_list <- Heatmap(
  pl_mat_X,
  name = "Expr \nZ-score",
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


pdf(file = paste0(out_dir,script_ind, "Module_activity_heatmap_by_sample.pdf"),
    width = 10, height = 10)
{
  draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
}
dev.off()



######################################################################
### calculate corr plot sample module activity vs covars
######################################################################

message("\n\n   *Generate module activity vs covariate corr plot \n")

mod_act_mat = bulk_data$wgcna$mod_activity_mat
m1 = t(mod_act_mat)

t1 = bulk_data$meta[rownames(m1),]

identical(rownames(t1), rownames(m1))

meta_mod_tab = cbind(t1, m1)


### check covariates for correlation/colinearity

#add modules to covariate analysis formula
form_mod = paste0(bulk_data$varPart_analysis$form_full, " + ", paste0(mods, collapse = " + "))

cor_mat <- canCorPairs(as.formula(form_mod), meta_mod_tab)

vars = c(bulk_data$varPart_analysis$candidate_covars, prim_vars)

m2 = cor_mat[vars, mods]

t1 = as_tibble(cbind(var = rownames(m2), m2))
write_csv(t1, file = paste0(out_dir, script_ind, "module_covar_CCA.csv"))


pdf(file = paste0(out_dir,script_ind, "module_covar_CCA.pdf"),
    width = 7, height = 7)
{
  p1 = pheatmap::pheatmap(m2,
                          cluster_rows = TRUE, cluster_cols = TRUE,
                          show_rownames = TRUE, show_colnames = TRUE,
                          color = colorRampPalette(c("white", "red"))(250),
                          breaks = seq(0, 1, length.out = 251),
                          border_color = NA, cellwidth = 10, cellheight = 10,
                          fontsize = 10,
                          main = "Covariates vs modules (canonical correlation annalysis")
}
dev.off()



######################################################################
### Fit variance partition model for module activity
######################################################################

### Fit variance partition model (module activity)
### NB: variancePartition::fitExtractVarPartModel is BROKEN in this environment
### (variancePartition 1.32.5 imports lme4::findbars, which lme4 2.0 removed/moved
### to reformulas). On the module matrix it fails for some modules and SILENTLY
### DROPS them (we saw only 4 of 8 returned) - the same incompatibility that hit
### H02. We therefore compute the per-module variance fractions directly: fit each
### module with lme4::lmer (its own parser works fine) and decompose with
### variancePartition::calcVarPart() on the FITTED model (post-fit, so it does not
### re-trigger findbars). All modules are reported; any genuine failure is flagged.

meta = bulk_data$meta[colnames(mod_act_mat),]
form_full = bulk_data$varPart_analysis$form_full
form_mod_fit = as.formula(paste("y", form_full))   # y ~ <covariates>

vp_list = list()
for (mod1 in rownames(mod_act_mat)){
  d = meta
  d$y = mod_act_mat[mod1, rownames(d)]
  vp_list[[mod1]] = tryCatch({
    fit = suppressWarnings(lme4::lmer(form_mod_fit, data = d, REML = TRUE))
    calcVarPart(fit)
  }, error = function(e){
    message("   !! module varPart failed for ", mod1, ": ", conditionMessage(e))
    NULL
  })
}

vp_list = vp_list[!sapply(vp_list, is.null)]
varPart = as.data.frame(do.call(rbind, vp_list))
rownames(varPart) = names(vp_list)

message("\n   module varPart computed for ", nrow(varPart), " of ",
        nrow(mod_act_mat), " modules\n")

### Summaries (keep module label as a column so it survives write_csv)

write_csv(cbind(module = rownames(varPart), as_tibble(varPart)),
          file = paste0(out_dir, script_ind, "varPart_variance_expl_by_module.csv"))

pdf(file = paste0(out_dir,script_ind, "varPart_variance_expl_by_module.pdf"),
    width = 10, height = 10)
{
  p1 = pheatmap::pheatmap(t(varPart[,colnames(varPart) != "Residuals"]),
                          cluster_rows = TRUE, cluster_cols = TRUE,
                          show_rownames = TRUE, show_colnames = TRUE,
                          color = colorRampPalette(c("white", "red"))(250),
                          breaks = seq(0, max(varPart[,colnames(varPart) != "Residuals"]), length.out = 251),
                          border_color = NA, cellwidth = 10, cellheight = 10,
                          fontsize = 10,
                          main = "Covariates fraction variance explained vs modules")
}
dev.off()



######################################################################
###plot module activity vs group scatterplot by sample and module
######################################################################

m1 = bulk_data$wgcna$mod_activity_mat

meta = bulk_data$meta[colnames(m1),]

gr = unique(meta$group)

pl = list()

for (mod1 in rownames(m1)){

  t1 = cbind(meta, module = mod1, activity = m1[mod1, rownames(meta)])

  t1_sum <- t1 %>%
    group_by(group) %>%
    summarise(
      mean_activity = mean(activity),
      sd_activity = sd(activity),
      ymin_activity = mean_activity - sd_activity,
      ymax_activity = mean_activity + sd_activity
    )

  pl[[paste0(mod1)]] = ggplot()+
    geom_point(data = t1, aes(x = group, y = activity, fill = group), shape = 21,
               position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.7),
               size = 0.5, stroke = 0.05, color = "grey50")+
    geom_errorbar(data = t1_sum, aes(x = group, ymin = ymin_activity, ymax = ymax_activity),
                  width = 0.2, color = "black") +
    geom_point(data = t1_sum, aes(x = group, y = mean_activity, fill = group), size = 2,
               shape = 21, color = "black") +
    geom_hline(yintercept = 0)+
    scale_x_discrete(limits = gr)+
    scale_color_manual(limits = gr, values = pal(gr))+
    scale_fill_manual(limits = gr, values = pal(gr))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    labs(title = paste0(mod1, " - Module activity by group"))

}



pdf(file = paste0(out_dir,script_ind, "Module_activity_vs_group_scatterplot.pdf"),
    width = 3, height = 3)
{
  lapply(pl, function(x){x})
}
dev.off()



######################################################################
###plot module activity vs group scatterplot by sample and module (combined)
######################################################################

m1 = bulk_data$wgcna$mod_activity_mat

meta = bulk_data$meta[colnames(m1),]

mods = unique(bulk_data$wgcna$mod_gene_tab$module)

gr = unique(meta$group)

pl_tab = NULL

for (mod1 in rownames(m1)){

  t1 = cbind(meta, module = mod1, activity = m1[mod1, rownames(meta)])

  t1_sum <- t1 %>%
    group_by(group) %>%
    summarise(
      mean_activity = mean(activity),
      sd_activity = sd(activity)
    )

  t1_sum$ymin_activity = t1_sum$mean_activity - t1_sum$sd_activity
  t1_sum$ymax_activity = t1_sum$mean_activity + t1_sum$sd_activity


  t1$mean_activity = t1_sum$mean_activity[match(t1$group, t1_sum$group)]
  t1$ymin_activity = t1_sum$ymin_activity[match(t1$group, t1_sum$group)]
  t1$ymax_activity = t1_sum$ymax_activity[match(t1$group, t1_sum$group)]

  pl_tab = rbind(pl_tab, t1)

}

###plot combined by cluster vs module

t1 = pl_tab
t1$module = factor(t1$module, levels = mods)

p1 = ggplot()+
  geom_point(data = t1, aes(x = group, y = activity, fill = group), shape = 21,
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.7),
             size = 0.3, stroke = 0.02, color = "grey50")+
  geom_errorbar(data = t1, aes(x = group, ymin = ymin_activity, ymax = ymax_activity),
                width = 0.2, color = "black") +
  geom_point(data = t1, aes(x = group, y = mean_activity, fill = group), size = 2,
             shape = 21, color = "black") +
  geom_hline(yintercept = 0)+
  facet_grid(rows = vars(module), scales = "free")+
  scale_x_discrete(limits = gr)+
  scale_color_manual(limits = gr, values = pal(gr))+
  scale_fill_manual(limits = gr, values = pal(gr))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(title = paste0("Module activity by group"))


pdf(file = paste0(out_dir,script_ind, "Module_activity_vs_group_scatterplot_by_module_combined.pdf"),
    width = 3, height = 6)
{
  plot(p1)
}
dev.off()



#########################################
#plot gene z-scores all DEGs ordered by module (uncorr)
#########################################

pl_genes = mod_gene_tab$gene

pl_mat_X = bulk_data$gene_Z_scores_uncorr

lims_X = 0.1*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pl_meta = bulk_data$meta[colnames(pl_mat_X), ]

pl_mat_X = pl_mat_X[pl_genes, rownames(pl_meta)]

#Create Column Annotation Data (Top Bars)
annot_vars = c(bulk_data$varPart_analysis$candidate_covars, prim_vars)

col_anno_df <- as.data.frame(pl_meta[,c(annot_vars)])

for (col1 in names(col_anno_df)){
    if (!(is.numeric(col_anno_df[[col1]]) | is.factor(col_anno_df[[col1]]))){
      col_anno_df[[col1]] = factor(col_anno_df[[col1]], levels = unique(col_anno_df[[col1]]))}
  }
rownames(col_anno_df) = pl_meta$cluster_sample

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


pdf(file = paste0(out_dir,script_ind, "Gene_z_score_uncorr_all_var_genes_by_module.pdf"),
    width = 10, height = 10)
{
  draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "right")
}
dev.off()



#########################################
#plot gene z-scores all DEGs ordered by module (corr)
#########################################

pl_genes = mod_gene_tab$gene

pl_mat_X = bulk_data$gene_Z_scores_corr

lims_X = 0.1*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pl_meta = bulk_data$meta[colnames(pl_mat_X), ]

pl_mat_X = pl_mat_X[pl_genes, rownames(pl_meta)]

#Create Column Annotation Data (Top Bars)
annot_vars = c(bulk_data$varPart_analysis$candidate_covars, prim_vars)

col_anno_df <- as.data.frame(pl_meta[,annot_vars])

for (col1 in names(col_anno_df)){
    if (!(is.numeric(col_anno_df[[col1]]) | is.factor(col_anno_df[[col1]]))){
      col_anno_df[[col1]] = factor(col_anno_df[[col1]], levels = unique(col_anno_df[[col1]]))}
  }
rownames(col_anno_df) = pl_meta$cluster_sample

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


pdf(file = paste0(out_dir,script_ind, "Gene_z_score_corr_all_var_genes_by_module.pdf"),
    width = 10, height = 10)
{
  draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "right")
}
dev.off()



#################################################
### GO-BP over-representation analysis of all DEGs combined => enrichment by module
#################################################

message("\n\n          *** GO analysis DEGs combined ... ", Sys.time(),"\n\n")

### extract GO terms

comp_genes = unique(unlist(bulk_data$wgcna$mod_genes))

ego = enrichGO(gene         = comp_genes,
               OrgDb         = org.Hs.eg.db,
               keyType       = 'SYMBOL',
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05)


GO_results_tab = ego@result[ego@result$p.adjust<=0.05,]

### identify overlap of GO genes with module genes

mod_gene_list = bulk_data$wgcna$mod_genes

v1 = GO_results_tab$geneID
l1 = str_split(v1, "/")
names(l1) = GO_results_tab$Description

GO_gene_list = l1

go_mod_mat = matrix(nrow = length(GO_gene_list), ncol = length(mod_gene_list),
                    dimnames = list(names(GO_gene_list), names(mod_gene_list)))

for (mod1 in names(mod_gene_list)){

  for (go in names(GO_gene_list)){
    mod1_genes = mod_gene_list[[mod1]]
    go_genes = GO_gene_list[[go]]
    genes_overlap = intersect(mod1_genes, go_genes)
    go_mod_mat[go,mod1] = length(genes_overlap)
  }
}


### save overlap of GO genes with module genes

t1 = cbind(GO_results_tab, go_mod_mat)
bulk_data$wgcna$GO_results_comb$GO_results_tab = t1
write_csv(t1, paste0(out_dir, script_ind, "GO_results_combined_with_N_genes_by_mod.csv"))


### plot gene overlap module vs GO term (number and fraction of all enriched genes; for all or top 20 terms)

go_mod_heatmap = function(m1, main = ""){
  pheatmap::pheatmap(m1, show_rownames=TRUE, show_colnames = TRUE,
                     cluster_rows = TRUE, cluster_cols = FALSE,
                     clustering_distance_rows = "euclidean",
                     clustering_method = "ward.D2",
                     treeheight_row = 10, treeheight_col = 10,
                     color = colorRampPalette(c("white", "blue"))(250),
                     breaks = seq(0, max(m1), length.out = 251),
                     border_color = NA, fontsize = 10,
                     cellwidth = 10, cellheight = 10,
                     main = main
  )
}


pdf(file = paste0(out_dir,script_ind, "GO_results_combined_with_N_genes_by_mod.pdf"),
    width = 18, height = 30)
{
  m1 = go_mod_mat

  go_mod_heatmap(m1, main = "N genes - GO terms vs modules")

  m2 = go_mod_mat/lengths(mod_gene_list)

  go_mod_heatmap(m2, main = "Fract module genes in GO terms")

  m1 = go_mod_mat
  if(nrow(m1)>50){m1 = m1[1:50,]}

  go_mod_heatmap(m1, main = "N genes - Top50 GO terms vs modules")

  m2 = go_mod_mat/lengths(mod_gene_list)
  if(nrow(m2)>50){m2 = m2[1:50,]}

  go_mod_heatmap(m2, main = "Fract module genes in Top50 GO terms")

  m1 = go_mod_mat
  if(nrow(m1)>20){m1 = m1[1:20,]}

  go_mod_heatmap(m1, main = "N genes - Top20 GO terms vs modules")

  m2 = go_mod_mat/lengths(mod_gene_list)
  if(nrow(m2)>20){m2 = m2[1:20,]}

  go_mod_heatmap(m2, main = "Fract module genes in Top20 GO terms")

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

bulk_data$wgcna$GO_results_by_mod$GO_list = GO_list
bulk_data$wgcna$GO_results_by_mod$GO_res = GO_results_tab

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
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        axis.text.y  = element_text(size = 12),
        axis.title   = element_text(size = 14),
        plot.title   = element_text(size = 15),
        legend.text  = element_text(size = 11),
        legend.title = element_text(size = 12))+
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
  t1 = t1[t1$p.adjust<0.01,]

  if (nrow(t1)>2){
    edo <- pairwise_termsim(edo)
    p1 = emapplot(edo, showCategory = 100)+labs(title = comp)
    pl[[comp]] <- p1
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

GO_results_tab_top20 = t2


#extract pl matrix and annotations

pl_genes = mod_gene_tab$gene

pl_mat_X = bulk_data$gene_Z_scores_uncorr

lims_X = 0.1*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pl_meta = bulk_data$meta[colnames(pl_mat_X),]

pl_mat_X = pl_mat_X[pl_genes, rownames(pl_meta)]

#Create Column Annotation Data (Top Bars)
v1 = bulk_data$varPart_analysis$candidate_covars
annot_vars = c(bulk_data$varPart_analysis$candidate_covars, prim_vars)

col_anno_df <- as.data.frame(pl_meta[,annot_vars])

for (col1 in names(col_anno_df)){
    if (!(is.numeric(col_anno_df[[col1]]) | is.factor(col_anno_df[[col1]]))){
      col_anno_df[[col1]] = factor(col_anno_df[[col1]], levels = unique(col_anno_df[[col1]]))}
  }
rownames(col_anno_df) = pl_meta$cluster_sample

col_annot = create_heatmap_annot(annot_df = col_anno_df)


###extract genes to plot
t2 = GO_results_tab_top20

go_genes_list = str_split(t2$geneID, "/")
names(go_genes_list) = paste0(t2$module, "_", t2$Description)


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
    name = "Expr \nZ-score",
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


pdf(file = paste0(out_dir,script_ind, "GO_results_Gene_expression_heatmaps_Top20_terms.pdf"),
    width = 10, height = 10)
{
  for (pl_set in names(go_genes_list)){
    draw(ht_list[[pl_set]], heatmap_legend_side = "right", annotation_legend_side = "right")
  }
}
dev.off()



############################################################################
### Plot heatmaps of GOI genes
############################################################################

pl_mat_X = bulk_data$gene_Z_scores_uncorr

lims_X = 0.1*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pl_meta = bulk_data$meta[colnames(pl_mat_X),]

###plot genes for all GO terms

ht_list = list()

for (pl_set in names(GOI)){

  pl_genes = intersect(GOI[[pl_set]], rownames(pl_mat_X))

  ht_list[[pl_set]] <- Heatmap(
    pl_mat_X[pl_genes,],
    name = "Expr \nZ-score",
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


pdf(file = paste0(out_dir,script_ind, "Gene_expression_heatmap_GOI_by_sample.pdf"),
    width = 10, height = 20)
{
  for (pl_set in names(GOI)){
    draw(ht_list[[pl_set]], heatmap_legend_side = "right", annotation_legend_side = "right")
  }
}
dev.off()



############################################################################
### Plot heatmaps of subtype markers
############################################################################

pl_gene_list = subtype_marker_list

pl_mat_X = bulk_data$gene_Z_scores_uncorr

lims_X = 0.1*c(-max(abs(pl_mat_X)), max(abs(pl_mat_X)))

pl_meta = bulk_data$meta[colnames(pl_mat_X),]


###plot genes for all GO terms

ht_list = list()

for (pl_set in names(pl_gene_list )){

  pl_genes = intersect(pl_gene_list[[pl_set]], rownames(pl_mat_X))

  if (length(pl_genes)>1){

    ht_list[[pl_set]] <- Heatmap(
      pl_mat_X[pl_genes,],
      name = "Expr \nZ-score",
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


pdf(file = paste0(out_dir,script_ind, "Gene_expression_heatmap_subtype_marker_list_by_sample.pdf"),
    width = 10, height = 10)
{
  for (pl_set in names(pl_gene_list )){
    draw(ht_list[[pl_set]], heatmap_legend_side = "right", annotation_legend_side = "right")
  }
}
dev.off()






#get info on version of R, used packages etc
sessionInfo()


message("\n\n##########################################################################\n",
        "# Completed LD_H03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")
