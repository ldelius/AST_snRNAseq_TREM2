message("\n\n##########################################################################\n",
        "# Start LD_H01: Pseudobulk variance partitioning ", Sys.time(),
        "\n##########################################################################\n\n")

# The final model excludes N_cells_scaled because it is collinear with
# cluster_name, and retains RNA_counts_scaled because pathology remains visible
# in the collinearity diagnostic. This affects only the H01 diagnostic.
library(qs)
library(tidyverse)
library(Seurat)
library(DESeq2)
library(colorRamps)
library(viridis)
library(pheatmap)
library(ggrepel)
library(variancePartition)
library(BiocParallel)


### define directories and script index

main_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/"
setwd(main_dir)

#specify output directory
out_dir = paste0(main_dir,"LD_H_VarPartition_output/")
if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}

#specify script/output index as prefix for file names
script_ind = "LD_H01_v02_"


### load dataset

bulk_data = qread(file = paste0(main_dir, "LD_F_DESeq_pseudobulk_WGCNA/LD_F01_v02_bulk_data.qs"))


###select cluster_samples to keep (match Michael: drop Controls, drop NA BrainRegion/APOEgroup)

t1 = as.data.frame(bulk_data$meta)
rownames(t1) = t1$cluster_sample
names(t1)

t1 = t1[t1$NeuropathologicalDiagnosis != "Control" & !is.na(t1$BrainRegion) &
          !is.na(t1$APOEgroup),]

# #keep clusters with at least 1 samples per cluster and genotype
# t2 = t1 %>% group_by(cluster_name, APOEgroup, TREM2Variant) %>% summarise(N = n())
# t3 = t2[t2$N>0,] %>% group_by(cluster_name) %>% summarise(N = n())
# t1 = t1[t1$cluster_name %in% t3$cluster_name[t3$N>1], ]

bulk_data$meta = t1

### adapt metadata for modelling, order by APOE vs TREM2

t1 = bulk_data$meta

t1$sample = factor(t1$sample, levels = unique(t1$sample))
t1$cluster_name = factor(t1$cluster_name, levels = unique(t1$cluster_name))
t1$Braak = factor(t1$Braak, levels = unique(t1$Braak))
t1$TREM2Variant = factor(t1$TREM2Variant, levels = c("CV", "R62H", "R47H"))
t1$APOEgroup = factor(t1$APOEgroup, levels = c("APOE4-neg", "APOE4-pos"))
t1$CD33Group = factor(t1$CD33Group, levels = c("CV","CD33var"))
t1$BrainRegion = factor(t1$BrainRegion, levels = c("MTG","SSC"))
t1$cohort = factor(t1$cohort, levels = unique(t1$cohort))

t1$group = paste0(t1$APOEgroup, "_", t1$TREM2Variant)
t1$group = factor(t1$group, levels = unique(t1$group))

#coerce numeric covariates (AST metadata stores some as character) before scaling
t1$PostMortemInterval = as.numeric(t1$PostMortemInterval)
t1$plaque_dens        = as.numeric(t1$plaque_dens)
t1$Age                = as.numeric(t1$Age)
t1$pctAT8PositiveArea  = as.numeric(t1$pctAT8PositiveArea)
t1$pctPHF1PositiveArea = as.numeric(t1$pctPHF1PositiveArea)
t1$pct4G8PositiveArea  = as.numeric(t1$pct4G8PositiveArea)

t1$PMI_scaled = scale(t1$PostMortemInterval, center = FALSE)[,1]
t1$plaque_dens_scaled = scale(t1$plaque_dens, center = FALSE)[,1]
t1$Age_scaled = scale(t1$Age, center = FALSE)[,1]
t1$RNA_counts_scaled =  scale(apply(bulk_data$counts[,t1$cluster_sample], 2, sum),
                              center = FALSE)[,1]

t1 = t1[order(t1$APOEgroup, t1$TREM2Variant, t1$cluster_name),]

bulk_data$meta = t1

m1 = bulk_data$counts
m2 = m1[,bulk_data$meta$cluster_sample]
bulk_data$counts = m2

identical(t1$cluster_sample, colnames(m2))

# REQUIRED: define your primary variable(s) of interest
primary_var <- c("TREM2Variant")

# Candidate covariates to consider
names(bulk_data$meta)

candidate_covars <- c("cluster_name", "cohort", "PMI_scaled","RNA_counts_scaled",
                      "BrainRegion","Age_scaled",
                      "Sex", "CD33Group", "APOEgroup",
                       "Braak", "plaque_dens_scaled",
                      "pctAT8PositiveArea" , "pctPHF1PositiveArea", "pct4G8PositiveArea"
                       )


### drop rows with NA in any modelled covariate (needed for AST data: CD33Group/Braak/PMI
### contain NAs that Michael's data did not; fitExtractVarPartModel cannot fit NA terms)

vars_model = c(candidate_covars, primary_var)

t1 = bulk_data$meta
keep_rows = complete.cases(t1[, vars_model])
message("\n   Dropping ", sum(!keep_rows), " of ", nrow(t1),
        " pseudobulks with NA in a modelled covariate\n")
t1 = t1[keep_rows,]

bulk_data$meta = t1
bulk_data$counts = bulk_data$counts[, t1$cluster_sample]

identical(t1$cluster_sample, colnames(bulk_data$counts))




###########################################################
# Preparing dataset for analysis
###########################################################

message("\n\n          *** Preparing dataset for variancePartition analysis... ", Sys.time(),"\n\n")

t1 = bulk_data$meta

#remove clusters present in <4 samples
t2 = t1 %>% dplyr::group_by(cluster_name) %>% dplyr::summarise(N_samples = n())

comp_clusters = unique(t2$cluster_name[t2$N_samples>3])

comp_meta = t1[t1$cluster_name %in% comp_clusters,]

bulk_data$meta = comp_meta

#keep genes with counts >0 in >20% of cluster_samples

m1 = bulk_data$counts

identical(comp_meta$cluster_sample, colnames(m1))

m2 = m1
m2[m2>0] = 1
v2 = rowSums(m2)

keep_genes = names(v2)[v2 > ncol(m2)*0.2]

comp_counts = m1[keep_genes, comp_meta$cluster_sample]

identical(comp_meta$cluster_sample, colnames(comp_counts))

bulk_data$counts = comp_counts



####################################################################################
# Test impact of covariates with variancePartition (categorical as random effects)
####################################################################################
#  as recommended: https://www.bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/variancePartition.html)

counts  <- bulk_data$counts
meta <- as.data.frame(bulk_data$meta)

dds = DESeqDataSetFromMatrix(counts, colData = meta,
                             design = ~1)
bulk_data$deseq_dataset_groups_clusters_combined = dds

vsd = assay(vst(dds))
bulk_data$vst_mat_uncorr = vsd

# Convert character columns to factor
for (nm in colnames(meta)) {
  if (!(is.numeric(meta[[nm]]) | is.factor(meta[[nm]])) ){meta[[nm]] <- factor(meta[[nm]])}
}

bulk_data$meta = meta

### Generate formula for candidate covariates

form_full = "~ "

vars_full = c(candidate_covars, primary_var)

for (var1 in vars_full){

  if (!is.numeric(meta[[var1]]) ){
    form_full = paste0(form_full, "(1 | ",var1,")")
    } else {
      form_full = paste0(form_full, var1)
    }
  if (var1 != last(vars_full)){form_full = paste0(form_full, " + ")}
}

message("\n   Model formula:\n   ", form_full, "\n")


### check covariates for correlation/colinearity

cor_mat <- canCorPairs(as.formula(form_full), meta)

t1 = as_tibble(cbind(var = rownames(cor_mat), cor_mat))
write_csv(t1, file = paste0(out_dir, script_ind, "varPart_mixed_model_covariate_correlation.csv"))


pdf(file = paste0(out_dir,script_ind, "varPart_mixed_model_covariate_correlation.pdf"),
    width = 7, height = 7)
{
  plotCorrMatrix(cor_mat)
}
dev.off()


### Fit variance partition model

param <- SnowParam(6, "SOCK", progressbar = TRUE)   # adjust cores

varPart <- fitExtractVarPartModel(vsd, as.formula(form_full), meta, BPPARAM = param)


### Summaries

vp_mean <- sort(colMeans(varPart), decreasing = TRUE)
print(vp_mean)
vp_mean_tab = tibble(covariate = names(vp_mean), Variance_explained_mean = vp_mean)

write_csv(varPart, file = paste0(out_dir, script_ind, "varPart_variance_expl_by_gene.csv"))
write_csv(vp_mean_tab, file = paste0(out_dir, script_ind, "varPart_variance_expl_mean_ord.csv"))


# Plot mean variance explained by covariate
p1 <- plotVarPart(varPart) + ggtitle("Variance explained by covariates")+
  scale_x_discrete(limits = names(vp_mean))

p2 <- p1 + scale_y_continuous(limits = c(0,20)) + ggtitle("Variance explained by covariates (max 20%)")

pdf(file = paste0(out_dir,script_ind, "varPart_variance_expl_mean_ord.pdf"),
    width = 10, height = 10)
{
  print(p1)
  print(p2)
}
dev.off()



###save bulk_dataset with cleaned metadata and analysed covars

bulk_data$varPart_analysis$primary_var = primary_var
bulk_data$varPart_analysis$candidate_covars = candidate_covars
bulk_data$varPart_analysis$form_full = form_full
bulk_data$varPart_analysis$varPart = varPart

qsave(bulk_data, file = paste0(out_dir,script_ind, "bulk_data.qs"))




#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed LD_H01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")
