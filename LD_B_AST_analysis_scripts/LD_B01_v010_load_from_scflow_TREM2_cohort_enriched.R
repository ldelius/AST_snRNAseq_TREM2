# this file takes the scflow pipwlinw outputs and assembles them into a single Seurat object. saved as B01_seur.qs

message("\n\n##########################################################################\n",
        "# Start B01: Load sce dataset from scflow into Seurat ", Sys.time(),
        "\n##########################################################################\n",
        "\n   Load metadata, create basic plot, \n",
        "\n##########################################################################\n\n")
# this just prints a message when the script starts including a timestamp


#set environment/load packages
library(qs)
library(tidyverse)
library(Seurat) # main single cell analysis framework
library(colorRamps)
library(viridis)
library(ComplexHeatmap)
library(patchwork)

### define directories and script index

main_dir = paste0("/rds/general/user/mlattke/projects/ukdrmultiomicsproject/live/Users/MichaelL/",
                  "P02E03_260121_TREM2_AD_tissue_snRNAseq_MIC_AST_scflow/")
setwd(main_dir)

#specify script/output index as prefix for file names
script_ind = "B01_"

#specify output directory
out_dir = paste0(main_dir,"B_load_from_scflow_subcluster/")
if (!dir.exists(out_dir)){dir.create(out_dir)}

#main scflow dataset directory
dataset_dir = paste0("/rds/general/user/mlattke/projects/ukdrmultiomicsproject/live/",
                        "MAP_analysis/TREM2_enriched_scflow/")


### load and curate group and file info for analysis dataset

t1 = read_tsv("A_input/TREM2_Samplesheet_snRNAseq_GliaEnriched_with_cohort_info.tsv")
names(t1)

t2 = t1[order(match(t1$NeuropathologicalDiagnosis, c("Control", "AD")),  # sorts the samples in an order
              t1$TREM2Variant, t1$pctPHF1PositiveArea),]

t3 = tibble(group = paste0(t2$NeuropathologicalDiagnosis, "_", t2$TREM2Variant))

gr_tab = cbind(t3, t2)

write_csv(gr_tab, file = "A_input/group_tab.csv")



### load counts and metadata, convert to seurat

#load counts matrix

m1 = ReadMtx(mtx = paste0(dataset_dir, "results/final/SCE/final_sce_added_metadata/full_sce/matrix.mtx.gz"),
             cells = paste0(dataset_dir, "results/final/SCE/final_sce_added_metadata/full_sce/barcodes.tsv.gz"),
             features = paste0(dataset_dir, "results/final/SCE/final_sce_added_metadata/full_sce/features.tsv.gz"),
             feature.column = 1)

seur = CreateSeuratObject(m1)

#add cell metadata

t1 = read_tsv(paste0(dataset_dir, "results/final/SCE/final_sce_added_metadata/full_sce/sce-coldata.tsv"))

t2 = seur@meta.data
t3 = cbind(t2, group = paste0(t1$NeuropathologicalDiagnosis, "_", t1$TREM2Variant), t1)

seur@meta.data = t3

#add gene names

t1 = read_tsv(paste0(dataset_dir, "results/final/SCE/final_sce_added_metadata/full_sce/sce-rowdata.tsv"))

t2 = seur@assays$RNA
identical(rownames(t2), t1$ensembl_gene_id)
t1$gene[duplicated(t1$gene)]

rownames(seur@assays$RNA) = t1$gene


# add dimensionality reductions Liger

t1 = read_tsv(paste0(dataset_dir, "results/final/SCE/final_sce_added_metadata/full_sce/ReducedDim_Liger.tsv"))
rownames(t1) = seur$barcode
colnames(t1) = paste0("liger", 1:ncol(t1))

seur@reductions$Liger = t1

t1 = read_tsv(paste0(dataset_dir, "results/final/SCE/final_sce_added_metadata/full_sce/ReducedDim_UMAP_Liger.tsv"))
rownames(t1) = seur$barcode
colnames(t1) = paste0("umapLiger", 1:ncol(t1))

seur@reductions$UMAP_Liger = t1

qsave(seur, file = paste0(out_dir, script_ind, "seur.qs"))




####################################
#Functions
####################################

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


#distinct scale (larger number of colours, colour vector gets shuffled)
pal_dist = function(v){
  v2 = length(unique(v))
  if (v2 < 6){p2 = matlab.like(6)[1:v2]} else {
    p2 = matlab.like(v2)
    set.seed(12)
    p2 = sample(p2)
  }
  return(p2)
}



#######################################
#plot umap with samples, cell types, pathology, TREM2 genotype, cohort
#######################################

#seur = qread(file = paste0(out_dir, script_ind, "seur.qs"))


### plot whole dataset as UMAP with named clusters/cell classes

#define grouping variables
gr = unique(gr_tab$group)
samples = unique(gr_tab$sample)
cell_types = unique(seur$cluster_celltype)
cell_types = cell_types[order(cell_types)]


#subsample seurat object for plotting for large datasets (else plots become too large (pdf with hundreds of MB))

seur0 = seur #seur0 saves a copy of the full object

cells = seur$barcode

if (length(cells)>100000){ # if more than 100,000 cells, randomly sample 100,000 for plotting
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


#umap plots for cell type (originally asigned cell type; manual plotting with ggplot (seurat plots do not work without dimensionality reductions))
pl = list()

t1 = seur@meta.data
t2 = seur@reductions$UMAP_Liger
t2 = t2[rownames(t1),]

pl_tab = cbind(t1, t2)
set.seed(1234)
pl_tab = pl_tab[sample(rownames(pl_tab), size =100000, replace=F),]


p1 = ggplot(pl_tab, aes(x = umapLiger1, y = umapLiger2))+
  theme_bw()

pl[["cell_type"]] = p1+
  geom_point(aes(colour = cluster_celltype), size = 0.5)+
  scale_color_manual(limits = cell_types, values = pal(cell_types))+
  labs(title = "UMAP by cell type")

pl[["group"]] = p1+
  geom_point(aes(colour = group), size = 0.5)+
  scale_color_manual(limits = gr, values = pal(gr))+
  labs(title = "UMAP by group")

pl[["sample"]] = p1+
  geom_point(aes(colour = sample), size = 0.5)+
  scale_color_manual(limits = samples, values = pal(samples))+
  labs(title = "UMAP by sample")

pl[["sample_no_legend"]] = p1+
  geom_point(aes(colour = sample), size = 0.5)+
  scale_color_manual(limits = samples, values = pal(samples))+
  labs(title = "UMAP by group")+NoLegend()

pl[["cohort"]] = p1+
  geom_point(aes(colour = cohort), size = 0.5)+
  scale_color_manual(limits = unique(gr_tab$cohort), values = pal(unique(gr_tab$cohort)))+
  labs(title = "UMAP by cohort")

pl[["diagnosis"]] = p1+
  geom_point(aes(colour = NeuropathologicalDiagnosis), size = 0.5)+
  scale_color_manual(limits = unique(gr_tab$NeuropathologicalDiagnosis), 
                     values = pal(unique(gr_tab$NeuropathologicalDiagnosis)))+
  labs(title = "UMAP by diagnosis")

pl[["TREM2Variant"]] = p1+
  geom_point(aes(colour = TREM2Variant), size = 0.5)+
  scale_color_manual(limits = unique(gr_tab$TREM2Variant), 
                     values = pal(unique(gr_tab$TREM2Variant)))+
  labs(title = "UMAP by TREM2Variant")


pdf(file = paste0(out_dir,script_ind, "UMAP_Liger_complete_dataset.pdf"), width = 6, height = 5)
lapply(pl, function(x){x})
dev.off()


message("\n\n##########################################################################\n",
        "# Completed B01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")

