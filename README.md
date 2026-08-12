# Single-nucleus transcriptomic analysis of astrocyte responses to Alzheimer’s disease in human TREM2 R47H and R62H carriers
### An MRes Thesis Project (Biomedical Research - Data Science, Imperial College London)

## Abstract

Alzheimer’s disease (AD) is the most common form of dementia. Among the genetic variants raising its risk are two rare coding variants in TREM2, R47H and R62H. TREM2 is a receptor expressed mainly by microglia in the brain. R47H and R62H alter microglial responses to AD pathology, with partly shared and partly variant-specific effects. Given that microglia signal to astrocytes, which support neuronal metabolism and synaptic function, TREM2 variants may also alter the astrocyte response to AD. While this could contribute substantially to the increased disease susceptibility, astrocyte responses in variant carriers have received little attention.
This study examines whether and how R47H and R62H influence astrocyte transcriptional responses to AD by analysing a glia-enriched single-nucleus RNA-sequencing dataset from human post-mortem cortex. The dataset comprises 225,292 astrocyte nuclei from 70 donors (43 AD, 27 control), including 10 R47H and 17 R62H carriers. Astrocytes were first subclustered and characterised, and their abundance was tested across diagnosis and genotype. Differential expression between these groups was then tested on pseudobulk profiles, and the associated biological processes were identified by Gene Ontology over-representation, gene set enrichment and co-expression network analysis. The study design and analysis workflow are summarised in Figure 1.

<p align="center">
  <img src="doc/images/optimise_pp_method.png?v=2" width="900">
</p>

**Figure 1.** Study design and analysis workflow.

Input data, intermediate analysis objects, and generated results are not included in this repository.

## Analyses and Visualisation Scripts

Thesis status is given after each analysis description.

#### Astrocyte Preprocessing and Characterisation

- **Load scFlow output into Seurat** *(used in thesis)*

  - [`LD_B01_v010_load_from_scflow_TREM2_cohort_enriched.R`](LD_B_AST_analysis_scripts/LD_B01_v010_load_from_scflow_TREM2_cohort_enriched.R)

- **Astrocyte subsetting, integration, and first-round clustering** *(used in thesis)*

  - [`LD_B02a_v010_subset_reintegrate_reclustering_tests_AST_round1.R`](LD_B_AST_analysis_scripts/LD_B02a_v010_subset_reintegrate_reclustering_tests_AST_round1.R)

- **Removal of non-astrocytic clusters and second-round subclustering** *(used in thesis)*

  - [`LD_B03a_v010_subset_reintegrate_reclustering_tests_AST_round2.R`](LD_B_AST_analysis_scripts/LD_B03a_v010_subset_reintegrate_reclustering_tests_AST_round2.R)

- **Subcluster annotation and final astrocyte object** *(used in thesis)*

  - [`LD_B04a_v036_subcluster_charact_abund_analysis_by_plaque_dens.R`](LD_B_AST_analysis_scripts/LD_B04a_v036_subcluster_charact_abund_analysis_by_plaque_dens.R)


#### Pseudobulk Differential Expression and Enrichment

- **Pseudobulk generation by astrocyte subcluster and sample** *(used in thesis)*

  - [`LD_E01_v051_seur_pseudobulk_generation_min_20cells_per_pb_by_cluster.R`](LD_E_DESeq_pseudobulk_scripts/LD_E01_v051_seur_pseudobulk_generation_min_20cells_per_pb_by_cluster.R)

- **Five-covariate differential expression** *(used in thesis)*

  - [`LD_E02c_v01_pseudobulk_DESeq2_LRT_by_clust_by_diagn_5covar_corr_cohort_APOE_CD33_BrainReg_Sex.R`](LD_E_DESeq_pseudobulk_scripts/LD_E02c_v01_pseudobulk_DESeq2_LRT_by_clust_by_diagn_5covar_corr_cohort_APOE_CD33_BrainReg_Sex.R)

- **Hallmark and Green et al. gene set enrichment** *(used in thesis)*

  - [`LD_E03c_GSEA_5covar.R`](LD_E_DESeq_pseudobulk_scripts/LD_E03c_GSEA_5covar.R)

- **Cumulative age, post-mortem interval, and Braak adjustment** *(used in thesis)*

  - [`LD_E04c_incremental_covar_5covar_base.R`](LD_E_DESeq_pseudobulk_scripts/LD_E04c_incremental_covar_5covar_base.R)

- **Earlier differential-expression and enrichment workflows** *(superseded)*

  - [`LD_E02a2_v061_pseudobulk_DESeq2_LRT_by_clust_by_diagn_TREM2_R47H_vs_R62H_corr_APOE_CD33_BrainReg.R`](LD_E_DESeq_pseudobulk_scripts/LD_E02a2_v061_pseudobulk_DESeq2_LRT_by_clust_by_diagn_TREM2_R47H_vs_R62H_corr_APOE_CD33_BrainReg.R)
  - [`LD_E03a2_v061_DEG_characterisation_pairwise_comps_GO.R`](LD_E_DESeq_pseudobulk_scripts/LD_E03a2_v061_DEG_characterisation_pairwise_comps_GO.R)
  - [`LD_E04a_v01_incremental_covar_adj_log2FC_corr_sel_comps.R`](LD_E_DESeq_pseudobulk_scripts/LD_E04a_v01_incremental_covar_adj_log2FC_corr_sel_comps.R)


#### Co-expression Network Analysis

- **Pseudobulk generation for WGCNA** *(used in thesis)*

  - [`LD_F01_v060_seur_pseudobulk_generation_min_20cells_per_pb_by_cluster.R`](LD_F_DESeq_pseudobulk_WGCNA_scripts/LD_F01_v060_seur_pseudobulk_generation_min_20cells_per_pb_by_cluster.R)

- **Seven-covariate selection of genes associated with either TREM2 variant** *(used in thesis)*

  - [`LD_F02c_v01_pseudobulk_DESeq2_LRT_comb_AD_both_variants_7covar_corr_cohort_APOE_CD33_Region_Sex_PMI_libsize.R`](LD_F_DESeq_pseudobulk_WGCNA_scripts/LD_F02c_v01_pseudobulk_DESeq2_LRT_comb_AD_both_variants_7covar_corr_cohort_APOE_CD33_Region_Sex_PMI_libsize.R)

- **Final WGCNA network, robustness analyses, hub genes, and Figure 7** *(used in thesis)*

  - [`LD_F03c_v02_WGCNA_AD_both_variants_7covar_corr_DEGseed.R`](LD_F_DESeq_pseudobulk_WGCNA_scripts/LD_F03c_v02_WGCNA_AD_both_variants_7covar_corr_DEGseed.R)
  - [`LD_F03d_v01_WGCNA_AD_DEGseed_pathology_robustness.R`](LD_F_DESeq_pseudobulk_WGCNA_scripts/LD_F03d_v01_WGCNA_AD_DEGseed_pathology_robustness.R)
  - [`LD_F03e_v01_WGCNA_kME_hub_genes_pathology_LMM.R`](LD_F_DESeq_pseudobulk_WGCNA_scripts/LD_F03e_v01_WGCNA_kME_hub_genes_pathology_LMM.R)
  - [`LD_F03e_v01_replot.R`](LD_F_DESeq_pseudobulk_WGCNA_scripts/LD_F03e_v01_replot.R)

- **Earlier R62H-only WGCNA workflow** *(superseded)*

  - [`LD_F02a1_v060_pseudobulk_DESeq2_LRT_comb_AD_TREM2_R62H_by_diagn_corr_APOE_CD33_Region_low_p_cutoff.R`](LD_F_DESeq_pseudobulk_WGCNA_scripts/LD_F02a1_v060_pseudobulk_DESeq2_LRT_comb_AD_TREM2_R62H_by_diagn_corr_APOE_CD33_Region_low_p_cutoff.R)
  - [`LD_F03a1_v060_DEG_charact_WGCNA_TREM2_R62H_by_diagn_corr_APOE_CD33_Region_low_p_cutoff.R`](LD_F_DESeq_pseudobulk_WGCNA_scripts/LD_F03a1_v060_DEG_charact_WGCNA_TREM2_R62H_by_diagn_corr_APOE_CD33_Region_low_p_cutoff.R)

- **TREM2 variance within WGCNA GO terms** *(not used in thesis / exploratory)*

  - [`LD_F04c_v01_TREM2varExpl_per_v02_WGCNA_GO_term.R`](LD_F_DESeq_pseudobulk_WGCNA_scripts/LD_F04c_v01_TREM2varExpl_per_v02_WGCNA_GO_term.R)


#### Variance Partitioning

- **Candidate-covariate variance partitioning and Supplementary Figure 3** *(used in thesis)*

  - [`LD_H01_pseudobulk_varPart_mixed_model.R`](LD_H_VarPartition_scripts/LD_H01_pseudobulk_varPart_mixed_model.R)

- **TREM2-associated variance and WGCNA follow-up** *(not used in thesis / exploratory)*

  - [`LD_H02_pseudobulk_varPart_var_genes_char_corr_tech_bio_vars.R`](LD_H_VarPartition_scripts/LD_H02_pseudobulk_varPart_var_genes_char_corr_tech_bio_vars.R)
  - [`LD_H03_varPart_char_WGCNA_corr_tech_bio_vars_TREM2_thr0.05.R`](LD_H_VarPartition_scripts/LD_H03_varPart_char_WGCNA_corr_tech_bio_vars_TREM2_thr0.05.R)


#### Microglia-Astrocyte Communication Analysis

- **Preparation of astrocyte and microglial inputs for communication analysis** *(not used in thesis / exploratory)*

  - [`LD_G01a_AST_MIC_pseudobulk_merge.R`](LD_G_MIC_AST_communication_analysis/LD_G01a_AST_MIC_pseudobulk_merge.R)
  - [`LD_G01b_v001_AST_MIC_seurat_merge.R`](LD_G_MIC_AST_communication_analysis/LD_G01b_v001_AST_MIC_seurat_merge.R)

- **CellChat inference with a DEG-filtered interaction database** *(not used in thesis / exploratory)*

  - [`LD_G02a_ungrouped_DEG_filter_cellchat_analysis_MIC_AST.R`](LD_G_MIC_AST_communication_analysis/LD_G02a_ungrouped_DEG_filter_cellchat_analysis_MIC_AST.R)
  - [`LD_G02b_by_TREM2_variant_DEG_filter_cellchat_analysis_MIC_AST.R`](LD_G_MIC_AST_communication_analysis/LD_G02b_by_TREM2_variant_DEG_filter_cellchat_analysis_MIC_AST.R)
  - [`LD_G02c_by_TREM2variantXdiagnosis_DEG_filter_cellchat_analysis_MIC_AST.R`](LD_G_MIC_AST_communication_analysis/LD_G02c_by_TREM2variantXdiagnosis_DEG_filter_cellchat_analysis_MIC_AST.R)

- **CellChat visualisation and pairwise comparisons** *(not used in thesis / exploratory)*

  - [`LD_G03a_visualisation_ungrouped_cellchat_analysis_MIC_AST.R`](LD_G_MIC_AST_communication_analysis/LD_G03a_visualisation_ungrouped_cellchat_analysis_MIC_AST.R)
  - [`LD_G03b_visualisation_byTREM2variant_cellchat_analysis_MIC_AST.R`](LD_G_MIC_AST_communication_analysis/LD_G03b_visualisation_byTREM2variant_cellchat_analysis_MIC_AST.R)
  - [`LD_G03c_visualisation_byTREM2variantXdiagnosis_cellchat_analysis_MIC_AST.R`](LD_G_MIC_AST_communication_analysis/LD_G03c_visualisation_byTREM2variantXdiagnosis_cellchat_analysis_MIC_AST.R)
  - [`LD_G03d_CellChat_comparison.R`](LD_G_MIC_AST_communication_analysis/LD_G03d_CellChat_comparison.R)

The G02 and G03 scripts use CellChat (1) with G01 preparing their inputs. Cell-cell communication analysis is not included in the submitted thesis. For future work and exploration will adapt LIANA (2), replacing CellChat. The G workflow also depends on microglial objects not included in this repository.

#### References
1: Jin S, Guerrero-Juarez CF, Zhang L, Chang I, Ramos R, Kuan C-H, et al. Inference and analysis of cell-cell communication using CellChat. *Nature Communications*. 2021;12:1088. doi:10.1038/s41467-021-21246-9.

2: Dimitrov D, Türei D, Garrido-Rodriguez M, Burmedi PL, Nagai JS, Boys C, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nature Communications*. 2022;13:3224. doi:10.1038/s41467-022-30755-0. 


#### Thesis Tables and Figures

- **Cohort summary and Supplementary Table 1** *(used in thesis)*

  - [`LD_X01_cohort_summary_tables.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X01_cohort_summary_tables.R)

- **Figure 2 and Supplementary Figure 2** *(used in thesis)*

  - [`LD_X02_overview_UMAP_figure.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X02_overview_UMAP_figure.R)
  - [`LD_X03a_supp_marker_dotplot_round1_res1p5.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X03a_supp_marker_dotplot_round1_res1p5.R)
  - [`LD_X12_Supp_Fig2_abc.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X12_Supp_Fig2_abc.R)

- **Figure 3, Supplementary Tables 3–5, and Supplementary Figure 5** *(used in thesis)*

  - [`LD_X04_B_characterisation_plots.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X04_B_characterisation_plots.R)
  - [`LD_X04c_green_state_LOO_sensitivity.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X04c_green_state_LOO_sensitivity.R)
  - [`LD_X11_Green24_astrocyte_states_table.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X11_Green24_astrocyte_states_table.R)
  - [`LD_X15_replot_uncurated_GO.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X15_replot_uncurated_GO.R)

- **Figure 4 and Supplementary Table 6** *(used in thesis)*

  - [`LD_X05_abundance_plots.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X05_abundance_plots.R)
  - [`LD_X05_v02_subtype_abundance_plot.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X05_v02_subtype_abundance_plot.R)
  - [`LD_X16_sccomp_supplementary_tables.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X16_sccomp_supplementary_tables.R)

- **Figure 5, Supplementary Table 7, and Supplementary Figure 6** *(used in thesis)*

  - [`LD_X07_v04_combined_DEG_pairwise.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X07_v04_combined_DEG_pairwise.R)
  - [`LD_X08_GO_quadrants_family_pooled_5covar.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X08_GO_quadrants_family_pooled_5covar.R)
  - [`LD_X17_covariate_robustness_table.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X17_covariate_robustness_table.R)
  - [`LD_X17b_recompute_robustness_measures.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X17b_recompute_robustness_measures.R)
  - [`LD_X22_GO_quadrants_pooled_supplementary.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X22_GO_quadrants_pooled_supplementary.R)

- **Figure 6 and Supplementary Figure 7** *(used in thesis)*

  - [`LD_X10b_GSEA_heatmap_combined_prep.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X10b_GSEA_heatmap_combined_prep.R)
  - [`LD_X19_GSEA_heatmaps_supplementary_all.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X19_GSEA_heatmaps_supplementary_all.R)

- **Software versions and Supplementary Table 2** *(used in thesis)*

  - [`LD_X13_package_versions_table.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X13_package_versions_table.R)

- **Variance-partition replot for Supplementary Figure 3** *(used in thesis)*

  - [`LD_X14_varPart_replot_bigger_title.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X14_varPart_replot_bigger_title.R)

- **WGCNA Supplementary Figure 8 and Supplementary Tables 8–9** *(used in thesis)*

  - [`LD_X18_replot_WGCNA_GO_dotplot.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X18_replot_WGCNA_GO_dotplot.R)
  - [`LD_X20_hub_genes_table.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X20_hub_genes_table.R)
  - [`LD_X21_WGCNA_pathology_robustness_table.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X21_WGCNA_pathology_robustness_table.R)

- **Earlier figure-generation workflows** *(superseded)*

  - [`LD_X03_supp_resolution_QC_UMAPs.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X03_supp_resolution_QC_UMAPs.R)
  - [`LD_X07_pseudobulk_DESeq2_figures_5covar.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X07_pseudobulk_DESeq2_figures_5covar.R)
  - [`LD_X07_v03_combined_DEG_pairwise.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X07_v03_combined_DEG_pairwise.R)

- **Additional validation and exploratory figures** *(not used in thesis / exploratory)*

  - [`LD_X05b_sccomp_contrast_verification.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X05b_sccomp_contrast_verification.R)
  - [`LD_X06_pathology_analysis.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X06_pathology_analysis.R)
  - [`LD_X09_IFN_hallmark_dotplots.R`](LD_X_Thesis_Presentatioon_Figures_scripts/LD_X09_IFN_hallmark_dotplots.R)


## Reproducibility

- Analyses were run in R on the Imperial College HPC service.
- The PBS submission wrapper is provided in [`qsub_cx3_R_R43_240426_1c32g8h.sh`](qsub_cx3_R_R43_240426_1c32g8h.sh). It loads Miniforge and activates the `R43_240426` environment.
- Package requirements are declared through `library()` calls within each script. No package lockfile is included.
- Random seeds are set within scripts where stochastic procedures are used.
- Most analysis scripts record `sessionInfo()` in their job output.
