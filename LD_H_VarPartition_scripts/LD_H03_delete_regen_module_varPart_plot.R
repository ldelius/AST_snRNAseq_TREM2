### TEMP / DELETE-ME: regenerate ONLY the module-level varPart CSV + PDF for H03 v03
### using the robust per-module lme4::lmer + calcVarPart approach
### (fitExtractVarPartModel is broken in this env: lme4 2.0 removed findbars -> it
### silently drops modules). Reads existing LD_H03_v03_bulk_data.qs; does NOT rerun
### WGCNA or anything else. Overwrites only:
###   LD_H03_v03_varPart_variance_expl_by_module.csv
###   LD_H03_v03_varPart_variance_expl_by_module.pdf
###
### Run on the HPC:  Rscript LD_H03_delete_regen_module_varPart_plot.R   (light, ~1-2 min)

library(qs)
library(tidyverse)
library(pheatmap)
library(variancePartition)   # for calcVarPart()
library(lme4)

out_dir = "/rds/general/user/lvd25/home/AST_scRNAseq_TREM2/LD_H_VarPartition_output/"
script_ind = "LD_H03_v03_"

bulk_data = qread(paste0(out_dir, script_ind, "bulk_data.qs"))

mod_act_mat = bulk_data$wgcna$mod_activity_mat
meta        = bulk_data$meta[colnames(mod_act_mat), ]
form_full   = bulk_data$varPart_analysis$form_full
form_mod_fit = as.formula(paste("y", form_full))   # y ~ <covariates>

### per-module variance decomposition (all modules; robust to fit failures)
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
        nrow(mod_act_mat), " modules: ", paste(rownames(varPart), collapse = ", "), "\n")

### write CSV (module label kept as a column)
write_csv(cbind(module = rownames(varPart), as_tibble(varPart)),
          file = paste0(out_dir, script_ind, "varPart_variance_expl_by_module.csv"))

### plot heatmap (same style as H03)
pdf(file = paste0(out_dir, script_ind, "varPart_variance_expl_by_module.pdf"),
    width = 10, height = 10)
{
  pheatmap::pheatmap(t(varPart[, colnames(varPart) != "Residuals"]),
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     show_rownames = TRUE, show_colnames = TRUE,
                     color = colorRampPalette(c("white", "red"))(250),
                     breaks = seq(0, max(varPart[, colnames(varPart) != "Residuals"]), length.out = 251),
                     border_color = NA, cellwidth = 10, cellheight = 10,
                     fontsize = 10,
                     main = "Covariates fraction variance explained vs modules")
}
dev.off()

message("\n   done - regenerated module varPart CSV + PDF for ", script_ind, "\n")
