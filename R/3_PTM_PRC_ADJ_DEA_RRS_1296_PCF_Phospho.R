####===================================
# Author: Dennis Friedel, PhD
# Date: 2026-09-07
# Bioinformatician,
# Department of Neuropahtology, University Clinic Heidelberg
####===================================

dataset  = '20260813_RRS_1296_PCF_Phospho' # - Name of the dataset that is going to be analysed.
analysis = paste0('PTM_PRC_ADJ_DEA/', Sys.Date()) # - Name of the analysis e.g Marker Identification

set.seed(2905)
library("dplyr")
library("patchwork")
library("ggplot2")
library("SummarizedExperiment")
source("./R/utils/rip_functions.R")
source("./R/utils/utils_module.R")
source("./R/utils/import_module.R")
source("./R/utils/preprocess_module.R")
source("./R/utils/compare_module.R")

#####------- 0 Load RDS data for addtional Normalization by Protein -------#####
ptm_se <- readRDS("./output/20260813_RRS_1296_PCF_Phospho/PTM_PRC_DEA/2026-09-08/ptm_se_prc.rds")$filt

prot_se <- readRDS("./output/20260817_RR1296_PCF/WP_PRC_DEA/2026-09-08/ms_se_prc.rds")$filt
prot_se@elementMetadata[["id"]] <- gsub(";.*", "", prot_se@elementMetadata[["id"]])

#####------- Step 1.: acquire long formated dataframes with ID annotation ####
ptm_df <- reshape2::melt((assay(ptm_se)))
ptm_df$ProteinID <- plyr::mapvalues(ptm_df$Var1,
                                    ptm_se@elementMetadata[["Name"]],
                                    ptm_se@elementMetadata[["id"]]) %>%
  gsub(";.*", "", .)
colnames(ptm_df) <- c("ResidueID", "SampleID", "Intensity", "ProteinID")

wp_df <- reshape2::melt((assay(prot_se)))
wp_df$ProteinID <- plyr::mapvalues(wp_df$Var1,
                                   prot_se@elementMetadata[["row_id"]],
                                   prot_se@elementMetadata[["id"]]) %>%
  gsub(";.*", "", .)
colnames(wp_df) <- c("Gene", "SampleID", "Intensity", "ProteinID")

#####-------  Step 2.: Join Phospho to Protein data and calculate the ratio ####
phospho_adj <- ptm_df %>%
  dplyr::left_join(wp_df, by = c("ProteinID", "SampleID")) %>%
  dplyr::mutate(Intensity_adj = Intensity.x - Intensity.y)

#####------- Step 3.: Create new SummarizedExperiment object with adjusted values ####
wide_df <- reshape(
  phospho_adj[, c("ResidueID", "SampleID", "Intensity_adj")],
  idvar = "ResidueID",
  timevar = "SampleID",
  direction = "wide"
)
rownames(wide_df) <- wide_df$ResidueID
wide_df$ResidueID <- NULL
colnames(wide_df) <- gsub(".*[.]", "", colnames(wide_df))

ptm_se_adj<- ptm_se <- add_assay(
  se_object = ptm_se,
  new_assay_name = "adjusted",
  new_assay  = as.matrix(wide_df)
) %>%
  set_primary_assay(., primary_assay_name = "adjusted")

#####------- Step 4.: Impute -------#####
assay_imp <-
  imputeLCMD::impute.MAR.MNAR(
    dataSet.mvs = SummarizedExperiment::assay(ptm_se_adj),
    method.MNAR = "MinProb",
    model.selector = rep(0, nrow(ptm_se_adj))
)

ptm_se_adj_imp <-
  add_assay(ptm_se_adj,
            new_assay = as.matrix(assay_imp),
            new_assay_name = "imputed") %>%
  set_primary_assay(., primary_assay_name = "imputed")

norm_dist <- plot_normalized_dist(se_object = ptm_se_adj_imp,
                                  assay_name = "adjusted",
                                  group = "group")

norm_imp_dist <- plot_normalized_dist(se_object = ptm_se_adj_imp,
                                  assay_name = "imputed",
                                  group = "group")

ggplot2::ggsave(
  filename = save_here(
    object_name = paste0("adjusted_norm_plot.pdf"),
    analysis = analysis
  ),
  plot = norm_dist,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300
)

saveRDS(ptm_se_adj_imp, save_here(object_name = "ptm_se_adjprc.rds"))

#####------- Step 5.: Differential Expression -------#####
ptm_se_adj_imp$group <- toupper(ptm_se_adj_imp$group)
dea_res_adj <- wrapper_dea_gsea(
  ms_se = ptm_se_adj_imp,
  wrp_group = "group",
  wrp_test = "AMP",
  wrp_contrast = list(c("WT")),
  wrp_mode = "MANUAL",
  wrp_pv_fil = 0.05,
  wrp_fc_fil = 0.58,
  wrp_use_padj = T,
  gene_set_catalouge = NULL
)

ggplot2::ggsave(
  filename = save_here(
    object_name = paste0("WT_AMP_dea_padj_plot.pdf"),
    analysis = analysis
  ),
  plot = dea_res_adj$enhanced_volcanos[[1]],
  width = 10,
  height = 10,
  units = "in",
  dpi = 300
)

dea_res <- wrapper_dea_gsea(
  ms_se = ptm_se_adj_imp,
  wrp_group = "group",
  wrp_test = "AMP",
  wrp_contrast = list(c("WT")),
  wrp_mode = "MANUAL",
  wrp_pv_fil = 0.05,
  wrp_fc_fil = 0.58,
  wrp_use_padj = F,
  gene_set_catalouge = NULL
)
ttresult <- dea_res$tt_combined
ttresult[ttresult$significant,]

dea_res$tt_combined <- NULL
purrr::map2(dea_res, names(dea_res), function(x, y) {
  ggplot2::ggsave(
    filename = save_here(
      object_name = paste0((y), "_dea_plot.pdf"),
      analysis = analysis
    ),
    plot = x,
    width = 10,
    height = 10,
    units = "in",
    dpi = 300
  )
})
openxlsx::write.xlsx(ttresult, save_here(object_name = "TopTable_protadjusted_WT_AMP.xlsx"))
saveRDS(ptm_se_adj_imp, save_here(object_name = "ptm_se_imp.rds"))

#####------- Step 6.: Compare with non-normalized  -------#####
ttunadjust<-openxlsx::read.xlsx("./output/20260813_RRS_1296_PCF_Phospho/PTM_PRC_DEA/2026-09-09/TopTable_WT_AMP.xlsx")
ttun<-ttunadjust$genes[ttunadjust$significant]
ttadj<-ttresult$genes[ttresult$significant]
intersect(ttadj,ttun)

vd<-ggVennDiagram::ggVennDiagram(x = list("Site Adjusted" = ttadj, "Unadjusted" =
                                        ttun)) + coord_flip()
ggsave(plot = vd,
       filename = 
       save_here(object_name = "Venndiagram_Sites.pdf"),
       width = 7.5,
       height = 5)
