####===================================####
# Author: Dennis Friedel, PhD
# Date: 2026-09-07
# Bioinformatician,
# Department of Neuropahtology, University Clinic Heidelberg
####===================================####

dataset  = '20260813_RRS_1296_PCF_Phospho_Protnorm' # - Name of the dataset that is going to be analysed.
analysis = paste0('PTM_PRC_ADJ_DEA/',Sys.Date()) # - Name of the analysis e.g Marker Identification

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

### Load data
mq_results<-
  load_ms_results(ms_result_dir = "./data//20260813_RRS_1296_PCF_Phospho/MaxQuant2.4.2.0",
                             file_names = c("summary.txt",
                             "peptides.txt",
                             "evidence.txt",
                             "proteinGroups.txt",
                             "Phospho_STY_Sites.txt"))

mq_list <- create_maxquant_se_list(
  mq_results = mq_results,
  ptm_probability_flt = 0.75,
  remove_low_quality = T
)

### Annoate
metadata<-openxlsx::read.xlsx("./data/Sample_annotation.xlsx",sheet = 1)
rownames(metadata)<-metadata$sample_id

ptm_se<-mq_list$ptm_se
ptm_se<-ptm_se[!duplicated(ptm_se@elementMetadata$Name),]
rownames(ptm_se)<-ptm_se@elementMetadata$Name

ptm_se<-ptm_se[,ptm_se$sample_id%in%metadata$original_id]
colnames(ptm_se)<-metadata$sample_id
SummarizedExperiment::colData(ptm_se)<-S4Vectors::DataFrame(metadata[match(ptm_se$sample_id,metadata$original_id),])

ptm_se_flt<-ptm_se[,ptm_se$group%in%c("Amp","WT")]
saveRDS(ptm_se_flt,save_here(object_name = "ptm_se_raw.rds"))

### Preprocess & Normalize by Protein
ptm_se_prc <- pre_processing_wrapper(
  se_object = ptm_se_flt,
  group_sel = "group",
  filter_fractioned = T,
  thr_sample = 3500,
  thr_feature = 0.75,
  normalization_method = "median_center",
  imputation_method  = "MinProb"
)
qc_plots <- qp_plots_wrapper(se_proc_ls = ptm_se_prc,
                             thr_sample = 4000,
                             group_sel = "group")
qc_plots
gc()

purrr::map2(qc_plots, names(qc_plots), function(x, y) {
  ggplot2::ggsave(
    filename = save_here(
      object_name = paste0((y), "_qc_plot.pdf"),
      analysis = analysis
    ),
    plot = x,
    width = 10,
    height = 10,
    units = "in",
    dpi = 300
  )
})

#### Addtional Normalization by Protein
prot_se<-readRDS("./output/20260817_RR1296_PCF/WP_PRC_DEA/2026-09-08/ms_se_prc.rds")$filt
prot_se@elementMetadata[["id"]]<-gsub(";.*","",prot_se@elementMetadata[["id"]])

# Step 1.: acquire long formated dataframes with ID annotation
ptm_df<-reshape2::melt((assay(ptm_se_prc$filt)))
ptm_df$ProteinID <- plyr::mapvalues(ptm_df$Var1,
                                    ptm_se_prc$filt@elementMetadata[["Name"]],
                                    ptm_se_prc$filt@elementMetadata[["id"]])%>%
  gsub(";.*","",.)
colnames(ptm_df)<-c("ResidueID","SampleID","Intensity","ProteinID")
wp_df<-reshape2::melt((assay(prot_se)))
wp_df$ProteinID <- plyr::mapvalues(wp_df$Var1,
                                   prot_se@elementMetadata[["row_id"]],
                                   prot_se@elementMetadata[["id"]])%>%
  gsub(";.*","",.)
colnames(wp_df)<-c("Gene","SampleID","Intensity","ProteinID")

# Step 2.: Join Phospho to Protein data and calculate the ratio
phospho_adj<-ptm_df%>%
  dplyr::left_join(wp_df,by=c("ProteinID","SampleID"))%>%
  dplyr::mutate(Intensity_adj=Intensity.x-Intensity.y)
rownames(phospho_adj)
# Step 3.: Create new SummarizedExperiment object with adjusted values
wide_df <- reshape(
  phospho_adj[, c("ResidueID", "SampleID", "Intensity_adj")],
  idvar = "ResidueID",
  timevar = "SampleID",
  direction = "wide"
)
rownames(wide_df)<-wide_df$ResidueID
wide_df$ResidueID<-NULL
colnames(wide_df)<-gsub(".*[.]","",colnames(wide_df))

ptm_se_prc$filt<-add_assay(se_object = ptm_se_prc$filt,new_assay_name = "adjusted",
          new_assay  = as.matrix(wide_df))
ptm_se_prc$filt<-set_primary_assay(ptm_se_prc$filt,primary_assay_name = "adjusted")

# Step 4.: Impute
se_object_adjust<-ptm_se_prc$filt
assay_imp <-
  imputeLCMD::impute.MAR.MNAR(
    dataSet.mvs = SummarizedExperiment::assay(se_object_adjust),
    method.MNAR = "MinProb",
    model.selector = rep(0, nrow(se_object_adjust))
  )
se_object_adjimp <-
  add_assay(se_object_adjust,
            new_assay = as.matrix(assay_imp),
            new_assay_name = "imputed") %>%
  set_primary_assay(., primary_assay_name = "imputed")

norm_dist <- plot_normalized_dist(
  se_object = se_object_adjimp,
  assay_name = "adjusted",
  group = "group"
)
norm_dist
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

saveRDS(se_object_adjimp, save_here(object_name = "ptm_se_adjprc.rds"))

### Differential Expression
ptm_se_imp <- se_object_adjimp
ptm_se_imp$group <- toupper(ptm_se_imp$group)

dea_res <- wrapper_dea_gsea(
  ms_se = ptm_se_imp,
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
sum(ttresult$significant)

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
openxlsx::write.xlsx(ttresult, save_here(object_name = "TopTable_WT_AMP.xlsx"))
saveRDS(ptm_se_imp, save_here(object_name = "ptm_se_imp.rds"))
ttresult
