####===================================####
# Author: Dennis Friedel, PhD
# Date: 2026-09-07
# Bioinformatician,
# Department of Neuropahtology, University Clinic Heidelberg
####===================================####

dataset  = '20260813_RRS_1296_PCF_Phospho' # - Name of the dataset that is going to be analysed.
analysis = paste0('PTM_PRC_DEA/',Sys.Date()) # - Name of the analysis e.g Marker Identification
"PTM_PRC_DEA_20260813_RRS_1296_PCF_Phospho.R"

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

### Preprocess
ptm_se_prc <- pre_processing_wrapper(
  se_object = ptm_se_flt,
  group_sel = "group",
  filter_fractioned = T,
  thr_sample = 3500,
  thr_feature = 0.75,
  normalization_method = "median_center",
  imputation_method  = "MinProb"
)

saveRDS(ptm_se_prc, save_here(object_name = "ptm_se_prc.rds"))

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

### Differential Expression
ptm_se_imp <- ptm_se_prc$imp
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
#ttresult[ttresult$AMP_vs_WT_adj_P_Val < 0.05, ]

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
