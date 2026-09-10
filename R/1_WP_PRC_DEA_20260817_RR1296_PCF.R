####===================================####
# Author: Dennis Friedel, PhD
# Date: 2026-09-07
# Bioinformatician,
# Department of Neuropahtology, University Clinic Heidelberg
####===================================####

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

dataset  = '20260817_RR1296_PCF' # - Name of the dataset that is going to be analysed.
analysis = paste0('WP_PRC_DEA/',Sys.Date()) # - Name of the analysis e.g Marker Identification
create_rip_envir(rip_dir = "./")

#### Import data #####
#### Note #####
# In order to run this code you need to download the Data Stored in PRIDE 
# For running DIANN analyis with create_diann_se_list 
# you will need the parquet file and the stats file  

diann_result_path<-"./data/20260817_RR1296_PCF/diann_2.3.1/Lib_free/"
remove_low_quality <- TRUE
only_proteotypic <- TRUE
number_cores = 2
verbose = TRUE

msdata <- load_ms_results(ms_result_dir = diann_result_path)
msdata$result_stats_tsv$file_name<-gsub(".*/","",msdata$result_stats_tsv$file_name)%>%gsub("[.]raw","",.)

msdata<-create_diann_se_list(
  result_tsv=msdata$result_parquet,
  stats_tsv=msdata$result_stats_tsv,
  remove_low_quality = TRUE,
  use_only_proteotypic = TRUE,
  threshold_precursor = 0.05,
  protein_theshold = 0.05,
  pq_q_threshold = 0.05,
  gg_q_threshold = 0.05,
  pg_quantity=c("pg_max_lfq"),#"pg_normalised"
  gg_quantity=c("genes_max_lfq"),#"genes_normalised"
  verbose=TRUE
)

saveRDS(msdata,save_here(object_name="diann_list.rds"))

### Annoate with Supplemental Table 1 
metadata<-openxlsx::read.xlsx("./data/Sample_annotation.xlsx",sheet = 1)
metadata$original_id<-gsub("_p_raw","",metadata$original_id)
rownames(metadata)<-metadata$sample_id

### Use Intensities of Proteins mapped to Genes 
ms_se<-msdata$ggu_se
ms_se$sample_id<-plyr::mapvalues(colnames(ms_se),metadata$original_id,metadata$sample_id)
new_cd<-plyr::join(as.data.frame(colData(ms_se)),metadata,by = "sample_id")%>%DataFrame()
rownames(new_cd)<-new_cd$sample_id
colnames(ms_se)<-plyr::mapvalues(colnames(ms_se),metadata$original_id,metadata$sample_id)
colData(ms_se)<-new_cd

### Drop samples which are not amp or WT
ms_se_flt<-ms_se[,ms_se$group%in%c("WT","Amp")]

### Run Preproc and DEA
### Preprocess
ms_se_prc <- pre_processing_wrapper(
  se_object = ms_se_flt,
  group_sel = "group",
  filter_fractioned = T,
  thr_sample = 5000,
  thr_feature = 0.75,
  normalization_method = "none",
  imputation_method  = "MinProb"
)

saveRDS(ms_se_prc, save_here(object_name = "ms_se_prc.rds"))

qc_plots <- qp_plots_wrapper(se_proc_ls = ms_se_prc,
                                        thr_sample = 4000,
                                        group_sel = "group")
qc_plots

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
ms_se_imp <- ms_se_prc$imp
ms_se_imp$group <- toupper(ms_se_imp$group)

dea_res <- wrapper_dea_gsea(
  ms_se = ms_se_imp,
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
saveRDS(ms_se_imp, save_here(object_name = "ms_se_imp.rds"))
