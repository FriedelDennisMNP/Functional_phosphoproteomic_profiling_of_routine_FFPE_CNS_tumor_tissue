####=============== Rscript: unsupervised_analysis===============####
# Author:Dennis Friedel
# Date: 2026-08-20

####===================================####
dataset  = '20260817_RR1296_PCF' # - Name of the dataset that is going to be analysed.
analysis = paste0('WP_PRC_DEA/',Sys.Date()) # - Name of the analysis e.g Marker Identification

library("dplyr")
library("patchwork")
library("ggplot2")
library("SummarizedExperiment")
source("./R/utils/rip_functions.R")
source("./R/utils/utils_module.R")
source("./R/utils/import_module.R")
source("./R/utils/run_pathview_in.R")

#### Import data #####
diann_result_path<-"/mnt/add50/PATHO-PROTEOMICS/bioinformatics/results/DIA/20260817_RR1296_PCF/diann_2.3.1/Lib_free/"
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

### Annoate
metadata<-openxlsx::read.xlsx("./projects/20260512_RRS_1296_PCF_Phospho2/sample_sheets/Annotated_sample_sheet3_20260707.xlsx",sheet = 1)
metadata$original_id<-gsub("_p_raw","",metadata$original_id)
rownames(metadata)<-metadata$sample_id

ms_se<-msdata$ggu_se
ms_se$sample_id<-plyr::mapvalues(colnames(ms_se),metadata$original_id,metadata$sample_id)
new_cd<-plyr::join(as.data.frame(colData(ms_se)),metadata,by = "sample_id")%>%DataFrame()
rownames(new_cd)<-new_cd$sample_id
colnames(ms_se)<-plyr::mapvalues(colnames(ms_se),metadata$original_id,metadata$sample_id)
colData(ms_se)<-new_cd

### Drop samples which are not amp or WT
ms_se_flt<-ms_se[,ms_se$group%in%c("WT","Amp")]

### Run Preproc and DEA
tt_comparison<-purrr::map(c(0.8,0.75),function(flt){
  
  ### Preprocess
  ms_se_prc <- proteoLab::pre_processing_wrapper(
    se_object = ms_se_flt,
    group_sel = "group",
    filter_fractioned = T,
    thr_sample = 5000,
    thr_feature = flt,
    normalization_method = "median_center",
    imputation_method  = "MinProb"
  )
  
  saveRDS(ms_se_prc,save_here(object_name = "ms_se_prc.rds"))
  
  qc_plots<-proteoLab::qp_plots_wrapper(se_proc_ls = ms_se_prc,
                                        thr_sample = 4000,
                                        group_sel = "group")
  qc_plots
  
  purrr::map2(qc_plots, names(qc_plots),function(x,y){
    ggplot2::ggsave(
      filename = save_here(object_name = paste0((y),"_qc_plot.pdf"),analysis = analysis),
      plot = x,
      width = 10,
      height = 10,
      units = "in",
      dpi = 300
    )
  })
  
  ### Differential Expression
  ms_se_imp<-ms_se_prc$imp
  ms_se_imp$group<-toupper(ms_se_imp$group)
  
  dea_res <- proteoLab::wrapper_dea_gsea(
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
  
  ttresult<-dea_res$tt_combined
  ttresult[ttresult$AMP_vs_WT_adj_P_Val<0.05,]
  
  dea_res$tt_combined<-NULL
  purrr::map2(dea_res, names(dea_res),function(x,y){
    ggplot2::ggsave(
      filename = save_here(object_name = paste0((y),"_dea_plot.pdf"),analysis = analysis),
      plot = x,
      width = 10,
      height = 10,
      units = "in",
      dpi = 300
    )
  })
  openxlsx::write.xlsx(ttresult,save_here(object_name = "TopTable_WT_AMP.xlsx"))
  saveRDS(ms_se_imp,save_here(object_name = "ms_se_imp.rds"))
  ttresult
})
tt80<-tt_comparison[[1]]
tt75<-tt_comparison[[2]]

dim(tt80);dim(tt75)

overlap_filter<-ggVennDiagram::ggVennDiagram(x = list("75%(3/5)"=tt75$genes[tt75$significant],
                                                      "80%(4/5)"=tt80$genes[tt80$significant]))+coord_flip()
ggsave(plot = overlap_filter,save_here(object_name = "Overlap_SigSites_75_vs_90.pdf"),width = 10,height = 10)

only_80<-tt80[which(!tt80$genes[tt80$significant]%in%tt75$genes[tt75$significant]),]
only_75<-tt75[which(!tt75$genes[tt75$significant]%in%tt80$genes[tt80$significant]),]
intersect_75_90<-tt75[which(tt75$genes[tt75$significant]%in%tt80$genes[tt80$significant]),]

tt75[c("EGFR","NES","VIM","TP53BP1","PSIP1","MAP2"),]

