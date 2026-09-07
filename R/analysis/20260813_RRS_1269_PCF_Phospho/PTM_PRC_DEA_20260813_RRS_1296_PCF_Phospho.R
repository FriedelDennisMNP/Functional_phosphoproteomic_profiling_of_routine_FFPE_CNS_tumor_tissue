####=============== Rscript: unsupervised_analysis===============####
# Author:Dennis Friedel
# Date: 2024-06-03
# Modification: 2024-06-04
# Description: Unsupervised analysis of beta values of the sample cohort.
# Major aim quality control prior protoemic analyis, idenitifcation of membership of samples and segregation 
# into right cns methylation classes. 
# Detail:
#'
####===================================####
dataset  = '20260813_RRS_1296_PCF_Phospho' # - Name of the dataset that is going to be analysed.
analysis = paste0('PTM_PRC_DEA/',Sys.Date()) # - Name of the analysis e.g Marker Identification
"PTM_PRC_DEA_20260813_RRS_1296_PCF_Phospho.R"
library("rip")
library("proteoLab")
library("dplyr")
library("patchwork")
library("ggplot2")
library("SummarizedExperiment")
source("./R/utils/run_pathview_in.R")

### Load data
mq_results<-proteoLab::load_ms_results(ms_result_dir = "/mnt/add50/PATHO-PROTEOMICS/bioinformatics/results/DDA/20260813_RRS_1296_PCF_Phospho/MaxQuant2.4.2.0/")
mq_list<-proteoLab::create_maxquant_se_list(mq_results = mq_results,ptm_probability_flt = 0.75,remove_low_quality = T)

### Annoate
metadata<-openxlsx::read.xlsx("./projects/20260512_RRS_1296_PCF_Phospho2/sample_sheets/Annotated_sample_sheet3_20260707.xlsx",sheet = 1)
rownames(metadata)<-metadata$sample_id

ptm_se<-mq_list$ptm_se
ptm_se<-ptm_se[!duplicated(ptm_se@elementMetadata$Name),]
rownames(ptm_se)<-ptm_se@elementMetadata$Name

ptm_se<-ptm_se[,ptm_se$sample_id%in%metadata$original_id]
colnames(ptm_se)<-metadata$sample_id
SummarizedExperiment::colData(ptm_se)<-S4Vectors::DataFrame(metadata[match(ptm_se$sample_id,metadata$original_id),])

ptm_se_flt<-ptm_se[,ptm_se$group%in%c("Amp","WT")]
saveRDS(ptm_se_flt,save_here(object_name = "ptm_se_raw.rds"))

test<-colSums(!is.na(assay(ptm_se_flt)))-colSums(!is.na(assay(readRDS("./output/20260512_RRS_1296_PCF_Phospho2/PTM_PRC_DEA/2026-08-17/ptm_se_raw.rds"))))
mean((test/colSums(!is.na(assay(ptm_se_flt))))*100)

ptm_se_flt$sample_id[ptm_se_flt$EGFR=="WT"]
#ptm_se_flt<-ptm_se_flt[,ptm_se_flt$sample_id!="P13"]

tt_comparison<-purrr::map(c(0.8,0.75),function(flt){
  
  ### Preprocess
  ptm_se_prc <- proteoLab::pre_processing_wrapper(
    se_object = ptm_se_flt,
    group_sel = "group",
    filter_fractioned = T,
    thr_sample = 3500,
    thr_feature = flt,
    normalization_method = "median_center",
    imputation_method  = "MinProb"
  )
  
  saveRDS(ptm_se_prc,save_here(object_name = "ptm_se_prc.rds"))
  
  qc_plots<-proteoLab::qp_plots_wrapper(se_proc_ls = ptm_se_prc,
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
  ptm_se_imp<-ptm_se_prc$imp
  ptm_se_imp$group<-toupper(ptm_se_imp$group)
  
  dea_res <- proteoLab::wrapper_dea_gsea(
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
  saveRDS(ptm_se_imp,save_here(object_name = "ptm_se_imp.rds"))
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

tt75[c("EGFR_Y1110","EGFR_Y1197","EGFR_T693","NES_S1492","VIM_S7","VIM_S72","VIM_S205","TP53BP1_S1101","PSIP1_T272","MAP2_S1795","MAP2_S725"),]
tt80[c("EGFR_Y1110","EGFR_Y1197","EGFR_T693","NES_S1492","VIM_S7","VIM_S72","VIM_S205","TP53BP1_S1101","PSIP1_T272","MAP2_S1795","MAP2_S725"),]
