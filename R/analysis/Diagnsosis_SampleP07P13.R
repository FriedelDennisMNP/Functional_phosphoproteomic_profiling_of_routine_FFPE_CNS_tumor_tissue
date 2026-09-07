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
analysis = paste0('/Diagnosis/',Sys.Date()) # - Name of the analysis e.g Marker Identification
dataset = paste0('/')
library("rip")
library("proteoLab")
library("SummarizedExperiment")
library("dplyr")
library("patchwork")
library("ggplot2")
source("/mnt/NAS4/user_data/np-dennis/projects/Functions/Helpful_functions/nphd_functions/nphd_get_latest_braintumor_classification.R")
classifier_path<-"/mnt/NAS4/methylation/results/"
calc_cv_fraction<-function(x,na.rm=T){
  sd(x,na.rm=T)/mean(x,na.rm=T)
}

### Set colors ###
layers <- c(
  'EGFR-amplified' = "#F77576",
  'Non-amplified' = "#91BFD1"
)


##### Load PTM data ###
proc_ms<-readRDS("./output/20260813_RRS_1296_PCF_Phospho//PTM_PRC_DEA/2026-08-19///ptm_se_prc.rds")
clin_anno<-openxlsx::read.xlsx("./data/clinanno.xlsx")
clin_anno$sample_id<-clin_anno$P_ID

ms_imp<-proc_ms$imp
colData(ms_imp)<-plyr::join(as.data.frame(colData(ms_imp)),as.data.frame(clin_anno),by = "sample_id")%>%DataFrame()
ids_flt<-get_latest_classifer_result(anno_table = as.data.frame(colData(ms_imp)),txt_idat = ms_imp$Sentrix.ID)

classifier_results  =  purrr::map(ids_flt$Sentrix.ID[ids_flt$sample_id%in%c("P07","P09","P13")], function(x) {
  classifer_file  =
    paste0(classifier_path,
           gsub("_.*", "", x),
           "/",
           x,
           "/",
           classifier_version,
           "/",
           x,
           "_scores_cal.csv")
  
  tmp_csv  =  data.table::fread(classifer_file) %>% as.data.frame()
  tmp_csv  =  tmp_csv[order(tmp_csv[, 2], decreasing = T), ]
  colnames(tmp_csv)  =  c(classifier_version, "scores_cal")
  head(tmp_csv)
})
classifier_results

###### get EGFR scores
ms_imp$cnv_gene_detail<-paste0("/mnt/NAS4/methylation/results/",gsub("_.*","",ms_imp$Sentrix.ID),
                                                "/",ms_imp$Sentrix.ID,"/cnvp_v5.4/",ms_imp$Sentrix.ID,".detail.txt")
gene_cnb_details<-purrr::map(ms_imp$sample_id,function(x){
  tmp_res<-data.table::fread(ms_imp$cnv_gene_detail[ms_imp$sample_id==x]) %>% as.data.frame()
  tmp_res<-data.frame("value"=tmp_res$value,row.names = tmp_res$name)
  colnames(tmp_res)<-x
  as.data.frame(t(tmp_res))
})%>%do.call(plyr::rbind.fill,.)
rownames(gene_cnb_details)<-ms_imp$sample_id

topanno<-ComplexHeatmap::HeatmapAnnotation("EGFR"=ms_imp$EGFR,"Class"=ids_flt$predictBrain_v12.8)
ComplexHeatmap::Heatmap(cor(t(gene_cnb_details)),top_annotation = topanno)
ComplexHeatmap::Heatmap(t(gene_cnb_details),
                        top_annotation = topanno,
                        column_split = ms_imp$EGFR)

