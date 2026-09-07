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
dataset  = '20260512_RRS_1296_PCF' # - Name of the dataset that is going to be analysed.
analysis = paste0('Correlation_WHOP_PTM/',Sys.Date()) # - Name of the analysis e.g Marker Identification

library("rip")
library("proteoLab")
library("dplyr")
library("patchwork")
library("ggplot2")
source("./R/utils/run_pathview_in.R")

ptm_se<-readRDS("./output/20260512_RRS_1296_PCF_Phospho2/PTM_PRC_DEA/2026-08-17/ptm_se_prc.rds")
wp_se<-readRDS("./output/20260512_RRS_1296_PCF/WHOP_PRC_DEA/2026-08-17/ms_se_prc.rds")

genes<-ptm_se$unfilt@elementMetadata$gene
wp_se$unfilt

ptm_prots<-gsub("_.*","",rownames(ptm_se$imp))
whop_prots<-rownames(wp_se$imp)
intprots<-intersect(whop_prots,ptm_prots)
x<-intprots[64]

corres<-purrr::map(intprots,function(x){
  
  x_prot<-assay(wp_se$imp)[rownames(wp_se$imp)==x,]
  x_ptm<-assay(ptm_se$imp)[grep(x,rownames(ptm_se$imp)),]
  
  if(is.null(nrow(x_ptm))){
    cors<-cor(x_ptm,x_prot,use="na.or.complete")  
    data.frame(protein=x,site=rownames(ptm_se$imp)[grep(x,rownames(ptm_se$imp))],correlation=cors)
  }else{
    cors<-purrr::map(1:nrow(x_ptm),function(y){
      cor(x_ptm[y,],x_prot,use="na.or.complete")  
    })%>%unlist()
    data.frame(protein=x,site=rownames(x_ptm),correlation=cors)
  }
},.progress = T)%>%do.call(rbind,.)

median(corres$correlation)
mean(corres$correlation)
openxlsx::write.xlsx(corres,save_here(object_name = "correraltion_df.xlsx"))

