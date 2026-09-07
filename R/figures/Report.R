####=============== Rscript: unsupervised_analysis===============####
# Author:Dennis Friedel
# Date: 2024-06-03
# Modification: 2026-08-11
# 
# Create a potential diagnostic report for a sample analyzed by LC-MS 
#
####==============================================================####
analysis = paste0('/',Sys.Date()) # - Name of the analysis e.g Marker Identification
dataset = paste0('/Report/')
library("rip")
library("proteoLab")
library("SummarizedExperiment")
library("dplyr")
library("patchwork")
library("ggplot2")
calc_cv_fraction<-function(x,na.rm=T){
  sd(x,na.rm=T)/mean(x,na.rm=T)
}

##### Load DB ressources #####
# enzsub<-OmnipathR::enzyme_substrate(resources = "ProtMapper",cache = F)
resources<-c("ProtMapper","SIGNOR","PhosphoSite","PhosphoNetworks","phosphoELM","KEA")
enzsub<-data.table::fread("./data/omnipath_webservice_enz_sub.tsv")%>%as.data.frame()
enzsub[intersect(grep("EGFR",enzsub$enzyme_genesymbol),grep("EGFR",enzsub$substrate_genesymbol)),]

enzsub_phos<-enzsub[enzsub$modification=="phosphorylation",]
enzsub_dbs<-purrr::map(resources,function(x){
  tmp_enzsub<-enzsub_phos[grep(x,enzsub_phos$sources),]%>%as.data.frame()
  if(nrow(tmp_enzsub)<5){return(NULL)}
  df<-data.frame("p_site" = paste0(tmp_enzsub$substrate_genesymbol,"_",tmp_enzsub$residue_type,tmp_enzsub$residue_offset), 
             "enzyme_genesymbol"= tmp_enzsub$enzyme_genesymbol, 
             "mor"= 1, 
             "likelihood"= 1,db=x)
  distinct(df)
})
names(enzsub_dbs)<-resources
enzsub_dbs<-plyr::compact(enzsub_dbs)

##### Load samples #####
ptm_se_raw<-readRDS(".//output/20260512_RRS_1296_PCF_Phospho2/PTM_PRC_DEA/2026-08-06//ptm_se_raw.rds")
ptm_sse_raw<-ptm_se_raw[,1]
ptm_sse_raw<-ptm_sse_raw[!is.na(assay(ptm_sse_raw)),]
#### Single Kinase Score 
dc_kinacts<-purrr::map(enzsub_dbs,function(x){
  dc_kin_act <- dc_kin_activities(
    phospho_differential_analysis = ptm_sse_raw,
    dc_ptm = x,
    net_idenitifiers = "genes",
    inference_method = "fgsea"
  )
})
names(dc_kinacts)<-resources
