####===================================####
# Author: Dennis Friedel, PhD
# Date: 2026-09-07
# Bioinformatician,
# Department of Neuropahtology, University Clinic Heidelberg
####===================================####

#######----------- 0.Load Packages and Functions -----------#####
set.seed(2905)
library("purrr")
library("dplyr")
library("ggplot2")
library("patchwork")
library("reshape2")
#library(pathview)
library("ComplexHeatmap")
library("SummarizedExperiment")
source("./R/utils/rip_functions.R")
source("./R/utils/utils_module.R")
source("./R/utils/import_module.R")
source("./R/utils/preprocess_module.R")
source("./R/utils/compare_module.R")
source("./R/utils/ptm_module.R")

# Set variables for save_here
dataset  = '20260813_RRS_1296_PCF_Phospho' # - Name of the dataset that is going to be analysed.
analysis = paste0('PTM_Decoupler_Kinase_Infernce/', Sys.Date()) # - Name of the analysis e.g Marker Identification


#######----------- 1.Load Data from previous results for Kinase Inference --#####
ptm_se_imp <- readRDS(
  "./output/20260813_RRS_1296_PCF_Phospho//PTM_PRC_ADJ_DEA/2026-09-09/ptm_se_adjprc.rds")

ttresult <- openxlsx::read.xlsx(
  "./output/20260813_RRS_1296_PCF_Phospho/PTM_PRC_ADJ_DEA/2026-09-09/TopTable_protadjusted_WT_AMP.xlsx"
)
rownames(ttresult) <- ttresult$genes

####----------- 2.Load Enzyme Substrate-databases -----------####
resources <- c("ProtMapper",
               "SIGNOR",
               "PhosphoSite",
               "PhosphoNetworks",
               "phosphoELM",
               "KEA")

### Use the latest omnipath which contains EGFR-EGFR interactions
enzsub <- data.table::fread("./data/omnipath_webservice_enz_sub.tsv") %>%
  as.data.frame()
enzsub[intersect(grep("EGFR", enzsub$enzyme_genesymbol),
                 grep("EGFR", enzsub$substrate_genesymbol)), ]
# Sanity check
subset(enzsub, enzyme == "P00533" & substrate == "P00533")

# Only use Phosphorylation as moidification
enzsub_phos <- enzsub[enzsub$modification == "phosphorylation", ]
enzsub_dbs <- purrr::map(resources, function(x) {
  tmp_enzsub <- enzsub_phos[grep(x, enzsub_phos$sources), ] %>% as.data.frame()
  if (nrow(tmp_enzsub) < 5) {
    return(NULL)
  }
  df <- data.frame(
    "p_site" = paste0(
      tmp_enzsub$substrate_genesymbol,
      "_",
      tmp_enzsub$residue_type,
      tmp_enzsub$residue_offset
    ),
    "enzyme_genesymbol" = tmp_enzsub$enzyme_genesymbol,
    "mor" = 1,
    "likelihood" = 1,
    db = x
  )
  distinct(df)
})
names(enzsub_dbs)<-resources
enzsub_dbs<-plyr::compact(enzsub_dbs)

####----------- 3.Differential active Kinases  -----------####
dea_dc_kinacts<-purrr::map(enzsub_dbs,function(x){
  sdc_kin_act <- dc_kin_activities(
    phospho_differential_analysis = ttresult,
    dc_ptm = x,
    net_idenitifiers = "genes",
    inference_method = "viper"
  )
})
names(dea_dc_kinacts)<-names(enzsub_dbs)

#### Signifikant Active Kinase Barplot 
dc_barplots <- purrr::map2(dea_dc_kinacts, names(dea_dc_kinacts), function(x, y) {
  sig_dc_acts<-x[x$p_value<0.05,]
  make_dc_barplot(
    contrast_acts = sig_dc_acts,
    n_kin = 30,
    color_gradient = c("#08306B", "#FFFFFF", "#67000D"),
    facet_by_condition = FALSE
  ) + ggtitle(y)
})

ggsave(
  dc_barplots[[1]],
  file = save_here(object_name = "Kinases_protmapper_diff_EGFR_amp.pdf"),
  width = 7,
  height = 7
)
openxlsx::write.xlsx(dc_barplots[[1]]@data,
                     file = save_here(object_name = "Kinases_diff_EGFR_amp.xlsx"))
ggplot2::ggsave(
  filename = save_here(object_name = paste0("AMP_vs_WT_dc_kinase_viper_inference.pdf"),analysis = analysis),
  plot = dc_barplots,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

#### Save tables 
openxlsx::write.xlsx("",save_here(object_name = paste0("Kinase_DEA_significant.xlsx")))
openxlsx::createWorkbook(save_here(object_name = paste0("Kinase_DEA_significant.xlsx")))
sig_kin_list<-purrr::map(names(dea_dc_kinacts),function(x){
  dc_kinact<-dea_dc_kinacts[[x]]
  wb<-openxlsx::loadWorkbook(save_here(object_name = paste0("Kinase_DEA_significant.xlsx")))
  sig_dc_acts<-dc_kinact[dc_kinact$p_value<0.05,]
  openxlsx::addWorksheet(wb = wb,x)
  openxlsx::writeData(wb,sheet=x,x=sig_dc_acts)
  openxlsx::saveWorkbook(wb,save_here(object_name = paste0("Kinase_DEA_significant.xlsx")),overwrite = T)
  dc_kinact$db<-x
  dc_kinact
})
names(sig_kin_list)<-names(dea_dc_kinacts)

####----------- 4. Single Sample Active Kinases  -----------------------####

#### Estimate activities of significant Kinases in single sample context 
dc_kinacts<-purrr::map(enzsub_dbs,function(x){
  dc_kin_act <- dc_kin_activities(
    phospho_differential_analysis = ptm_se_imp,
    dc_ptm = x,
    net_idenitifiers = "genes",
    inference_method = "viper"
  )
})
names(dc_kinacts)<-resources

### Make Heatmaps of significant kinases 
hms<-purrr::map(names(dc_kinacts),function(y){
  
  score_matrix <- acast(dc_kinacts[[y]], condition ~ source, value.var = "score")%>%t()
  kinase<-sig_kin_list[[y]]$source[sig_kin_list[[y]]$p_value<0.05]
  
  score_matrix<-score_matrix[kinase,]
  hm_lim <- max(abs(range(score_matrix, na.rm = TRUE)))
  hm_col <- circlize::colorRamp2(
    c(-hm_lim, 0, hm_lim),
    c("#08306B", "#FFFFFF", "#67000D")
  )
  
  ### make large legend
  topanno<-ComplexHeatmap::HeatmapAnnotation(
    df = data.frame(Condition = plyr::mapvalues(colnames(score_matrix),ptm_se_imp$sample_id,ptm_se_imp$group)),
    col = list(Condition = c("Amp" = "forestgreen", "WT" = "#377EB8")), 
    annotation_legend_param = list(title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
                                   labels_gp = grid::gpar(fontsize = 12),
                                   legend_height = grid::unit(20, "cm")))
  if(nrow(score_matrix)>5){
    nrow_km = 2  
  }else{
    nrow_km=NULL
  }
  
  
  hm<-ComplexHeatmap::Heatmap(
    score_matrix,
    column_title = y,
    column_split = plyr::mapvalues(colnames(score_matrix),ptm_se_imp$sample_id,ptm_se_imp$group),
    show_row_dend = F,show_column_dend = F,
    top_annotation = topanno,
    #clustering_method_columns = "ward.D2",
    #clustering_distance_columns =  "euclidean",
    name = "Kinase Score",
    row_names_side =  "left",
    row_names_gp = grid::gpar(fontsize = 10, fontface = "bold",angle=45),
    cluster_columns = T,cluster_rows=T,
    column_title_side =  "top",
    column_names_gp = grid::gpar(fontsize = 7, fontface = "bold"),
    col = hm_col,
    border = T,
    heatmap_width = grid::unit(12, "cm"),
    rect_gp = grid::gpar(col = "grey60", lwd = 0.5),
    heatmap_legend_param = list(
      title = "Kinase Score",
      title_position = "leftcenter-rot",
      title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 12),
      legend_height = grid::unit(20, "cm"))
  )
  complexheatmap_to_ggplot(hm)
})

ggsave(
  plot = hms[[1]],
  save_here(
    object_name = "Figure3c_Protmapper_Significant_Kinases_Activities.pdf",
    analysis = paste0('/'),
    dataset_name = paste0('/Figures/')
  ),
  width = 12,
  height = 10
)

####----- 6. Single Sample Protmapper result (Supplemental Table 3) --------####
y<-names(dc_kinacts)[[1]]
message("Use ",y)

dc_kin_act<-dc_kinacts[[y]]
dea_dc_act<-dea_dc_kinacts[[y]]

pm<-enzsub_dbs$ProtMapper[enzsub_dbs$ProtMapper$p_site%in%rownames(ptm_se_imp)&
                            enzsub_dbs$ProtMapper$enzyme_genesymbol%in%dea_dc_act$source,]
dea_dc_act$supportive_sites<-plyr::mapvalues(dea_dc_act$source,names(table(pm$enzyme_genesymbol)),table(pm$enzyme_genesymbol))

xlsx_path<-save_here(object_name = "Supplemental_table3.xlsx",dataset_name = paste0('/Supplemental_table/'),analysis_name = "")

openxlsx::write.xlsx("",xlsx_path);openxlsx::createWorkbook(xlsx_path)
wb<-openxlsx::loadWorkbook(xlsx_path)

openxlsx::addWorksheet(wb = wb,sheetName = "Differnital_active_kinases")
openxlsx::writeData(wb,sheet="Differnital_active_kinases",x=dea_dc_act)

openxlsx::addWorksheet(wb = wb,sheetName = "Single_sample_kinact")
openxlsx::writeData(wb,sheet="Single_sample_kinact",x=dc_kin_act)

openxlsx::saveWorkbook(wb,xlsx_path,overwrite = T)

