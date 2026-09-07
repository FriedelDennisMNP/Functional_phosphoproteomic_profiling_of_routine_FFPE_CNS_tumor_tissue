dataset  = '20260512_RRS_1296_PCF_Phospho2' # - Name of the dataset that is going to be analysed.
analysis = paste0('PTM_Decoupler_Kinase_Infernce/',Sys.Date()) # - Name of the analysis e.g Marker Identification

library("rip")
library("proteoLab")
library("dplyr")
library("patchwork")
library("ggplot2")
library("ComplexHeatmap")
source("./R/utils/run_pathview_in.R")

###### Kinase Inference #####
ptm_se_imp<-readRDS("./output/20260512_RRS_1296_PCF_Phospho2/PTM_PRC_DEA/2026-07-31/ptm_se_imp.rds")
ttresult<-openxlsx::read.xlsx("./output/20260512_RRS_1296_PCF_Phospho2/PTM_PRC_DEA/2026-07-31/TopTable_WT_AMP.xlsx")
rownames(ttresult)<-ttresult$genes
# 
# ### Load Enzyme Substarte-databases ####
resources<-OmnipathR::enzsub_resources()
resources<-c("ProtMapper","SIGNOR","PhosphoSite","PhosphoNetworks","phosphoELM","KEA")

### Instead of using OmnipathR use the latest omnipath_webservice_enz_sub which contains EGFR-EGFR interactions
# enzsub<-OmnipathR::enzyme_substrate(resources = "ProtMapper",cache = F)
enzsub<-data.table::fread("./data/omnipath_webservice_enz_sub.tsv")%>%as.data.frame()
enzsub[intersect(grep("EGFR",enzsub$enzyme_genesymbol),grep("EGFR",enzsub$substrate_genesymbol)),]

subset(enzsub,enzyme=="P00533"&substrate=="P00533")
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

###### DEA evaluation ######
# Get the kinases that use to phosphorylate the EGFR Site 
EGFR_candidates<-ttresult[grep("EGFR",rownames(ttresult)),] %>%
  .[.$significant,]

EGFR_candidates_enzymes<-purrr::map(rownames(EGFR_candidates),function(x){
  db_match<-purrr::map(enzsub_dbs,function(y){
    tmp<-y[y$p_site %in% x, ]$enzyme_genesymbol
  })%>%unlist()%>%table()%>%as.matrix()%>%t()
  db_match<-as.data.frame(db_match)
  db_match$residue<-x
  db_match
})%>%do.call(plyr::rbind.fill,.)%>%cbind(EGFR_candidates,.)
EGFR_candidates_enzymes$residue<-NULL
EGFR_candidates_enzymes
openxlsx::write.xlsx(EGFR_candidates_enzymes,save_here(object_name = "EGFR_candidates_enzymes.xlsx"))


###### EGFR - Kinase Substrate Heatmap ######
cands<-colnames(EGFR_candidates_enzymes)[11:ncol(EGFR_candidates_enzymes)]
EGFR_candidates_enzymes_mat<-EGFR_candidates_enzymes[,cands]
EGFR_candidates_enzymes_mat[is.na(EGFR_candidates_enzymes_mat)]<-0
id_col <- circlize::colorRamp2(
  c(0, 7),
  c("#FFFFFF", "#67000D")
)
EGFR_kin_sub_heatmap<-ComplexHeatmap::Heatmap(
  t(EGFR_candidates_enzymes_mat),
  col = id_col,name = "ocurrence",
  border = T,row_names_side = "left",
  show_row_dend = F,show_column_dend = F,
  row_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
  column_names_gp = grid::gpar(fontsize = 8,fontface = "bold"),
  rect_gp = grid::gpar(col = "grey20", lwd = 0.5),
  heatmap_width = grid::unit(10,"cm"),
  heatmap_legend_param = list(
    title = "Count",
    title_position = "leftcenter-rot",
    title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 12),
    legend_height = grid::unit(10, "cm"))
)%>%complexheatmap_to_ggplot(.)
EGFR_kin_sub_heatmap
ggsave(plot = EGFR_kin_sub_heatmap,
       save_here(object_name = "EGFR_kin_sub_heatmap.pdf"),
       width = 10,
       height = 10)

#### Run DEA Decoupler ####
dea_dc_kinacts<-purrr::map(enzsub_dbs,function(x){
  sdc_kin_act <- dc_kin_activities(
    phospho_differential_analysis = ttresult,
    dc_ptm = x,
    net_idenitifiers = "genes",
    inference_method = "viper"
  )
})
names(dea_dc_kinacts)<-names(enzsub_dbs)
decoupleR::run_viper

#### Signifikant Active Kinase Barplot ####
dc_barplots <- purrr::map2(dea_dc_kinacts, names(dea_dc_kinacts), function(x, y) {
  sig_dc_acts<-x[x$p_value<0.05,]
  make_dc_barplot(
    contrast_acts = sig_dc_acts,
    n_kin = 30,
    color_gradient = c("#08306B", "#FFFFFF", "#67000D"),
    facet_by_condition = FALSE
  ) + ggtitle(y)
})
openxlsx::write.xlsx(dc_barplots[[1]]@data,file = save_here(object_name = "Kinases_diff_EGFR_amp.xlsx"))

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

#### Make Conensus ####
sig_kins<-do.call(rbind,sig_kin_list)
sig_kins$significant<-sig_kins$p_value<0.05

kinsensus<-t(as.matrix(table(sig_kins$source,sig_kins$significant)/length(dea_dc_kinacts)))*100
kinsensus<-kinsensus[2,]
kinsensus<-kinsensus[kinsensus>0]

hm_col <- circlize::colorRamp2(
  c(0, 100),
  c("#FFFFFF", "#67000D")
)
cmp <- ComplexHeatmap::Heatmap((kinsensus),
                               name = paste0("Occurence in percent"),
                               heatmap_width = grid::unit(6, "cm"),
                               col = hm_col,
                               border = T,
                               show_row_dend = F,
                               show_column_dend = F,
                               show_column_names = F,
                               rect_gp = grid::gpar(col = "grey30", lwd = 0.5)
)
ggcmp<-complexheatmap_to_ggplot(heatmap = cmp)
ggsave(plot = ggcmp,save_here(object_name = "Sig_kinase_consensus.pdf"),width = 5,height = 7)

####´Check with Pathwview
library(pathview)
gene_vector<-sig_dc_acts$source

kinase_scores<-dea_dc_kinacts$ProtMapper$score
names(kinase_scores)<-dea_dc_kinacts$ProtMapper$source

pathids<-c("05214",# GLIOMA 
           "04012",# ERB Signaling
           "04630",#"JAK-STAT-KEGG pathway"
           "04151",# PI3K-AKT
           "04510",# MTOR
           "04110",#CellCycle
           "04510"#Focal Adhesion (SRC)
)
purrr::map(pathids,function(x){
  run_pathview_in(save_here(object_name = "pathview/ProtMapper"),
                  gene.data = kinase_scores,
                  pathway.id= x, # JAK-STAT-KEGG pathway
                  species="hsa",
                  gene.idtype = "SYMBOL",
                  limit=list(gene=max(abs(kinase_scores)))
  )
  
})

#### Single Sample Evaluation #####
#### Decoupler
dc_kinacts<-purrr::map(enzsub_dbs,function(x){
  dc_kin_act <- dc_kin_activities(
    phospho_differential_analysis = ptm_se_imp,
    dc_ptm = x,
    net_idenitifiers = "genes",
    inference_method = "viper"
  )
})
names(dc_kinacts)<-resources


# dc_kin_act<-dc_kinacts[[1]]
# y<-names(dc_kinacts)[[1]]
#### With all Kinases found
y<-names(dc_kinacts)[1]
dc_kin_act<-dc_kinacts$ProtMapper

####---------Supplemental table 3 ------------###
openxlsx::write.xlsx("",save_here(object_name = "Supplemental_table3.xlsx"))
openxlsx::createWorkbook(save_here(object_name = paste0("Supplemental_table3.xlsx")))
wb<-openxlsx::loadWorkbook(save_here(object_name = paste0("Supplemental_table3.xlsx")))

openxlsx::addWorksheet(wb = wb,sheetName = "Differnital_active_kinases")
openxlsx::writeData(wb,sheet="Differnital_active_kinases",x=dea_dc_kinacts$ProtMapper)

openxlsx::addWorksheet(wb = wb,sheetName = "Single_sample_kinact")
openxlsx::writeData(wb,sheet="Single_sample_kinact",x=dc_kin_act)

openxlsx::saveWorkbook(wb,save_here(object_name = paste0("Supplemental_table3.xlsx")),overwrite = T)

####---------Plot single sample score ------------###
hms<-purrr::map2(dc_kinacts,names(dc_kinacts),function(dc_kin_act,y){
  library(reshape2)
  score_matrix <- acast(dc_kin_act, condition ~ source, value.var = "score")%>%t()
  pvalue_matrix <- acast(dc_kin_act, condition ~ source, value.var = "p_value")%>%t()
  
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
  
  size_legend <- ComplexHeatmap::Legend(
    title = "-log10(p-value)",
    type = "points",
    at = round(seq(0, max(-log10(pvalue_matrix), na.rm = TRUE), length.out = 4), 2),
    legend_gp = gpar(fill = "gray50"),
    size = unit(seq(0.05, 0.35, length.out = 4), "cm"),
    pch = 21,
    border = "black"
  )
  
  hm<-ComplexHeatmap::Heatmap(
    score_matrix,
    column_title = y,
    column_split = plyr::mapvalues(colnames(score_matrix),ptm_se_imp$sample_id,ptm_se_imp$group),
    show_row_dend = F,show_column_dend = F,
    top_annotation = topanno,
    clustering_method_columns = "ward.D2",
    clustering_distance_columns =  "euclidean",
    name = "GSEA Score",
    row_names_side =  "left",
    row_names_gp = grid::gpar(fontsize = 10, fontface = "bold",angle=45),
    cluster_columns = T,cluster_rows=F,
    column_title_side =  "top",
    column_names_gp = grid::gpar(fontsize = 7, fontface = "bold"),
    col = hm_col,
    border = T,
    heatmap_width = grid::unit(12, "cm"),
    rect_gp = grid::gpar(col = "grey60", lwd = 0.5),
    heatmap_legend_param = list(
      title = "GSEA Score",
      title_position = "leftcenter-rot",
      title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 12),
      legend_height = grid::unit(20, "cm"))
  )
  final_hm<-draw(hm)
  final_hm
})

pdf(save_here(object_name = "Decoupler_Viper_Single_Sample_Kinase_Activities.pdf",analysis = analysis),
    width = 12,
    height = 15)
hms
dev.off()

names(sig_kin_list)<-names(dea_dc_kinacts)

### Only the signifikant kinases
hms<-purrr::map(names(dc_kinacts),function(y){
  library(reshape2)
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
  final_hm<-draw(hm)
  complexheatmap_to_ggplot(final_hm)
})

ggsave(
  plot = hms,
  save_here(object_name = "Decoupler_Viper_Significant_Kinases_Activities.pdf", analysis = analysis),
  width = 12,
  height = 10
)

ggsave(
  plot = hms[[1]],
  save_here(object_name = "Protmapper_Significant_Kinases_Activities.pdf", analysis = analysis),
  width = 12,
  height = 10
)


