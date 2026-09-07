####=============== Rscript: unsupervised_analysis===============####
# Author:Dennis Friedel
# Date: 2024-06-03
# Modification: 2026-08-11
# 
# Create a potential diagnostic report for a sample analyzed by LC-MS 
#
####==============================================================####
analysis = paste0('/',Sys.Date()) # - Name of the analysis e.g Marker Identification
dataset = paste0('/SingleSampleReport/')
library("rip")
library("proteoLab")
library("SummarizedExperiment")
library("dplyr")
library("patchwork")
library("ggplot2")
calc_cv_fraction<-function(x,na.rm=T){
  sd(x,na.rm=T)/mean(x,na.rm=T)
}


##### Load samples #####
ptm_se_raw<-readRDS(".//output/20260512_RRS_1296_PCF_Phospho2/PTM_PRC_DEA/2026-08-06//ptm_se_raw.rds")
ptm_sse_raw<-ptm_se_raw[,1]
ptm_sse_raw<-ptm_sse_raw[!is.na(assay(ptm_sse_raw)),]
ptm_sse_raw@elementMetadata$log2intensity<-assay(ptm_sse_raw)[,1]
ptmEWM<-ptm_sse_raw@elementMetadata

##### References #####
tt_WT_AMP<-openxlsx::read.xlsx(".//output/20260512_RRS_1296_PCF_Phospho2/PTM_PRC_DEA/2026-08-11//TopTable_WT_AMP.xlsx")
tt_WT_AMP<-tt_WT_AMP[tt_WT_AMP$AMP_vs_WT_significant,]
WT_candidates<-tt_WT_AMP[order(tt_WT_AMP$AMP_vs_WT_logFC),]%>%head(.,10)%>%.$genes
AMP_candidates<-tt_WT_AMP[order(tt_WT_AMP$AMP_vs_WT_logFC),]%>%tail(.,10)%>%.$genes

###### 1. Quality of Single Sample #####

###### 1.1 Summary of identified sites & their intensity #####
residue<-ptmEWM$amino_acid%>%table(.)%>%as.data.frame()
colnames(residue)<-c("AminoAcid","Freq")

phosphosite_occurence_summary<-ggpubr::ggdonutchart(residue,x = "Freq",fill = "AminoAcid",color="grey50",alpha=.7)+
  scale_fill_manual(name = paste0("Number and type of identified \nphosphorylated sites (n=",sum(residue$Freq),")"),
                    values = c("S"="dodgerblue4","T"="coral4","Y"="#B78ECB"),
                    labels = c("S" = paste0("S (n=",residue$Freq[residue$AminoAcid=="S"],")"),
                               "T" = paste0("T (n=",residue$Freq[residue$AminoAcid=="T"],")"),
                               "Y" = paste0("Y (n=",residue$Freq[residue$AminoAcid=="Y"],")")))+
  ggtitle(label = "Phosphosite occurence summary")+
  theme(legend.position = "right",plot.title = element_text(color="darkblue",face = "bold"))


intensity_distribution <-
  ggplot2::ggplot(ptmEWM, ggplot2::aes(x = amino_acid,y = log2intensity, fill = amino_acid)) +
  ggplot2::scale_fill_manual(name = "Amino Acid",values= c("S"="dodgerblue4","T"="coral4","Y"="#B78ECB"))+
  ggplot2::geom_violin(trim = F,alpha=.7) +
  ggplot2::geom_boxplot(width=0.1,outlier.shape = NA,fill="white") +
  labs(y = "log2 Intensity", x = "Amnio Acid") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text =  element_text(face = "bold"),
    axis.text =  element_text(angle = 45, hjust = 1)
  )
phosphosite_occurence_summary+intensity_distribution

##### 2. Association with characterized groups ####

###### 2.1 Rank of identified sites in toptable of groups #####
library(ggpubr)
ptmEWM$zscore<-scale(ptmEWM$log2intensity,center = T,scale = T)
ptmEWM_ordered<-ptmEWM[order(ptmEWM$zscore,decreasing=T),]
ptmEWM_ordered$rank<-1:nrow(ptmEWM_ordered)

groupptable<-ptmEWM_ordered[ptmEWM_ordered$Name%in%c(AMP_candidates,WT_candidates),c("gene","Name","zscore","rank")]%>%as.data.frame()
groupptable<-groupptable[order(groupptable$rank),]
colnames(groupptable)<-c("Protein","Residue","Score","Rank")
groupptable$Site<-gsub(".*_","",groupptable$Residue)
groupptable$Score<-round(groupptable$Score,digits = 2)

groupptable$EGFR_status<-ifelse(groupptable$Residue%in%AMP_candidates,"EGFR-amplified","Non_amplified")
RAmp<-which(groupptable$EGFR_status=="EGFR-amplified")+1
RWT<-which(groupptable$EGFR_status!="EGFR-amplified")+1


txttbl<-ggpubr::ggtexttable(
  groupptable[,c("Protein","Site","Score","Rank")],
  rows = NULL,
  theme = ttheme(
    tbody.style = tbody_style(fill = "white"),
    colnames.style = colnames_style(
      face = "bold",
      color = "darkblue",
      fill = "white"
    )
  ),
) %>% tab_add_hline(at.row =   1)
para<-ggparagraph("Rank of group specific Phosphosites (EGFR Amplified (red) vs Non-amplified (blue) )",face = "bold",color = "darkblue",size=14)
txttbl_colored<-txttbl %>% table_cell_font(row=RAmp,column = 2,color = "firebrick4")%>%
  table_cell_font(row=RWT,column = 2,color = "dodgerblue4")
toptable<-ggarrange(para,txttbl_colored,ncol=1,heights=c(0.2,1))

###### 2.2 Summarized Score ######
WT_candidates<-tt_WT_AMP[order(tt_WT_AMP$AMP_vs_WT_logFC),]%>%head(.,50)%>%.$genes
AMP_candidates<-tt_WT_AMP[order(tt_WT_AMP$AMP_vs_WT_logFC),]%>%tail(.,50)%>%.$genes

ptmEWM_topfilt<-ptmEWM[ptmEWM$Name%in%c(WT_candidates,AMP_candidates),]
perc_WT<-sum(WT_candidates%in%ptmEWM_topfilt$Name)/50
perc_AMp<-sum(AMP_candidates%in%ptmEWM_topfilt$Name)/50

ptmEWM_topfilt$group<-ifelse(ptmEWM_topfilt$Name%in%WT_candidates,"Non-amplified","EGFR amplified")

ptmEWM_topfilt$log2intensity[ptmEWM_topfilt$group=="EGFR amplified"]<-ptmEWM_topfilt$log2intensity[ptmEWM_topfilt$group=="EGFR amplified"]*perc_AMp
ptmEWM_topfilt$log2intensity[ptmEWM_topfilt$group!="EGFR amplified"]<-ptmEWM_topfilt$log2intensity[ptmEWM_topfilt$group!="EGFR amplified"]*perc_WT

candidate_distribution <-
  ggplot2::ggplot(ptmEWM_topfilt, ggplot2::aes(x = group,y = log2intensity, fill = group)) +
  ggplot2::scale_fill_manual(name = "EGFR Status",values= c("Non-amplified"="dodgerblue4","EGFR amplified"="coral4"),
                             )+
  ggplot2::geom_violin(trim = F,alpha=.7) +
  stat_compare_means(method="wilcox.test")+
  ggplot2::geom_boxplot(width=0.1,outlier.shape = NA,fill="white") +
  labs(y = "signature weighted log2 Intensity", x = "EGFR Status") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text =  element_text(face = "bold"),
    axis.text =  element_text(angle = 45, hjust = 1)
  )+ggtitle(label = "Intensity distribution by group sepcific sites")+
  theme(legend.position = "right",plot.title = element_text(color="darkblue",face = "bold"))

temporary_out<-(phosphosite_occurence_summary+intensity_distribution)/
  (candidate_distribution+toptable)
ggsave(plot=temporary_out,save_here(object_name = "Test_Sample_report.pdf"),height = 15,width = 12)


##### 3. Kinase Inferecne & Potential active pathways ######
top_kinases<-openxlsx::read.xlsx(".//output/20260512_RRS_1296_PCF_Phospho2/PTM_Decoupler_Kinase_Infernce/2026-08-11//Kinases_diff_EGFR_amp.xlsx")
top_kinases$group<-ifelse(top_kinases$score<0,"EGFR-amplified","non-amplified")

##### Load DB ressources #####
# enzsub<-OmnipathR::enzyme_substrate(resources = "ProtMapper",cache = F)
resources<-c("ProtMapper")
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

#### Single Kinase Score 
dc_kinacts<-purrr::map(names(enzsub_dbs),function(x){
  dc_kin_act <- dc_kin_activities(
    phospho_differential_analysis = ptm_sse_raw,
    dc_ptm = enzsub_dbs[[x]],
    net_idenitifiers = "genes",
    inference_method = "viper"
  )
  dc_kin_act$SKDB<-x
  dc_kin_act<-dc_kin_act[order(dc_kin_act$p_value,abs(dc_kin_act$score)),]
  dc_kin_act$rnk<-1:nrow(dc_kin_act)
  dc_kin_act
})%>%do.call(rbind,.)

### use the Kinases which where relevant in the comparison!
# evaluate top 20 Kinases from conensus of all databases
dc_kinacts_pval<-dc_kinacts
dc_kinacts_pval$log10p_value<--log10(dc_kinacts_pval$p_value)
dc_kinacts_pval$group<-NA
dc_kinacts_pval$group[dc_kinacts_pval$source%in%top_kinases$source[top_kinases$group=="EGFR-amplified"]]<-"EGFR-amplified"
dc_kinacts_pval$group[dc_kinacts_pval$source%in%top_kinases$source[top_kinases$group!="EGFR-amplified"]]<-"Non-amplified"
dc_kinacts_pval<-dc_kinacts_pval[!is.na(dc_kinacts_pval$group),]
dc_kinacts_pval<-dc_kinacts_pval[dc_kinacts_pval$p_value<0.05,]

active_kinases <- ggpubr::ggballoonplot(
  dc_kinacts_pval,
  x = "group",
  "source",
  fill = "score",
  size = "log10p_value"
) +
  ggplot2::scale_fill_gradientn(name = "Kinase Score", colours = c("#FFFFFF", "#67000D")) +
  xlab(label = "Kinase Substrate Database") +
  ylab(label = "Kinase") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text =  element_text(face = "bold"),
    axis.title =  element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
active_kinases

dc_kinacts_pval<-dc_kinacts

##### 4. Run pathway analysis ######
library(enrichplot)
data("biological_list", package = "proteoLab")
kegg<-msigdbr::msigdbr(species="Homo sapiens",category="C2",subcategory = "CP:KEGG_LEGACY")%>%
  dplyr::select(gs_name,gene_symbol)

ora_result <- clusterProfiler::enricher(
  gene = (dc_kinacts$source[dc_kinacts$p_value < 0.05]),
  pAdjustMethod = "BH",
  universe = unique(enzsub_dbs[[1]]$enzyme_genesymbol),
  TERM2GENE = kegg,
  pvalueCutoff = 0.05
)
ora_result@result<-ora_result@result[grep("SIGNALING",rownames(ora_result@result)),]
ora_result@result$Description<-gsub("KEGG_","",(ora_result@result$Description))
ora_result@result$Description<-gsub("_PATHWAY","",(ora_result@result$Description))
ora_barplot<-barplot(ora_result,showCategory = 10)

temporary_out2<-(phosphosite_occurence_summary+intensity_distribution)/
  (candidate_distribution+toptable)/
  (active_kinases+ora_barplot)

ggsave(plot=temporary_out2,save_here(object_name = "Test_Sample_report.pdf"),height = 17,width = 15)


##### 4. Carnival Inferecne & Potential active pathways ######
# library(decoupleR)
# library(OmnipathR)
# library(CARNIVAL)
# library(dplyr)
# 
# # ---- 1. Get the prior knowledge network from OmniPath ---- # Signed, directed interactions (signaling network) 
# omnipath_pkn <- import_omnipath_interactions(
#   organism = 9606,# human
#   filter_databases = c("SIGNOR", "PhosphoSite", "InnateDB", "SPIKE"))%>%# curated, signed sources) %>%
#   filter(consensus_direction == 1) %>%  # keep only directed
#   filter(consensus_stimulation == 1 |
#            consensus_inhibition == 1)  # keep only signed
# 
# # Convert to CARNIVAL's expected PKN format: source, interaction (1/-1), 
# pkn <- omnipath_pkn %>%
#   transmute(
#     source = source_genesymbol,
#     interaction = if_else(consensus_stimulation == 1, 1, -1),
#     target = target_genesymbol
#   ) %>%
#   distinct()
# 
# # ---- 2. Format kinase activity scores as CARNIVAL measurements ---- # 'measurements' = your decoupler kinase activity scores (named vector, kinase -> score) # CARNIVAL expects normalized values roughly in [-1, 1] or z-scores work fine too
# kinase_scores <- dc_kinacts[dc_kinacts$p_value<0.05,] %>%
#   dplyr::select(source, score) %>%
#   tibble::deframe()  # named vector: kinase symbol -> activity score
# 
# # Only keep kinases that are actually present in the PKN, otherwise CARNIVAL can't place them
# kinase_scores <- kinase_scores[names(kinase_scores) %in% c(pkn$source, pkn$target)]
# # ---- 3. (Optional) Define perturbation/input nodes ---- # If you know what was perturbed (e.g., a growth factor treatment), set it as the input # If unknown, CARNIVAL can run "invertedRun" mode which doesn't require a fixed input
# 
# # Example if you have a known perturbation:
# # input_nodes <- c("EGF" = 1)
# 
# # ---- 4. Run CARNIVAL ----
# # Requires a solver - CPLEX or CBC (free) via lpSolve/rcplex backend 
# carnival_result <- runCARNIVAL(
#   measObj = kinase_scores,
#   netObj = pkn[pkn$target%in%names(kinase_scores)|pkn$source%in%names(kinase_scores),],
#   solver = 'cbc'
#   # free solver; use "cplex" if you have a license (faster, better for large networks)
# )
# 
# # ---- 5. Inspect results ----
# # carnival_result$weightedSIF = the inferred subnetwork (source, sign, target, weight) # carnival_result$nodesAttributes = predicted activity/sign per node
# 
# head(carnival_result$weightedSIF)
# head(carnival_result$nodesAttributes)
# 
# # Check specifically for PI3K-AKT-mTOR axis presence/activity pi3k_axis <- c("PIK3CA", "PIK3R1", "AKT1", "AKT2", "MTOR", "PDK1", "PTEN") carnival_result$nodesAttributes %>%
# filter(Node %in% pi3k_axis)
