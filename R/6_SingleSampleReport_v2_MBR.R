####===================================####
# Author: Dennis Friedel, PhD
# Date: 2026-09-07
# Bioinformatician,
# Department of Neuropahtology, University Clinic Heidelberg
####===================================####

##### Load libraries #####
set.seed(2905)
library("dplyr")
library("patchwork")
library("ggplot2")
library("ggpubr")
library("enrichplot")
library("SummarizedExperiment")
source("./R/utils/rip_functions.R")
source("./R/utils/utils_module.R")
source("./R/utils/import_module.R")
source("./R/utils/preprocess_module.R")
source("./R/utils/compare_module.R")
source("./R/utils/ptm_module.R")

analysis = paste0('/',Sys.Date()) # - Name of the analysis e.g Marker Identification
dataset = paste0('/Figure/')

##### Load cohort #####
ptm_se_raw<-readRDS(".//output/20260813_RRS_1296_PCF_Phospho//PTM_PRC_DEA/2026-09-08/ptm_se_raw.rds")
phospho_samples<-colnames(ptm_se_raw)
ptm_se_raw$EGFR<-gsub("Amp","EGFR amplification",ptm_se_raw$EGFR)%>%gsub("WT","EGFR unamplified",.)

##### References #####
library(enrichplot)
kegg<-msigdbr::msigdbr(species="Homo sapiens",category="C2",subcategory = "CP:KEGG_LEGACY")%>%
  dplyr::select(gs_name,gene_symbol)
reactome<-msigdbr::msigdbr(species="Homo sapiens",category="C2",subcategory = "CP:REACTOME")%>%
  dplyr::select(gs_name,gene_symbol)

###### Load Therapeutic Pathways defined by Philip and selected fro Reactome by Claude
therapeutic_pathways<-openxlsx::read.xlsx("./data/reactome_dbs_filtered.xlsx")
relevant_kin<-openxlsx::read.xlsx("./data/cancer_kinases_by_category.xlsx")

relevant_kin$Kinase<-gsub(" .*","",relevant_kin$Kinase)
relevant_kin$Category_names<-plyr::mapvalues(relevant_kin$Category,therapeutic_pathways$Category,therapeutic_pathways$Category_names)

# Phosphosite-Kinase relationship
resources<-"ProtMapper"
enzsub<-data.table::fread("./data/omnipath_webservice_enz_sub.tsv")%>%as.data.frame()
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
enzsub_dbs_protmapper<-enzsub_dbs$ProtMapper

##### Iterate through samples ######
purrr::map(phospho_samples,function(sample){
  message(sample)
  
  ptm_sse_raw<-ptm_se_raw[,sample]
  peptide_sequences_identified<-ptm_sse_raw$peptide_sequences_identified
  ptm_sse_raw<-ptm_sse_raw[!is.na(assay(ptm_sse_raw)),]
  ptm_sse_raw@elementMetadata$log2intensity<-assay(ptm_sse_raw)[,1]
  ptm_sse_raw@elementMetadata$log2intensity_1<-assay(ptm_sse_raw,assayNames(ptm_sse_raw)[2])[,1]
  ptm_sse_raw@elementMetadata$log2intensity_2<-assay(ptm_sse_raw,assayNames(ptm_sse_raw)[3])[,1]
  ptm_sse_raw@elementMetadata$log2intensity_3<-assay(ptm_sse_raw,assayNames(ptm_sse_raw)[4])[,1]
  
  ptmEWM<-ptm_sse_raw@elementMetadata
  
  #####---------  A. Sample Info ------------------#####
  infotable<-data.frame(
    "sample" = c(
      "Sample ID:  ",
      "Diagnosis:  ",
      "Molecular Context:  ",
      "Methylation Class:  "
    ),
    "values" =
      c(
        ptm_sse_raw$sample_id,
        "Glioblastoma WT",
        ptm_sse_raw$EGFR,
        ptm_sse_raw$Methylation_Class
      )
  )
  
  gginfo<-ggpubr::ggtexttable(
    infotable,
    rows = NULL,cols = NULL,
    theme = ttheme(
      tbody.style = tbody_style(fill = "white"),
    )
  )
  gginfo
  
  #####--------- B. Quality Control ------------------#####
  thesholds <- c(5000,"65-70%","90:10:0.1", "10-20%","200-500")
  quality_table <- data.frame(
    "QC_Metric" = c(
      "Phosphosites Identified",
      "Enrichment efficiency",
      "Ratio S/T/Y-Distribution",
      "Protmapper Coverage",
      "Mapped Sites"
    ),
    "Value" = c(
      nrow(ptmEWM),
      paste0(round((nrow(ptmEWM)/peptide_sequences_identified)*100,digits = 1),"%"),
      paste0(round(table(ptm_sse_raw@elementMetadata$amino_acid)/nrow(ptmEWM)*100,1),collapse = ":"),
      paste0(round(sum(ptmEWM$Name%in%enzsub_dbs_protmapper$p_site)/nrow(ptmEWM)*100),"%"),
      sum(ptmEWM$Name%in%enzsub_dbs_protmapper$p_site)
    ),
    "Threshold" = thesholds)
  quality_table$Status<-c(
    as.numeric(quality_table[1, 2]) > quality_table[1, 3],
    as.numeric(gsub("%", "", quality_table[2, 2])) > 65,
    all(abs((
      table(ptm_sse_raw@elementMetadata$amino_acid) / nrow(ptmEWM) * 100
    ) - c(90, 10, 0.1)) < c(5, 2, 1)),
    round(
      sum(ptmEWM$Name %in% enzsub_dbs_protmapper$p_site) / nrow(ptmEWM) * 100
    ) > 10,
    sum(ptmEWM$Name %in% enzsub_dbs_protmapper$p_site) > 200
  )
  quality_table$Status<-ifelse(quality_table$Status,"Pass","Failed")
  
  ggqctable<-ggpubr::ggtexttable(
    quality_table,
    rows = NULL,
    theme = ttheme(
      tbody.style = tbody_style(fill = "white"),
      colnames.style = colnames_style(
        face = "bold",
        color = "white",
        fill = "dodgerblue4"
      )
    ),
  ) %>% tab_add_hline(at.row =   1) %>%
    tab_add_vline(at.column = c(2,3,4))
  if(any(quality_table$Status=="Pass")){
    ggqctable<-ggqctable%>%table_cell_font(row=which(quality_table$Status=="Pass")+1,column = 4,color = "forestgreen")
  }
  if(any(quality_table$Status!="Pass")){
    ggqctable<-ggqctable%>%table_cell_font(row=which(quality_table$Status!="Pass")+1,column = 4,color = "firebrick")
  }
  
  #####--------- C. Pathway activity summary ------------------#####
  
  ######--------- Kinase activity inference ------ 
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
  
  dc_kinacts_pval<-dc_kinacts[dc_kinacts$p_value<0.05,]
  dc_kinacts_pval$log10p_value<--log10(dc_kinacts_pval$p_value)
  
  active_kinases <- ggpubr::ggballoonplot(
    dc_kinacts_pval,
    x = "condition",
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
    )+coord_flip()
  active_kinases
  
  ######--------- ORA ------ 
  ora_result <- clusterProfiler::enricher(
    gene = (dc_kinacts$source[dc_kinacts$p_value < 0.05]),
    pAdjustMethod = "BH",
    universe = unique(enzsub_dbs[[1]]$enzyme_genesymbol),
    TERM2GENE = reactome[reactome$gs_name%in%therapeutic_pathways$Pathway,],
    pvalueCutoff = 1,qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 300
  )
  
  res<-ora_result@result
  res$category<-plyr::mapvalues(res$Description,therapeutic_pathways$Pathway,therapeutic_pathways$Category_names)
  res$Description<-gsub("REACTOME_","",(res$Description))%>%
    gsub("_PATHWAY","",.)%>%
    gsub("_"," ",.)
  res$Description<-factor(res$Description,levels = res$Description[order(res$Count,decreasing = F)])
  
  ####### Plot Active Pathways #####
  ## --- Hit count / hit rate per category ---
  ## Normalized by number of pathways *tested* in that category, since
  ## categories with more pathways (e.g. CDK_CELL_CYCLE, n=57) will always
  ## rack up more raw hits than a sparse category (e.g. SRC_FAMILY, n=2).
  ora_df <- as.data.frame(res) %>%
    mutate(
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      is_sig = p.adjust < 0.05        # your significance threshold
    )
  ora_df$Pathway<-ora_df$Description
  
  category_hit_summary <- ora_df %>%
    group_by(category) %>%
    summarise(
      n_pathways_tested   = n(),
      n_pathways_sig      = sum(is_sig, na.rm = TRUE),
      hit_rate            = n_pathways_sig / n_pathways_tested,
      median_padj         = median(p.adjust, na.rm = TRUE),
      min_padj            = min(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(hit_rate))
  
  
  ## --- 3b. Best representative pathway per category ---
  ## The single most significant pathway per category = your reportable
  ## "headline" hit for that signaling axis.
  category_top_pathway <- ora_df %>%
    group_by(category) %>%
    slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
    dplyr::select(category, Pathway, p.adjust, GeneRatio, geneID) %>%
    arrange(p.adjust)
  
  ## --- 3c. Union of driving kinases per category ---
  ## Pools the overlapping/leading-edge genes across ALL significant
  ## pathways within a category, so you see how much of your active-kinase
  ## list is actually explained by that axis (independent of pathway count).
  category_driving_kinases <- ora_df %>%
    filter(is_sig) %>%
    group_by(category) %>%
    summarise(
      driving_kinases = paste(sort(unique(unlist(strsplit(geneID, "/")))), collapse = "/"),
      n_unique_kinases = length(unique(unlist(strsplit(geneID, "/")))),
      .groups = "drop"
    ) %>%
    arrange(desc(n_unique_kinases))
  
  ## --- combine into one summary table ---
  category_summary_full <- category_hit_summary %>%
    left_join(category_top_pathway,
              by = "category") %>%
    left_join(category_driving_kinases, by = "category")
  
  plot_df <- ora_df %>%
    filter(is_sig) %>%             # or drop this filter to show all tested pathways
    mutate(
      Pathway_short = stringr::str_replace(Pathway, "^REACTOME_", ""),
      Pathway_short = stringr::str_trunc(Pathway_short, 40),
      neglog10padj  = -log10(p.adjust)
    )
  
  plot_df$category<-factor(plot_df$category)
  plot_df$Pathway_short<-factor(plot_df$Pathway_short,levels=plot_df$Pathway_short[order(plot_df$category)])
  
  active_pathways_plot<-ggplot(plot_df, aes(x = category, y = Pathway_short)) +
    geom_point(aes(size = Count, fill = neglog10padj),color="black",shape=21,stroke=0.6) +
    scale_fill_gradient(name = expression(-log[10]~adj.~p),low = "steelblue", high = "firebrick") +
    scale_size(name = "Kinase Count",range = c(6,16))+
    theme_bw(base_size = 14) +
    theme(text = ggplot2::element_text(size = 12),
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5,color="darkblue"),
          strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
          strip.text =  ggplot2::element_text(face = "bold"),
          axis.text =  element_text(color = "grey25"),
          axis.text.x= element_text(angle = 45, hjust = 1,size = 13)
    ) +
    labs(
      y = NULL,x=NULL,
      size = "Gene count",
      color = expression(-log[10]~adj.~p),
      title = "Significant Reactome pathways by signaling category"
    )
  
  #####--------- D. Supporting Phosphosites ------------------#####
  ptmEWM$zscore<-scale(ptmEWM$log2intensity,center = T,scale = T)
  ptmEWM_ordered<-ptmEWM[order(ptmEWM$zscore,decreasing=T),]
  ptmEWM_ordered$rank<-1:nrow(ptmEWM_ordered)
  
  #Get driving Kinases
  driving_kinases<-strsplit(category_summary_full$driving_kinases,"/")%>%
    unlist()%>%
    table()%>%
    sort()%>%
    names()
  
  clinkin<-relevant_kin$Kinase
  clin_p_sites<-enzsub_dbs$ProtMapper$p_site[enzsub_dbs$ProtMapper$enzyme_genesymbol%in%clinkin]
  
  groupptable<-ptmEWM_ordered[ptmEWM_ordered$Name%in%clin_p_sites,c("gene","Name","zscore","rank")]%>%as.data.frame()
  colnames(groupptable)<-c("Protein","Residue","Score","Rank")
  
  pmdb<-enzsub_dbs$ProtMapper[enzsub_dbs$ProtMapper$p_site%in%groupptable$Residue&
                                enzsub_dbs$ProtMapper$enzyme_genesymbol%in%clinkin,c("p_site","enzyme_genesymbol")]
  pmdb$score<-plyr::mapvalues(pmdb$p_site,groupptable$Residue,groupptable$Score)%>%as.numeric()
  pmdb$rank<-plyr::mapvalues(pmdb$p_site,groupptable$Residue,groupptable$Rank)%>%as.numeric()
  pmdb$protein<-plyr::mapvalues(pmdb$p_site,groupptable$Residue,groupptable$Protein)
  pmdb$category<-plyr::mapvalues(pmdb$enzyme_genesymbol,relevant_kin$Kinase,relevant_kin$Category_names)
  
  topkin<-split(pmdb,pmdb$category)%>%
    purrr::map(.,function(x){
      x[order(x$rank),]%>%head(.,3)
    })%>%do.call(rbind,.)
  rownames(topkin)<-NULL
  
  if(is.null(topkin)){
    topcatphosphosites<-ggplot() + annotate("text", 0, 0, label = "Not enough evidence found")+theme_void()
  }else{
    colnames(topkin)<-c("Site","Kinase","Score","Rank","Protein","Category")
    topkin$Category<-factor(topkin$Category)
    
    catsites<-purrr::map(levels(topkin$Category),function(x){
      which(topkin$Category==x)
    })
    x<-1:length(catsites)
    
    
    topcatphosphosites<-ggpubr::ggtexttable(
      topkin[,c("Protein","Site","Kinase","Score","Rank","Category")],
      rows = NULL,
      theme = ttheme(
        tbody.style = tbody_style(fill = "white"),
        colnames.style = colnames_style(
          face = "bold",
          color = "darkblue",
          fill = "white"
        )
      ),
    ) %>% tab_add_hline(at.row =   1)%>%
      tab_add_vline(at.column = c(2,3,4,5,6))
    
    if(length(x)>1){
      color_pos<-unlist(catsites[x[x%%2==0]])+1
      topcatphosphosites<-topcatphosphosites%>%table_cell_bg(row=color_pos,
                                                             column = 1:ncol(topkin),
                                                             fill = "grey90",alpha=0.8)
      color_pos<-unlist(catsites[x[x%%2!=0]])+1
      topcatphosphosites<-topcatphosphosites%>%table_cell_bg(row=color_pos,
                                                             column = 1:ncol(topkin),
                                                             fill = "grey99",alpha=0.8)
    }
  }
  
  #####--------- Report ------------------#####
  report_layout  =  c(
    patchwork::area(1, 1),
    patchwork::area(1, 2),
    patchwork::area(1, 4),
    patchwork::area(2,3),
    patchwork::area(2, 2),
    patchwork::area(2, 4),
    patchwork::area(2, 5)
  )
  plot(report_layout)
  report <-
    plot_spacer()+
    gginfo +
    ggqctable+
    plot_spacer()+
    active_pathways_plot+
    topcatphosphosites+
    plot_spacer()+
    plot_layout(
      design = report_layout,
      ncol = 5,
      widths = c(0.1,10,4, 20,1),
      heights = c(1.0, 2.0)
    )+
    plot_annotation(
      title = "Phosphoproteomic Profiling",
      subtitle = "Research use only"
    )
  report
  ggsave(plot = report, 
         save_here(object_name = paste0("Report_", sample, "_", ptm_sse_raw$EGFR, ".pdf"),dataset_name = "Figures",analysis_name = ""),width = 20,height = 12)
  
  if(sample=="P02"){
    xlsx_path<-save_here(object_name = "Supplemental_table4.xlsx",dataset_name = paste0('/Supplemental_table/'),analysis_name = "")
    openxlsx::write.xlsx("",xlsx_path)
    openxlsx::createWorkbook(xlsx_path)
    wb<-openxlsx::loadWorkbook(xlsx_path)
    
    openxlsx::addWorksheet(wb = wb,sheetName = "Reactome_ORA")
    openxlsx::writeData(wb,sheet="Reactome_ORA",x=ora_df)
    
    openxlsx::addWorksheet(wb = wb,sheetName = "Representative_Sites")
    openxlsx::writeData(wb,sheet="Representative_Sites",x=topkin)
    openxlsx::saveWorkbook(wb,xlsx_path,overwrite = T)
  }
  return(report)
})
