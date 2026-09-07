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
analysis = paste0('/',Sys.Date()) # - Name of the analysis e.g Marker Identification
dataset = paste0('/Figure3/')
library("rip")
library("proteoLab")
library("SummarizedExperiment")
library("dplyr")
library("patchwork")
library("ggplot2")
calc_cv_fraction<-function(x,na.rm=T){
  sd(x,na.rm=T)/mean(x,na.rm=T)
}

### Set colors ###
layers <- c(
  'EGFR-amplified' = "#F77576",
  'Non-amplified' = "#91BFD1"
)

ptm_se_prc<-readRDS(".//output/20260512_RRS_1296_PCF_Phospho2/PTM_PRC_DEA/2026-08-06///ptm_se_prc.rds")


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
limma_res<-dea_res$tt_combined

#### Plot Custom enhanced Voclano
alpha = 0.05
lfc = .58
use_padj = F
pointsize=5
labsize=5
dot_colors=
  c(
    "#D3D3D3",# not significant
    "#A1C9F4",# FC only
    "#BFE3A1",# p-value only 
    "#FFB3BA")# both significant


# Convert to toptable to format for enhanced volcano
limma_res_input<-limma_res
list_limma_results<-convert_toptable_to_envo_input(limma_res = limma_res_input,
                                                   tests_end_with="logFC",
                                                   colname_adjP = "adj_P_Val",
                                                   colname_LogFC = "logFC",
                                                   colname_pvalue = "P_Value")
tmp_toptable<-list_limma_results[[1]]
tmp_toptable$site<-rownames(tmp_toptable)
openxlsx::write.xlsx(tmp_toptable,save_here(object_name = "Supplemental_table2.xlsx"))

A<-gsub("vs_.*"," ",tmp_toptable$group)%>%gsub("_"," ",.)%>%unique()
B<-gsub(".*vs_"," ",tmp_toptable$group)%>%gsub("_"," ",.)%>%unique()
tmp_toptable$significant <-
  tmp_toptable$`p-value` < alpha &
  abs(tmp_toptable$log2FoldChange) > lfc
use_this_pvalue <- "p-value"
y_axis_label <- bquote( ~ -Log[10] ~ italic(P))


top_candidates<-rownames(tmp_toptable)[tmp_toptable$significant]
genes_of_interest<-paste(c("EGFR_Y1110","EGFR_Y1197","EGFR_T693","PSIP1","VIM","MAP2","NES","TP53BP"),collapse = "|")
show_top_candidates<-top_candidates[grepl(genes_of_interest,top_candidates)]

## Plot Enhanced Volcano
enVo <-
  EnhancedVolcano::EnhancedVolcano(
    toptable = tmp_toptable,
    title = "Comparison of GBM with EGFR Status Amplified vs. Non-Amplified",
    subtitle = paste0(
      "Number significant phosphosites: ",
      sum(tmp_toptable$significant, na.rm = T)
    ),
    lab = rownames(tmp_toptable),
    selectLab = show_top_candidates,
    x = 'log2FoldChange',
    y = use_this_pvalue,
    col = dot_colors,
    ylim = c(0, max(tmp_toptable$'log2FoldChange') + 0.5),
    pCutoff = alpha,
    FCcutoff = lfc,
    pointSize = 4,
    colAlpha = 0.7,
    legendPosition = "top",
    legendLabSize = 12,
    legendIconSize = 5.0,
    boxedLabels = T,
    drawConnectors = T,
    widthConnectors = 0.75,lengthConnectors = 1,
    colConnectors = "grey10",
    max.overlaps = 15,
    arrowheads = FALSE,
    endsConnectors = "first",
    typeConnectors = "closed",
    maxoverlapsConnectors = 10,
    min.segment.length = 20,labSize = 3
  ) +
  ggplot2::annotate(
    "text",
    x = max(tmp_toptable$'log2FoldChange'),
    y = 0,
    label = "Amplified",
    size = ggplot2::unit(6, "pt"),
    parse = F
  ) +
  ggplot2::annotate(
    "text",
    x = min(tmp_toptable$'log2FoldChange'),
    y = 0,
    label = "Non-Amplified",
    size = ggplot2::unit(6, "pt"),
    parse = F
  ) +
  ggplot2::ylab(y_axis_label)
enVo

ggsave(
  plot = enVo,
  save_here(object_name = "Figrue3_Volcano.pdf"),
  width = 10,
  height = 10
)
  