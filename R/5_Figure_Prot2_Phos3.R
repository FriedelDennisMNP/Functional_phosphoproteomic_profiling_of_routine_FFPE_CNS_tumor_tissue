####===================================####
# Author: Dennis Friedel, PhD
# Date: 2026-09-07
# Bioinformatician,
# Department of Neuropahtology, University Clinic Heidelberg
####===================================####

######------- 0. Load Packages and Functions -----------####
set.seed(2905)
library("dplyr")
library("patchwork")
library("ggplot2")
library("SummarizedExperiment")
source("./R/utils/rip_functions.R")
source("./R/utils/utils_module.R")
source("./R/utils/import_module.R")
source("./R/utils/preprocess_module.R")
source("./R/utils/compare_module.R")

# Set variables for save_here
analysis = paste0('/')

### Set colors 
layers <- c(
  'EGFR-amplified' = "#F77576",
  'Non-amplified' = "#91BFD1"
)

######------- 1. Load Data -----------####
ptm<-
  readRDS("./output/20260813_RRS_1296_PCF_Phospho/PTM_PRC_ADJ_DEA/2026-09-09/ptm_se_adjprc.rds")
wp<-
  readRDS("./output/20260817_RR1296_PCF/WP_PRC_DEA/2026-09-08/ms_se_prc.rds")

datasets<-list("PTM"=ptm,
               "WP"=wp$imp)

######------- 2. Figure 2-3 PCA  -----------####
figure2_3_A<-purrr::map(datasets,function(ms_imp){
  ms_imp$group<-ifelse(ms_imp$group=="WT",'Non-amplified','EGFR-amplified')
  ms_imp$Group<-factor(ms_imp$group,levels = c('Non-amplified','EGFR-amplified'))
  
  ms_pca<-runPCA_se(ms_imp,scale = T,center=T)
  
  egfstat_annova <-pca_anova(
    variable = "group",
    mat = t(SummarizedExperiment::assay(ms_pca)),
    data_df = as.data.frame(SummarizedExperiment::colData(ms_pca)))
  subtit<-paste0("PC1 ",round(egfstat_annova$PoV_sum,3)*100,"% contribution to dataset variability ")
  
  pca_x <- as.data.frame(get_reduction_se(ms_pca, "PCA")$x)
  pca_x$EGFR<-ms_pca$Group
  Figure2a <-
    ggplot2::ggplot(pca_x, ggplot2::aes(x = PC1, y = PC2, color = EGFR)) +
    ggplot2::geom_point(size = 5) +
    ggplot2::xlab(label = paste0(colnames(pca_x)[1])) +
    ggplot2::ylab(label = paste0(colnames(pca_x)[2])) +
    scale_color_manual(name="EGFR status",values = layers)+
    ggplot2::stat_ellipse(type = "norm", linetype = 1)+
    ggtitle(label = "Distribution of samples in PC1 and PC2 ")+
    theme_minimal()+
    theme(
      panel.grid.major = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text =  element_text(face = "bold"),
      axis.text =  element_text(angle = 45, hjust = 1),
      axis.title =  element_text(size=14)
    )+
    geom_text(
      data = aggregate(PC1 ~ EGFR,pca_x,median),
      aes(x=PC1+0.5,y=PC1-5,label = EGFR),
      hjust=0,fontface="bold"
    )
})

#### Figure 2A PCA Whole Proteome
ggsave(
  plot = figure2_3_A[[2]],
  save_here(dataset_name = "Figures",
              object_name = "Figure2A.pdf"),
  width = 6,
  height = 6
)

#### Figure 3A PCA Phospho Proteome
ggsave(
  plot = figure2_3_A[[1]],
  save_here(dataset_name = "Figures",
            object_name = "Figure3A.pdf"),
  width = 6,
  height = 6
)

######------- 3. Figure 2 B Protein Volcano  -----------####
wp_se_imp<-wp$imp
wp_se_imp$group<-toupper(wp_se_imp$group)

### Get Hallmarks from MSIGDB
msigdbr_coll<-msigdbr::msigdbr_collections()%>%as.data.frame()
categories<-paste0(msigdbr_coll$gs_collection,".",msigdbr_coll$gs_subcollection)
misg_df<-msigdbr::msigdbr(species = "Homo sapiens",collection = "H")  
misg_df <- split(x = misg_df$gene_symbol, f = misg_df$gs_name)
hallmarks<-misg_df

dea_res <- wrapper_dea_gsea(
  ms_se = wp_se_imp,
  wrp_group = "group",
  wrp_test = "AMP",
  wrp_contrast = list(c("WT")),
  wrp_mode = "MANUAL",
  wrp_pv_fil = 0.05,
  wrp_fc_fil = 0.58,
  wrp_use_padj = F,
  wrp_fgsea = T,
  wrp_minSize = 5,wrp_maxSize = 300,
  gene_set_catalouge = list("Hallmarks"=hallmarks)
)

#### Set parameters for enhanced Voclano
alpha = 0.05
lfc = .58
use_padj = F
pointsize=5
labsize=5
use_this_pvalue <- "p-value"
dot_colors=
  c(
    "#D3D3D3",# not significant
    "#A1C9F4",# FC only
    "#BFE3A1",# p-value only 
    "#FFB3BA")# both significant


# Convert to toptable to format for enhanced volcano
tmp_toptable<-convert_toptable_to_envo_input(limma_res = dea_res$tt_combined,
                                                   tests_end_with="logFC",
                                                   colname_adjP = "adj_P_Val",
                                                   colname_LogFC = "logFC",
                                                   colname_pvalue = "P_Value")[[1]]
A<-gsub("vs_.*"," ",tmp_toptable$group)%>%gsub("_"," ",.)%>%unique()
B<-gsub(".*vs_"," ",tmp_toptable$group)%>%gsub("_"," ",.)%>%unique()

tmp_toptable$significant <-
  tmp_toptable$`p-value` < alpha &
  abs(tmp_toptable$log2FoldChange) > lfc
y_axis_label <- bquote( ~ -Log[10] ~ italic(P))

top_candidates<-rownames(tmp_toptable)[tmp_toptable$significant]

## Plot Enhanced Volcano
enVo <-
  EnhancedVolcano::EnhancedVolcano(
    toptable = tmp_toptable,
    title = "Comparison of GBM with EGFR Status Amplified vs. Non-Amplified",
    subtitle = paste0(
      "Number significant proteins: ",
      sum(tmp_toptable$significant, na.rm = T)
    ),
    lab = rownames(tmp_toptable),
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
  save_here(dataset_name = "Figures",
            object_name = "Figrue2B.pdf"),
  width = 10,
  height = 10
)

######------- 4. Figure 3 C Protein Enrichment  -----------####
reactome_result<-dea_res$fgsea_result$Hallmarks@data
Figure3C <- ggplot2::ggplot(reactome_result, ggplot2::aes(x = reorder(pathway, NES), y = NES, fill = logAPV)) +
  ggplot2::geom_col(alpha = 0.9) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_gradient(name = expression(-log[10]~adj.~P),
                               low = "#08306B", high = "#67000D") +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Pathway", y = "Normalized Enrichment Score")+
  theme(text = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5,color="darkblue"),
        strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
        strip.text =  ggplot2::element_text(face = "bold"),
        axis.text =  element_text(color = "grey25"),
        axis.text.x= element_text(angle = 45, hjust = 1,size = 13)
  )
Figure3C
ggsave(
  plot = Figure3C,
  save_here(dataset_name = "Figures",
            object_name = "Figure3C.pdf"),
  width = 10,
  height = 10
)

######------- 5. Figure 3 B PTM Volcano  -----------####
ptm_se_imp<-ptm
ptm_se_imp$group<-toupper(ptm_se_imp$group)
dea_res <- wrapper_dea_gsea(
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
openxlsx::write.xlsx(tmp_toptable,save_here(dataset_name = "Supplemental_table",
                                            object_name = "Supplemental_table2.xlsx"))

A<-gsub("vs_.*"," ",tmp_toptable$group)%>%gsub("_"," ",.)%>%unique()
B<-gsub(".*vs_"," ",tmp_toptable$group)%>%gsub("_"," ",.)%>%unique()
tmp_toptable$significant <-
  tmp_toptable$`p-value` < alpha &
  abs(tmp_toptable$log2FoldChange) > lfc
use_this_pvalue <- "p-value"
y_axis_label <- bquote( ~ -Log[10] ~ italic(P))

top_candidates <- rownames(tmp_toptable)[tmp_toptable$significant]
genes_of_interest <- paste(
  c(
    "EGFR_Y1110",
    "EGFR_Y1197",
    "EGFR_T693",
    "PSIP1",
    "VIM",
    "MAP2",
    "NES",
    "TP53BP"
  ),
  collapse = "|"
)
show_top_candidates <- top_candidates[grepl(genes_of_interest, top_candidates)]

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
  save_here(dataset_name = "Figures",
            object_name = "Figure3B.pdf"),
  width = 10,
  height = 10
)
