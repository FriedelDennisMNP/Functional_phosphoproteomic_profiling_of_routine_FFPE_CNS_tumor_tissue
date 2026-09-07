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
dataset = paste0('/Figure2/')
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


##### Load PTM data ###
proc_ms<-readRDS("./output/20260512_RRS_1296_PCF_Phospho2/PTM_PRC_DEA/2026-08-06//ptm_se_prc.rds")

ms_imp<-proc_ms$imp
ms_imp$group<-ifelse(ms_imp$group=="WT",'Non-amplified','EGFR-amplified')
ms_imp$Group<-factor(ms_imp$group,levels = c('Non-amplified','EGFR-amplified'))

#### Figure 2A PCA #####
ms_pca<-runPCA_se(ms_imp,scale = T,center=T)

egfstat_annova <-pca_anova(
  variable = "Group",
  mat = t(SummarizedExperiment::assay(ms_pca)),
  data_df = as.data.frame(SummarizedExperiment::colData(ms_pca)))
subtit<-paste0(round(egfstat_annova$PoV_sum,2)*100,"% of EGFR Status contribute to sample variability in PC1")

pca_x <- as.data.frame(get_reduction_se(ms_pca, "PCA")$x)
pca_x$EGFR<-ms_pca$Group
Figure2a <-
  ggplot2::ggplot(pca_x, ggplot2::aes(x = PC1, y = PC2, color = EGFR)) +
  ggplot2::geom_point(size = 5) +
  ggplot2::xlab(label = paste0(colnames(pca_x)[1])) +
  ggplot2::ylab(label = paste0(colnames(pca_x)[2])) +
  scale_color_manual(name="EGFR status",values = layers)+
  ggplot2::stat_ellipse(type = "norm", linetype = 1)+
  ggtitle(label = "Distribution of samples in PC1 and PC2 ",subtitle = subtit)+
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
ggsave(plot = Figure2a,save_here(object_name = "Figure2.pdf"),width = 6,height = 6)

#### Figure 2B Correlation #####
cormat<-stats::cor(SummarizedExperiment::assay(ms_imp)[, ], method = "spear")
hm_lim_plus <- max((range(cormat, na.rm = TRUE)))
hm_lim_min <- max((range(cormat, na.rm = TRUE)))
hm_col <- circlize::colorRamp2(
  c(hm_lim_min, 0, hm_lim_plus),
  c("#08306B", "#FFFFFF", "#67000D")
)
library(ComplexHeatmap)

top_anno<-HeatmapAnnotation("Group"=ms_imp$Group,col = list("Group"=layers))
use_linkage <- "ward.D2"
corhm <- ComplexHeatmap::Heatmap(
  matrix = cormat,
  top_annotation = top_anno,
  name = "Correlation",
  border = T,
  clustering_distance_rows = function(x)
    stats::as.dist(1 - x),
  clustering_distance_columns = function(x)
    stats::as.dist(1 - x),
  clustering_method_columns = use_linkage,row_names_side = "left",
  clustering_method_rows = use_linkage,show_column_dend = F,
  heatmap_height = grid::unit(15, "cm"),
  heatmap_width = grid::unit(15, "cm"),
  rect_gp = grid::gpar(col = "grey60", lwd = 1),
  heatmap_legend_param = list(
    title = "Spearman Correlation",title_position = "leftcenter-rot",
    title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 12),
    legend_height = grid::unit(14, "cm")))

ggcorhm <- complexheatmap_to_ggplot(corhm)
ggcorhm


#### Figure 2C CV #####
# raw_ms<-proc_ms$unfilt
# raw_ms$group<-ifelse(raw_ms$group=="WT",'Non-amplified','EGFR-amplified')
# raw_ms$Group<-factor(raw_ms$group,levels = c('Non-amplified','EGFR-amplified'))
# group<-levels(raw_ms$Group)
# 
# assay_name<-"intensity"
# mat <-
#   return_matrix_from_se(se_object = raw_ms, assay_name = assay_name)
# linear_mat<-2^mat
# cv_values<-apply(linear_mat,1,calc_cv_fraction)
# identical(names(cv_values),rownames(raw_ms))
# 
# assay_long_cv<-purrr::map(group,function(layers_tmp){
#     tmpsamples<-raw_ms$sample_id[raw_ms$Group==layers_tmp]
#     cv_values<-apply(linear_mat[,tmpsamples],1,calc_cv_fraction)%>%as.data.frame()
#     cv_values$group<-layers_tmp
#     cv_values$gene<-rownames(cv_values)
#     colnames(cv_values)[1]<-"value"
#     rownames(cv_values)<-NULL
#     cv_values
#   })%>%do.call(rbind,.)
# 
# assay_long_cv$group<-factor(assay_long_cv$group,levels = c('Non-amplified','EGFR-amplified'))
# 
# Figure2B <-
#   ggplot2::ggplot(assay_long_cv, ggplot2::aes(x = value,y = group, fill = group)) +
#   coord_flip()+
#   ggplot2::scale_fill_manual(name="EGFR status",values=layers)+
#   ggplot2::geom_violin(trim = F,alpha=.5) +
#   ggplot2::geom_boxplot(width=0.1,outlier.shape = NA,fill="white") +
#   labs(y = "Layer", x = "Coefficient of variation") +
#   theme_minimal(base_size = 14) +
#   theme(
#     panel.grid.major = element_blank(),
#     strip.background = element_rect(fill = "grey90", color = NA),
#     strip.text =  element_text(face = "bold"),
#     axis.text =  element_text(angle = 45, hjust = 1)
#   )

# ggsave(plot = ggbp,rip::save_here(object_name = "CV_all_residue.pdf"),width = 5,height = 5)

