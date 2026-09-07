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
dataset = paste0('/Figure1/')
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

##### Load data ###
proc_ms<-readRDS("./output/20260512_RRS_1296_PCF_Phospho2/PTM_PRC_DEA/2026-08-06/ptm_se_prc.rds")
raw_ms<-proc_ms$unfilt
clindat<-colData(raw_ms)[, c("experiment",
                    "Rhaissa_sample_id",
                    "experiment",
                    "EGFR",
                    "Class",
                    "Methylation_Class")]
openxlsx::write.xlsx(clindat,save_here(object_name = "Supplemental_table1.xlsx"))

raw_ms$group<-ifelse(raw_ms$group=="WT",'Non-amplified','EGFR-amplified')
raw_ms$Group<-factor(raw_ms$group,levels = c('Non-amplified','EGFR-amplified'))

df_long <- count_idenitfied_features(
  se_object = raw_ms,
  assay_name = "intensity",
  group = "Group",
  column_wise = T,
  in_percent = F
)
df_long$group<-ifelse(df_long$group==1,'Non-amplified','EGFR-amplified')
df_long$xpos<-c(1:5,7:11)
cohort_median<-median(df_long$value)

##### Figure 1 Panel A - Total PTM Identifications ######
fig1_panA <-
  ggplot2::ggplot(df_long, ggplot2::aes(x = xpos,y = value,fill=group))+
  ggplot2::geom_col(width=0.8)+
  ggplot2::scale_fill_manual(name="EGFR status",
                             values = layers,
                             labels = c("Non-amplified" = "Non-amplified (n=5)", "EGFR-amplified" = "EGFR-amplified (n=5)"))+
  ggplot2::scale_x_continuous(breaks=df_long$xpos,
                              labels = df_long$sample_id)+
  ggplot2::geom_hline(yintercept = cohort_median,linetype="dashed",linewidth=0.8)+
  ggplot2::annotate(
    "text",
    y = cohort_median,
    label = paste0("Cohort median = ",cohort_median),
    x = 0,
    hjust = 1.1,
    vjust = 0
  )+
  ggplot2::labs(y = "Number of confidently localized phosphosites", x = "Sample") +
  ggplot2::theme_minimal(base_size = 14) +
  ggtitle("Identified phosphosites per sample")+
  ggplot2::theme(
    plot.title = ggplot2::element_text(face="bold",hjust=0.5),
    panel.grid.major = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
    strip.text =  ggplot2::element_text(face = "bold"),
    axis.text =  ggplot2::element_text(angle = 45, hjust = 1))+
  coord_flip()
fig1_panA
ggsave(plot = fig1_panA,filename = save_here(object_name = "Figure1A.pdf"),width = 7,height = 7)

##### Figure 1 Panel B - PTM Coverage ######
binary_assay<-!is.na(assay(proc_ms$unfilt))
binary_assay<-as.data.frame(binary_assay)
binary_assay$occurence<-rowSums(binary_assay)
binary_assay<-binary_assay[binary_assay$occurence!=0,]

histdata<-as.data.frame(table(binary_assay$occurence))
colnames(histdata)<-c("Samples","Freq")
histdata$Samples

bluescale<-c(blues9,"darkblue")
names(bluescale)<-levels(histdata$Samples)
fig1_panelB<-ggpubr::ggbarplot(
  histdata,
  x = "Samples",
  y = "Freq",
  fill = "Samples",
  width = 0.7,
  color = NA
) +
  scale_fill_manual(values = bluescale, name = "Detected in \n x samples") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  geom_text(aes(x = Samples, y = Freq, label = Freq), vjust = -0.5) +
  theme(legend.position = "right") +
  ggplot2::labs(y = "Number of phosphosites", x = "Detected in x samples") +
  ggtitle("Phosphosite detection frequency across the Cohort (n=10)") +
  ggplot2::theme(
    text = ggplot2::element_text(size = 14),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    panel.grid.major = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
    strip.text =  ggplot2::element_text(face = "bold"),
    axis.text =  ggplot2::element_text(angle = 45, hjust = 1)
  )
ggsave(plot = fig1_panelB,filename = save_here(object_name = "Figure1B.pdf"),width = 10,height = 7)

##### Figure 1 Panel C - Distribution ######
long_df_int<-proteoLab:::get_intensity_long(se_object = proc_ms$filt,assay_name = "intensity_norm",
                                            group = "group")
long_df_int$EGFR<-ifelse(long_df_int$group=="WT",'Non-amplified','EGFR-amplified')
long_df_int$group<-NULL
fig1_panelC <-
  ggplot2::ggplot(long_df_int, ggplot2::aes(x = sample_id,y = value, fill = EGFR)) +
  ggplot2::geom_violin(trim = F,alpha=.5) +
  ggplot2::geom_jitter(aes(color=EGFR),width=0.1,size=1.5,alpha=0.1) +
  ggplot2::geom_boxplot(width=0.1,outlier.shape = NA,fill="white") +
  labs(y = expression(Log[2] ~ "normalized Intensities (Median Center)"), x = "") +
  ggplot2::scale_fill_manual(name="EGFR status",values=layers)+
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text =  element_text(face = "bold"),
    axis.text =  element_text(angle = 45, hjust = 1)
  )

ggsave(plot = fig1_panelC,filename = save_here(object_name = "Figure1C.pdf"),width = 10,height = 7)

figure1<-(fig1_panA+fig1_panelB)/fig1_panelC+plot_annotation(tag_levels = "A")
figure1
ggsave(plot = figure1,filename = save_here(object_name = "Figure1.pdf"),width = 15,height = 12)

##### B) Coefficient of variation
# assay_name<-"intensity"
# mat <-
#   return_matrix_from_se(se_object = raw_ms, assay_name = assay_name)
# linear_mat<-2^mat
# cv_values<-apply(linear_mat,1,calc_cv_fraction)
# identical(names(cv_values),rownames(raw_ms))
# summary(cv_values)



# ### B) Calculate CV in fraction ######
# group<-levels(raw_ms$Group)
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
# assay_long_cv$group<-factor(assay_long_cv$group,levels = c("WT","AMP"))
# assay_long_cv
# 
# coefficient_of_var_bp <-
#   ggplot2::ggplot(assay_long_cv, ggplot2::aes(x = value,y = group, fill = group)) +
#   coord_flip()+
#   ggplot2::scale_fill_manual(values=layers)+
#   ggplot2::geom_violin(trim = F,alpha=.7) +
#   ggplot2::geom_boxplot(width=0.1,outlier.shape = NA,fill="white") +
#   labs(y = "Layer", x = "Coefficient of variation") +
#   theme_minimal(base_size = 14) +
#   theme(
#     panel.grid.major = element_blank(),
#     strip.background = element_rect(fill = "grey90", color = NA),
#     strip.text =  element_text(face = "bold"),
#     axis.text =  element_text(angle = 45, hjust = 1)
#   )
# coefficient_of_var_bp
# ggsave(plot = ggbp,rip::save_here(object_name = "CV_all_residue.pdf"),width = 5,height = 5)
# 
# ### C) Calculate coverage over replicates ######
# binary_assay<-!is.na(assay(proc_ms$filt))
# binary_assay<-as.data.frame(binary_assay)
# binary_assay$protein<-rownames(proc_ms$filt)
# 
# long_df <-reshape(
#   binary_assay,
#   varying = setdiff(colnames(binary_assay), "protein"),
#   v.names = "detected",
#   times = setdiff(colnames(binary_assay), "protein"),
#   direction = "long"
# )
# 
# rownames(long_df)<-NULL
# long_df$time<-gsub("x","",long_df$time)
# 
# ### Add group and Patient Info 
# long_df$group<-proc_ms$filt$group[match(long_df$time,proc_ms$filt$sample_id)]%>%toupper()
# agg1<-aggregate(detected ~ protein +group,data=long_df,FUN=sum)
# names(agg1)[names(agg1)=="detected"]<-"n_detected"
# 
# agg1<-agg1[agg1$n_detected>0,]
# 
# bar_df<-as.data.frame(table(agg1$group,agg1$n_detected))
# names(bar_df)<-c("group","n_detected","n_residues")
# bar_df$n_detected<-factor(bar_df$n_detected,levels = c(1:5))
# 
# id_per_group<-ggplot(bar_df,aes(x=group,y=n_residues,fill=n_detected))+
#   geom_col(color="black",width=0.7)+
#   scale_fill_manual(values = c(
#     "5" = "grey5",
#     "4" = "grey30",
#     "3" = "grey60",
#     "2" = "grey80",
#     "1" = "grey100"),
#     name = "Detected in",labels = c("5"= "5 replicates",
#                                     "4"= "4 replicates",
#                                     "3"= "3 replicates",
#                                     "2"= "2 replicates",
#                                     "1"= "1 replicates"
#     )
#   )+
#   labs(x=NULL,y="Identified phosphorylated Sites")+
#   theme_minimal(base_size=14)+
#   theme(
#     panel.grid.major = element_blank(),
#     strip.background = element_rect(fill="grey90",color=NA),
#     strip.text =  element_text(face="bold"),
#     axis.text =  element_text(angle=45,hjust=1))
# ggsave(plot = id_per_group,rip::save_here(object_name = "Residue_identification_per_group.pdf"),width = 5,height = 5)
# 
# plot_id/(group_numbers+coefficient_of_var_bp+id_per_group)
# 
# #### D ######
