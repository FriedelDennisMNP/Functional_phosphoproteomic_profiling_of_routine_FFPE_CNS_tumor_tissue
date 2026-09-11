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

### Set colors ###
layers <- c(
  'EGFR-amplified' = "#F77576",
  'Non-amplified' = "#91BFD1"
)

######------- 1. Load Data -----------####
proc_ms<-
  readRDS("./output/20260813_RRS_1296_PCF_Phospho/PTM_PRC_DEA/2026-09-11/ptm_se_prc.rds")
proc_wp<-
  readRDS("./output/20260817_RR1296_PCF/WP_PRC_DEA/2026-09-11//ms_se_prc.rds")

clindat<-colData(proc_ms$unfilt)[, c("experiment",
                             "Rhaissa_sample_id",
                             "experiment",
                             "EGFR",
                             "Class",
                             "Methylation_Class")]

openxlsx::write.xlsx(clindat,save_here(
  dataset_name = "Supplemental_table",
  object_name = "Supplemental_table1.xlsx"))

######------- 2. Plot Figure/SuppFigure 1 A PTM/WP data Quality  -----------####
proc_list<-list("PTM" = proc_ms$unfilt, "WP" = proc_wp$unfilt)
figure_1A<-purrr::map(names(proc_list),function(raw_ms_name){
  raw_ms<-proc_list[[raw_ms_name]]
  
  if(raw_ms_name=="WP"){
    plottext<-"proteins"
  }else{
    plottext<-"phosphosites"
  }
  
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
  

  fig1_panA <-
    ggplot2::ggplot(df_long, ggplot2::aes(x = xpos, y = value, fill = group)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_fill_manual(
      name = "EGFR status",
      values = layers,
      labels = c("Non-amplified" = "Non-amplified (n=5)", 
                 "EGFR-amplified" = "EGFR-amplified (n=5)")
    ) +
    ggplot2::scale_x_continuous(breaks = df_long$xpos, 
                                labels = df_long$sample_id) +
    ggplot2::geom_hline(yintercept = cohort_median,
                        linetype = "dashed",
                        linewidth = 0.8) +
    ggplot2::annotate(
      "text",
      y = cohort_median,
      label = paste0("Cohort median = ", cohort_median),
      x = 0,
      hjust = 1.1,
      vjust = 0
    ) +
    ggplot2::labs(y = paste0("Number of confidently idenitfied ", plottext),
                  x = "Sample") +
    ggplot2::theme_minimal(base_size = 14) +
    ggtitle(paste0("Identified ", plottext, " per sample")) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.major = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
      strip.text =  ggplot2::element_text(face = "bold"),
      axis.text =  ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    coord_flip()
  fig1_panA
})

ggsave(plot = figure_1A[[1]],
       filename = save_here(dataset_name = "Figures",
                            object_name = "Figure1A.pdf"),
       width = 7,height = 7)

ggsave(plot = figure_1A[[2]],
       filename = save_here(dataset_name = "Figures",
                            object_name = "Supplemental_Figure1A.pdf"),
       width = 7,height = 7)

######------- 3. Plot Figure/SuppFigure 1B PTM/WP data Quality  -----------####
figure_1B <- purrr::map(names(proc_list), function(raw_ms_name) {
  raw_ms <- proc_list[[raw_ms_name]]
  
  if (raw_ms_name == "WP") {
    plottext <- "Protein"
  } else{
    plottext <- "Phosphosite"
  }
  binary_assay <- !is.na(assay(raw_ms))
  binary_assay <- as.data.frame(binary_assay)
  binary_assay$occurence <- rowSums(binary_assay)
  binary_assay <- binary_assay[binary_assay$occurence != 0, ]
  
  histdata <- as.data.frame(table(binary_assay$occurence))
  colnames(histdata) <- c("Samples", "Freq")
  histdata$Samples
  
  bluescale <- c(blues9, "darkblue")
  names(bluescale) <- levels(histdata$Samples)
  fig1_panelB <- ggpubr::ggbarplot(
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
    ggplot2::labs(y = paste0("Number of ", plottext), x = "Detected in x samples") +
    ggtitle(paste0(plottext, " detection frequency across the Cohort (n=10)")) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.major = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
      strip.text =  ggplot2::element_text(face = "bold"),
      axis.text =  ggplot2::element_text(angle = 45, hjust = 1)
    )
})
ggsave(
  plot = figure_1B[[1]],
  filename = save_here(dataset_name = "Figures",
                       object_name = "Figure1B.pdf"),
  width = 10,
  height = 7
)
ggsave(
  plot = figure_1B[[2]],
  filename = save_here(
    dataset_name = "Figures",
    object_name = "Supplemental_Figure1B.pdf"),
  width = 10,
  height = 7
)

######------- 4. Plot Figure/SuppFigure 1C PTM/WP data Quality  -----------####
proc_ms_adj<-
  readRDS("./output/20260813_RRS_1296_PCF_Phospho/PTM_PRC_ADJ_DEA//2026-09-11/ptm_se_adjprc.rds")

proc_list_flt<-list("PTM" = proc_ms_adj, "WP" = proc_wp$filt)
figure_1C<-purrr::map(names(proc_list_flt),function(raw_ms_name){
  
  if(raw_ms_name=="WP"){
    long_df_int<-get_intensity_long(se_object = proc_list_flt[[raw_ms_name]],assay_name = "intensity",
                                                group = "group")
  }else{
    long_df_int<-get_intensity_long(se_object = proc_list_flt[[raw_ms_name]],assay_name = "adjusted",
                                                group = "group")
  }
    
  long_df_int$EGFR<-ifelse(long_df_int$group=="WT",'Non-amplified','EGFR-amplified')
  long_df_int$group<-NULL
  fig1_panelC <-
    ggplot2::ggplot(long_df_int, ggplot2::aes(x = sample_id,y = value, fill = EGFR)) +
    ggplot2::geom_violin(trim = F,alpha=.5) +
    ggplot2::geom_jitter(aes(color=EGFR),width=0.1,size=1.5,alpha=0.1) +
    ggplot2::geom_boxplot(width=0.1,outlier.shape = NA,fill="white") +
    labs(y = expression(Log[2] ~ "normalized Intensities (Median Center+ proteom adjusted)"), x = "") +
    ggplot2::scale_fill_manual(name="EGFR status",values=layers)+
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text =  element_text(face = "bold"),
      axis.text =  element_text(angle = 45, hjust = 1)
    )
  
})

ggsave(
  plot = figure_1C[[1]],
  filename = save_here(dataset_name = "Figures",
                       object_name = "Figure1D.pdf"),
  width = 15,
  height = 7
)
ggsave(
  plot = figure_1C[[2]],
  filename = save_here(dataset_name = "Figures",
                       object_name = "Supplemental_Figure1C.pdf"),
  width = 15,
  height = 7
)


figure1 <- (figure_1A[[1]] + figure_1B[[1]]) / figure_1C[[1]] + 
  plot_annotation(tag_levels = "A")
suppfigure1 <- (figure_1A[[2]] + figure_1B[[2]]) / figure_1C[[2]] + 
  plot_annotation(tag_levels = "A")

ggsave(
  plot = figure1,
  filename = save_here(dataset_name = "Figures",
                       object_name = "Figure1.pdf"),
  width = 15,
  height = 12
)
ggsave(
  plot = suppfigure1,
  filename = save_here(dataset_name = "Figures",
                       object_name = "Supplemental_Figure1.pdf"),
  width = 15,
  height = 12
)

