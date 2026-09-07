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
dataset  = '20260512_RRS_1296_PCF_Phospho2' # - Name of the dataset that is going to be analysed.
analysis = paste0('enhanced_volcano_chromsome/',Sys.Date()) # - Name of the analysis e.g Marker Identification
library("proteoLab")
library("dplyr")
library("patchwork")

ptm_se<-readRDS("./projects/20260512_RRS_1296_PCF_Phospho2/rds/ptm_se_id_1_preproc.rds")

resiude<-gsub(".*;","",ptm_se@elementMetadata$id_ptmgsea)%>%gsub("-p","",.)
gene_residue<-paste0(gsub(";.*","",ptm_se@elementMetadata$gene),"_",paste0(gsub("[0-9]","",resiude),gsub("[STY]","",resiude)))
rownames(ptm_se)<-gene_residue
dea_res<-wrapper_dea_gsea(
  ms_se = ptm_se,
  wrp_group = "group",
  wrp_test = "WT",
  wrp_contrast = list(c("Amp")),
  wrp_mode = "MANUAL",
  wrp_pv_fil = 0.05,
  wrp_fc_fil = 0.58,wrp_use_padj = T,
  gene_set_catalouge = NULL
)

tt<-dea_res$tt_combined
tt$gene<-gsub("_.*","",rownames(tt))

enhanced_voclanos <- plot_enhanced_volcano(
  limma_res = dea_res$tt_combined,
  tests_end_with = "logFC",
  colname_adjP = "adj_P_Val",
  colname_LogFC = "logFC",
  colname_pvalue = "P_Value",
  pointsize = 3,
  use_padj = F,
  alpha = 0.05,
  lfc =  0.58
)
c1_set<-msigdbr::msigdbr(species = "Homo sapiens",collection="C1")
c1_set%>% filter(grepl("chr7[p]",gs_name))%>%select(gene_symbol)%>%distinct%>%as.data.frame()->chr7p_set
c1_set%>% filter(grepl("chr7[q]",gs_name))%>%select(gene_symbol)%>%distinct%>%as.data.frame()->chr7q_set

### Run gene set test for chr7 ###
indexp<-tt$gene %in% chr7p_set$gene_symbol
resp<-limma::geneSetTest(
  indexp,
  tt$WT_vs_Amp_t,
  alternative = "down",
  ranks.only = F,
  nsim = 100000
)
resp

indexq<-tt$gene %in% chr7q_set$gene_symbol
resq<-limma::geneSetTest(
  indexq,
  tt$WT_vs_Amp_t,
  alternative = "down",
  ranks.only = F,
  nsim = 100000
)

### Barcodeplots
pdf(rip::save_here(object_name = "chr7_barcode_plots.pdf"),width = 10,height = 5)
limma::barcodeplot(
  statistics = tt$WT_vs_Amp_t,
  index = indexp,
  main= paste0(
  "Chr7p genes among all t-statisitcs (Single gene set \n test p-val:",
  round(resp, 3),
  ")"
  )
)
limma::barcodeplot(
  tt$WT_vs_Amp_t, 
  indexq, 
  main = paste0(
  "Chr7q genes among all t-statisitcs (Single gene set \n test p-val:", 
  round(resq,2),
  ")"
  )
)
dev.off()

### Volcano
chr7p_set <- tt$gene[tt$gene %in% chr7p_set$gene_symbol] %>%unique()
chr7q_set <- tt$gene[tt$gene %in% chr7q_set$gene_symbol] %>%unique()

p<-enhanced_voclanos[[1]]
volcano_7p <- p + 
  ggnewscale::new_scale_color() +
  ggplot2::geom_point(
    data = tt[tt$gene %in% chr7p_set, ],
    ggplot2::aes(
      x = WT_vs_Amp_logFC,
      y = -log10(WT_vs_Amp_P_Value),
      color = "Chr7p"
    ),
    size = 3
  ) + ggplot2::scale_color_manual(name = "Chr7p", values = ("Chr7p" = "green4"))


volcano_7q <- p +
  ggnewscale::new_scale_color() +
  ggplot2::geom_point(
    data = tt[tt$gene %in% chr7q_set, ],
    ggplot2::aes(
      x = WT_vs_Amp_logFC,
      y = -log10(WT_vs_Amp_P_Value),
      color = "Chr7q"
    ),
    size = 3
  ) + ggplot2::scale_color_manual(name = "Chr7q", values = ("chr7q" = "wheat"))

p_chr7<-list(p,volcano_7p,volcano_7q)
ggplot2::ggsave(plot = p_chr7,rip::save_here(object_name = "chr7_volcano.pdf"),height = 10,width = 12)
