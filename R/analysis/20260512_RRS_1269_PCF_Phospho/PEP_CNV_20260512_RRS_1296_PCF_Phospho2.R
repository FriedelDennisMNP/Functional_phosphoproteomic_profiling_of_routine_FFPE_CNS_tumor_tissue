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
analysis = paste0('PEP_CNV/',Sys.Date()) # - Name of the analysis e.g Marker Identification

library("rip")
library("proteoLab")
library("dplyr")
library("patchwork")
library("ggplot2")
library("SummarizedExperiment")
source("./R/utils/run_pathview_in.R")
source("./R/utils/copro.R")

### Load data
mq_results<-proteoLab::load_ms_results(ms_result_dir = "/mnt/add50/PATHO-PROTEOMICS/bioinformatics/results/DDA/20260512_RRS_1296_PCF_Phospho2///MaxQuant2.4.2.0/")
mq_list<-proteoLab::create_maxquant_se_list(mq_results = mq_results,ptm_probability_flt = 0.75,remove_low_quality = T)

### Annoate
metadata<-openxlsx::read.xlsx("./projects/20260512_RRS_1296_PCF_Phospho2/sample_sheets/Annotated_sample_sheet3_20260707.xlsx",sheet = 1)
rownames(metadata)<-metadata$sample_id

pep_se<-mq_list$pep_se

pep_se<-pep_se[,pep_se$sample_id%in%metadata$original_id]
colnames(pep_se)<-metadata$sample_id
SummarizedExperiment::colData(pep_se)<-S4Vectors::DataFrame(metadata[match(pep_se$sample_id,metadata$original_id),])

pep_se<-pep_se[,pep_se$group%in%c("Amp","WT")]
saveRDS(pep_se,save_here(object_name = "pep_se_raw.rds"))

### Run per sample CNV analysis

## Get chromosomal location
library(biomaRt)
genes<-unique(pep_se@elementMetadata$gene)
mart<-useMart("ensembl",dataset = "hsapiens_gene_ensembl")

locations<-getBM(
  attributes = c(
    "external_gene_name",
    "chromosome_name",
    "start_position",
    "end_position","band"
  ),
  filters = "external_gene_name",
  values = genes,
  mart = mart
)

locations_flt<-locations[locations$chromosome_name%in%c(1:22,"X"),]
locations_flt<-locations_flt[locations_flt$external_gene_name%in%pep_se@elementMetadata$gene,]

# Position of Peptides in locus 
# Unprocessed data 
fasta_path<-"/mnt/add50/PATHO-PROTEOMICS/bioinformatics/analysis/fasta_files/uniprotkb_Human_AND_reviewed_true_AND_m_2025_01_16.fasta"
fasta<-proteoLab::read_fasta(fasta_path)
names(fasta)<-gsub(".*GN=","",names(fasta))%>%gsub(" .*","",.)
fasta_flt<-fasta[names(fasta)%in%pep_se@elementMetadata$gene_names]

## Annotate ###
pep_anno<-pep_se
pep_anno<-pep_anno[pep_anno@elementMetadata$is_proteotypic,]
pep_anno<-pep_anno[pep_anno@elementMetadata$gene_names%in%names(fasta),]

### get position of sequence in gene to define order 
AnnoElementMetadata<-purrr::map(unique(pep_anno@elementMetadata$gene),function(gene_tmp){
  
  fasta_sel<-fasta_flt[names(fasta_flt)%in%gene_tmp]
  fasta_sel<-fasta_sel[[1]]
  
  sc_gene<-pep_anno@elementMetadata[pep_anno@elementMetadata$gene==gene_tmp,]%>%as.data.frame()
  sc_gene$seq_postion<-purrr::map(sc_gene$sequence,function(x)regexpr(x,fasta_sel)[[1]])%>%unlist()%>%factor()
  sc_gene$seq_order<-order(sc_gene$seq_postion)
  
  locations_gene<-locations_flt[locations_flt$external_gene_name==gene_tmp,]
  locations_gene<-locations_gene[1,]
  
  sc_gene$chromosome<-locations_gene$chromosome_name
  sc_gene$start_position<-locations_gene$start_position
  sc_gene$end_position<-locations_gene$end_position
  sc_gene$band<-locations_gene$band
  sc_gene$arm<-gsub("[0-9].*","",sc_gene$band)
  sc_gene
})%>%do.call(rbind,.)
rownames(AnnoElementMetadata)<-AnnoElementMetadata$row_id

pep_anno_flt<-pep_anno[rownames(AnnoElementMetadata),]
pep_anno_flt@elementMetadata<-S4Vectors::DataFrame(AnnoElementMetadata)

##### order by Chrom 
pep_ordered<-pep_anno_flt
pep_ordered@elementMetadata$chromosome<-factor(pep_ordered@elementMetadata$chromosome,levels = as.character(c(1:22,"X","Y")))
pep_ordered<-pep_ordered[order(
  pep_ordered@elementMetadata$chromosome,
  pep_ordered@elementMetadata$start_position,
  pep_ordered@elementMetadata$seq_order
),]
pep_ordered<-pep_ordered[!is.na(pep_ordered@elementMetadata$chromosome),]
pep_ordered@elementMetadata$chromosome_arm<-paste0(pep_ordered@elementMetadata$chromosome,pep_ordered@elementMetadata$arm)
chromosome_arm<-unique(pep_ordered@elementMetadata$chromosome_arm)

#### Iterate througb single sample
chromplot_list<-purrr::map(pep_ordered$sample_id,function(sample_name){
  
  ### ChromPlot settings
  col_high = "orange3"
  col_mid = "grey99"
  col_low = "dodgerblue4"
  dot_alpha=1
  cumu_chr_pos = get("chr_line_pos")
  cumu_pq_pos = get("pq_line_pos")
  xtick_pos = get("anno_xtick_pos")
  genome_size = get("anno_genome_size")
  genome_anno = get("anno_genome_chr")
  
  ##### Plot preproc
  scapegoat<-(pep_ordered)[,sample_name]
  EGFRgroup<-scapegoat$group
  scapegoat@elementMetadata$value<-assay(scapegoat)[,1]
  chromprot_input<-scapegoat@elementMetadata[,c("chromosome","start_position","seq_order","value","chromosome_arm","gene")]
  chromprot_input<-chromprot_input[!is.na(chromprot_input$value),]
  chromprot_input<-as.data.frame(chromprot_input)
  genome_center<-mean(chromprot_input$value,na.rm = T)
  
  #### Update Chromosome Position
  chromprot_input_cumsum<-purrr::map(levels(chromprot_input$chromosome),function(x){
    tmpdf<-chromprot_input[chromprot_input$chromosome==x,]
    if(!x%in%names(cumu_chr_pos)){
      return(tmpdf[tmpdf$chromosome==x,])
    }else{
      tmpdf$start_position<-tmpdf$start_position+
        as.numeric(cumu_chr_pos[names(cumu_chr_pos)==x])  
      return(tmpdf[tmpdf$chromosome==x,])
    }
  })%>%do.call(rbind,.)
  
  #### Get Chromosome Arm Center / Segments
  chromprot_input_med<-purrr::map(unique(chromprot_input_cumsum$chromosome_arm),function(x){
    tmpdf<-chromprot_input_cumsum[chromprot_input_cumsum$chromosome_arm==x,]
    tmpdf$arm_value_med<-mean(tmpdf$value)
    return(tmpdf)
  })%>%do.call(rbind,.)
  
  segments<-chromprot_input_med[,c("chromosome","chromosome_arm","arm_value_med")]
  rownames(segments)<-NULL
  segments<-distinct(segments)
  
  cumu_chr_p_pos<-cumu_chr_pos
  names(cumu_chr_p_pos)<-paste0(names(cumu_chr_p_pos),"p")
  
  names(cumu_pq_pos)<-c(seq(1:22),"X","Y")%>%paste0(.,"q")
  chompos_mapping<-c(c("1p"=0),sort(c(cumu_chr_p_pos,cumu_pq_pos)))
  
  segments$start_pos<-plyr::mapvalues(segments$chromosome_arm,names(chompos_mapping),chompos_mapping)%>%as.numeric()
  segments$end_pos<-c(segments$start_pos[2:length(segments$start_pos)],cumu_chr_pos["Y"])
  segments$arm_value_med<-segments$arm_value_med-genome_center
  segments$arm_col<-ifelse(segments$arm_value_med>0,"orange4","dodgerblue4")
  
  #### Smooth to identify trends
  chromprot_input_med$value_mid<-chromprot_input_med$value-genome_center
  chromprot_input_smooth<-purrr::map(unique(chromprot_input_med$chromosome_arm),function(x){
    tmpdf<-chromprot_input_med[chromprot_input_med$chromosome_arm==x,]
    tmpdf$value_roll<-smooth_roller(vals = tmpdf$value_mid,
                                    window_size = 5,
                                    iteration = 1,
                                    "mean")
    tmpdf
  })%>%do.call(rbind,.)
  
  egfr_df<-chromprot_input_smooth[chromprot_input_smooth$gene=="EGFR",]
  
  #### Plot Segment Plots
  chrom_plot <-
    ggplot(chromprot_input_smooth,
           aes(x = start_position,
               y = value_mid,
               color = value_mid)) +
    geom_point() +
    scale_color_gradient2(
      high = col_high,
      mid = col_mid,
      low = col_low,
      name = "CNV score"
      ,
      na.value = NA,
      midpoint = 0,
      guide = "colourbar") +
    geom_point(data = egfr_df,aes(y = value_mid,fill="EGFR"),color="coral")+
    scale_fill_manual(values = c("EGFR"="coral"),name=NULL)+
    geom_hline(yintercept = 0,
               color = "black",
               size = .5)+
    geom_vline(xintercept = cumu_chr_pos,
               color = "black",
               size = .3)+
    geom_segment(data = segments,
                 aes(y = arm_value_med, 
                     x = start_pos, 
                     xend = end_pos),
                 color=segments$arm_col,linewidth = 2)+
    geom_vline(
      xintercept = cumu_pq_pos,
      color = "black",
      size = .3,
      linetype = "dotted"
    )+
    geom_vline(xintercept = 0,
               color = "black",
               size = 1) +
    geom_vline(
      xintercept = sum(as.numeric(genome_size$X1)),
      color = "black",
      size = 1
    ) +
    xlab(label = "Chromosome") +
    ylab(label = "CNV score") +
    scale_x_continuous(breaks=xtick_pos,labels = genome_anno$X1) +
    theme(axis.text.x= element_text(size=9,angle = 45),
          panel.background = element_rect(fill = "white", colour = "white"))+
    ggtitle(paste0(sample_name,"- EGFR ", EGFRgroup))
  chrom_plot
})
chromplot_list
ggsave(
  plot = chromplot_list,
  filename = save_here(object_name = "CNV_plots_peptides_20260512_RRS_1296_PCF_Phopsho.pdf"),
  width = 12,
  height = 6
)
