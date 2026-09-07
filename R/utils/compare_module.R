### Describe the function of this script here
# Author: Dennis Friedel
# Date: 2025-12-29
# Modification: 2026-01-13
# Description of functions in this script:
# - dep_gsea_wrp_wrapper: Wrapper function that performs differential expression analysis using limma for all against one comparisons, generates summary plots (upset plot, volcano plot, heatmap), and conducts gene set enrichment analysis using fgsea for a given SummarizedExperiment object and a specified group column. The function returns a list containing the results of the analyses and the generated plots.
# - limma_analysis: Function that performs differential expression analysis using limma for a given SummarizedExperiment object, with options for control group comparison or all group comparisons, and returns the results in the rowData of the SummarizedExperiment object.
# - run_limma_group_comparison: Function that performs differential expression analysis using limma for a specified group comparison in a SummarizedExperiment object, and returns a data frame containing the results of the analysis, including log fold changes, p-values, adjusted p-values, and significance status for each protein.
# - convert_toptable_to_envo_input: Function that converts the output of a limma toptable into a format suitable for enhanced volcano plotting, allowing for customization of column names and test identifiers.
# - plot_enhanced_volcano: Function that creates an enhanced volcano plot from the output of convert_toptable_to_envo_input, allowing for customization of p-value and log fold change thresholds, as well as point size and color gradients.
# - plot_wrp_upset_plot: Function that generates an upset plot summarizing the significant differentially expressed proteins across all contrasts in an AAGO analysis.
# - plot_wrp_lfc_overview_heatmap: Function that creates a heatmap overview of log fold changes for the significant differentially expressed proteins across all contrasts in an AAGO analysis.
# - plot_correlation_heatmap: Function that generates a correlation heatmap of the significant differentially expressed proteins across samples in a SummarizedExperiment object, with options for hierarchical clustering methods and distance metrics.
# - run_fgsea_analysis: Function that performs fast gene set enrichment analysis (fgsea) for a given list of gene sets and ranked lists of genes from differential expression analysis, returning the significant enriched gene sets based on specified thresholds.
# - plot_fgsea_ballon: Function that creates a balloon plot to visualize the results of fgsea analysis, showing the top enriched gene sets for each group in the analysis.
# - rank_dea_result: Function that ranks the results of differential expression analysis based on log fold changes, returning a named vector suitable for input into fgsea analysis.
# - get_misg_test_db: Function that retrieves a test database of gene sets for use in gene set enrichment analysis, returning a list of gene sets categorized by database type (e.g., GO, KEGG, Reactome).
# - download_ependymoma_preproc_diann_rds: Function that downloads a preprocessed and normalized SummarizedExperiment object for ependymoma data from a specified source, returning the object for use in differential expression analysis and gene set enrichment analysis.
# =============================================================#

"%>%" <- magrittr::"%>%"

###### Wrapper #####

#' wrapper_dea_gsea
#'
#' @param ms_se SummarizedExperiment object; preprocessed and normalized data ready for differential expression analysis
#' @param wrp_group character; name of the column in the colData of ms_se that contains the group information for the AAGO analysis
#' @param wrp_pv_fil numeric; p-value threshold for filtering significant differentially expressed proteins in the AAGO analysis
#' @param wrp_fc_fil numeric; log2 fold change threshold for filtering significant differentially expressed proteins in the AAGO analysis
#' @param corhm_linkage character; method for hierarchical clustering in the correlation heatmap (e.g., "ward.D2", "complete", "average")
#' @param corhm_dist character; method for calculating distance in the correlation heatmap (e.g., "spearman", "pearson", "euclidean")
#' @param gene_se_catalouge list; list of gene sets to be used for the gene set enrichment analysis, where each element is a list of gene sets for a specific database (e.g., GO, KEGG, Reactome)
#' @param wrp_use_padj logical; whether to use adjusted p-values for filtering significant differentially expressed proteins in the AAGO analysis
#'
#' @returns List; list containing the results of the differential expression analysis and the gene set enrichment analysis, as well as the generated plots.
#' @export 
#'
#' @examples NULL
#'   
#' biological_list <- get_misg_test_db()
#' se_proc<-download_ependymoma_preproc_diann_rds()
#' aago_res<-wrapper_dea_gsea(
#'   ms_se = se_proc,
#'   wrp_mode = "AAGO",
#'   wrp_test = NULL,
#'   wrp_contrast = NULL,
#'   wrp_group = "group",
#'   gene_set_catalouge = biological_list,wrp_fgsea = F
#' )
## biological_list <- get_misg_test_db()
## se_proc<-readRDS("./projects/testMQ/rds/ptm_se_id_3_preproc.rds")
## manual_res<-wrapper_dea_gsea(
##   ms_se = se_proc,
##   wrp_mode = "MANUAL",
##  wrp_test = c("EGF"),
##   wrp_contrast = list(c("woEGF")),
##   wrp_group = "group",wrp_use_padj = F,
##  gene_set_catalouge = biological_list,wrp_fgsea = F
## )
wrapper_dea_gsea <- function(ms_se,
                             wrp_mode = c("AAGO", "MANUAL"),
                             wrp_test = NULL,
                             wrp_contrast = NULL,
                             wrp_group = "group",
                             wrp_pv_fil = 0.05,
                             wrp_fc_fil = 0.58,
                             corhm_linkage = "ward.D2",
                             corhm_dist = "spearman",
                             wrp_minSize  = 10,
                             wrp_maxSize  = 500,
                             wrp_nPermSimple = 5000,
                             wrp_use_padj = T,
                             wrp_fgsea = F,
                             gene_set_catalouge=NULL) {
  
  #=============================================================#
  #### 0.Assert that input is correct  #######
  #=============================================================#
  wrp_list <- list()
  
  assertthat::assert_that(inherits(ms_se, "SummarizedExperiment"))
  assertthat::assert_that(is.character(wrp_group), length(wrp_group) == 1)
  assertthat::assert_that(is.numeric(wrp_pv_fil), length(wrp_pv_fil) == 1)
  assertthat::assert_that(is.numeric(wrp_fc_fil), length(wrp_fc_fil) == 1)
  assertthat::assert_that(is.character(corhm_linkage), length(corhm_linkage) == 1)
  assertthat::assert_that(is.character(corhm_dist), length(corhm_dist) == 1)
  
  if(is.null(gene_set_catalouge)){
    message("Deactivate FGSEA")
    wrp_fgsea<- F
  }
  
  if(wrp_fgsea){
    assertthat::assert_that(is.list(gene_set_catalouge))
    assertthat::assert_that(all(sapply(gene_set_catalouge, is.list)))
    assertthat::assert_that(all(sapply(gene_set_catalouge, function(x)
      all(sapply(x, is.character)))))
  }
  
  
  #=============================================================#
  #### 1.Run DEA (Differential Expression Analysis)  #######
  #=============================================================#
  if(wrp_mode=="AAGO"){message("Run All Against One (AAGO) Moderated T-test using Bayesian statistics (Limma)")}
  if(wrp_mode!="AAGO"){message("Run Moderated T-test using Bayesian statistics (Limma) on manually provided groups")}
  tts_combined_new <- run_limma_se(
    limma_se = ms_se,
    group = wrp_group,
    pv_fil = wrp_pv_fil,
    fc_fil = wrp_fc_fil,
    use_padj = wrp_use_padj,
    mode = wrp_mode,
    lm_contrast = wrp_contrast,
    lm_test = wrp_test
  )
  
  # Stop when there is no significant protein and return empty list with message
  if (sum(tts_combined_new$significant) == 0) {
    message(
      "No significant differentially expressed proteins found with the given thresholds. Please adjust the p-value and log fold change thresholds and try again."
    )
    wrp_list <- list(tts_combined_new)
    return(wrp_list)
  }
  
  #=============================================================#
  #### 2.DEA plots  #######
  #=============================================================#
  
  message("Create Summary Plots ")
  #####  2.1. Upset plot #####
  wrp_list$tt_combined <- tts_combined_new
  wrp_list$upset_plot <- plot_wrp_upset_plot(wrp_res = tts_combined_new)
  
  ##### 2.2. Plot Volcanoplot showing the top expressed proteins ####
  message("Volcano Plot ")
  enhanced_voclanos <- plot_enhanced_volcano(
    limma_res = tts_combined_new,
    tests_end_with = "logFC",
    colname_adjP = "adj_P_Val",
    colname_LogFC = "logFC",
    colname_pvalue = "P_Value",
    pointsize = 3,
    use_padj = wrp_use_padj,
    alpha = wrp_pv_fil,
    lfc =  wrp_fc_fil
  )
  wrp_list$enhanced_volcanos <- enhanced_voclanos
  
  ##### 2.3. Plot Heatmap showing the top expressed proteins ####
  message("LFC Overview ")
  lfc_overview_heatmap <- plot_wrp_lfc_overview_heatmap(wrp_res = tts_combined_new)
  wrp_list$lfc_overview_hm <- lfc_overview_heatmap
  
  ##### 2.4. Correlation Heatmap ####
  message("Correlation Heatmap ")
  corhm_plot <- plot_correlation_heatmap(
    se_object = ms_se,
    group = wrp_group,
    cmLink = corhm_linkage,
    cmDist = corhm_dist,
    goi = tts_combined_new$genes[tts_combined_new$significant]
  )
  wrp_list$correlation_heatmap <- corhm_plot
  
  if(is.null(wrp_fgsea)){
    return(wrp_list)
  }
  
  #=============================================================#
  #### 3. Gene set enrichment analysis  #######
  #=============================================================#
  if(wrp_fgsea&(!is.null(gene_set_catalouge))){
    message("Run Fast Gene Set enrichment Analysis ")

    tests <-colnames(tts_combined_new)[grep("_logFC", colnames(tts_combined_new))] %>%
      gsub("_logFC","",.)
    ranked_res_list <- purrr::map(tests, function(tn) {
      rank_dea_result(
        tmp_gs = tts_combined_new,
        colname_lfc = paste0(tn, "_logFC"),
        colname_feature = "genes"
      )
      
    })
    names(ranked_res_list) <- tests
    
    fgsea_db_results <- purrr::map(gene_set_catalouge, function(gene_set_list) {
      fgsea_res<-run_fgsea_analysis(
        gene_list = gene_set_list,
        ranks_list = ranked_res_list,
        minSize  = wrp_minSize,
        maxSize  = wrp_maxSize,
        nPermSimple = wrp_nPermSimple,
        padj_threshold = wrp_pv_fil
      )
    })
    
    ballon_plots <- purrr::map(fgsea_db_results, function(x)
      plot_fgsea_ballon(fgsea_res = x, terms_show = 10))%>%purrr::compact()
    wrp_list$fgsea_result <- ballon_plots
  }
  message("Done return results as list")
  return(wrp_list)
}

###### DEA Analysis Functions ######


#' design_matrix_from_se
#'
#' @param se_object SummarizedExperimentObject
#' @param group character
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#' design_matrix_from_se()
#' se_proc<-download_ependymoma_preproc_diann_rds()
#' design_matrix_from_se(se_proc,group = "group2")
design_matrix_from_se<-function(se_object,group){
  assertthat::assert_that(inherits(group,"character"),group%in%colnames(SummarizedExperiment::colData(se_object)))
  # Define design matrix and contrast matrix for all against one comparison 
  des_mat <- model.matrix(formula(paste0("~0+", group)), data = as.data.frame(SummarizedExperiment::colData(se_object)))
  colnames(des_mat) <- gsub(group, "", colnames(des_mat))
  des_mat
}

#' contrast_all_against_one
#'
#' @param group character; vector which contains  
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#' se_proc<-download_ependymoma_preproc_diann_rds()
#' contrast_all_against_one(group = se_proc$group2)
contrast_all_against_one<-function(group){
  
  assertthat::assert_that(inherits(group,"character"),
                          length(group)>1,
                          length(unique(group))>1)
  
  # Create contrast matrix for all against one comparison
  all_against_one <- purrr::map(unique(group), function(target) {
    others <- setdiff(group, target)
    contrast <- paste0(target,
                       "-(",
                       paste0(others, collapse = "+"),
                       ")/",
                       length(others))
  }) %>% unlist()
  names(all_against_one)<-paste0(unique(group),"_vs_all")
  all_against_one
}

#' contrast_manual
#'
#' @param group 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#'se_proc<-download_ependymoma_preproc_diann_rds()
#'
#' cm<-contrast_manual(group = se_proc$group2,test = c("D","D"),list(c("Myxo"),c("Myxo","SPN")))
#' cont_mat <- limma::makeContrasts(contrasts = cm, levels = design_matrix_from_se(se_object = se_proc,group = "group2"))
#' 
#' fit_all <- limma::lmFit(as.matrix(assay(se_proc)), design = des_mat)
#' fit_cont <- limma::contrasts.fit(fit_all, contrasts = cont_mat)
#' fit_bay <- limma::eBayes(fit_cont)
#' limma::topTable(fit_bay)
#'
contrast_manual<-function(group,test,contrast){
  
  assertthat::assert_that(inherits(group, "character"),
                          length(group) > 1,
                          length(unique(group)) > 1)
  assertthat::assert_that(inherits(contrast, "list"), length(test) == length(contrast))
  assertthat::assert_that(all(test%in%group))
  assertthat::assert_that(all(unlist(contrast)%in%group))
  
  # Create contrast matrix for all against one comparison
  limma_formula <- purrr::map2(test,contrast,function(target,contrast_vec) {
    contrast_forumla <- paste0(target,
                       "-(",
                       paste0(contrast_vec, collapse = "+"),
                       ")/",
                       length(contrast_vec))
    names(contrast_forumla)<-paste0(target,
                                    "_vs_",
                                    paste0(contrast_vec, collapse = "&"))
    contrast_forumla
  }) %>% unlist()
  limma_formula
}

#' run_limma_se
#' @author Dennis Friedel
#' @param limma_se SummarizedExperiment object; preprocessed and normalized data ready for differential expression analysis
#' @param group character; name of the column in the colData of limma_se that contains the group information for the comparison
#' @param pv_fil numeric; p-value threshold for filtering significant differentially expressed proteins in the comparison
#' @param fc_fil numeric; log2 fold change threshold for filtering significant differentially expressed proteins in the comparison
#' @param use_padj logical; whether to use adjusted p-values for filtering significant differentially expressed proteins in the comparison
#'
#' @returns data.frame; a data frame containing the results of the differential expression analysis, including log fold changes, p-values, adjusted p-values, and significance status for each protein.
#'
#'
#' @examples NULL
#' se_proc<-download_ependymoma_preproc_diann_rds()
#'wrp_res<-run_limma_se(limma_se = se_proc,group = "group",mode = "AAGO")
#'man_res<-run_limma_se(limma_se = se_proc,group = "group",mode = "MANUAL",lm_test = c("EPN_MPE","EPN_SPINE"),lm_contrast = list(c("EPN_SPINE"),c("EPN_MPE")))
run_limma_se <- function(limma_se,
                           group = "group",
                           mode=c("AAGO","MANUAL"),
                           lm_test=NULL,
                           lm_contrast=NULL,
                           pv_fil = 0.05,
                           fc_fil = 0.58,
                           use_padj = T) {
  
  # Assert that input is correct
  assertthat::assert_that(inherits(limma_se, "SummarizedExperiment"))
  assertthat::assert_that(is.character(group),
                          length(group) == 1,
                          group %in% colnames(limma_se@colData))
  assertthat::assert_that(any(mode%in%c("AAGO","MANUAL")),length(mode) == 1)
  assertthat::assert_that(is.numeric(pv_fil), length(pv_fil) == 1)
  assertthat::assert_that(is.numeric(fc_fil), length(fc_fil) == 1)
  assertthat::assert_that(is.logical(use_padj), length(use_padj) == 1)
  
  # Set name and id in rowdata to the rownames of the assay for merging with toptable later on
  limma_se@elementMetadata$name <- limma_se@elementMetadata$row_id
  limma_se@elementMetadata$id <- limma_se@elementMetadata$id
  
  # Define design matrix 
  des_mat<-design_matrix_from_se(se_object = limma_se,group = group)
  
  # Define contrast matrix 
  if(mode=="AAGO") {
    cont_mat <- limma::makeContrasts(contrasts = 
                    contrast_all_against_one(group = limma_se[[group]]),
                    levels = des_mat)
  }
  if (mode == "MANUAL") {
    assertthat::assert_that(!is.null(lm_test),!is.null(lm_contrast))
    cont_mat <- limma::makeContrasts(contrasts = contrast_manual(
      group = limma_se[[group]],
      test = lm_test,
      contrast = lm_contrast),levels = des_mat)
  }
  
  # Run Limma 
  fit_all <- limma::lmFit(as.matrix(SummarizedExperiment::assay(limma_se)), design = des_mat)
  fit_cont <- limma::contrasts.fit(fit_all, contrasts = cont_mat)
  fit_bay <- limma::eBayes(fit_cont)
  
  # Retrieve toptable for all contrasts
  tt_all <- limma::topTable(fit_bay, number = Inf)
  tt_all <- data.frame(genes = rownames(tt_all), tt_all)
  tt_all_list <- purrr::map(colnames(cont_mat), function(c) {
    df <- limma::topTable(fit_bay, coef = c, number = Inf)
    df <- data.frame(genes = rownames(df), df)
  })
  names(tt_all_list) <- colnames(cont_mat)
  
  # Combine toptable for all contrasts into one table and add significant column based on p-value and logFC thresholds
  if(use_padj){
    pv_col <- "adj.P.Val"
  } else{
    pv_col <- "P.Value"
  }
  
  tts_combined <- purrr::map(names(tt_all_list), function(b) {
    tt <- tt_all_list[[b]]
    tt[["significant"]] <- (tt[[pv_col]] <= pv_fil) &
      (abs(tt$logFC) >= (fc_fil))
    tt[["Padj_Rank"]] <- order(tt$P.Value, decreasing = F)
    colnames(tt)[-c(1)] <- paste(b, colnames(tt)[-c(1)], sep = "_")
    colnames(tt) <- gsub("\\.", "_", colnames(tt))
    tt
  }, .progress = T) %>% purrr::reduce(., dplyr::full_join)
  rownames(tts_combined) <- tts_combined$genes
  tts_combined_table <- tts_combined
  
  if(length(grep("significant", colnames(tts_combined_table)))>1){
    tts_combined_table$significant <- rowSums(tts_combined_table[, grep("significant", colnames(tts_combined_table))]) >
      0  
  }else{
    tts_combined_table$significant <- tts_combined_table[, grep("significant", colnames(tts_combined_table))]
  }

  return(tts_combined_table)
}

###### FGSEA Analysis Functions ######

#' rank_dea_result
#'
#' @param tmp_gs data.frame; Limma toptable 
#' @param colname_feature Character; Name of the column in tmp_gs containing features (e.g. Genes)
#' @param colname_lfc Character; Name of the column in tmp_gs containing logFold chnages (e.g. logFC)
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
rank_dea_result<-function(tmp_gs,colname_feature="genes",colname_lfc="logFC"){
  tmp_gs <- tmp_gs[, c(colname_feature, colname_lfc)]
  tmp_gs <- tmp_gs[order(tmp_gs[[colname_lfc]], decreasing = T), ]
  gsls <- tmp_gs[[paste0(colname_lfc)]]
  names(gsls) <- tmp_gs$genes
  gsls
}

#' run_fgsea_analysis
#'
#' @param pathways List of gene sets (named list of character vectors).
#' @param stats Named numeric vector of gene-level statistics.
#' @param minSize Integer; minimum gene set size.
#' @param maxSize Integer; maximum gene set size.
#' @param nPermSimple Integer; number of permutations for simple p-value estimation.
#' @param padj Numeric; adjusted p-value threshold.
#' @param padj_threshold Numeric; adjusted p-value cut-off for filtering results.
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
run_fgsea_analysis <- function(gene_list,
                               ranks_list,
                               minSize  = 10,
                               maxSize  = 500,
                               nPermSimple = 5000,
                               padj_threshold = 0.05) {
  assertthat::assert_that(!is.null(names(ranks_list)))
  
  ### Run FGSEA analysis for each ranked list
  fgsea_results <- purrr::map2(ranks_list, names(ranks_list), function(ranks, groups) {
    fgseaRes <- fgsea::fgsea(
      pathways = gene_list,
      stats    = ranks,
      minSize  = minSize,
      maxSize  = maxSize,
      nPermSimple = nPermSimple
    )
    fgseaRes$group <- groups
    as.data.frame(fgseaRes)
  }, .progress = T) %>% purrr::reduce(., dplyr::full_join)
  
  fgsea_results <- fgsea_results[fgsea_results$padj < padj_threshold, ]
  fgsea_results$logAPV <- -log10(fgsea_results$padj)
  fgsea_results <- fgsea_results[!is.na(fgsea_results$group), ]
  
  return(fgsea_results)
}


###### Plotting Functions ######
pastel_color_sets<-list(
  "skyblue_whiter_rawsalmon"=c("skyblue","white","salmon"),
  "softblue_whiter_salmon50degree10min"=c("#BFD7EA","#FFFDF7","#F7C6C7"),
  "pastelpurple_cream_teal"=c("#CDB4DB","#FAD3DD","#A8DADC"),
  "nature_style"=c("#A6CEE3","#F7F7F7","#FB9A99")
)

#' convert_toptable_to_envo_input
#'
#' @param limma_res 
#' @param tests_end_with 
#' @param colname_adjP 
#' @param colname_LogFC 
#' @param colname_pvalue 
#' @param ... 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
# tests_end_with="logFC"
# colname_adjP = "adj_P_Val"
# colname_LogFC = "logFC"
# colname_pvalue = "P_Value"
# convert_toptable_to_envo_input(limma_res = limma_res_input,
#                                colname_adjP = "adj_P_Val",
#                                colname_LogFC = "logFC",
#                                colname_pvalue = "P_Value")
convert_toptable_to_envo_input<-function(limma_res,
                                         tests_end_with="logFC",
                                         colname_adjP = "adjusted P-value",
                                         colname_LogFC = "logFC",
                                         colname_pvalue = "p.val",...){
  assertthat::assert_that(class(limma_res)=="data.frame")
  assertthat::assert_that(length(grep(colname_adjP,colnames(limma_res)))>0)
  assertthat::assert_that(length(grep(colname_LogFC,colnames(limma_res)))>0)
  assertthat::assert_that(length(grep(colname_pvalue,colnames(limma_res)))>0)
  
  test_names <- colnames(limma_res)[tidyselect::ends_with(var = colnames(limma_res), tests_end_with)] %>%
    gsub(tests_end_with, "", .)
  
  purrr::map(test_names,function(y){
    tmp_toptable <-
      limma_res[, dplyr::starts_with(vars = colnames(limma_res), y)]
    
    colnames(tmp_toptable) <-
      gsub(paste0(y), "", colnames(tmp_toptable))
    
    tmp_toptable <-
      dplyr::rename(
        tmp_toptable,
        c(
          "adjusted P-value" = colname_adjP,
          "log2FoldChange" = colname_LogFC,
          "p-value" = colname_pvalue
        )
      )
    tmp_toptable$group<-y
    tmp_toptable
  })
}

#' plot_enhanced_volcano
#'
#' @param limma_res data.frame; toptable in the format of the output of convert_toptable_to_envo_input
#' @param alpha Numeric; Define threshold for p-value in Volcano (Y-Axis)
#' @param lfc Numeric; Define threshold for Log2fold chnage in Volcano (X-Axis)
#' @param use_padj Boolean;use P.adjusted Value
#' @param pointsize Numeric; Define size of dots in the volcano
#' @param labsize Numeric; Define size of labs in the volcano
#' @param tests_end_with convert_toptable_to_envo_input;
#' @param colname_adjP convert_toptable_to_envo_input;
#' @param colname_LogFC convert_toptable_to_envo_input;
#' @param colname_pvalue convert_toptable_to_envo_input;
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
# limma_res = tts_combined_table_sig_tmp
# tests_end_with="logFC"
# colname_adjP = "adj_P_Val"
# colname_LogFC = "logFC"
# colname_pvalue = "P_Value"
# alpha = 0.05
# lfc = .58
# use_padj = T
# pointsize=5
# labsize=5
plot_enhanced_volcano<-function(limma_res,
                                alpha = 0.05,
                                lfc = .58,
                                use_padj = T,
                                pointsize=5,
                                labsize=5,
                                dot_colors=
                                  c(
                                    "#D3D3D3",# not significant
                                    "#A1C9F4",# FC only
                                    "#BFE3A1",# p-value only 
                                    "#FFB3BA"# both significant
                                  ),
                                ...){
  
  
  
  # Assert that input is correct
  assertthat::assert_that(class(limma_res)=="data.frame")
  assertthat::assert_that(length(grep("logFC",colnames(limma_res)))>0)
  assertthat::assert_that(length(grep("adj_P_Val",colnames(limma_res)))>0)
  assertthat::assert_that(length(grep("P_Value",colnames(limma_res)))>0)
  assertthat::assert_that(is.numeric(alpha), length(alpha) == 1)
  assertthat::assert_that(is.numeric(lfc), length(lfc) == 1)
  assertthat::assert_that(is.logical(use_padj), length(use_padj) == 1)
  assertthat::assert_that(is.numeric(pointsize), length(pointsize) == 1)
  assertthat::assert_that(is.numeric(labsize), length(labsize) == 1)
  
  # Convert to toptable to format for enhanced volcano
  limma_res_input<-limma_res
  list_limma_results<-convert_toptable_to_envo_input(limma_res = limma_res_input,...)

  # Create enhanced volcano plot for each contrast
  purrr::map(list_limma_results,function(tmp_toptable){
    
    A<-gsub("vs_.*"," ",tmp_toptable$group)%>%gsub("_"," ",.)%>%unique()
    B<-gsub(".*vs_"," ",tmp_toptable$group)%>%gsub("_"," ",.)%>%unique()
    if(A==B){
      B<-""
    }
    if (use_padj) {
      tmp_toptable$significant <-
        tmp_toptable$`adjusted P-value` < alpha &
        abs(tmp_toptable$log2FoldChange) > lfc
      use_this_pvalue <- "adjusted P-value"
      y_axis_label <- bquote(~-Log[10] ~ italic(adj.P))
    } else{
      tmp_toptable$significant <-
        tmp_toptable$`p-value` < alpha &
        abs(tmp_toptable$log2FoldChange) > lfc
      use_this_pvalue <- "p-value"
      y_axis_label <- bquote(~-Log[10] ~ italic(P))
    }
    
    ## Plot Enhanced Volcano
    enVo <-
      EnhancedVolcano::EnhancedVolcano(
        toptable = tmp_toptable,
        title = gsub("_"," ",unique(tmp_toptable$group)),
        subtitle = paste0(
          "Number significant proteins: ",
          sum(tmp_toptable$significant, na.rm = T)
        ),
        lab = rownames(tmp_toptable),
        selectLab = rownames(tmp_toptable)[tmp_toptable$significant],
        x = 'log2FoldChange',
        y = use_this_pvalue,
        col= dot_colors,
        ylim = c(0,max(tmp_toptable$'log2FoldChange')+1.5),
        pCutoff = alpha,
        FCcutoff = lfc,
        pointSize = pointsize,
        colAlpha = 0.7,
        legendPosition = "top",
        legendLabSize = 12,
        legendIconSize = 5.0,
        boxedLabels = F,
        drawConnectors = T,
        widthConnectors = 0.5,
        endsConnectors = "first",
        colConnectors = "grey10",
        max.overlaps = 15,
        maxoverlapsConnectors = NULL,
        min.segment.length = 20
      ) +
      ggplot2::annotate(
        "text",
        x = max(tmp_toptable$'log2FoldChange'),
        y = 0,
        label = paste0(paste0(A)),
        size = ggplot2::unit(10, "pt"),
        parse = F
      ) +
      ggplot2::annotate(
        "text",
        x = min(tmp_toptable$'log2FoldChange'),
        y = 0,
        label = paste0(B),
        size = ggplot2::unit(10, "pt"),
        parse = F
      )+
      ggplot2::ylab(y_axis_label)
    enVo
  })
}

#' plot_wrp_upset_plot
#' @author Dennis Friedel
#' @param wrp_res 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#' se_proc<-download_ependymoma_preproc_diann_rds()
#' wrp_result<-run_limma_se(limma_se = se_proc,group="group",pv_fil=0.05,fc_fil=0.58,use_padj=T)
#' plot_wrp_upset_plot(wrp_result)
plot_wrp_upset_plot<-function(wrp_res){
  
  suffix <- "_significant"
  if(length(grep(suffix, colnames(wrp_res)))<=1){
    return()
  }
  
  if (length(grep(suffix, colnames(wrp_res)))>=2){
    message("Upset Plot ")
    
    tt_binary <- wrp_res[wrp_res$significant, grep(suffix, colnames(wrp_res))]
    colnames(tt_binary) <- gsub(suffix, "", colnames(tt_binary))
    tib_tt_binary <- dplyr::tibble(tt_binary)
    
    upset_dep <- ComplexUpset::upset(
      data = tib_tt_binary,
      intersect = colnames(tib_tt_binary),
      name = "Common differential Expressed Proteins among layers",
      min_size = 5
    )
  }
}

#' plot_wrp_lfc_overview_heatmap
#' @author Dennis Friedel
#' @param wrp_res 
#' @examples 
#'se_proc<-download_ependymoma_preproc_diann_rds()
#'wrp_res<-run_limma_se(limma_se = se_proc,group = "group",mode = "AAGO")
#'plot_wrp_lfc_overview_heatmap(wrp_res)
#'
plot_wrp_lfc_overview_heatmap<-function(wrp_res,color_set=pastel_color_sets$nature_style){
  
  ### Assert that it is a AAGO table
  assertthat::assert_that(inherits(wrp_res, "data.frame"),
                          "significant"%in%colnames(wrp_res)
                          )
  
  if(length(grep("logFC", colnames(wrp_res)))==1){
    return()
  }
  
  heatmap_column_title <- paste0(
    "LogFC of differnital expressed Proteins (N = ",
    sum(wrp_res$significant),
    " / ",
    nrow(wrp_res),
    ")\n"
  )
  
  wrp_sig<-wrp_res[wrp_res$significant,]
  lfc_mat <- as.matrix(wrp_sig[, grep("logFC", colnames(wrp_sig))])
  colnames(lfc_mat)<-gsub("_"," ",colnames(lfc_mat))
  
  lim<-stats::quantile(abs(lfc_mat),0.99,na.rm=T)
  col_fun<-circlize::colorRamp2(
    c(-lim,0,lim),
    color_set
  )
  
  lfc_overview_hm <-
    ComplexHeatmap::Heatmap(
      lfc_mat,
      cluster_column_slices = F,
      column_title = heatmap_column_title,
      column_names_side = "bottom",
      cluster_columns = F,
      column_names_rot = 45,
      row_names_gp = grid::gpar(fontsize = 6),
      show_column_names = TRUE,
      show_column_dend = F,
      show_row_names = F,
      width = ggplot2::unit(4, "cm"),
      name = "LogFC",
      border = T,
      col = col_fun
    )
  complexheatmap_to_ggplot(lfc_overview_hm)
}


#' plot_correlation_heatmap
#' @author Dennis Friedel
#' @param se_object SummarizedExperiment object; preprocessed and normalized data ready for differential expression analysis
#' @param group Character; Column name of the group that will be shown as top annotation of the Heatmap
#' @param goi Vector: Feature names of the SummarizedExperiment object that will be used for subsetting 
#' @param cmDist Character; Distance to use for estimation of correaltion
#' @param cmLink Character; Linkage algorithm to use for estimation of clustering
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#se_proc<-download_ependymoma_preproc_diann_rds()
# wrp_res<-run_limma_se(limma_se = se_proc,group = "group",mode = "AAGO")
# plot_correlation_heatmap(
#  se_object = se_proc,
#  group = "group",
#  cmLink = "ward.D2",
#  cmDist = "spearman",
#  goi = wrp_res$genes[wrp_res$significant]
# )
plot_correlation_heatmap<-function(
    se_object,
    group,
    goi,
    cmDist="spearman",
    cmLink="ward.D2",
    color_set=pastel_color_sets$nature_style
){
  assertthat::assert_that(inherits(se_object, "SummarizedExperiment"))
  assertthat::assert_that(is.character(group), length(group) == 1)
  assertthat::assert_that(is.vector(group))
  assertthat::assert_that(is.character(cmDist), length(cmDist) == 1,any(cmDist%in%c("pearson", "kendall", "spearman")))
  assertthat::assert_that(is.character(cmLink), length(cmLink) == 1,any(cmLink%in%c("ward.D", "single", "complete", "average", 
                                                                          "mcquitty", "median", "centroid", "ward.D2")))
  
  ## Return when goi is null (makes only sense in App)
  if(sum(goi%in%rownames(se_object))<=1){return(NULL)}
  
  corassay <- cor(SummarizedExperiment::assay(se_object)[goi, ], method = cmDist)
  
  col_fun<-circlize::colorRamp2(
    c(min(corassay),median(corassay),max((corassay))),
    color_set
  )
  
  cor_hm <- ComplexHeatmap::Heatmap(
    name = paste0(cmDist," Correlation"),
    matrix = corassay,show_column_names = F,
    top_annotation = ComplexHeatmap::HeatmapAnnotation("Condition"=gsub("_"," ",SummarizedExperiment::colData(se_object)[[group]])),
    clustering_distance_rows = function(x)
      as.dist(1 - x),
    clustering_distance_columns = function(x)
      as.dist(1 - x),
    clustering_method_columns = cmLink,
    clustering_method_rows = cmLink,
    col=col_fun
  )
  ggpubr::as_ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(cor_hm,
                               heatmap_legend_side="left",
                               annotation_legend_side="left",
                               legend_grouping = "original")))  
}

#' make_longnames_nice
#'
#' @param x Character; Vector of character strings that will be separated into two lines if they are longer than 30 characters
#'
#' @returns Character; Vector of character strings that have been separated into two lines if they are longer than 30 characters
#' @export
#'
#' @examples NULL
make_longnames_nice<-function(x,nc=25){
  # x<-"GOBP MYELOID CELL DIFFERENTIATION"
  # nc=25
  # get the length of each word separated by " "
  lenght_words<-purrr::map(strsplit(x," "),function(y)nchar(y))%>%unlist()+1
  
  # estimate the cumulative sum of the length of the words and separate the names at the point where the cumulative sum is longer than 30 characters
  lenght_words_cumsum<-cumsum(lenght_words)
  
  position<-lenght_words_cumsum[which(lenght_words_cumsum<=nc)[length(which(lenght_words_cumsum<=nc))]]
  
  paste0(substr(x, 1, position), "\n", substr(x, position+1, nchar(x)))
  
}

#' plot_fgsea_ballon
#'
#' @param fgsea_res 
#' @param terms_show 
#' @param up_color 
#' @param mid_color 
#' @param down_color 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
# se_proc<-download_ependymoma_preproc_diann_rds()
# se_proc<-gg_se_id_2_preproc
# biological_list <- get_misg_test_db()
# wrp_result<-run_limma_se(limma_se = se_proc,mode=c("AAGO"),group="group",pv_fil=0.05,fc_fil=0.58,use_padj=T)
# 
# tests <-colnames(wrp_result)[grep("_logFC", colnames(wrp_result))] %>%
#   gsub("_logFC","",.)
# ranked_res_list <- map(tests, function(tn) {
#   rank_dea_result(
#     tmp_gs = wrp_result,
#     colname_lfc = paste0(tn, "_logFC"),
#     colname_feature = "genes"
#   )
# })
# names(ranked_res_list)<-tests
# fgsea_res_test<-run_fgsea_analysis(
# gene_list = biological_list$`C2.CP:KEGG_MEDICUS`,
# ranks_list = ranked_res_list,minSize = 10,maxSize = 50
# )
# fgsea_res_test
# #ggsave(plot=plot_fgsea_ballon(fgsea_res = fgsea_res_test,terms_show = 5),"test.pdf")
plot_fgsea_ballon<-function(fgsea_res,terms_show,up_color="#FFB3BA",mid_color="#D3D3D3",down_color="#A1C9F4"){
  # fgsea_res = fgsea_res_test
  # terms_show = 5
  
  assertthat::assert_that(inherits(fgsea_res, "data.frame"))
  assertthat::assert_that(is.character(up_color), length(up_color) == 1)
  assertthat::assert_that(is.character(mid_color), length(mid_color) == 1)
  assertthat::assert_that(is.character(down_color), length(down_color) == 1)
  assertthat::assert_that(is.numeric(terms_show))
  
  
  # 'get top term per cluster'
  top_terms <- purrr::map(unique(fgsea_res$group), function(x) {
    tmp_df <- fgsea_res[fgsea_res$group %in% x, ]
    head(tmp_df$pathway[order(tmp_df$logAPV, decreasing = T)], terms_show)
  })%>%unlist()%>%unique()
  
  fgsea_res_final <- fgsea_res[fgsea_res$pathway %in% top_terms, ]
  
  # Make Names short for better visualization by substituting the frist word behind a "_" with an empty string
  fgsea_res_final$pathway <- gsub("_", " ", fgsea_res_final$pathway)
  fgsea_res_final$group <- gsub("_", " ", fgsea_res_final$group)
  
  # Separate pathway names with names taken togehter are longer than 30 characters into two lines for better visualization
  nc_n<-25
  fgsea_res_final$pathway<-purrr::map(fgsea_res_final$pathway,function(x)make_longnames_nice(x,nc=nc_n))%>%unlist()
  
  signficiant_ballon_plots <- ggpubr::ggballoonplot(fgsea_res_final,
                                                    "pathway",
                                                    "group",
                                                    fill = "NES",
                                                    size = "logAPV") + 
    ggplot2::scale_fill_gradientn(colours = c(down_color, mid_color, up_color)) +
    ggplot2::coord_flip()
  
  return(signficiant_ballon_plots)
}

###### Load comparison module test data function ######

#' Download function to obtain test data for testing functions of the  `proteoR` package.
#' 
#' This dataset is larger DIA (Data-Independent Acquisition) dataset of 
#' Ependymomas which has been analyzed by DIA-NN.
#' The data has been publihsed and demonstrates the use of proteomics in 
#' identifiying novel IHC (Immunohistochemsitry) markers for improving brain cancers
#' diagnostics. The dataset can here be used demo data for testing the functions in the Package 
#' It contains information about proteins and peptides their corresponding proteins, 
#' and their quantification across different samples.
#'
#' @param path 
#' 
download_ependymoma_preproc_diann_rds<-function(path=paste0(tempdir(),"/diann_example/")){
  dir.create(path)
  
  url1<-"https://heibox.uni-heidelberg.de/f/2be0fb7e9e744ed68b69/?dl=1"
  dest1 <- file.path(path,"se_proc.rds")
  if(!file.exists(dest1)){
    download.file(url = url1,destfile = dest1)
  }
  
  se_object_proc<-readRDS(dest1)
  
  return(se_object_proc)
}

#' Download function to obtain msigdb data for testing functions of the  `proteoR` package.
#' 
#' This list conatins genesets from the MsigDB Databases 
#'
#' @param path 
#' 
get_misg_test_db<-function(subcategories=c("C5.GO:BP","C5.GO:CC","C2.CP:REACTOME","C2.CP:KEGG_LEGACY","C2.CP:KEGG_MEDICUS","H.")){
  
  msigdbr_coll<-msigdbr::msigdbr_collections()%>%as.data.frame()
  categories<-paste0(msigdbr_coll$gs_collection,".",msigdbr_coll$gs_subcollection)
  assertthat::assert_that(all(subcategories%in%categories))
  
  get_all_misg<-purrr::map(subcategories,function(x){
    
    msig_cat<-gsub("[.].*","",x)
    
    msig_subcat<-gsub(paste0(msig_cat,"[.]"),"",x)
    if(msig_subcat==""){msig_subcat<-NULL}
    
    misg_df<-msigdbr::msigdbr(species = "Homo sapiens",
                              category = msig_cat,
                              subcategory = msig_subcat)  
    
    misg_df <- split(x = misg_df$gene_symbol, f = misg_df$gs_name)
  },.progress = T)
  
  names(get_all_misg)<-subcategories
  
  return(get_all_misg)
}

