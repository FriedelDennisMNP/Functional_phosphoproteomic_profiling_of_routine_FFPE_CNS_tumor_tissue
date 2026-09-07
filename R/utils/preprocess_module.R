### Describe the function of this script here
# Author: Dennis Friedel
# Date: 2025-12-29
# Modification: 2026-01-13
# Description of functions in this script:
# - pre_processing_wrapper: Wrapper function for filtering and imputation of ms data.
# - qp_plots_wrapper: Wrapper function for quality control plots of ms data.
# - filter_features: Filter features based on percentage of identified values.
# - filter_samples: Filter samples based on percentage of identified values.
# - filter_ms_data: Wrapper function for filtering samples and features of ms data.
# - normalize_data: Normalize ms data using different methods.
# - correct_mse_batch_effect: Correct batch effects in ms data using ComBat.
# - proteom_cohort_summary: Summarize cohort statistics of ms data.
# - calculate_cv: Calculate coefficient of variation for ms data.
# - count_idenitfied_features: Count number of identified features per sample or feature.
# - get_intensity_long: Get intensity values in long format for plotting.
# - get_coefficient_variation: Get coefficient of variation per sample or group.
# - summarize_qc_metrics: Summarize quality control metrics of ms data.
# - plot_ggbarplot: Plot barplot of identified features per sample or group.
# - plot_horizontal_ggboxplot: Plot horizontal boxplot of coefficient of variation or intensity distribution.
# - plot_feature_overlap: Plot feature overlap between samples or groups.
# - plot_density_missing_features: Plot density of missing features per sample or group.
# - plot_normalized_dist: Plot distribution of normalized intensity values per sample or group.
# - plot_correlation_heatmap: Plot correlation heatmap of samples based on intensity values.
# - plot_pca: Plot PCA of samples based on intensity values.
# - plot_umap: Plot UMAP of samples based on intensity values.
# - plot_tsne: Plot t-SNE of samples based on intensity values.
# =============================================================#

#### Required libraries
### Load required packages and Functions

"%>%" <- magrittr::"%>%"

#### Wrappers #####

#' pre_processing_wrapper
#'
#' @param se_object SummarizedExperimentObject. Containing Protein/peptide Intensities n Assay, Sample Information in ColData and Protein infromaitn in RowData
#' @param thr_sample Numeric. 
#' @param thr_feature Numeric.
#' @param group_sel Character. Name of the Column in ColData to use for Fractionized Missing Value Filtering  
#' @param filter_fractioned Boolean. TRUE for Fractionized Missing Value Filtering default = FASLE
#' @param imputation_method Character. Name of the imputation method to use 
#' @param normalization_method Character. Name of the normalization method to use (none, median_center, vsn, quantile)
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#' 
## se_object = readRDS("/mnt/NAS4/user_data/np-dennis/projects/AG_Sahm_projects/Spatial_Proteomics_RRS/projects/20260702_RRSSpatialProteomics_samples_0726/rds/ggu_se_raw.rds")
## anno<-openxlsx::read.xlsx("/mnt/NAS4/user_data/np-dennis/projects/AG_Sahm_projects/Spatial_Proteomics_RRS/projects/20260702_RRSSpatialProteomics_samples_0726/sample_sheets/sample_sheet_anno4.xlsx")
## 
## se_object<-se_object[,colnames(se_object)%in%anno$original_id]
## rownames(anno)<-anno$original_id
## SummarizedExperiment::colData(se_object)<-S4Vectors::DataFrame(anno[colnames(se_object),])
## colnames(se_object)<-se_object$sample_id
## 
## qc_test<-pre_processing_wrapper(se_object = se_object,
##                        thr_sample = 0.4,
##                        thr_feature = 0.5,
##                        group_sel = "group",
##                        filter_fractioned = TRUE,
##                        imputation_method = "MinDet",
##                        normalization_method = "none")
pre_processing_wrapper<-function(se_object,
                                 thr_sample,
                                 thr_feature,
                                 group_sel,
                                 filter_fractioned,
                                 imputation_method,
                                 normalization_method = "none"
){
  ##### Filter #####
  se_object_flt1 <- se_object[, colSums(!is.na(SummarizedExperiment::assay(se_object))) > thr_sample]
  if(is.null(validate_prc_object(se_object = se_object_flt1))){return(NULL)}
  
  se_object_flt2 <- filter_ms_data(
    se_object = se_object_flt1,
    assay_name = "intensity",
    threshold_sample = 0,
    threshold_feature = thr_feature,
    group = group_sel,
    filter_by_group = filter_fractioned,
    verbose = TRUE
  )
  if(is.null(validate_prc_object(se_object = se_object_flt2))){return(NULL)}

  ##### Normalize #####
  se_object_proc <- se_object_flt2
  if (!identical(normalization_method, "none")) {
    se_object_proc <- normalize_data(
      se_object = se_object_flt2,
      assay_name = "intensity",
      by_method = normalization_method
    )
  }
  
  ##### Impute #####
  if(!anyNA(SummarizedExperiment::assay(se_object_proc))){
    se_proc_list <- list("unfilt" = se_object,
                         "filt" = se_object_proc)
    attr(se_proc_list, "normalization_method") <- normalization_method
    return(se_proc_list)
  }
  
  ###Zero
  if (imputation_method == "zero") {
    tmp_assay <- SummarizedExperiment::assay(se_object_proc)
    tmp_assay[is.na(tmp_assay)] <- NA
    se_object_imp <-
      add_assay(se_object_proc,
                new_assay = as.matrix(tmp_assay),
                new_assay_name = "imputed") %>%
      set_primary_assay(., primary_assay_name = "imputed")
  }
  
  ###MAR
  if (imputation_method == "SVD" |
      imputation_method == "KNN") {
    assay_imp <-
      imputeLCMD::impute.MAR.MNAR(
        dataSet.mvs = SummarizedExperiment::assay(se_object_proc),
        method.MAR = imputation_method,
        model.selector = rep(0, nrow(se_object_proc))
      )
    se_object_imp <-
      add_assay(se_object_proc,
                new_assay = as.matrix(assay_imp),
                new_assay_name = "imputed") %>%
      set_primary_assay(., primary_assay_name = "imputed")
  }
  
  ###MNAR
  if (imputation_method == "MinDet" |
      imputation_method == "MinProb" |
      imputation_method == "QRILC") {
    assay_imp <-
      imputeLCMD::impute.MAR.MNAR(
        dataSet.mvs = SummarizedExperiment::assay(se_object_proc),
        method.MNAR = imputation_method,
        model.selector = rep(0, nrow(se_object_proc))
      )
    
    se_object_imp <-
      add_assay(se_object_proc,
                new_assay = as.matrix(assay_imp),
                new_assay_name = "imputed") %>%
      set_primary_assay(., primary_assay_name = "imputed")
  }
  
  ##### Summarize #####
  se_proc_list <- list("unfilt" = se_object,
                       "filt" = se_object_proc,
                       "imp" = se_object_imp)
  attr(se_proc_list, "normalization_method") <- normalization_method
  return(se_proc_list)
}

#' qp_plots_wrapper
#'
#' @param se_proc_ls 
#' @param group_sel 
#' @param thr_sample 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
## se_proc_ls = qc_test
## group_sel = "group"
## thr_sample = 5000
##qp_plots_wrapper(se_proc_ls = qc_test,group_sel = "group",thr_sample = 5000)
qp_plots_wrapper<-function(se_proc_ls,group_sel,thr_sample){
  
  plotlist<-list()
  normalization_method <- attr(se_proc_ls, "normalization_method")
  if (is.null(normalization_method)) {
    normalization_method <- "none"
  }

  filt_assay_name <- if ("intensity_norm" %in% SummarizedExperiment::assayNames(se_proc_ls$filt)) {
    "intensity_norm"
  } else {
    "intensity"
  }

  pl_norm2_title <- if (identical(normalization_method, "none")) {
    "Filtered data"
  } else {
    paste0("Filtered data (", normalization_method, " normalization)")
  }
  
  # Sepcial case for shiny app
  if(is.null(se_proc_ls)){return(NULL)}
  
  # Check if filtering left NAs
  if(!anyNA(SummarizedExperiment::assay(se_proc_ls$filt))){
    noNA_matrix<-T
  }else{
    noNA_matrix<-F
  }
  
  
  ### Before filtering
  idplot_list <- purrr::map(se_proc_ls[names(se_proc_ls)%in%c("unfilt","filt")], function(se_tmp) {
    
    cohort_stats <- proteom_cohort_summary(se_object = se_tmp) %>% ggpubr::ggtexttable()
    df_long <- count_idenitfied_features(
      se_object = se_tmp,
      assay_name = "intensity",
      group = group_sel,
      column_wise = T,
      in_percent = F
    )
    
    plot_id2 <-
      plot_ggbarplot(
        df_long = df_long,
        by_group = TRUE,
        feature_name = "Identified features ",
        threshold_value = thr_sample
      )
    prot_id <- patchwork::wrap_plots(cohort_stats,plot_id2,ncol = 2)
  })
  
  plotlist$prc_id_unfilt<-idplot_list[[1]]
  plotlist$prc_id_filt<-idplot_list[[2]]
  
  ##### Intensity Distribution ####
  norm_dist <- plot_normalized_dist(
    se_object = se_proc_ls$unfilt,
    assay_name = "intensity",
    group = group_sel
  )
  norm_dist2 <- plot_normalized_dist(
    se_object = se_proc_ls$filt,
    assay_name = filt_assay_name,
    group = group_sel
  ) + ggplot2::ggtitle(pl_norm2_title)
  plotlist$pl_norm <- norm_dist
  plotlist$pl_norm2 <- norm_dist2
  
  
  ##### Missing value distributuion ####
  if(noNA_matrix){
    plotlist$pl_feature_overlap <- NULL
    plotlist$pl_density_missing_features <- NULL
    se_tmp<-se_proc_ls$filt
  }else{
    plotlist$pl_feature_overlap <- (plot_feature_overlap(se_object = se_proc_ls$filt, assay_name = "intensity"))
    plotlist$pl_density_missing_features <- (plot_density_missing_features(se_object = se_proc_ls$filt, assay_name = "intensity"))
    se_tmp<-se_proc_ls$imp
  }
  
  ##### Correlation Heatmap ####
  use_linkage <- "ward.D2"
  corhm <- ComplexHeatmap::Heatmap(
    matrix = stats::cor(SummarizedExperiment::assay(se_tmp)[, ], method = "spearman"),
    name = "Correlation",
    clustering_distance_rows = function(x)
      stats::as.dist(1 - x),
    clustering_distance_columns = function(x)
      stats::as.dist(1 - x),
    clustering_method_columns = use_linkage,
    clustering_method_rows = use_linkage
  )
  corhm <- complexheatmap_to_ggplot(corhm)
  plotlist$corhm <- (corhm)
  return(plotlist)
}

#### Preprocessing functions ######

##### Data Pre-processing #####

#' filter_features
#' @author Dennis Friedel
#' @param mat matrix; assay or matrix which contains NAs
#' @param grps character; vector indicating the group for groupwise filtering
#' @param percent numeric; Percentage of proteins without missing values
#' @param n integer; Number of conditions which have more than Percentage proteins
#' @param assay Character; name of the assay to operate on.
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#filter_features(se_object = pg_se,assay_name = "intensity",threshold = 0.8,group = "group",by_group = TRUE,verbose = TRUE)
filter_features <-
  function (se_object,
            assay_name = NULL,
            threshold = 0.5,
            group = NULL,
            by_group = TRUE,
            verbose = TRUE) {
    ## Select assay from assay_list
    mat <-
      return_matrix_from_se(se_object = se_object, assay_name = assay_name)
    assertthat::assert_that(is.numeric(mat), anyNA(mat))
    assertthat::assert_that(
      is.numeric(threshold),
      is.null(group) | is.character(group),
      is.logical(by_group)
    )
    
    ## Count measured values across all samples and filter features with to many NA
    mat_long <- reshape2::melt(mat)
    colnames(mat_long) <- c("feature", "sample_id", "value")
    
    identified_feature <- mat_long %>%
      dplyr::group_by(feature) %>%
      dplyr::summarise("percent_identified" = (sum(!is.na(value)) / length(feature))) %>%
      as.data.frame()
    
    if (!is.null(group) & by_group) {
      assertthat::assert_that("sample_id" %in% colnames(SummarizedExperiment::colData(se_object)))
      mat_long$group <-
        plyr::mapvalues(mat_long[["sample_id"]], se_object[["sample_id"]], se_object[[group]]) %>%
        as.factor(.)
      identified_feature <- mat_long %>%
        dplyr::group_by(group, feature) %>%
        dplyr::summarise("percent_identified" = (sum(!is.na(value)) / length(group))) %>%
        as.data.frame()
    }
    # Keep features that are x% identified in the dedicated group
    keep_feature <-
      identified_feature$feature[identified_feature$percent_identified > threshold] %>%
      unique() %>%
      as.character()
    
    if (verbose)
      message("Keep features: ", length(keep_feature))
    
    return(se_object[keep_feature,])
  }

#' filter_samples
#' @author Dennis Friedel
#' @param se_object SummarizedExperiment; input data object.
#' @param assay_name Character; name of the assay to use.
#' @param threshold numeric; number between 0 and 1 indicating the maxium tolerated percentage, default is 0.4 (40%)
#' of missing values within some sample
#'
#' @return SummarizedExperiment object with columns filtered according to identified features
#' @export
#'
#' @examples NULL
#'filter_samples(se_object = pg_se,assay_name = "intensity",threshold = 0.4)
#'filter_samples(se_object = pg_se,assay_name = "intensity",threshold = 0.5)
filter_samples <-
  function (se_object,
            assay_name = NULL,
            threshold = 0.4,
            verbose = TRUE) {
    ## Select assay from assay_list
    mat <-
      return_matrix_from_se(se_object = se_object, assay_name = assay_name)
    assertthat::assert_that(is.numeric(mat), anyNA(mat))
    
    ## Count idenitfied values across all samples and filter samples
    na_mat <- !is.na(mat)
    keep_samples <-
      na_mat[, colSums(na_mat) / nrow(na_mat) >= threshold] %>% colnames()
    removed_samples <-
      colnames(na_mat)[!(colnames(na_mat) %in% keep_samples)]
    
    if (verbose & length(removed_samples) > 0)
      message(
        "removed ",
        length(removed_samples),
        " samples : ",
        paste0(removed_samples, collapse = "\n -")
      )
    
    ## return Summarized experiment
    se_object[, keep_samples]
  }

#' filter_ms_data
#' @author Dennis Friedel
#' @param se_object SummarizedExperiment object with columns filtered according to identified features
#' @param assay_name character; name of the assay to use for filtering 
#' @param threshold_sample numeric; number between 0 and 1 indicating the percent of at identified values per samples
#' @param threshold_feature numeric; number between 0 and 1 indicating the percent of at identified values per feature 
#' @param verbose logical; Show warnings and messages
#' @param filter_by_group logical; filter features per idenitfied features in biological group
#' @param group character; Indicating the group to use for feature filtering
#' @description
#' Wrapper function combining sample and feature filtering of ms data. 
#' 
#' @return SummarizedExperiment with features and samples filtered
#' @export
#'
#' @examples NULL
filter_ms_data <- function(se_object,
                           assay_name = NULL,
                           threshold_sample = 0.4,
                           threshold_feature = 0.5,
                           filter_by_group = FALSE,
                           group =NULL,
                           verbose = TRUE) {
  assertthat::assert_that(isClass(se_object,"SummarizedExperiment"),
                          is.null(group)|(is.character(group)&!is.null(se_object[[group]])),
                          is.numeric(threshold_sample),
                          is.numeric(threshold_feature),
                          is.logical(filter_by_group),
                          is.logical(verbose)
  )
  ## Filter samples
  se_proc_sample_object <-
    filter_samples(
      se_object = se_object,
      assay_name = assay_name,
      threshold = threshold_sample,
      verbose = verbose
    )
  ## Filter features
  se_proc_feature_object <-
    filter_features(
      se_object = se_proc_sample_object,
      assay_name = assay_name,
      threshold = threshold_feature,
      verbose = verbose,
      group = group,
      by_group = filter_by_group
    )
  se_proc_feature_object
}

#' normalize_data
#'
#' @param se_object SummarizedExperiment object with columns filtered according to identified features
#' @param assay_name character; name of assay to use for normalization
#' @param by_method character; normalization method to use, one of "median_center", "vsn", or "quantile"
#'
#' @returns SummarizedExperiment object with normalized assay added and set as primary assay
#' @export
#'
#' @examples NULL
normalize_data <- function(se_object, assay_name, by_method = c("median_center", "vsn", "quantile")) {
  # Simple, safe normalization helper. By default provides median-centering across samples.
  by_method <- match.arg(by_method)
  mat <- return_matrix_from_se(se_object = se_object, assay_name = assay_name)

  if (!is.numeric(mat))
    stop("Assay matrix must be numeric")

  if (by_method == "median_center") {
    # subtract column medians, add global median to preserve overall scale
    col_medians <- apply(mat, 2, stats::median, na.rm = TRUE)
    global_median <- median(col_medians, na.rm = TRUE)
    mat_norm <- sweep(mat, 2, col_medians, FUN = "-")
    mat_norm <- mat_norm + global_median
  } else if (by_method == "vsn") {
    if (!requireNamespace("vsn", quietly = TRUE))
      stop("Package 'vsn' required for vsn normalization")
    vsn_fit <- vsn::vsnMatrix(2^mat)
    mat_norm <- vsn::predict(vsn_fit, 2^mat)
  } else if (by_method == "quantile") {
    mat_norm <- limma::normalizeBetweenArrays(mat, method = "quantile")
  }

  # write back normalized assay and return se
  se_out <- se_object
  se_out <- add_assay(se_out, new_assay = as.matrix(mat_norm), new_assay_name = paste0(assay_name, "_norm"))
  se_out <- set_primary_assay(se_out, primary_assay_name = paste0(assay_name, "_norm"))
  return(se_out)
}

#' correct_mse_batch_effect
#'
#' @param se_object 
#' @param batch_col 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#'test_se<-readRDS("./projects/20260629_RRSSpatialProteomics_samples_0626/rds/ggu_se_id_3_preproc.rds")
#'test_correct<-correct_mse_batch_effect(se_object = test_se,batch_col = "patient",use_assay = "imputed")
correct_mse_batch_effect<-function(se_object,batch_col,use_assay="imputation"){
  dm <-
    return_matrix_from_se(se_object = se_object,assay_name = use_assay)
  se_anno<-as.data.frame(SummarizedExperiment::colData(se_object))
  
  dm_correct <-sva::ComBat(dat = dm, batch = se_anno[[batch_col]], par.prior = TRUE, prior.plots = FALSE)
  
  se_correct<-se_object
  se_correct <-
    add_assay(se_object = se_object,
              new_assay = dm_correct[rownames(se_object),colnames(se_object)],
              new_assay_name = "intensity_correct")
  se_correct<-set_primary_assay(se_correct,primary_assay_name = "intensity_correct")
  return(se_correct)
}


##### Data Quality assessment #####

#' proteom_cohort_summary
#'
#' @param se_object ; SummarizedExperiment Object obtained from import msdata
#'
#' @returns DataFrame Object which summarizes Mean, Median, SD and Percentage Misinngness in the assay Object 
#' @export
#'
#' @examples NULL
#' proteom_cohort_summary(se_object = ms_se)%>%ggtexttable()
proteom_cohort_summary <- function(se_object) {
  assertthat::assert_that(class(se_object)=="SummarizedExperiment")
  df_tmp <- data.frame(
    round(mean(colSums(!is.na(
      SummarizedExperiment::assay(se_object)
    ))), digits = 3),
    round(median(colSums(!is.na(
      SummarizedExperiment::assay(se_object)
    ))), digits = 3),
    round(sd(colSums(!is.na(
      SummarizedExperiment::assay(se_object)
    ))), digits = 3),
    round(sum(is.na(
      SummarizedExperiment::assay(se_object)
      )) / (nrow(SummarizedExperiment::assay(
      se_object
    )) * ncol((
      SummarizedExperiment::assay(se_object)
    ))) * 100, digits = 3)) %>% t()
  df_names <- c("Protein Identification (Mean)",
                "Protein Identification (Median)",
                "Protein Identification (SD)",
                "Missingness [%]")
  rownames(df_tmp) <- df_names
  colnames(df_tmp) <- paste0("N = ", ncol(se_object), "")
  df_tmp
}

#' calculate_cv
#' @author Dennis Friedel
#' @param x matrix;
#' @param na.rm boolean; remove missing values if any
#' @description
#' Behold. If you use this function will calculate the coeficient of variation sample wise or group wise.
#'
#' @description
#' Calculate coefficient of variation within dataset. Dispersion of data aroudn mean.
#'
#' @return A dataframe with columns sample_id, value (which is the coefficient of variation)
#' and group (if you like colors)
#' @export
#'
#' @examples NULL
calculate_cv <- function(x, na.rm = TRUE) {
  sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
}

#' count_idenitfied_features
#' @author Dennis Friedel
#' @param se_object SummarizedExperimentObject
#'
#' @return dataframe; sample_id and value which is the number of !NA per sample/feature
#' @export
#'
#' @examples NULL
#' count_idenitfied_features(se_object = pg_se,assay_name = "intensity",group = "group",column_wise = FALSE,in_percent = T)
#' count_idenitfied_features(se_object = pg_se,assay_name = "intensity",group = "group",column_wise = TRUE)
#' count_missing_features(se_object = pg_se,assay_name = "intensity",group = "group",column_wise = TRUE)
count_idenitfied_features = function(se_object,
                                     assay_name,
                                     group = NULL,
                                     column_wise = TRUE,
                                     in_percent = FALSE) {
  assay_eval <-
    return_matrix_from_se(se_object = se_object, assay_name = assay_name)
  
  assertthat::assert_that(is.numeric(assay_eval))
  if (column_wise) {
    na_df <- colSums(!is.na(assay_eval))
    n_size <- nrow(assay_eval)
  } else{
    na_df <- rowSums(!is.na(assay_eval))
    n_size <- ncol(assay_eval)
  }
  na_df <- na_df %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_id")
  colnames(na_df)[2]<-"value"
  if (in_percent) {
    na_df$value <- (na_df$value / n_size) * 100
  }
  
  if (!is.null(group) & column_wise) {
    assertthat::assert_that("sample_id" %in% colnames(SummarizedExperiment::colData(se_object)))
    message("add group")
    na_df$group <-
      plyr::mapvalues(na_df[["sample_id"]], se_object[["sample_id"]], se_object[[group]]) %>%
      as.factor(.)
    na_df <- dplyr::arrange(na_df, group)
    na_df$sample_id <-
      factor(na_df$sample_id, levels = unique(na_df$sample_id))
  }
  na_df<-na_df[order(na_df$group,na_df$value),]
  na_df$sample_id<-factor(na_df$sample_id,levels = na_df$sample_id)
  return(na_df)
}

#' get_intensity_long
#' @author Dennis Friedel
#' @param se_object SummarizedExperiment;
#' @param assay_name character; name of the assay
#' @param group character; name of the column in coldata containing group information
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#' se_object = pg_se
#' assay_name = "intensity"
#' group = "group"
#' long_df_test<-get_intensity_long(se_object = pg_se,assay_name = "intensity",group = "group")
get_intensity_long <-
  function(se_object, assay_name, group = NULL) {
    mat <-
      return_matrix_from_se(se_object = se_object, assay_name = assay_name)
    assay_long <- reshape2::melt(mat)
    colnames(assay_long) <- c("feature","sample_id", "value")
    if (!is.null(group)) {
      assertthat::assert_that("sample_id" %in% colnames(SummarizedExperiment::colData(se_object)))
      assay_long$group <-
        plyr::mapvalues(assay_long[["sample_id"]], se_object[["sample_id"]], se_object[[group]]) %>%
        as.factor(.)
      
      assay_long <- dplyr::arrange(assay_long, group)
      assay_long$sample_id <-
        factor(assay_long$sample_id, levels = unique(assay_long$sample_id))
    }
    assay_long
  }

#' get_cv_group
#'
#' @param se_obj Summarized experiment object
#' @param group character; Indicating the group for cv calculcation
#'
#' @return long dataframe shoing the CV values of proteins within  group 
#' @export
#' @description
#' Get percentage of dispersion of data aroudn mean in each group.
#'
#' @examples NULL
## se_object = pg_se
## assay_name = "intensity"
## group = "group"
## test_cv<-get_coefficient_variation(se_object = mse,assay_name = "intensity",group = "group")
get_coefficient_variation <-
  function(se_object,
           assay_name,
           group = "group",
           islog2 = TRUE) {
    assertthat::assert_that(group %in% colnames(SummarizedExperiment::colData(se_object)))
    mat <-
      return_matrix_from_se(se_object = se_object, assay_name = assay_name)
    linear_mat<-2^mat
    if (islog2) {
      linear_mat<-2^mat
    }else{
      linear_mat<-mat
    }
    ## Calculate CV per fraction
    groups<-unique(se_object$group)
    assay_long_cv<-purrr::map(groups,function(layers_tmp){
        tmpsamples<-se_object$sample_id[se_object$group==layers_tmp]
        cv_values<-apply(linear_mat[,tmpsamples],1,calculate_cv)%>%as.data.frame()
        cv_values$group<-layers_tmp
        cv_values$gene<-rownames(cv_values)
        colnames(cv_values)[1]<-"value"
        rownames(cv_values)<-NULL
        cv_values
    })%>%do.call(rbind,.)
    assay_long_cv
}

#' summarize_qc_metrics
#'
#' @param se SummarizedExperimentObject
#' @param number_proteins number rows of summarized experiment
#' @param percentage_na percentage missing values summarized experiment
#' @param mean_proteins_per_sample mean of proteins identified per samples
#' @param sd_variance SD of identified proteins per sample
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#' summarize_qc_metrics(se = se_object)
summarize_qc_metrics = function(se,
                                number_proteins = nrow(se),
                                percentage_na = get_percentage_na_se(se),
                                mean_proteins_per_sample = mean(get_features_per_sample(se)),
                                sd_variance = sd(colSums(!is.na(SummarizedExperiment::assay(se)))),
                                get_cv_condi = get_cv_group(se_obj = se))
{
  QC_table = data.frame(
    "number_proteins" = number_proteins,
    "percentage_na" = percentage_na,
    "Mean Samplewise features datastet" = round(mean_proteins_per_sample, digits = 2),
    "SD Samplewise features" = round(sd_variance, digits = 2)
  ) %>% t()
  return(rbind(QC_table, get_cv_condi))
}

##### Data Quality Visualization #####

#' plot_ggbarplot
#' @author Dennis Friedel
#'
#' @param df_long 
#' @param by_group 
#' @param return_data 
#' @param statbar 
#'
#' @examples 
#' df_long<-count_idenitfied_features(se_object = pg_se,assay_name = "intensity",group = "group",column_wise = T,in_percent = T)
plot_ggbarplot <-
  function(df_long,
           by_group = FALSE,
           return_data = FALSE,
           statbar = c("identity"),
           feature_name="value",
           threshold_value=5000
  ) {
    assertthat::assert_that("sample_id" %in% colnames(df_long))
    assertthat::assert_that("value" %in% colnames(df_long))
    assertthat::assert_that(is.logical(by_group),
                            is.logical(return_data))
    assertthat::assert_that(is.numeric(threshold_value))
    #identical(levels(df_long$sample_id),df_long$sample_id)
    if (by_group) {
      ggbp <-
        ggplot2::ggplot(df_long, ggplot2::aes(
          x = sample_id,
          y = value, 
          fill = group))
    } else{
      ggbp <-
        ggplot2::ggplot(df_long, ggplot2::aes(x = sample_id,
                                              y = value))
    }
    
    res <- ggbp +
      ggplot2::geom_bar(stat = statbar) +
      ggplot2::xlab(label = feature_name) +
      ggplot2::geom_hline(yintercept = threshold_value)+
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", colour = "grey50"),
        axis.text.x = ggplot2::element_text(
          angle = 80,
          vjust = 1,
          hjust = 1,
        ))
    res
    if (return_data) {
      return(res$data)
    } else{
      return(res)
    }
  }

#' plot_horizontal_ggboxplot
#' @author Dennis Friedel
#' @param df_long dataframe; object
#' @param group character; name of the column in coldata containing group information
#' @description
#' Generic function to plot horizontal boxplots for coefficient of variation,
#' intensity distribution
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#' # Intensity distribution
#' # Show coefficient of variation
#'plot_horizontal_ggboxplot(df_long = test_cv,by_group = TRUE,return_data = FALSE,feature_name = "coefficient of variation")
plot_horizontal_ggboxplot <-
  function(df_long,
           by_group = FALSE,
           return_data = FALSE,
           feature_name="value") {
    assertthat::assert_that("sample_id" %in% colnames(df_long))
    assertthat::assert_that("value" %in% colnames(df_long))
    assertthat::assert_that(is.logical(by_group),
                            is.logical(return_data))
    if (by_group) {
      data_median <- dplyr::summarise(group_by(df_long, group), MD = round(median(value),digits = 3))
      data_median$offset<-max(df_long$value)
      ggbp <-
        ggplot2::ggplot(df_long, ggplot2::aes(x = value,y = group, fill = group)) +
        ggplot2::geom_text(
          data = data_median,
          aes(x = offset,y = group,label = MD),
          position = position_dodge(width = 0.8),
          size = 3,
          vjust = 1
        )
    } else{
      ggbp <-
        ggplot2::ggplot(df_long, ggplot2::aes(x = value, y = sample_id, fill = group))
    }
    res <- ggbp +
      ggplot2::geom_boxplot(color = "black", outlier.color = "black") +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", colour = "grey50"),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          vjust = 1,
          hjust = 1
        ))+
      ggplot2::xlab(label = feature_name)
    
    if (return_data) {
      return(res$data)
    } else{
      return(res)
    }
  }

#' plot_feature_overlap
#'
#' @param se_object SummarizedExperiment; input data object.
#' @param assay_name Character; name of the assay to use.
#'
#' @return NULL
#' @export
#'
#' @examples NULL
plot_feature_overlap <- function(se_object,
                                 assay_name = NULL) {
  ## Select assay from assay_list
  mat <-
    return_matrix_from_se(se_object = se_object, assay_name = assay_name)
  assertthat::assert_that(is.numeric(mat), 
                          anyNA(mat))
  ## Count measured values across all samples and filter features with to many NA
  mat_long <- reshape2::melt(mat)
  colnames(mat_long) <- c("feature", "sample_id", "value")
  
  identified_feature <- mat_long %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise("value" = sum(!is.na(value))) %>%
    as.data.frame()
  identified_feature[identified_feature$value == "0", ]
  overlap_df <- as.data.frame(table(identified_feature$value))
  colnames(overlap_df) <- c("n_samples", "count")
  
  ggplot2::ggplot(overlap_df,
                  ggplot2::aes(x = n_samples, y = count, fill = n_samples)) +
    ggplot2::geom_bar(color = "black", stat = "identity") +
    ggplot2::scale_fill_grey() +
    ggplot2::labs(title = "Protein identifications overlap", x = "Identified in number of samples") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 45,
      vjust = 0.5,
      hjust = 1
    ))
}

#' plot_density_missing_features
#'
#' @param se_object SummarizedExperiment; input data object.
#' @param assay_name Character; name of the assay to use.
#'
#' @return NULL
#' @export
#'
#' @examples NULL
plot_density_missing_features <- function(se_object,
                                          assay_name = NULL) {
  ## Select assay from assay_list
  mat <-
    return_matrix_from_se(se_object = se_object, assay_name = assay_name)
  assertthat::assert_that(is.numeric(mat), anyNA(mat))
  
  ## Count measured values across all samples and filter features with to many NA
  mat_long <- reshape2::melt(mat)
  colnames(mat_long) <- c("feature", "sample_id", "value")
  
  identified_feature <- mat_long %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise("value" = sum(!is.na(value)) >= ncol(mat)) %>%
    as.data.frame()
  
  mat_long$group <-
    ifelse(
      mat_long$feature %in% identified_feature$feature[identified_feature$value],
      "no missing values",
      "contains missing values"
    )
  ggplot2::ggplot(mat_long, ggplot2::aes(x = value, color = group)) +
    ggplot2::geom_density() + ggplot2::scale_color_manual(values = c("gold","purple"))+
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = "grey50")
    )+
    ggplot2::labs(title = "Intensity distribution of features w/wo missing values", x = "Identified in number of samples")
}

#' missval_heatmap
#' @author Dennis Friedel
#' @param se SummarizedExperimentObject
#' @param groups_to_show character; vector with names of the columns of interest
#' @param hm_title character; title of the Heatmap
#' @description
#' Plot missing values in samples an their biological group
#'
#' @return NULL
#' @export
#'
#' @examples NULL
# show_biological_groups<-c("group","codel1p_19q","gewebe_von_primartumor_oder_rezidiv")
# show_technial_groups = c("batch_protein")
plot_missval_heatmap <-
  function (se_object,
            assay_name,
            show_biological_groups = c("group"),
            show_technial_groups = NULL
  ) {
    
    mat <-
      return_matrix_from_se(se_object = se_object, assay_name = assay_name)
    ## Check input
    assertthat::assert_that(is.numeric(mat), 
                            anyNA(mat),
                            all(show_biological_groups%in%colnames(SummarizedExperiment::colData(se_object))))
    
    missval <- ifelse(is.na(mat), 0, 1)
    ### Show only rows which contain features with missing values 
    missval<-missval[!(rowSums(missval)==ncol(missval)),]
    
    dfanno <-
      as.data.frame(SummarizedExperiment::colData(se_object)[, show_biological_groups])
    colnames(dfanno)<-show_biological_groups
    ## Make head annotation for missing value heatmap
    ht_top_anno <-
      ComplexHeatmap::HeatmapAnnotation(df = dfanno,
                                        col = ht_cat_color(se_object = se_object, groups_to_show = show_biological_groups))
    
    if(!is.null(show_technial_groups)){
      dfanno <-
        as.data.frame(SummarizedExperiment::colData(se_object)[, show_technial_groups])
      colnames(dfanno)<-show_technial_groups
      
      ht_bottom_anno <-
        ComplexHeatmap::HeatmapAnnotation(df = dfanno,
                                          col = ht_cat_color(se_object = se_object, groups_to_show = show_technial_groups))
      missing_heatmap <-
        ComplexHeatmap::Heatmap(
          missval,
          col = c("white", "black"),
          top_annotation = ht_top_anno,
          bottom_annotation = ht_bottom_anno,
          column_names_side = "bottom",column_names_rot = 45,
          column_title = paste0("Missing values pattern (n(Features)=",nrow(missval),")"),
          show_row_names = FALSE,
          show_column_names = F,
          name = "Missing values pattern",
          heatmap_legend_param = list(
            at = c(0, 1),
            labels = c("Missing value", "Valid value")
          )
        )
      missing_heatmap
      missing_heatmap<-complexheatmap_to_ggplot(missing_heatmap)
    }else{
      missing_heatmap <-
        ComplexHeatmap::Heatmap(
          missval,
          col = c("forestgreen", "purple4"),
          top_annotation = ht_top_anno,
          column_names_side = "bottom",column_names_rot = 45,
          column_title = paste0("Missing values pattern (n(Features)=",nrow(missval),")"),
          show_row_names = FALSE,
          show_column_names = F,
          name = "Missing values pattern",
          heatmap_legend_param = list(
            at = c(0, 1),
            labels = c("Missing value", "Valid value")
          )
        )
      missing_heatmap<-complexheatmap_to_ggplot(missing_heatmap)
    }
    return(missing_heatmap)
  }

#' plot_pca_overview_single_group
#' @author Dennis Friedel
#' @param se 
#' @param group_of_interest 
#' @param group_shows 
#' @param PC1 
#' @param PC2 
#' @param n_features 
#' @details
#' the main advantage of this plot is that it provides a convient overview of variation across the samples of on category. 
#' the aim of this function is to provide a fast overview of biological/technical effects present within a likely homogenous group.
#' @return list of ggarangment showing pca of single categories from group_of_interest
#' 
plot_pca_overview_single_group <- function(se,
                                           group_of_interest = "group",
                                           group_shows = technical_groups,
                                           PC1 = 1,
                                           PC2 = 2,
                                           n_features = 100) {
  single_group = unique(se[[group_of_interest]])
  subset_pca_plots = purrr::map(single_group, function(cat) {
    message(cat)
    se_subset = se[, se[[group_of_interest]] == cat]
    ## Plot Summary of subset 
    cat_pca_plots<-plot_pca_summary(
      se = se_subset,
      show_groups = group_shows,
      PC1 = 1,
      PC2 = 2,
      n_features = n_features
    )%>% ggpubr::annotate_figure(., top = paste0("PCA ",n_features," of ",cat," from ",group_of_interest))
  })
  subset_pca_plots
}

#' plot_normalized_dist
#'
#' @param se_object 
#' @param assay_name 
#' @param group 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
plot_normalized_dist<-function(se_object, assay_name, group = NULL){
  long_df_int<-get_intensity_long(se_object = se_object,assay_name = assay_name,group = group)
  int_boxplot<-plot_horizontal_ggboxplot(df_long = long_df_int,by_group = FALSE,return_data = FALSE,feature_name=assay_name)
  int_boxplot
}


#' plot_normalization_comparison
#'
#' @param se_obj SummarizedExperimentObject
#' @param show_groups character; groups to show in normalization graph
#' @description
#' compare distribution between   if there are shifts in samples means in context of technical groups
#'
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#' plot_normalization_comparison(se_obj = pg_se_norm,assay_name1 ="intensity" ,assay_name2 = "intensity_medianCentering",group = "group")
plot_normalization_comparison = function(se_obj,assay_name1,assay_name2,group){
  assay_names<-c(assay_name1,assay_name2)
  purrr::map(assay_names,function(mat){
    long_df_int<-get_intensity_long(se_object = se_obj,assay_name = mat,group = group)
    int_boxplot<-plot_horizontal_ggboxplot(df_long = long_df_int,by_group = FALSE,return_data = FALSE,feature_name=mat)
    int_boxplot
  })%>%ggpubr::ggarrange(plotlist = .,ncol = 1,nrow = 2)%>%
    ggpubr::annotate_figure(., top = paste0("Distribution ",assay_names[1]," vs. ",assay_names[2]))
}

#' plot_coefficent_variation_comparison
#'
#' @param se_obj SummarizedExperimentObject
#' @param show_groups character; groups to show in normalization graph
#' @description
#' compare distribution between   if there are shifts in samples means in context of technical groups
#'
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#' plot_normalization_comparison(se_obj = pg_se_norm,assay_name1 ="intensity" ,assay_name2 = "intensity_medianCentering",group = "group")
plot_coefficent_variation <- function(se_obj,assay_name="intensity",group="group",use_colors){
  long_var<-get_coefficient_variation(se_object = se_object,assay_name = assay_name,group = group)
  ggbp <-
    ggplot2::ggplot(long_var, ggplot2::aes(x = value,y = group, fill = group)) +
    ggplot2::geom_violin(trim = F,alpha=.7) +
    ggplot2::geom_boxplot(width=0.1,outlier.shape = NA,fill="white") +
    ggplot2::labs(y = "Layer", x = "Coefficient of variation") +
    ggplot2::coord_flip()+
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text =  element_text(face = "bold"),
      axis.text =  element_text(angle = 45, hjust = 1)
    )
  if(!is.null(use_colors)){
    ggbp+ggplot2::scale_fill_manual(values=use_colors)
  }
  ggbp
}

#### Proteolab report function ####

#' validate_prc_object
#'
#' @param se_object 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
validate_prc_object<-function(se_object){
  if(nrow(se_object)<=1){
    return(NULL)
  }
  if(ncol(se_object)<=1){
    return(NULL)
  }
  return("Passed")
}


#' readme_preproc_param
#'
#' @param pjname Character; Name of projects 
#' @param pjdes Character; Text Description of projects 
#' @param input_type Character; Either DIANN or MaxQuant
#' @param input_path Character; Folder containing DIA-NN/MaxQuant input
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
readme_preproc_param<-function(
    id,
    se_name,
    se_object,
    thr_sample,
    thr_feature,
    group_sel,
    filter_fractioned,
  normalization_method,
    imputation_method,
    batch_correction,
    batch_select
){
  preproc_config<-data.frame(
    "ID"=id,
    "Input_type"=se_name,
    "N Features"=nrow(se_object),
    "N Samples"=ncol(se_object),
    "Created"=as.character(Sys.time()),
    "Threshold Samples" = thr_sample,
    "Threshold Features" = thr_feature,
    "Fractionized Filtering" = filter_fractioned,
    "Fractionized Filtering Group" = group_sel,
    "Normalization Method" = normalization_method,
    "Imputation Method" = imputation_method,
    "Correct Batch Effect"=batch_correction,
    "Group Batch Effect"=batch_select
  )
  return(preproc_config)
}
