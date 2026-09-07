##### Utils modules #####

# Required libraries
#library("SummarizedExperiment")
#library("Rtsne")
#library("umap")
# library("ggpubr")
# library("ggrepel")

## For test purposes ##
# nrows <- 500
# ncols <- 10
# counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
# colData <- DataFrame(group = rep(c("ChIP", "Input"), 5),
#                      row.names = LETTERS[1:10])
# se0 <- SummarizedExperiment(assays = SimpleList(counts = counts),
#                             colData = colData)

##### Generic functions #####


#' set_primary_assay
#' @author Dennis Friedel
#' @param se_object SummarizedExperiment
#' @param new_assay matrix with same rownames and colnames as se_object
#' @param new_assay_name character; new name of the assay 
#'
#' @return NULL
#' @export
#'
#' @examples NULL
set_primary_assay <- function(se_object, primary_assay_name) {
  assertthat::assert_that(
    isClass(se_object, SummarizedExperiment),
    is.character(primary_assay_name)
  )
  new_order<-c(grep(primary_assay_name,names(SummarizedExperiment::assays(se_object))),grep(primary_assay_name,names(SummarizedExperiment::assays(se_object)),invert = T))
  SummarizedExperiment::assays(se_object)<-SummarizedExperiment::assays(se_object)[new_order]
  se_object
}


#' add_assay
#' @author Dennis Friedel
#' @param se_object SummarizedExperiment
#' @param new_assay matrix with same rownames and colnames as se_object
#' @param new_assay_name character; new name of the assay 
#'
#' @return NULL
#' @export
#'
#' @examples NULL
add_assay <- function(se_object, new_assay, new_assay_name) {
  assertthat::assert_that(
    isClass(se_object, "SummarizedExperiment"),
    is.matrix(new_assay),
    identical(colnames(new_assay), colnames(se_object)),
    all(rownames(new_assay) %in% rownames(se_object))
  )
  new_assay <- new_assay[rownames(se_object),]
  SummarizedExperiment::assays(se_object)[[new_assay_name]] <- new_assay
  se_object
}


#' convert_to_df_long
#' @author Dennis friedel
#' @param dfdata Dataframe/Summarized Experiment
#' @param sample_id character; name of the column containing the sample id 
#' @param value character; name of the column containing the value
#' @param group character; name of the column containing biological/technical group infromation by default NULL
#'
#' @return long data frame which can be used for plot_ggbarplot()
#' @export
#'
#' @examples NULL
convert_to_df_long <-
  function(dfdata,
           sample_id = "sample_id",
           value,
           group = NULL) {
    if (methods::is(dfdata, "SummarizedExperiment")) {
      df <- SummarizedExperiment::colData(dfdata) %>%
        as.data.frame()
    } else if (is.data.frame(dfdata)) {
      df <- dfdata
    } else {
      stop("dfdata must be a SummarizedExperiment or data.frame")
    }

    if (!is.null(group)) {
      df_long <- df[, c(sample_id, value, group), drop = FALSE]
      colnames(df_long)[colnames(df_long) == value] <- "value"
      colnames(df_long)[colnames(df_long) == group] <- "group"
    } else {
      df_long <- df[, c(sample_id, value), drop = FALSE]
      colnames(df_long)[colnames(df_long) == value] <- "value"
    }

    return(df_long)
  }

#' return_matrix_from_se
#' @author Dennis Friedel
#' @param se_object SummarizedExperiment
#' @param assay_name character
#'
#' @return matrix
#' @export
#'
#' @examples NULL
return_matrix_from_se <- function(se_object, assay_name = NULL) {
  assertthat::assert_that(
    is.character(assay_name) | is.null(assay_name),
    methods::is(se_object, "SummarizedExperiment")
  )
  if (!is.null(assay_name)) {
    se_assay <- names(SummarizedExperiment::assays(se_object))
    if (!assay_name %in% se_assay)
      stop(paste0(
        "Assay name not listed in SummarizedExperiment. Available are :",
        paste0(se_assay, collapse = "\t")
      ))
    mat <-
      as.matrix(SummarizedExperiment::assay(se_object, assay_name))
  } else{
    mat <- as.matrix(SummarizedExperiment::assay(se_object, 1))
  }
  return(mat)
}

#' ht_cat_color
#'
#' @param se_object SummarizedExperimentObject; input data object.
#' @param groups_to_show character; vector with names of the columns of interest
#' @param is_dataframe Boolean; if TRUE se_object is a dataframe, if FALSE se_object is a SummarizedExperimentObject
#' @param ... 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
ht_cat_color <- function(se_object, groups_to_show,is_dataframe=F, ...) {
  if(is_dataframe){
    se_data <-se_object
  }else{
    se_data <- as.data.frame(SummarizedExperiment::colData(se_object))
  }

  dress_code <-
    c(ggsci::pal_simpsons, ggsci::pal_futurama, ggsci::pal_aaas, ggsci::pal_cosmic, ggsci::pal_gsea)
  dataframe_to_col <- se_data[, groups_to_show]
  dataframe_to_col[is.na(dataframe_to_col)] <- "No information"
  if (is(class(dataframe_to_col), "character"))
    dataframe_to_col <- as.data.frame(dataframe_to_col)
  
  cat_colors <- c()
  for (i in 1:length(groups_to_show)) {
    colors <-
      rip::cat_color_code(ggsci_colorpalette = dress_code[[i]](),
                          categorial_vector = dataframe_to_col[, i])
    cat_colors[[i]] <- colors
  }
  names(cat_colors) <- groups_to_show
  cat_colors
}


#' cat_color_code
#'
#' @param ggsci_colorpalette 
#' @param categorial_vector 
#' @param use_ggsci 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
cat_color_code<-function(ggsci_colorpalette,
                         categorial_vector,
                         use_ggsci = T) {
  cat <- unique(categorial_vector)
  if (use_ggsci) {
    colors <- ggsci_colorpalette(length(cat))
  } else{
    colors <- ggsci_colorpalette[1:length(cat)]
  }
  names(colors) <- (cat)
  colors
}

#' se_select_topvaribale_features
#'
#' @param se Summarized Experiment object
#' @param n integer number of top variable features
#'
#' @return Returns filtered SummarizedExperiment object
#' @export
#'
#' @examples NULL
#'
#' se_flt = se_select_topvaribale_features(se0,n=100)
#' se_flt
se_select_topvaribale_features <- function(se, n = 1000) {
  var <- apply(SummarizedExperiment::assay(se), 1,FUN =  var)
  se_object_filt <- se[order(var, decreasing = TRUE)[seq_len(n)]]
  return(se_object_filt)
}

#' addReduction
#'
#' @param se SummarizedExperiment; input object.
#'
#' @return returns Summarized experiment object with list reduction
#' @export
#'
#' @examples NULL
#' addReduction(se = se)
addReduction = function(se) {
  assertthat::assert_that(
    isClass(se, "SummarizedExperiment")
  )
  # Check if reduction is available in metadata of Summarized experiment
  if (is.null(S4Vectors::metadata(se)[["reduction"]])) {
    S4Vectors::metadata(se)[["reduction"]] = list()
    se
  }
  else{
    se
  }
}

#' get_reduction_se
#'
#' @param se Summarized experiment object with reduction in metdadata and result from the reduction method
#' @param red.method chracter indicating the result. Either PCA,TSNE or UMAP
#' @description
#' This function returns the dimensional reduction result of the selected method PCA,TSNE,UMAP
#' stored in metadata of the SummarizedExperiment object
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#' get_reduction_se(se0,red.method="PCA")
get_reduction_se = function(se, red.method = c("PCA", "TSNE", "UMAP")) {
  if (is.null(S4Vectors::metadata(se)$reduction[[red.method]])) {
    message("please execute run",
            red.method,
            "_se to get the desired output")
    return()
  }
  res = S4Vectors::metadata(se)[["reduction"]][[red.method]]
  res
}


### Load required packages and Functions
#library("paran")
#library("cluster")
#library("factoextra")

#####----------------------------------------------------------------------------------#
# Inspect contribution of principal components to total variance in the data
#' se_fviz_plot_overview
#' @author Dennis Friedel
#' @param se_object SummarizedExperiment; input data object.
#' @param nfeatures Integer; number of top variable features to retain.
#'
#' @return NULL
#' @export
#'
#' @examples NULL
se_fviz_plot_overview <-
  function(se_object,
           nfeatures = c(500, 1000, 1500)){
    se_flt_objects <-
      purrr::map(nfeatures, function(x){
        se_select_topvaribale_features(se_object, n = x)
      })
    
    se_pca_objects <-
      purrr::map(se_flt_objects, function(x){
        runPCA_se(
          se = x,
          center = T,
          scale = T
        )
      })
    fviz_plots <-
      purrr::map(se_pca_objects, function(x) {
        factoextra::fviz_eig(get_reduction_se(x, red.method = "PCA"), ncp = 30) +
          ggplot2::ggtitle(label = paste0("Percentage variance over Dimensions n_features ", nrow(x)))
      })
    fviz_plots_arrangement <- ggpubr::ggarrange(plotlist = fviz_plots)
    fviz_plots_arrangement
  }

#####----------------------------------------------------------------------------------#
#' se_paran_plot_overview
#' @author Dennis Friedel
#' @param se_object SummarizedExperiment;
#' @param nfeatures Integer; number features
#' @param n_perm Integer; number permuations
#' @description
#' Inspect contribution of principal components to total variance in the data, uses paran to performs Horn's parallel analysis
#' for a principal component or common factor analysis, so as to adjust for finite sample bias in the retention of components.
#'
#' @return NULL
#' @export
#'
#' @examples NULL
# se_object = beta_genescore_pca
# nfeatures = 500
# n_perm<-100
# se_paran_plot_overview <-
#   function(se_object,
#            nfeatures = c(500, 1000, 1500),
#            n_perm = 100,
#            plot_result = FALSE,
#            result_path = "") {
#     assertthat::assert_that(
#       is.numeric(nfeatures),
#       is.numeric(n_perm),
#       is.logical(plot_result),
#       is.character(result_path)
#     )
#     # Select features
#     se_flt_objects <-
#       purrr:::map(nfeatures, function(x)
#         se_select_topvaribale_features(se_object, n = x))
#     paran_result_vector<-c()
#     for (i in 1:length(se_flt_objects)) {
#       out <-
#         rip::save_here(
#           output_dir = result_path,
#           object_name =  paste0("paran_n", nrow(se_flt_objects[[i]]), "_perm", n_perm, ".pdf")
#         )
#       out <- gsub("pdf", "png", out)
#       png(out, width = 1000, height = 1000)
#       paran_result <- print(
#         paran::paran(
#           SummarizedExperiment::assay(se_flt_objects[[i]]),
#           iterations = n_perm,
#           quietly = FALSE,
#           status = TRUE,
#           all = TRUE,
#           cfa = FALSE,
#           graph = plot_result,
#           color = TRUE,
#           col = c("black", "red", "blue"),
#           lty = c(1, 2, 3),
#           lwd = 1,
#           legend = TRUE,
#           width = 1000,
#           height = 1000,
#           grdevice = "png",
#           seed = 0,
#           mat = NA,
#           n = NA
#         )
#       )
#       dev.off()
#       paran_result_vector <-
#         c(paran_result_vector, paran_result$Retained)
#     }
#     names(paran_result_vector) <- paste0("nfeatures", nfeatures)
#     paran_result_vector
#   }



#=====================================#
##### Dimension reduction methods #####
#=====================================#

# PCA
#' runPCA_se
#' @author Dennis Friedel
#' @param se Summarized Experiment object
#' @param center boolean if values should be centered
#' @param scale boolean if values shoudl be scaled
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#'  se_pca = runPCA_se(se=se0,center = F, scale = F)
#'  str(get_reduction_se(se_pca,red.method = "PCA"))
#se = runPCA_se(se=se0,center = F, scale = F)
runPCA_se <- function(se, center = F, scale = F) {
  assertthat::assert_that(methods::is(se, "SummarizedExperiment"))
  se = addReduction(se = se)
  pca_data <- stats::prcomp(t(SummarizedExperiment::assay(se)), center = center, scale = scale)
  S4Vectors::metadata(se)$reduction[["PCA"]] =  pca_data
  se
}

#' runTSNE_se
#' @author Dennis Friedel
#' @param se SummarizedExperiment object
#' @param perpl numeric; Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details RTSNE for interpretation
#' @param iter 	integer; Number of iterations
#' @param use.seed integer; seed to use
#' @param theta_se 	numeric; Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
#' @param pca_se 	logical; Whether an initial PCA step should be performed (default: TRUE)
#' @param pca_scale_se 	logical; Should data be scaled before pca is applied? (default: FALSE)
#' @param pca_center_se 	logical; Should data be centered before pca is applied? (default: TRUE)
#' @param check_duplicates_se 	logical; Checks whether duplicates are present. It is best to make sure there are no duplicates present and set this option to FALSE, especially for large datasets (default: TRUE)
#' @param normalize_se logical; Should data be normalized internally prior to distance calculations with normalize_input? (default: FALSE)
#' @param initial_dims_se integer; the number of dimensions that should be retained in the initial PCA step (default: 50)
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#' tsne_test_se = runTSNE_se(se = se_pca, perpl = 2)
#' res = get_reduction_se(se = tsne_test_se, red.method = "TSNE")
runTSNE_se <- function(se,
                       perpl = 10,
                       iter = 100,
                       theta_se = 0,
                       pca_se = TRUE,
                       pca_scale_se = FALSE,
                       pca_center_se = FALSE,
                       check_duplicates_se = FALSE,
                       normalize_se = FALSE,
                       initial_dims_se = 50,
                       use.seed = 2905) {
  
  set.seed(use.seed)
  assertthat::assert_that(methods::is(se, "SummarizedExperiment"))
  se = addReduction(se = se)
  
  tsne_result <-
    Rtsne::Rtsne(
      X = t(SummarizedExperiment::assay(se)),
      perplexity = perpl,
      max_iter = iter,
      theta = theta_se,
      pca = pca_se,
      pca_scale = pca_scale_se,
      pca_center = pca_center_se,
      normalize = normalize_se,
      initial_dims = initial_dims_se,
      check_duplicates = check_duplicates_se
    )
  S4Vectors::metadata(se)$reduction[["TSNE"]] =  tsne_result
  se
}

#' runUMAP_se
#' @author Dennis Friedel
#' @param se SummarizedExperiment object
#' @param method_umap character; either naive or umpa learn
#' @param umap_metric character; either euclidean or cosine
#' @param neighbors_se integer; number of neigbours similar to TSNE
#' @param min_dist_se numeric;
#' @param initial_dims_se integer how many initial PCA components to use, only required when use_PCA = T
#' @param seed.use integer
#' @param use_PCA boolean; to use the PCA matrix as input or not
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#' testumap=runUMAP_se(se = se_pca,method_umap = "naive",use_PCA = F,umap_metric = "cosine",initial_dims_se = 5)
#' res = get_reduction_se(testumap,red.method = "UMAP")
#'
#' testumap=runUMAP_se(se = se_pca,method_umap = "naive",use_PCA = T,umap_metric = "cosine",initial_dims_se = 5)
#' res = get_reduction_se(testumap,red.method = "UMAP")
runUMAP_se <-
  function(se,
           method_umap = c("naive", "umap-learn"),
           use_PCA = F,
           umap_metric = c("euclidean", "cosine"),
           neighbors_se = 10,
           min_dist_se = 0.1,
           initial_dims_se = 50,
           seed.use = 2905) {
    
    set.seed(seed.use)
            assertthat::assert_that(methods::is(se, "SummarizedExperiment"))
    
    if (use_PCA) {
      if(is.null(get_reduction_se(se, red.method = "PCA"))){
        se<-runPCA_se(se = se,center = T,scale = T)
        d_matrix <-
          get_reduction_se(se, red.method = "PCA")$x[, 1:initial_dims_se]
      }else{
        d_matrix <-
          get_reduction_se(se, red.method = "PCA")$x[, 1:initial_dims_se]
      }
    } else{
      d_matrix <- t((assay(se)))
    }
    
    umap_res <-
      umap::umap(
        d = d_matrix,
        method = method_umap,
        metric = umap_metric,
        preserve.seed = T,
        n_neighbors = neighbors_se,
        min_dist = min_dist_se
      )
    S4Vectors::metadata(se)$reduction[["UMAP"]] =  umap_res
    se
  }


# require("Rtsne")
# require("umap")
#===========================================================#
##### Plotting functions for Dimension reduction results #####
#===========================================================#

#' ggplot_pca_se
#' @author Dennis Friedel
#' @param se SummarizedExperiment with PCA in metadata reduction
#' @param PCX integer; number of PC in pccromp$x which shoudl be shown in X Axis
#' @param PCY integer; number of PC in pccromp$x which shoudl be shown in Y Axis
#' @param col character; character or numeric vector by which the dots shoudl be colorized
#' @param size integer; size of dots
#' @param colorvalues character vector; colors to use
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#' PCX = 1
#' PCY = 2
#' size = 4
ggplot_pca_se <-
  function(se,
           PCX = 1,
           PCY = 2,
           col = "",
           size = 4,
           colorvalues = NULL) {
    methods::is(se, "SummarizedExperiment")
    
    group_values = se[[col]]
    pca_x <- as.data.frame(get_reduction_se(se, "PCA")$x)
    
    p <-
      ggplot2::ggplot(pca_x, ggplot2::aes(x = pca_x[, PCX], y = pca_x[, PCY], color = group_values)) +
      ggplot2::geom_point(size = size) +
      ggplot2::xlab(label = paste0(colnames(pca_x)[PCX])) +
      ggplot2::ylab(label = paste0(colnames(pca_x)[PCY])) +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", colour = "grey50"))
    
    if (is.numeric(group_values)) {
      p <- p +
        ggplot2::scale_colour_gradient2(
          high = "red",
          mid = "grey",
          low = "blue",
          midpoint = median(group_values,na.rm = T)
        )
    }
    if (!is.null(colorvalues)) {
      if (!is.factor(group_values))
        group_values <- factor(group_values)
      names(colorvalues) <- levels(group_values)
      p <- p + ggplot2::scale_colour_manual(name = group, values = colorvalues)
    }
    return(p)
  }

#TSNE
#' ggplot_tsne_se
#' @author Dennis Friedel
#' @param se SummarizedExperiment with PCA in metadata reduction
#' @param col character; character or numeric vector by which the dots shoudl be colorized
#' @param size integer; size of dots
#'
#' @return NULL
#' @export
#'
#' @examples NULL
ggplot_tsne_se <- function(se,
                           col = "",
                           size = 4,
                           colorvalues = NULL) {
  methods::is(se, "SummarizedExperiment")
  group_values = se[[col]]
  tsne.coords = as.data.frame(get_reduction_se(se, "TSNE")$Y)
  
  p <-
    ggplot2::ggplot(tsne.coords, ggplot2::aes(x = V1, y = V2, color = group_values)) +
    ggplot2::geom_point(size = size) +
    ggplot2::xlab(label = paste0("X")) +
    ggplot2::ylab(label = paste0("Y")) + 
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = "grey50"))
  
  if (is.numeric(group_values)) {
    p <- 
      p + 
      ggplot2::scale_colour_gradient2(
        high = "red",
        mid = "grey",
        low = "blue",
        midpoint = median(group_values,na.rm = T)
      )
    
  }
  if (!is.null(colorvalues)) {
    if (!is.factor(group_values))
      group_values <- factor(group_values)
    names(colorvalues) <- levels(group_values)
    p<-p + ggplot2::scale_colour_manual(name = group, values = colorvalues)
  }
  return(p)
}

#UMAP
#' ggplot_umap_se
#' @author Dennis Friedel
#' @param se SummarizedExperiment with PCA in metadata reduction
#' @param col vector; character or numeric vector by which the dots shoudl be colorized
#' @param size integer; size of dots
#'
#' @return NULL
#' @export
#'
#' @examples NULL
ggplot_umap_se <- function(se,
                           col = "",
                           size = 4,
                           colorvalues = NULL) {
  methods::is(se, "SummarizedExperiment")
  umap.coords = as.data.frame(get_reduction_se(se, "UMAP")$layout)
  group_values = se[[col]]
  p <-
    ggplot2::ggplot(umap.coords, ggplot2::aes(x = V1, y = V2, color = group_values)) + 
    ggplot2::geom_point(size =size) +
    ggplot2::xlab(label = paste0("UMAP_1")) +
    ggplot2::ylab(label = paste0("UMAP_2")) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = "grey50"))
  if (is.numeric(group_values)) {
    p <- 
      p + 
      ggplot2::scale_colour_gradient2(
        high = "red",
        mid = "grey",
        low = "blue",
        midpoint = median(group_values,na.rm = T)
      )
  }
  if (!is.null(colorvalues)) {
    if (!is.factor(group_values))
      group_values <- factor(group_values)
    names(colorvalues) <- levels(group_values)
    p<-p + ggplot2::scale_colour_manual(name = col, values = colorvalues)
  }
  return(p)
}

#' plotPCA_summary
#' @author Dennis Friedel
#' @param se SummarizedExperimentObject
#' @param show_groups character; vector with names of the columns of interest
#' @param PC1 integer; PC to show  in X
#' @param PC2 integer; PC to show  in Y
#' @param n_features integer; number of top variable proteins to use for PCA
#' @param dotsize integer; size of dots in ggplot
#' 
#' @return NULL
#' @export
#'
#' @examples NULL
plot_pca_summary = function(se,
                            show_groups,
                            PC1=1,
                            PC2=2,
                            n_features,
                            dotsize = 4) {
  pca_plot_groups <- purrr::map(show_groups, function(sg) {
    test <- se_select_topvaribale_features(se, n = n_features) %>%
      runPCA_se(., center = T, scale = T) %>%
      ggplot_pca_se(
        se = .,
        PCX = PC1,
        PCY = PC2,
        col = sg,
        size = dotsize,
        colorvalues = NULL
      ) + 
      ggtitle(label = sg)
  }) %>%
    ggpubr::ggarrange(plotlist = .) %>% 
    ggpubr::annotate_figure(., top = paste0("PCA top ", n_features, " features"))
}

#====================================#
##### Advanced functions #####
#====================================#

#' Identify variation using pca_anova
#'
#' @param variable 
#' @param mat 
#' @param data_df 
#'
#' @return NULL
#' @export
#'
#' @examples NULL
pca_anova <- function(variable, mat, data_df){
  pca_new <- stats::prcomp(mat)
  pca_res <- summary(pca_new)
  #data_df[[variable]][is.na(data_df[[variable]])]<-"unknown"
  factor_ids = factor(data_df[[variable]])
  anova = aov(pca_new$x ~ factor_ids)
  anova_summary = summary(anova)
  df_anova <- do.call("rbind", anova_summary) %>% na.omit(.) 
  rownames(df_anova) <- colnames(pca_new$x)
  df_anova$PoV <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  PoV_PCs <- data.frame('PoV_sum' = sum(df_anova[df_anova$`Pr(>F)` < 0.05,]$PoV), 
                        'PCs_sum' = toString(rownames(df_anova[df_anova$`Pr(>F)` < 0.05,])))
  rownames(PoV_PCs) = variable
  return(PoV_PCs)
}

#' plot_components
#'
#' @param se_object 
#' @param show_group 
#' @param components 
#' @param n_features 
#' @param dotsize 
#' @description
#' Plots several Principal components to visualize their contribution in context of a selected group
#' 
#' @returns NULL
#' @export
#'
#' @examples NULL
plot_components<-function(se_object,
                          show_group = "group",
                          components = 1:4,
                          n_features = NULL,
                          dotsize = 4){
  assertthat::assert_that(class(se_object)=="SummarizedExperiment")
  assertthat::assert_that(is.numeric(components))
  assertthat::assert_that(is.numeric(dotsize))
  assertthat::assert_that(length(components)>1)
  assertthat::assert_that(show_group%in%colnames(SummarizedExperiment::colData(se_object)))
  if(is.null(n_features)){
    n_features<-nrow(se_object)
  }
  if (ncol(se_object) < 3) {
    return(NULL)
  }
  pca_x <- se_select_topvaribale_features(se_object, n = n_features) %>%
    runPCA_se(., center = T, scale = T) %>%
    get_reduction_se(., red.method = "PCA") %>%
    .$x %>%
    .[, components]%>%
    as.data.frame()
  pca_x$sample<-rownames(pca_x)
  pca_x$group<-as.character(se_object[[show_group]])
  
  comps<-paste0("PC", components)
  ### Prepare legend 
  p <- ggplot2::ggplot(pca_x, ggplot2::aes_string(x = comps[1], y = comps[2], color = "group")) +
    ggplot2::geom_point()+
    ggplot2::theme(
      legend.key = ggplot2::element_rect(fill="white"),
      legend.text = ggplot2::element_text(size=ggplot2::rel(1.2)),
      legend.title = ggplot2::element_text(size=ggplot2::rel(1.2))
    )
  legend<-GGally::grab_legend(p)
  p <- GGally::ggpairs(
    pca_x,
    columns = paste0("PC", components),
    lower = list(continuous = GGally::wrap("points", size =
                                             dotsize)),
    diag = list(continuous = 'densityDiag'),
    upper = list(continuous = GGally::wrap("points", size =
                                             dotsize)),
    mapping = ggplot2::aes_string(color = "group"),
    title = show_group,legend=legend
  )
  p<-p+ggplot2::theme_bw()+ggplot2::theme(
    panel.grid.major  = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text  = ggplot2::element_blank()
  )
  return(p)
}

#==================================================#
###### Principle component gene set enrichment #####
#==================================================#

#' create.feature.set.matrix
#' @author Dennis Friedel
#' @param term_list list; list with pathway-name as name and and gene names as vectors
#' @description
#' Takes a list of genesets as input and creates a binary matrix object with Pathways/genesets as names and genes/ENSG as columns
#'
#' @return Binary matrix object with Pathway names as names and genes/ENSG as columns
#' @export
#'
#' @examples NULL
create.feature.set.matrix <- function(term_list) {
  #v1 using description from MOFA manual
  row_anno <- names(term_list)
  col_anno <- unique(unlist((term_list)))
  factors.matrix <-
    matrix(
      data = 0,
      ncol = length(col_anno),
      nrow = length(row_anno),
      dimnames = list(row_anno, col_anno)
    )
  factors.matrix <-
    purrr::map(names(term_list), function(term_description) {
      factors.matrix[term_description,] <-
        names(factors.matrix[term_description,]) %in% term_list[[term_description]]
      factors.matrix[term_description,]
    }) %>% do.call(rbind, .)
  rownames(factors.matrix) <- names(Hmsig)
  factors.matrix
}

#' Convert Complex Heatmap object into GGplot object
#'
#' @param heatmap S4Object; ComplexHeatmap saved in variable
#'
#' @returns It does what it says. You input a Complex Heatmap object and the 
#' function returns a ggplot object. You are welcome.
#' @export
#'
#' @examples NULL
complexheatmap_to_ggplot<-function(heatmap){
  ggpubr::as_ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(heatmap)))  
}



#====================================#
##### Under Development #####
#====================================#
#' msr.pca.enrichment
#'
#' @param object SummarizedExperimentObject with PCA result in metadata
#' @param feature.sets Matrix; binary feature-set membership matrix (features x sets).
#' @param factors integer vector
#' @param set.statistic character
#' @param statistical.test character; Which test to use "parametric" or "cor.adj.parametric"
#' @param sign character; which direction of enrichement "all", "positive", "negative"
#' @param min.size integer; min size of a geneset
#' @param nperm integer; number of permutations
#' @param p.adj.method character; Which method should be selected for p-value adjustement BH,FDR,BY
#' @param alpha numeric; threshold for
#' @param verbose Logical; if TRUE print progress messages.
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#' test = msr.pca.enrichment(
#'   object = MS_PCA,
#'   feature.sets = feature.sets,
#'   factors = 1:9,
#'   set.statistic = "mean.diff",
#'   statistical.test = "parametric",
#'   sign = "all",
#'   min.size = 10,
#'   p.adj.method = "BH",
#'   alpha = 0.5,
#'   verbose = TRUE)
#'
# msr.pca.enrichment <- function(object,
#                                feature.sets,
#                                factors = 1:9,
#                                set.statistic = c("mean.diff", "rank.sum"),
#                                statistical.test = c("parametric", "cor.adj.parametric"),
#                                sign = c("all", "positive", "negative"),
#                                min.size = 10,
#                                nperm = 1000,
#                                p.adj.method = "BH",
#                                alpha = 0.1,
#                                verbose = TRUE) {
#   library("PCGSE")
#   if (!is(object, "SummarizedExperiment"))
#     stop("'object' has to be an SummarizedExperiment")
#   if (!is(object@metadata$PCA, "prcomp"))
#     stop(
#       "please run pca_summarized_experiment on the SummarizedExperiment proir this analysis"
#     )
#   if (!(is(feature.sets, "matrix") & all(feature.sets %in%
#                                          c(0, 1))))
#     
#     sign <-
#       match.arg(sign, choices = c("all", "positive", "negative"))
#   set.statistic <-
#     match.arg(set.statistic, c("mean.diff", "rank.sum"))
#   statistical.test <-
#     match.arg(statistical.test, c("parametric", "cor.adj.parametric"))
#   
#   if (sign == "positive") {
#     pca.rotation[pca.rotation < 0] <- 0
#   }  else{
#     if (sign == "negative") {
#       pca.rotation[pca.rotation > 0] <- 0
#       pca.rotation <- abs(pca.rotation)
#     }
#   }
#   
#   idx <-
#     apply(pca.rotation, 2, function(x)
#       var(x, na.rm = TRUE)) == 0
#   if (sum(idx) >= 1) {
#     warning(sprintf(
#       "%d features were removed because they had no variance in the data.\n",
#       sum(idx)
#     ))
#     pca.rotation <- pca.rotation[!idx, ]
#   }
#   features <-
#     intersect(rownames(pca.rotation), colnames(feature.sets))
#   if (length(features) == 0)
#     stop("Feature names in feature.sets do not match feature names in model.")
#   if (verbose) {
#     message(
#       sprintf(
#         "Intersecting features names in the model and the gene set annotation results in a total of %d features.",
#         length(features)
#       )
#     )
#   }
#   pca.rotation <- object@metadata$PCA$rotation
#   feature.sets <- feature.sets[, features]
#   feature.sets <- feature.sets[rowSums(feature.sets) >= min.size,]
#   
#   if (verbose) {
#     message(
#       "\nRunning feature set Enrichment Analysis with the following options...\n",
#       sprintf("Number of feature sets: %d \n",
#               nrow(feature.sets)),
#       sprintf("Set statistic: %s \n",
#               set.statistic),
#       sprintf("Statistical test: %s \n",
#               statistical.test)
#     )
#     if (sign %in% c("positive", "negative"))
#       message(sprintf("Subsetting weights with %s sign",
#                       sign))
#     if (statistical.test == "permutation") {
#       message(sprintf("Number of permutations: %d", nperm))
#     }
#     message("\n")
#   }
#   
#   results <- PCGSE::pcgse(
#     data = t(assay(object)),
#     pc.indexes = factors,
#     gene.sets = feature.sets,
#     gene.set.statistic = set.statistic,
#     gene.set.test = statistical.test
#   )
#   
#   
#   
#   pathways <- rownames(feature.sets)
#   colnames(results$p.values) <-
#     colnames(results$statistics) <- factors
#   rownames(results$p.values) <-
#     rownames(results$statistics) <- pathways
#   
#   
#   if (!p.adj.method %in% p.adjust.methods)
#     stop("p.adj.method needs to be an element of p.adjust.methods")
#   adj.p.values <-
#     apply(results$p.values, 2, function(lfw)
#       p.adjust(lfw,
#                method = p.adj.method))
#   
#   if (sign %in% c("positive", "negative")) {
#     results$p.values[results$statistics < 0] <- 1
#     adj.p.values[results$statistics < 0] <- 1
#     results$statistics[results$statistics < 0] <- 0
#   }
#   if (sign %in% c("positive", "negative")) {
#     results$p.values[results$statistics < 0] <- 1
#     adj.p.values[results$statistics < 0] <- 1
#     results$statistics[results$statistics < 0] <- 0
#   }
#   sigPathways <-
#     lapply(factors, function(j)
#       rownames(adj.p.values)[adj.p.values[, j] <= alpha])
#   output <-
#     list(
#       feature.sets = feature.sets,
#       pval = results$p.values,
#       pval.adj = adj.p.values,
#       feature.statistics = results$feature.statistics,
#       set.statistics = results$statistics,
#       sigPathways = sigPathways
#     )
#   
#   return(output)
# }

##### Deprecarted ####
### PLot several plots to data
#' Title
#' #'
#' #' @param se_object SummarizedExperiment; input data object.
#' #' @param dimx
#' #' @param dimy
#' #' @param pca_center
#' #' @param pca_scale
#' #' @param show_categories
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' se_inspect_w_pca <-
#'   function(se_object,
#'            dimx = 1,
#'            dimy = 2,
#'            pca_center = T,
#'            pca_scale = T,
#'            show_categories = c("condition")) {
#'     se_object <-runPCA_se(se = se_object, center = pca_center, scale = pca_scale)
#'     PCA_plot <- map(show_categories, function(catg) {
#'       ggplot_pca(
#'         se_meta = get_reduction_se(se_object,red.method = "PCA"),
#'         color =  MS_DEP_dimred_topvar[[catg]],
#'         size = 4,
#'         PCX = dimx,
#'         PCY = dimy,
#'         colorvalues = cat_color_code(
#'           ggsci_colorpalette = pal_simpsons(),
#'           categorial_vector = factor(MS_DEP_dimred_topvar[[catg]])
#'         )
#'       ) +
#'         ggtitle(label = paste0(
#'           "Category ",
#'           catg,
#'           " (PCA Top ",
#'           nrow(se_object),
#'           " varibale proteins) "
#'         ))
#'     })
#'     PCA_result <- ggarrange(plotlist = list(PCA_plot))
#'   }
#'
#'
#'
#' ggplot_pca_enrichment <-
#'   function(se_object,
#'            GSEA_meta,
#'            PCX = 1,
#'            PCY = 2,
#'            show = c("Description", "Cluster")) {
#'     pca_y <- as.data.frame(se_object@metadata$PCA$rotation)
#'     GSEA_comp <- GSEA_meta
#'     GSEA_comp <- GSEA_comp[!duplicated(GSEA_comp$Description),]
#'     res_ggo_all <-
#'       map(1:nrow(GSEA_comp), function(x)
#'         data.frame(
#'           "Description" = GSEA_comp$Description[x],
#'           "Enriched_gene" = unlist(strsplit(GSEA_comp$core_enrichment[x], "/"))
#'         )) %>% do.call(rbind, .)
#'     if (show == "Cluster") {
#'       res_ggo_all$Cluster <-
#'         mapvalues(res_ggo_all$Description,
#'                   GSEA_comp$Description,
#'                   GSEA_comp$Cluster)
#'     }
#'
#'     #remove genes that occur in lower enriched results
#'     res_ggo_all <-
#'       res_ggo_all[!duplicated(res_ggo_all$Enriched_gene),]
#'     pca_y$genename <-
#'       ifelse(
#'         rownames(pca_y) %in% res_ggo_all$Enriched_gene,
#'         yes = rownames(pca_y),
#'         no = NA
#'       )
#'
#'     pca_y$class <-
#'       mapvalues(pca_y$genename, res_ggo_all$Enriched_gene, res_ggo_all[, show])
#'
#'     p <-
#'       ggplot(as.data.frame(pca_y),
#'              aes(x = pca_y[, PCX], y = pca_y[, PCY], label = genename)) +
#'       geom_point(aes(color = class), size = 3) + geom_text() +
#'       xlab(label = paste0(colnames(pca_y)[PCX])) +
#'       ylab(label = paste0(colnames(pca_y)[PCY])) +
#'       theme_classic()
#'     plot(p)
#'   }
#'
#' #Visulaization of ROATION
#' Title
#'
#' @param se_object SummarizedExperiment; input data object.
#' @param PCX
#' @param PCY
#' @param ont
#'
#' @return NULL
#' @export
#'
#' @examples NULL
# se_object = pg_se_imputed
# PCX = 1
# PCY = 2
# ont = c("MF")
# enrich_pca_components <-
#   function(se_object,
#            PCX = 1,
#            PCY = 2,
#            ont = c("MF", "BP", "CC")) {
#     library("clusterProfiler")
#     library("org.Hs.eg.db")
#     pca_y <- as.data.frame(get_reduction_se(se_object, "PCA")$rotation)
#     if (ont %in% c("MF", "BP", "CC")) {
#       GSEA_comp <- map(c(PCX, PCY), function(x) {
#         original_gene_list <- pca_y[, x]
#         names(original_gene_list) <- rownames(pca_y)
#         gene_list <- na.omit(original_gene_list)  # remove NAs
#         gene_list = sort(gene_list, decreasing = TRUE)#Rank by ordering
#         gene_list <-
#           gene_list[abs(gene_list) > 0]# Filter gene list by deg
#         
#         enrichedGO <- gseGO(
#           geneList = gene_list,
#           OrgDb = org.Hs.eg.db,
#           ont = ont,
#           keyType = "SYMBOL",
#           minGSSize = 10,
#           maxGSSize = 100,
#           pvalueCutoff = 0.05,
#           pAdjustMethod = "BH",
#           verbose = TRUE
#         )
#         if(nrow(enrichedGO@result)==0){
#           return(NULL)
#         }else{
#           res_ggo <-
#             enrichedGO@result[enrichedGO@result$p.adjust < 0.05, ]
#           res_ggo$PC<-x  
#           res_ggo
#         }
#       }) %>% compact()%>%do.call(rbind, .)
#     }
#   }
#'

