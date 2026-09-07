"%>%" <- magrittr::"%>%"

.proteolab_get_omnipath_version <- function() {
  if (!requireNamespace("OmnipathR", quietly = TRUE)) {
    return("unknown")
  }

  pkg_version <- tryCatch(
    as.character(utils::packageVersion("OmnipathR")),
    error = function(e) "unknown"
  )

  pkg_version
}

#' create_ptm_references
#'
#' @param ptm_databases 
#' @param database_dir 
#' @details
#' This Function performs an update of the database objects (rds for decoupler and gmt for ssGSEA-PTM) stored in data/ptm_databases. 
#' The new objects will be saved letting the user now the current Omnipath version and the PTM databases. 
#' Kinase substrate information is stored as chracater : 'gene_name_modified position' e.g.EEF2K_S366   
#' 
#' @returns NULL
#' @export
#'
#' @examples NULL
#' # Just run 
#create_ptm_references()
create_ptm_references<-function(ptm_databases=c("ProtMapper","PhosphoNetworks","PhosphoSite","SIGNOR","KEA"),
                                database_dir=here::here("data","ptm_databases")){
  op_version <- .proteolab_get_omnipath_version()
  if(!dir.exists(database_dir)){
    dir.create(database_dir)
  }
  purrr::map(ptm_databases,function(db_search_string){
    prot_mapper_source<-get_gmt_references(ptm_source = db_search_string)  
    net<-prot_mapper_source
    saveRDS(net,here::here(database_dir,paste0("ptm.",tolower(db_search_string),".omnipathversion",op_version,".rds")))
    
    message("Skipping GMT export; PTM-GSEA now uses .rds PTM databases directly.")
    message(paste0("Done PTM Db Update ",Sys.time()))
    return()
  })
}

#' convert_se_to_gct
#'
#' @param se_object SummarizedExperimentObject holding Phopphopeptide Intensisites
#' @param data_type 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
# gct_object<-convert_se_to_gct(se_object = phos_se_raw,"MaxQuant")
convert_se_to_gct <- function(se_object,
                              data_type = c("DIA-NN", "MaxQuant"),
                              fasta_path = NULL) {
  # Example.: LHPPPQLSPFLQPHG-p
  if (data_type == "MaxQuant") {
    
    site_info<-se_object@elementMetadata[,c("leading_proteins","gene","sequence_window","position")]
    peptide <- purrr::map(se_object@elementMetadata$phospho_sty_probabilities, function(x)
      extract_highest_prob_sequence(x),.progress = T) %>% do.call(rbind, .)
    site_info2 <- cbind(site_info, peptide)
    
    PTM_df <- as.data.frame(SummarizedExperiment::assay(se_object))
    PTM_df$Name <- paste0(gsub(";.*","",site_info2$gene),"_",site_info2$AminoAcid,site_info2$position)#paste0(site_info$gene, "_", site_info$residue)
  }
  if (data_type == "DIA-NN") {
    site_info <- get_residue_from_DIANN(diann_se = se_object, fasta_file_path  = fasta_path)
    PTM_df <- as.data.frame(SummarizedExperiment::assay(se_object)[rownames(se_object) %in% rownames(site_info), ])
    PTM_df$Name <- plyr::mapvalues(rownames(PTM_df), rownames(site_info), site_info$residue)
    
  }
  
  if (any(duplicated(PTM_df$Name))) {
    PTM_df_collpased <- PTM_df %>%
      group_by(Name) %>%
      summarise(across(everything(), ~ mean(.x, na.rm = T))) %>%
      as.data.frame()
    rownames(PTM_df_collpased) <- PTM_df_collpased$Name
    PTM_df_collpased$Name <- NULL
    PTM_df_collpased <- as.matrix(PTM_df_collpased)
    PTM_df_collpased[is.nan(PTM_df_collpased)] <- NA
    
  } else{
    
    PTM_df_collpased <- PTM_df
    rownames(PTM_df_collpased) <- PTM_df_collpased$Name
    PTM_df_collpased$Name <- NULL
    PTM_df_collpased <- as.matrix(PTM_df_collpased)
    
  }
  PTM_mat <- as.matrix(PTM_df_collpased)

  return(PTM_mat)
}

#' convert_se_to_gct
#'
#' @param se_object SummarizedExperimentObject holding Phopphopeptide Intensisites
#' @param data_type 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
get_gmt_references<-function(ptm_source=c("ProtMapper","PhosphoNetworks","PhosphoSite","SIGNOR","KEA")){
  uniprot_kinases <- OmnipathR::import_omnipath_annotations(resources = "UniProt_keyword") %>%
    dplyr::filter(value == "Kinase" & !grepl("COMPLEX", uniprot)) %>%
    dplyr::distinct() %>%
    dplyr::pull(genesymbol) %>%
    unique()
  
  omnipath_ptm <- OmnipathR::signed_ptms() %>%
    dplyr::filter(modification %in% c("dephosphorylation","phosphorylation")) %>%
    dplyr::filter(!(stringr::str_detect(sources, ptm_source) & n_resources == 1)) %>%
    dplyr::mutate(p_site = paste0(substrate_genesymbol, "_", residue_type, residue_offset),
                  mor = ifelse(modification == "phosphorylation", 1, -1)) %>%
    dplyr::transmute(p_site, enzyme_genesymbol, mor) %>%
    dplyr::filter(enzyme_genesymbol %in% uniprot_kinases)
  
  omnipath_ptm$likelihood <- 1
  omnipath_ptm$id <- paste(omnipath_ptm$p_site,omnipath_ptm$enzyme_genesymbol, sep ="")
  omnipath_ptm <- omnipath_ptm[!duplicated(omnipath_ptm$id),]
  omnipath_ptm <- omnipath_ptm[,-5]
  omnipath_ptm
}

##### PTM-SSGSEA #####

#' Normalize PTM Scores
#'
#' @param x vector of PTM scores
#' @param sample.norm.type character; one of rank, log, log.rank, none
#'
#' @returns vector with normalized PTM scores
#' @export
#'
#' @examples NULL
## x<-values
## sample.norm.type<-"rank"
.normalize_ptm_scores <- function(x, sample.norm.type = c("rank", "log", "log.rank", "none")) {
  sample.norm.type <- match.arg(sample.norm.type)
  ng <- length(x)
  
  if (sample.norm.type == "none") {
    return(x)
  }
  
  if (sample.norm.type == "log") {
    y <- x
    y[y < 1] <- 1
    return(log(y + exp(1)))
  }
  
  rank_vec <- rep(NA_real_, ng)
  keep_idx <- which(!is.na(x))
  
  x_rank<-rank(x[keep_idx], ties.method = "average")
  rank_vec[keep_idx] <- x_rank
  names(rank_vec)<-names(x)
  
  rank_vec <- 10000 * rank_vec / ng
  
  if (sample.norm.type == "rank") {
    return(rank_vec)
  }
  
  log(rank_vec + exp(1))
}

#' Calculate PTM Enrichment Scores
#'
#' @param values vector of PTM scores for a single sample
#' @param gene_set vector of PTM identifiers for a single kinase set
#' @param statistic character; one of "area.under.RES" or "Kolmogorov-Smirnov"
#' @param weight numeric; score exponent for enrichment weighting
#'
#' @returns numeric; enrichment score for the given sample and kinase set
#' @export
#'
#' @examples NULL 
.calc_ptm_enrichment_score <- function(values,
                                       gene_set,
                                       statistic = c("area.under.RES", "Kolmogorov-Smirnov"),
                                       weight = 0.75) {
  statistic <- match.arg(statistic)
  
  keep <- !is.na(values)
  values <- values[keep]
  if (length(values) == 0) {
    return(NA_real_)
  }
  
  ord <- order(values, decreasing = TRUE)
  values_ord <- values[ord]
  genes_ord <- names(values_ord)
  
  in_set <- genes_ord %in% gene_set
  n_hits <- sum(in_set)
  n_miss <- sum(!in_set)
  if (n_hits == 0 || n_miss == 0) {
    return(NA_real_)
  }
  
  alpha <- max(weight, 0)
  hit_weights <- abs(values_ord) ^ alpha
  hit_weights[!in_set] <- 0
  sum_hit <- sum(hit_weights)
  if (!is.finite(sum_hit) || sum_hit <= 0) {
    hit_weights[in_set] <- 1
    sum_hit <- sum(hit_weights)
  }
  
  p_hit <- cumsum(hit_weights / sum_hit)
  p_miss <- cumsum((!in_set) / n_miss)
  res <- p_hit - p_miss
  
  if (statistic == "area.under.RES") {
    return(mean(res, na.rm = TRUE))
  }
  
  max_pos <- max(res, na.rm = TRUE)
  min_neg <- min(res, na.rm = TRUE)
  if (abs(max_pos) >= abs(min_neg)) {
    return(max_pos)
  }
  min_neg
}

#' Resolve PTM database RDS 
#'
#' @param ptm_db character ; PTM database key (e.g. "ProtMapper") or path to .rds file.
#' @param ptm_database_dir character; directory containing PTM database .rds files.
#'
#' @returns character; path to the resolved .rds file for the specified PTM database.
#' @export
#'
#' @examples NULL
.resolve_ptm_db_rds <- function(ptm_db,
                                ptm_database_dir = here::here("data", "ptm_databases")) {
  if (!is.character(ptm_db) || length(ptm_db) != 1 || !nzchar(ptm_db)) {
    stop("ptm_db must be a single non-empty character string.")
  }
  
  if (file.exists(ptm_db) && grepl("\\.rds$", ptm_db, ignore.case = TRUE)) {
    return(ptm_db)
  }
  
  db_files <- list.files(ptm_database_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(db_files) == 0) {
    stop("No PTM database .rds files found in: ", ptm_database_dir)
  }
  
  hits <- db_files[grepl(tolower(ptm_db), tolower(basename(db_files)), fixed = TRUE)]
  if (length(hits) == 0) {
    stop("No PTM database .rds matching '", ptm_db, "' found in: ", ptm_database_dir)
  }
  
  if (length(hits) > 1) {
    hits <- sort(hits)
    warning("Multiple PTM databases matched '", ptm_db, "'. Using: ", basename(hits[[1]]))
  }
  
  hits[[1]]
}

#' Obtain PTM Matrix from Summarized Experiment
#'
#' @param se_object SummarizedExperiment object containing phosphosite-level values.
#' @param data_type character; one of "auto", "DIA-NN", or "MaxQuant". If "auto", the function will attempt to infer the data type based on the row names and rowData of the SummarizedExperiment.
#' @param fasta_path character; path to FASTA file required for DIA-NN conversion. Only used if data_type is "DIA-NN".
#'
#' @returns matrix; rows are PTM identifiers in the format GENE_RESIDUEPOSITION, columns are sample IDs, and values are the corresponding phosphosite intensities.
#' @export
#'
#' @examples NULL
## .ptm_matrix_from_se(se_object = se_object,
##                     data_type = "auto",
##                     fasta_path = fasta_path)
.ptm_matrix_from_se <- function(se_object,
                                data_type = c("auto", "DIA-NN", "MaxQuant"),
                                fasta_path = NULL) {
  if (!methods::is(se_object, "SummarizedExperiment")) {
    stop("se_object must be a SummarizedExperiment.")
  }
  
  data_type <- match.arg(data_type)
  ptm_df <- as.data.frame(SummarizedExperiment::assay(se_object))
  
  row_ids <- rownames(se_object)
  already_formatted <- grepl("^[^_]+_[A-Z][0-9]+$", row_ids)
  
  if (data_type == "auto") {
    if (all(already_formatted)) {
      data_type <- "formatted"
    } else if ("phospho_sty_probabilities" %in% colnames(SummarizedExperiment::rowData(se_object))) {
      data_type <- "MaxQuant"
    } else {
      data_type <- "DIA-NN"
    }
  }
  
  if (data_type == "formatted") {
    ptm_df$Name <- row_ids
  }
  
  if (data_type == "MaxQuant") {
    row_data <- as.data.frame(SummarizedExperiment::rowData(se_object))
    required_cols <- c("gene", "position", "phospho_sty_probabilities")
    if (!all(required_cols %in% colnames(row_data))) {
      stop("MaxQuant mode requires rowData columns: gene, position, phospho_sty_probabilities")
    }
    peptide <- purrr::map(row_data$phospho_sty_probabilities, function(x) {
      extract_highest_prob_sequence(x)
    }) %>% do.call(rbind, .)
    ptm_df$Name <- paste0(gsub(";.*", "", row_data$gene), "_", peptide$AminoAcid, row_data$position)
  }
  
  if (data_type == "DIA-NN") {
    if (is.null(fasta_path) || !nzchar(fasta_path)) {
      stop("DIA-NN mode requires fasta_path.")
    }
    site_info <- get_residue_from_DIANN(diann_se = se_object, fasta_file_path = fasta_path)
    keep <- rownames(se_object) %in% rownames(site_info)
    ptm_df <- ptm_df[keep, , drop = FALSE]
    ptm_df$Name <- plyr::mapvalues(rownames(ptm_df), rownames(site_info), site_info$residue)
  }
  
  if (!"Name" %in% colnames(ptm_df)) {
    stop("Failed to derive PTM identifiers in format GENE_RESIDUEPOSITION.")
  }
  
  if (any(duplicated(ptm_df$Name))) {
    ptm_df <- ptm_df %>%
      dplyr::group_by(Name) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")%>%
      as.data.frame()
  }
  rownames(ptm_df) <- ptm_df$Name
  ptm_df$Name <- NULL
  
  out <- as.matrix(ptm_df)
  out[is.nan(out)] <- NA
  out
}

#' run_ptm_gsea
#'
#' Native PTM-GSEA implementation using SummarizedExperiment input and PTM
#' databases stored as .rds files in data/ptm_databases.
#'
#' @param se_object SummarizedExperiment object containing phosphosite-level values.
#' @param ptm_db character; PTM database key (e.g. "ProtMapper") or path to .rds file.
#' @param ptm_database_dir character; directory containing PTM database .rds files.
#' @param output.prefix character; optional tag for output metadata.
#' @param output.directory character; optional output directory for parameter log.
#' @param sample.norm.type character; one of rank, log, log.rank, none.
#' @param weight numeric; score exponent for enrichment weighting.
#' @param statistic character; one of area.under.RES or Kolmogorov-Smirnov.
#' @param output.score.type character; one of NES or ES.
#' @param nperm integer; number of permutations for p-values and NES.
#' @param min.overlap integer; minimum overlap required between sample and gene set.
#' @param global.fdr logical; if TRUE adjust p-values globally, otherwise per sample.
#' @param extended.output logical; kept for API compatibility.
#' @param param.file logical; write a parameter file into output.directory.
#' @param log.file character; path to run log.
#' @param data_type character; auto, DIA-NN, or MaxQuant.
#' @param fasta_path character; FASTA path required for DIA-NN conversion.
#' @param verbose boolean; if TRUE, print progress messages to console.
#'
#' @returns data.frame in long format with columns: score, p_value, fdr, kinase, sample_id.
#' @export
#'
#' @examples NULL
## ptm_database_dir = here::here("data", "ptm_databases")
## output.prefix = "ProteoLab"
## output.directory = "."
## sample.norm.type = c("rank")
## weight = 0.75
## statistic = c("area.under.RES")
## output.score.type = c("NES")
## nperm = 1000
## min.overlap = 10
## global.fdr = FALSE
## extended.output = TRUE
## par = TRUE
## param.file = TRUE
## log.file = "run.log"
## data_type = c("auto", "DIA-NN", "MaxQuant")
## fasta_path = NULL
## verbose=FALSE
run_ptm_gsea <- function(se_object,
                         ptm_db,
                         ptm_database_dir = here::here("data", "ptm_databases"),
                         output.prefix = "ProteoLab",
                         output.directory = ".",
                         sample.norm.type = c("rank", "log", "log.rank", "none"),
                         weight = 0.75,
                         statistic = c("area.under.RES", "Kolmogorov-Smirnov"),
                         output.score.type = c("NES", "ES"),
                         nperm = 1000,
                         min.overlap = 10,
                         global.fdr = FALSE,
                         extended.output = TRUE,
                         par = FALSE,
                         param.file = TRUE,
                         log.file = "run.log",
                         data_type = c("auto", "DIA-NN", "MaxQuant"),
                         fasta_path = NULL,
                         verbose=FALSE
) {
  
  sample.norm.type <- match.arg(sample.norm.type)
  statistic <- match.arg(statistic)
  output.score.type <- match.arg(output.score.type)
  data_type <- match.arg(data_type)
  
  if (!methods::is(se_object, "SummarizedExperiment")) {
    stop("se_object must be a SummarizedExperiment.")
  }
  
  ptm_db_rds <- .resolve_ptm_db_rds(ptm_db = ptm_db, ptm_database_dir = ptm_database_dir)
  ptm_db_df <- readRDS(ptm_db_rds)
  
  req_cols <- c("enzyme_genesymbol", "p_site")
  if (!all(req_cols %in% colnames(ptm_db_df))) {
    stop("PTM database .rds must contain columns: enzyme_genesymbol, p_site")
  }
  
  gene_sets <- split(ptm_db_df$p_site, ptm_db_df$enzyme_genesymbol)
  gene_sets <- lapply(gene_sets, unique)
  
  score_mat <- .ptm_matrix_from_se(
    se_object = se_object,
    data_type = data_type,
    fasta_path = fasta_path
  )
  
  if (nrow(score_mat) == 0 || ncol(score_mat) == 0) {
    stop("No analyzable PTM values found in se_object.")
  }
  
  gene_names <- rownames(score_mat)
  keep_sets <- vapply(gene_sets, function(gs) {
    sum(gs %in% gene_names) >= min.overlap
  }, logical(1))
  gene_sets <- gene_sets[keep_sets]
  if (length(gene_sets) == 0) {
    stop("No PTM sets pass min.overlap against the provided se_object.")
  }
  
  if (!dir.exists(output.directory)) {
    dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (isTRUE(param.file)) {
    param_lines <- c(
      paste("##", Sys.time()),
      paste("output.prefix =", output.prefix, sep = "\t"),
      paste("ptm_db =", basename(ptm_db_rds), sep = "\t"),
      paste("sample.norm.type =", sample.norm.type, sep = "\t"),
      paste("weight =", weight, sep = "\t"),
      paste("statistic =", statistic, sep = "\t"),
      paste("output.score.type =", output.score.type, sep = "\t"),
      paste("nperm =", nperm, sep = "\t"),
      paste("min.overlap =", min.overlap, sep = "\t"),
      paste("global.fdr =", global.fdr, sep = "\t"),
      paste("extended.output =", extended.output, sep = "\t")
    )
    writeLines(param_lines, con = file.path(output.directory, paste0(output.prefix, "_parameters.txt")))
  }
  
  cat("##", format(Sys.time()), "\n", file = log.file)
  cat("Running native PTM-GSEA on", ncol(score_mat), "samples and", length(gene_sets), "kinase sets\n", file = log.file, append = TRUE)
  
  sample_ids <- colnames(score_mat)
  kinase_ids <- names(gene_sets)
  alpha <- max(weight, 0)
  do_perm <- is.finite(nperm) && nperm > 0
  
  # Fast ES calculator matching ssGSEA2::gsea_score + ssGSEA2::score
  .es_from_order <- function(order_idx, set_mask, abs_vals_pow, statistic) {
    n_total <- length(order_idx)
    hit_pos <- which(set_mask[order_idx])
    n_hits <- length(hit_pos)
    n_miss <- n_total - n_hits
    if (n_hits == 0L || n_miss == 0L) {
      return(NA_real_)
    }

    hit_w <- abs_vals_pow[order_idx[hit_pos]]
    sum_hit <- sum(hit_w)
    if (!is.finite(sum_hit) || sum_hit <= 0) {
      hit_w <- rep(1, n_hits)
      sum_hit <- n_hits
    }

    up <- hit_w / sum_hit
    gaps <- (c(hit_pos - 1L, n_total) - c(0L, hit_pos))
    down <- gaps / n_miss
    res <- cumsum(c(up, up[n_hits]) - down)
    valleys <- res[seq_len(n_hits)] - up
    max_es <- max(res)
    min_es <- min(valleys)

    if (statistic == "Kolmogorov-Smirnov") {
      if (max_es > -min_es) {
        return(signif(max_es, digits = 5))
      }
      return(signif(min_es, digits = 5))
    }

    gaps <- gaps + 1
    auc <- c(valleys, 0) * gaps + 0.5 * (c(0, res) - c(valleys, 0)) * gaps
    sum(auc)
  }
  
  process_sample <- function(sid) {
    if (verbose) cat("Processing sample:", sid, "\n")
    
    values <- score_mat[, sid]
    names(values) <- rownames(score_mat)
    values <- .normalize_ptm_scores(values, sample.norm.type = sample.norm.type)
    
    keep    <- !is.na(values)
    vals_ok <- values[keep]
    genes_ok <- names(vals_ok)
    if (length(vals_ok) < 2L) {
      return(list())
    }

    # Keep original gene index space and cache absolute weights once per sample.
    n_total <- length(vals_ok)
    abs_vals_pow <- abs(vals_ok)^alpha
    obs_order <- order(vals_ok, decreasing = TRUE)
    
    sample_res <- vector("list", length(kinase_ids))
    k_idx <- 1L
    
    for (kin in kinase_ids) {
      gs <- intersect(gene_sets[[kin]], genes_ok)
      if (length(gs) < min.overlap) next
      
      set_mask <- genes_ok %in% gs
      n_set <- sum(set_mask)

      es <- .es_from_order(
        order_idx = obs_order,
        set_mask = set_mask,
        abs_vals_pow = abs_vals_pow,
        statistic = statistic
      )
      
      p_val     <- NA_real_
      score_val <- es
      
      if (do_perm) {
        # Match ssGSEA2 null: randomize ordered gene list (sample(1:n.rows)).
        perm_es <- vapply(seq_len(nperm), function(i) {
          .es_from_order(
            order_idx = sample.int(n_total),
            set_mask = set_mask,
            abs_vals_pow = abs_vals_pow,
            statistic = statistic
          )
        }, numeric(1))
        
        if (output.score.type == "NES") {
          if (es >= 0) {
            pos_phi <- perm_es[perm_es >= 0]
            if (length(pos_phi) == 0) pos_phi <- 0.5
            pos_m     <- mean(pos_phi)
            score_val <- es / pos_m
            s         <- sum(pos_phi >= es) / length(pos_phi)
          } else {
            neg_phi <- perm_es[perm_es < 0]
            if (length(neg_phi) == 0) neg_phi <- 0.5
            neg_m     <- mean(neg_phi)
            score_val <- es / abs(neg_m)
            s         <- sum(neg_phi <= es) / length(neg_phi)
          }
        } else {
          if (es >= 0) {
            pos_phi <- perm_es[perm_es >= 0]
            if (length(pos_phi) == 0) pos_phi <- 0.5
            s <- sum(pos_phi >= es) / length(pos_phi)
          } else {
            neg_phi <- perm_es[perm_es < 0]
            if (length(neg_phi) == 0) neg_phi <- 0.5
            s <- sum(neg_phi <= es) / length(neg_phi)
          }
        }
        p_val <- ifelse(s == 0, 1 / nperm, s)
      }
      
      sample_res[[k_idx]] <- data.frame(
        score     = score_val,
        p_value   = p_val,
        kinase    = kin,
        sample_id = sid,
        n_overlap = n_set,
        stringsAsFactors = FALSE
      )
      k_idx <- k_idx + 1L
    }
    sample_res[!vapply(sample_res, is.null, logical(1))]
  }
  
  if (isTRUE(par)) {
    n_cores <- max(1L, parallel::detectCores() - 1L)
    res_nested <- parallel::mclapply(sample_ids, process_sample, mc.cores = n_cores)
  } else {
    res_nested <- lapply(sample_ids, process_sample)
  }
  
  res_list <- unlist(res_nested, recursive = FALSE)
  if (length(res_list) == 0) {
    stop("No kinase scores could be computed. Check overlap and identifier format.")
  }
  
  res <- do.call(rbind, res_list)
  if (isTRUE(global.fdr)) {
    res$fdr <- stats::p.adjust(res$p_value, method = "BH")
  } else {
    res <- res %>%
      dplyr::group_by(sample_id) %>%
      dplyr::mutate(fdr = stats::p.adjust(p_value, method = "BH")) %>%
      dplyr::ungroup()
  }
  res
}


#' make_ssGSEA2_output_long
#'
#' @param gct_obj 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
# signature_gct_path<-"/mnt/NAS4/user_data/np-dennis/projects/Proteomics/ProteoLab/projects/60_DDA_Phospho_Jurkat_MCF7_120min/results/id_1/ssGSEA_PTM/SingleSample/2026-06-10_16-54-20/"
# gct_obj<-read_ssGSEA2_output(signature_gct_path,file_pattern="ProteoLab-combined")
# make_ssGSEA2_output_long(gct_obj = gct_obj)
# openxlsx::write.xlsx(x = make_ssGSEA2_output_long(gct_obj = gct_obj),
#                      file = paste0( "/mnt/NAS4/user_data/np-dennis/projects/Proteomics/ProteoLab/projects/60_DDA_Phospho_Jurkat_MCF7_120min/results/id_1/ssGSEA_PTM/SingleSample/2026-06-10_16-54-20//ssgsea_res.xlsx"))
make_ssGSEA2_output_long<-function(gct_obj){
  if (is.data.frame(gct_obj)) {
    req_cols <- c("score", "p_value", "fdr", "kinase", "sample_id")
    if (all(req_cols %in% colnames(gct_obj))) {
      return(gct_obj)
    }
  }

  if (!methods::is(gct_obj, "GCT")) {
    stop("gct_obj must be a cmapR::GCT object or a long-format data.frame.")
  }

  enri_scores_mat<-as.data.frame(gct_obj@mat)
  enri_scores_mat<-janitor::clean_names(enri_scores_mat)
  samples<-colnames(enri_scores_mat)
  
  enri_scores_mat$id<-rownames(enri_scores_mat)
  enri_scores_stat<-as.data.frame(gct_obj@rdesc)%>%janitor::clean_names(.)
  ssgsea_combined<-cbind(enri_scores_mat,enri_scores_stat)
  
  ### Make long format
  long_ssgsea<-purrr::map(samples,function(tmp_sample){
    searchtag<-c(tmp_sample,paste0("pvalue_",tmp_sample),paste0("fdr_pvalue_",tmp_sample))
    
    ssgsea_combined_tmp<-ssgsea_combined[,searchtag]
    ssgsea_combined_tmp
    colnames(ssgsea_combined_tmp)<-c("score","p_value","fdr")
    ssgsea_combined_tmp$kinase<-rownames(ssgsea_combined_tmp)
    ssgsea_combined_tmp$sample_id<-tmp_sample
    rownames(ssgsea_combined_tmp)<-NULL
    ssgsea_combined_tmp
  })%>%do.call(rbind,.)
  return(long_ssgsea)
}

#' read_ssGSEA2_output
#'
#' @param signature_gct_path 
#' @param file_pattern 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#signature_gct_path<-"/mnt/NAS4/user_data/np-dennis/projects/Proteomics/ProteoLab/projects/60_DDA_Phospho_Jurkat_MCF7_120min/results/id_1/ssGSEA_PTM/SingleSample/2026-06-10_16-54-20/"
# test<-read_ssGSEA2_output(signature_gct_path,file_pattern="ProteoLab-combined")
read_ssGSEA2_output <- function(signature_gct_path,file_pattern="ProteoLab-combined") {
  ssgsea_xlsx <- file.path(signature_gct_path, "ssgsea_res.xlsx")
  if (file.exists(ssgsea_xlsx)) {
    return(openxlsx::read.xlsx(ssgsea_xlsx))
  }

  ssgsea_results_combined <- list.files(signature_gct_path,
                                        pattern = file_pattern,
                                        full.names = T)
  if (length(ssgsea_results_combined) == 0) {
    stop("No ssGSEA output found in directory: ", signature_gct_path)
  }

  if (!requireNamespace("cmapR", quietly = TRUE)) {
    stop("No ssgsea_res.xlsx found and cmapR is not available for parsing legacy gctx output.")
  }

  cmapR::parse_gctx(ssgsea_results_combined)
}

# Helper to compute hierarchical sample ordering based on score similarity.
#' get_hclust_order_by_score
#'
#' @param df dataframe; lng-format data frame containing features, samples, and scores
#' @param feature_col character; name of column containing featurees
#' @param sample_col character; name of column containing samples 
#' @param score_col character; name of column containing scores
#' @param cluster_distance character; method for distance calculation
#' @param cluster_method character; method for hierarchical clustering
#' @param order_dimension character; either "sample" or "feature" to specify which dimension to order
#' 
#' @returns ordered rownames of df 
#' @export
#'
#' @examples NULL
get_hclust_order_by_score <- function(df,
                                      feature_col,
                                      sample_col,
                                      score_col,
                                      cluster_distance = "euclidean",
                                      cluster_method = "complete",
                                      order_dimension = c("sample", "feature")) {
  order_dimension <- match.arg(order_dimension)
  df_tmp <- df[, c(feature_col, sample_col, score_col), drop = FALSE]
  colnames(df_tmp) <- c("feature", "sample", "score")
  df_tmp <- df_tmp[!is.na(df_tmp$feature) & !is.na(df_tmp$sample), , drop = FALSE]

  sample_order <- unique(as.character(df_tmp$sample))
  feature_order <- unique(as.character(df_tmp$feature))
  if (nrow(df_tmp) == 0) {
    if (order_dimension == "sample") {
      return(sample_order)
    }
    return(feature_order)
  }

  score_mat <- reshape2::acast(
    df_tmp,
    feature ~ sample,
    value.var = "score",
    fun.aggregate = mean,
    fill = 0
  )

  if (is.null(dim(score_mat))) {
    if (order_dimension == "sample") {
      return(sample_order)
    }
    return(feature_order)
  }

  if (order_dimension == "sample") {
    if (ncol(score_mat) < 2) {
      return(sample_order)
    }
    hc <- stats::hclust(
      stats::dist(t(score_mat), method = cluster_distance),
      method = cluster_method
    )
    return(colnames(score_mat)[hc$order])
  }

  if (nrow(score_mat) < 2) {
    return(feature_order)
  }
  hc <- stats::hclust(
    stats::dist(score_mat, method = cluster_distance),
    method = cluster_method
  )
  rownames(score_mat)[hc$order]
}

#' get_hclust_object_by_score
#'
#' @param df dataframe; lng-format data frame containing features, samples, and scores
#' @param feature_col character; name of column containing featurees
#' @param sample_col character; name of column containing samples
#' @param score_col character; name of column containing scores
#' @param order_dimension character either sample or feature
#' @param cluster_distance character; method for distance calculation
#' @param cluster_method character; method for hierarchical clustering
#'
#' @returns hierarchical clustering object
#' @export
#'
#' @examples NULL
get_hclust_object_by_score <- function(df,
                                       feature_col,
                                       sample_col,
                                       score_col,
                                       order_dimension = c("sample", "feature"),
                                       cluster_distance = "euclidean",
                                       cluster_method = "complete") {
  order_dimension <- match.arg(order_dimension)
  df_tmp <- df[, c(feature_col, sample_col, score_col), drop = FALSE]
  colnames(df_tmp) <- c("feature", "sample", "score")
  df_tmp <- df_tmp[!is.na(df_tmp$feature) & !is.na(df_tmp$sample), , drop = FALSE]

  if (nrow(df_tmp) == 0) {
    return(NULL)
  }

  score_mat <- reshape2::acast(
    df_tmp,
    feature ~ sample,
    value.var = "score",
    fun.aggregate = mean,
    fill = 0
  )

  if (is.null(dim(score_mat))) {
    return(NULL)
  }

  if (order_dimension == "sample") {
    if (ncol(score_mat) < 2) {
      return(NULL)
    }
    return(stats::hclust(
      stats::dist(t(score_mat), method = cluster_distance),
      method = cluster_method
    ))
  }

  if (nrow(score_mat) < 2) {
    return(NULL)
  }
  stats::hclust(
    stats::dist(score_mat, method = cluster_distance),
    method = cluster_method
  )
}

#' make_hclust_tree_plot
#'
#' @param hc hierachical clustering object
#' @param orientation character; either "top" or "left" to specify orientation of the dendrogram
#'
#' @returns ggplot object
#' @export
#'
#' @examples NULL
make_hclust_tree_plot <- function(hc, orientation = c("top", "left")) {
  orientation <- match.arg(orientation)
  if (is.null(hc)) {
    return(NULL)
  }

  n_leaf <- length(hc$order)
  if (n_leaf < 2) {
    return(NULL)
  }

  leaf_pos <- numeric(n_leaf)
  leaf_pos[hc$order] <- seq_len(n_leaf)
  node_x <- numeric(n_leaf - 1)
  node_h <- hc$height

  segs <- data.frame("x" = numeric(0), "y" = numeric(0), "xend" = numeric(0), "yend" = numeric(0))
  
  for (i in seq_len(n_leaf - 1)) {
    left <- hc$merge[i, 1]
    right <- hc$merge[i, 2]

    left_x <- if (left < 0) leaf_pos[-left] else node_x[left]
    right_x <- if (right < 0) leaf_pos[-right] else node_x[right]
    left_h <- if (left < 0) 0 else node_h[left]
    right_h <- if (right < 0) 0 else node_h[right]
    parent_h <- node_h[i]

    segs <- rbind(
      segs,
      data.frame("x" = left_x, "y" = left_h, "xend" = left_x, yend = parent_h),
      data.frame("x" = right_x, "y" = right_h, "xend" = right_x, yend = parent_h),
      data.frame("x" = left_x, "y" = parent_h, "xend" = right_x, yend = parent_h)
    )

    node_x[i] <- mean(c(left_x, right_x))
  }

  if (orientation == "top") {
    return(
      ggplot2::ggplot(segs) +
        ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.4, color = "grey20") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(0, 0, 0, 0)
        )
    )
  }

  segs_left <- data.frame(
    x = segs$y,
    y = segs$x,
    xend = segs$yend,
    yend = segs$xend
  )

  ggplot2::ggplot(segs_left) +
    ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.4, color = "grey20") +
    # Reverse x so leaf tips are drawn on the right, toward the dotplot panel.
    ggplot2::scale_x_reverse(expand = expansion(mult = c(0, 0))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
}

#' make_ssgsea_dotplot
#'
#' @param ssgsea_long 
#' @param by_group 
#' @param group 
#' @param by_group_color 
#' @param nkin 
#' @param color_gradient 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
# signature_gct_path<-"/mnt/NAS4/user_data/np-dennis/projects/Proteomics/ProteoLab/projects/masterStudents26/results/id_1/ssGSEA_PTM/2026-05-22_13-32-02/"
# gct_obj_tmp <- read_ssGSEA2_output(signature_gct_path = signature_gct_path)
# ssgsea_res <- make_ssGSEA2_output_long(gct_obj = gct_obj_tmp)
# ssgsea_long<-ssgsea_res
# sample_test<-unique(ssgsea_long$sample_id)
# names(sample_test)<-c(rep("A",7),rep("B",8))
# ssgsea_long$group<-plyr::mapvalues(ssgsea_res$sample_id,
#                                    sample_test,
#                                    names(sample_test))
# make_ssgsea_dotplot()
# make_ssgsea_dotplot(ssgsea_long = ssgsea_long,show_significance = F)
make_ssgsea_dotplot<-function(ssgsea_long,
                              group="none",
                              by_group_color=NULL,
                              show_significance=T,
                              cluster_samples = FALSE,
                              cluster_distance = "euclidean",
                              cluster_method = "complete",
                              nkin = 50,xsize=10,ysize=10,
                              legend_size = 10,
                              color_gradient = c("darkblue", "whitesmoke", "indianred")){
  
  ssgsea_long<-ssgsea_long[order(ssgsea_long$p_value,decreasing = F),]
  if (!is.null(nkin) && is.numeric(nkin) && nkin > 0) {
    top_kinases <- ssgsea_long %>%
      dplyr::group_by(kinase) %>%
      dplyr::summarize(best_p = min(p_value, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(best_p) %>%
      head(nkin) %>%
      dplyr::pull(kinase)
    ssgsea_long <- ssgsea_long[ssgsea_long$kinase %in% top_kinases, ]
  }
  ssgsea_long$log10_p_value<-(-log10(ssgsea_long$p_value))
  ssgsea_long$log10_fdr<-(-log10(ssgsea_long$fdr))
  
  ssgsea_long$sample_id <- factor(ssgsea_long$sample_id, levels = unique(ssgsea_long$sample_id))
  sample_order <- levels(ssgsea_long$sample_id)
  kinase_order <- unique(as.character(ssgsea_long$kinase))

  if (isTRUE(cluster_samples)) {
    clustered_order <- get_hclust_order_by_score(
      df = ssgsea_long,
      feature_col = "kinase",
      sample_col = "sample_id",
      score_col = "score",
      cluster_distance = cluster_distance,
      cluster_method = cluster_method,
      order_dimension = "sample"
    )
    if (length(clustered_order) > 0) {
      sample_order <- clustered_order
    }

    clustered_kinase_order <- get_hclust_order_by_score(
      df = ssgsea_long,
      feature_col = "kinase",
      sample_col = "sample_id",
      score_col = "score",
      cluster_distance = cluster_distance,
      cluster_method = cluster_method,
      order_dimension = "feature"
    )
    if (length(clustered_kinase_order) > 0) {
      kinase_order <- clustered_kinase_order
    }
  }
  
  if(!("group" %in% colnames(ssgsea_long))){
    group<-"none"
  }
  
  if (group!="none") {
    
    group_tmp <- distinct(ssgsea_long[, c("sample_id", "group")])
    if (!isTRUE(cluster_samples)) {
      group_tmp <- group_tmp[order(group_tmp$group), ]
    }
    group_tmp$class <- ""
    group_tmp$sample_id <- factor(group_tmp$sample_id, levels = sample_order)
    sample_order <- levels(group_tmp$sample_id)
    
    group_ballons <- ggplot2::ggplot(group_tmp, ggplot2::aes(x = sample_id, y = class, fill = group)) +
      ggplot2::geom_tile() +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "right",
        plot.background = element_rect(fill = "white", colour = NA)
      )
    if (!is.null(by_group_color)) {
      group_ballons <- group_ballons + ggplot2::scale_fill_manual(values = by_group_color)
    }
  }
  ### Adjust Size
  if (length(unique(ssgsea_long$Kinase)) > 80) {
    ysize <- 8
  }
  if (length(unique(ssgsea_long$Kinase)) > 120) {
    ysize <- 5
    xsize <- 7
  }
  
  ssgsea_long$Significant <- ifelse(ssgsea_long$fdr < 0.05, "FPR<0.05", "FPR>0.05")
  levels(ssgsea_long$sample_id) <- sample_order
  ssgsea_long$kinase <- factor(ssgsea_long$kinase, levels = kinase_order)

  # Dendrogram rendering disabled by request.
  
  if(show_significance){
    ssGSEA_ballons <- ggplot2::ggplot(ssgsea_long, ggplot2::aes(x = sample_id, y = kinase, color = Significant, size = log10_fdr)) +
      ggplot2::geom_point() +
      cowplot::theme_cowplot() +
      ggplot2::theme(axis.line = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = xsize, angle = 90, vjust = 0.5, hjust = 1)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = ysize)) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = legend_size),
            legend.title = ggplot2::element_text(size = legend_size)) +
      ggplot2::ylab("") +
      ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
      ggplot2::labs(x = "Sample", y = "Kinase") +
      ggplot2::scale_color_manual(values = c(
        "FPR>0.05" = "#D3D3D3",
        "FPR<0.05" = "#FFB3BA")) +
      ggplot2::scale_size_continuous(range = c(1, 7), name = "-log10(FDR)")
  }else{
    ssGSEA_ballons <- ggplot2::ggplot(ssgsea_long, ggplot2::aes(x = sample_id, y = kinase, color = score, size = log10_fdr)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_viridis_c(name = 'Score') +
      cowplot::theme_cowplot() +
      ggplot2::theme(axis.line = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = xsize, angle = 90, vjust = 0.5, hjust = 1)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = ysize)) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = legend_size),
            legend.title = ggplot2::element_text(size = legend_size)) +
      ggplot2::ylab("") +
      ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
      ggplot2::labs(x = "Sample", y = "Kinase") +
      ggplot2::scale_size_continuous(range = c(1, 7), name = "-log10(FDR)")
  }
  if ((group!="none")) {
    output <- cowplot::plot_grid(group_ballons, ssGSEA_ballons, ncol = 1, align = "v", rel_heights = c(0.45, 6.55))
  } else{
    output <- ssGSEA_ballons
  }
  return(output)
}

#' make_ssgsea_barplot
#'
#' @param ssgsea_long data.frame; long-format ssGSEA output
#' @param nkin integer; number of kinases to include
#' @param show_significance logical; if TRUE plot -log10(p-value), else plot score
#' @param color_gradient character vector of length 3; low/mid/high colors
#'
#' @returns ggplot object
#' @export
#'
#' @examples NULL
#' # make_ssgsea_barplot(ssgsea_long = ssgsea_res, nkin = 30)
make_ssgsea_barplot <- function(ssgsea_long,
                                nkin = 30,
                                show_significance = FALSE,
                                color_gradient = c("darkblue", "whitesmoke", "indianred")) {
  ssgsea_long <- ssgsea_long[order(ssgsea_long$p_value, decreasing = FALSE), ]
  top_kinases <- ssgsea_long %>%
    dplyr::group_by(kinase) %>%
    dplyr::summarize(best_p = min(p_value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(best_p) %>%
    head(nkin) %>%
    dplyr::pull(kinase)

  ssgsea_plot_df <- ssgsea_long %>%
    dplyr::filter(kinase %in% top_kinases)

  if (show_significance) {
    p <- ggplot2::ggplot(ssgsea_plot_df, ggplot2::aes(x = reorder(kinase, p_value), y = -log10(p_value), fill = -log10(p_value))) +
      ggplot2::geom_col(alpha = 0.9) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~sample_id, scales = "free_x") +
      ggplot2::scale_fill_gradient2(low = color_gradient[1], mid = color_gradient[2], high = color_gradient[3]) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Kinase", y = "-log10(p-value)")
  } else {
    p <- ggplot2::ggplot(ssgsea_plot_df, ggplot2::aes(x = reorder(kinase, score), y = score, fill = score)) +
      ggplot2::geom_col(alpha = 0.9) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~sample_id, scales = "free_x") +
      ggplot2::scale_fill_gradient2(low = color_gradient[1], mid = color_gradient[2], high = color_gradient[3]) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Kinase", y = "ssGSEA Score")
  }

  return(p)
}

##### KSTAR #####

#' make_KSTAR_input
#' @param se_object SummarizedExperimentObject holding Phopphopeptide Intensisites
#' @param data_type character; Indicate wheter DIANN or MaxQuant was used for primary analysis
#' @param show_significance logical; if TRUE plot -log10(p-value), else plot score
#' @examples NULL
make_KSTAR_input <- function(se_object, data_type=c("DIA-NN", "MaxQuant"), fasta_file = NULL) {

  # Example of first three mandatory, a different PTM_DF is needed
  # query_accession	mod_sites	peptide
  # Q6DHY5	S112	GMPMNIRGPMWsVLLNIEEMK
  
  if (data_type == "MaxQuant") {
    site_info<-se_object@elementMetadata[,c("leading_proteins","sequence_window","position")]
    peptide <- purrr::map(se_object@elementMetadata$phospho_sty_probabilities, function(x)
      extract_highest_prob_sequence(x),.progress = T) %>% do.call(rbind, .)
    site_info2 <- cbind(site_info, peptide)
    
    #The modified Site has to be in lower case
    KSTAR_mandatory <- data.frame(
      "query_accession" = gsub(";.*","",site_info2$leading_proteins),
      "mod_sites" = paste0(site_info2$AminoAcid,site_info2$position),
      "peptide" = site_info2$peptide
    )
    # To compare between samples scale 
    KSTAR_DF <- as.data.frame((SummarizedExperiment::assay(se_object)))
    colnames(KSTAR_DF)<-paste0("data:",colnames(KSTAR_DF))
    KSTAR_input <- as.data.frame(cbind(KSTAR_mandatory, KSTAR_DF))
    KSTAR_input[is.na(KSTAR_input)] <- 0
    
  }
  if (data_type == "DIA-NN") {
    if (is.null(fasta_file) || !nzchar(fasta_file) || !file.exists(fasta_file)) {
      stop("DIA-NN KSTAR input requires a valid fasta_file path.")
    }
    site_info <- get_residue_from_DIANN(diann_se = se_object,fasta_file_path = fasta_file)
    
    peptdide <- gsub("[(].*[)]", "_", rownames(site_info)) %>% strsplit(., "_") %>%
      purrr::map(., function(x) {
        mod_seq <- strsplit(x[1], "") %>% unlist()
        mod_seq[length(mod_seq)] <- tolower(mod_seq[length(mod_seq)])
        paste0(paste0(mod_seq, collapse = ""), x[2])
      }) %>% unlist()
    site_info$peptdide <- peptdide
    
    #The modified Site has to be in lower case
    KSTAR_mandatory <- data.frame(
      "query_accession" = site_info$uniprot_id,
      "mod_sites" = gsub(".*_", "", site_info$residue),
      "peptide" = site_info$peptdide
    )
    KSTAR_DF <- as.data.frame(SummarizedExperiment::assay(se_object))
    colnames(KSTAR_DF)<-paste0("data:",colnames(KSTAR_DF))
    KSTAR_input <- as.data.frame(cbind(KSTAR_mandatory, KSTAR_DF))
    KSTAR_input[is.na(KSTAR_input)] <- 0
    KSTAR_input[(KSTAR_input$query_accession!=0),]
  }
  return(KSTAR_input)
}
#' run_kstar_mapping
#'
#' @param data_path 
#' @param output_dir 
#' @param python_env 
#' @param pyscript_dir 
#' @param work_dir 
#' @param network_dir 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#' run_kstar_mapping( )
# run_kstar_mapping(
#   data_path = "/mnt/NAS4/user_data/np-dennis/projects/Proteomics/ProteoLab/projects/Sievers_GBMPhopsho_test/results/id_9/KSTAR_input/DEA/2026-07-10_10-01-08/KSTAR_input.xlsx",
#   project = "Sievers_GBMPhopsho_test",
#   output_dir = "./projects/Sievers_GBMPhopsho_test/results/id_9/KSTAR_input/DEA/test2/",
#   python_env = python_path1,
#   pyscript_dir = here::here("python_scripts", "Run_KSTAR_mapping.py"),
#   network_dir = here::here("data", "NETWORKS"),
#   work_dir = here::here()
# )
run_kstar_mapping<-function(data_path,
                            project,
                            output_dir,
                            python_env,
                            pyscript_dir,
                            work_dir,
                            network_dir){
  
  res<-system2(
    command = python_env,
    args = (
      c(
        pyscript_dir,
        "-p",
        work_dir,
        "-n",
        network_dir,
        "-d",
        data_path,
        "-o",
        output_dir,
        "-i",
        project
      )
    ),
    stderr = T,
    stdout = T
  )
  return(res)
}

#' run_kstar_mapping
#'
#' @param data_path 
#' @param output_dir 
#' @param python_env 
#' @param pyscript_dir 
#' @param work_dir 
#' @param network_dir 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#'run_kstar()
run_kstar<-function(kstarinput,#test "./projects/Sievers_GBMPhopsho_test/results/ptm_se_raw/KSTAR_input/2026-05-29_11-36-37/MAPPED_DATA/Sievers_GBMPhopsho_test_mapped.csv"
                            name,#test "Sievers_GBMPhopsho_test"
                            output,#test "./projects/Sievers_GBMPhopsho_test/results/ptm_se_raw/KSTAR_input/2026-05-27_11-16-32/"
                            python_env,#default "/home/users/np/dennis/.local/share/mamba/envs/kstar311//bin/python3.11"
                            pyscript_dir, #default "./python_scripts/Run_KSTAR.py"
                            workpath, #default "/mnt/NAS4/user_data/np-dennis/projects/Proteomics/ProteoLab/"
                            network, #default "./data/NETWORKS"
                            processes=6,
                            threshold=0.2){
  
  res<-system2(
    command = python_env,
    args = (
      c(
        pyscript_dir,
        "-p",
        workpath,
        "-n",
        network,
        "-d",
        kstarinput,
        "-o",
        output,
        "-i",
        name,
        "-t",
        threshold,
        "-c",
        processes
      )
    ),
    stderr = T,
    stdout = T
  )
  return(res)
}

#' load_KSTAR_res_from_dir
#'
#' @param kstar_res_dir 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#' kstar_res_dir<-"./projects//Sievers_GBMPhopsho_test/results/ptm_se_raw/KSTAR_input/2026-05-29_15-59-24/RESULTS/"
load_KSTAR_res_from_dir<-function(kstar_res_dir){
  kinase<-c("ST","Y")
  purrr::map(kinase,function(kin){
    res_dir<-here::here(kstar_res_dir,kin)
    assertthat::assert_that(dir.exists(res_dir))
    
    tmp_files<-list.files(res_dir)
    fpr<-data.table::fread(here::here(res_dir,tmp_files[grep("fpr.tsv",tmp_files)]))%>%as.data.frame()%>%reshape2::melt(.)%>%rename("value"="fpr")%>%.[,3]
    act<-data.table::fread(here::here(res_dir,tmp_files[grep("activities.tsv",tmp_files)]))%>%as.data.frame()%>%reshape2::melt(.)%>%rename("value"="score","variable"="sampleid","V1"="Kinase")
    kstar_res<-cbind(act,fpr)
    kstar_res$sampleid<-gsub("data[:]","",kstar_res$sampleid)
    kstar_res$residue<-kin 
    kstar_res
  })%>%do.call(rbind,.)
}

#' ballon_kstar
#'
#' @param kstar_res dataframe;
#' @param by_group bool; wheter to order by group or not 
#' @param by_group_color character vector
#' @param residue character; either Y or ST
#' @param man_colors character
#' @param ysize integer 
#' @param xsize integer
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#' testrds<-"./projects//Sievers_GBMPhopsho_test/rds/ptm_se_id_1_preproc.rds"
#'objectptm<-readRDS(testrds)
#'
# test<-"./projects//Sievers_GBMPhopsho_test/results/id_1/KSTAR_input/DEA/2026-06-09_13-19-12/RESULTS/"
# kstar_res<-load_KSTAR_res_from_dir(kstar_res_dir = test)

#'
## Add group
#'kstar_res<-kstar_res[kstar_res$sampleid%in%objectptm$sample_id,]
#'kstar_res$group<-plyr::mapvalues(kstar_res$sampleid,objectptm$sample_id,objectptm$group)
#'by_group_color <- c(
#'  "EGFR_AMP" = "#CDB4DB",
#'  "EGFR_AMP_Variant" = "#FAD3DD",
#'  "EGFR_WT" = "#A8DADC",
#'  "Variant" = "#FB9A99"
#')
#'
#'ballon_kstar(kstar_res = kstar_res,residue = "Y",by_group = F)
#'ballon_kstar(kstar_res = kstar_res,residue = "Y",by_group = T)
#'ballon_kstar(kstar_res = kstar_res,residue = "Y",by_group = T,by_group_color = by_group_color)
make_kstar_dotplot <- function(kstar_res,
                         by_group = F,
                         by_group_color = NULL,
                         show_significance=T,
                         cluster_samples = FALSE,
                         cluster_distance = "euclidean",
                         cluster_method = "complete",
                         residue = c("Y", "ST"),
                         nkin = 50,
                         man_colors = c("FPR>0.05" = "#D3D3D3",
                                        "FPR>0.05" = "#FFB3BA"),
                         color_gradient = c("darkblue", "whitesmoke", "indianred"),
                         ysize = 10,
                         xsize = 8,
                         legend_size = 10) {
  kstar_res$fpr[kstar_res$fpr == 0] <- 0.0001
  kstar_res$log10fpr <- (-log10(kstar_res$fpr))
  kstar_res_tmp <- kstar_res[kstar_res$residue %in% residue, ]
  
  ##remove redundant kinases
  removekinases <- purrr::map(unique(kstar_res_tmp$Kinase), function(kin) {
    tmp <- kstar_res_tmp[kstar_res_tmp$Kinase == kin, ]
    all(tmp$log10fpr == 0)
  }) %>% unlist()
  names(removekinases) <- unique(kstar_res_tmp$Kinase)
  
  kstar_res_tmp <- kstar_res_tmp[!kstar_res_tmp$Kinase %in% names(removekinases[removekinases]), ]
  if (!is.null(nkin) && is.numeric(nkin) && nkin > 0) {
    top_kinases <- kstar_res_tmp %>%
      dplyr::group_by(Kinase) %>%
      dplyr::summarize(best_fpr = min(fpr, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(best_fpr) %>%
      head(nkin) %>%
      dplyr::pull(Kinase)
    kstar_res_tmp <- kstar_res_tmp[kstar_res_tmp$Kinase %in% top_kinases, ]
  }
  kstar_res_tmp$sampleid <- factor(kstar_res_tmp$sampleid, levels = unique(kstar_res_tmp$sampleid))
  sample_order <- levels(kstar_res_tmp$sampleid)
  kinase_order <- unique(as.character(kstar_res_tmp$Kinase))

  if (isTRUE(cluster_samples)) {
    clustered_order <- get_hclust_order_by_score(
      df = kstar_res_tmp,
      feature_col = "Kinase",
      sample_col = "sampleid",
      score_col = "score",
      cluster_distance = cluster_distance,
      cluster_method = cluster_method,
      order_dimension = "sample"
    )
    if (length(clustered_order) > 0) {
      sample_order <- clustered_order
    }

    clustered_kinase_order <- get_hclust_order_by_score(
      df = kstar_res_tmp,
      feature_col = "Kinase",
      sample_col = "sampleid",
      score_col = "score",
      cluster_distance = cluster_distance,
      cluster_method = cluster_method,
      order_dimension = "feature"
    )
    if (length(clustered_kinase_order) > 0) {
      kinase_order <- clustered_kinase_order
    }
  }
  
  if (by_group) {
    assertthat::assert_that("group" %in% colnames(kstar_res_tmp))
    kstar_group_tmp <- distinct(kstar_res_tmp[, c("sampleid", "group")])
    if (!isTRUE(cluster_samples)) {
      kstar_group_tmp <- kstar_group_tmp[order(kstar_group_tmp$group), ]
    }
    kstar_group_tmp$class <- ""
    kstar_group_tmp$sampleid <- factor(kstar_group_tmp$sampleid, levels = sample_order)
    sample_order <- levels(kstar_group_tmp$sampleid)
    
    group_ballons <- ggplot2::ggplot(kstar_group_tmp, ggplot2::aes(x = sampleid, y = class, fill = group)) +
      ggplot2::geom_tile() +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "right",
        plot.background = element_rect(fill = "white", colour = NA)
      )
    if (!is.null(by_group_color)) {
      group_ballons <- group_ballons + scale_fill_manual(values = by_group_color)
    }
  }
  
  ### Adjust Size
  if (length(unique(kstar_res_tmp$Kinase)) > 80) {
    ysize <- 8
  }
  if (length(unique(kstar_res_tmp$Kinase)) > 120) {
    ysize <- 5
    xsize <- 7
  }
  
  kstar_res_tmp$Significant <- ifelse(kstar_res_tmp$fpr < 0.05, "FPR<0.05", "FPR>0.05")
  levels(kstar_res_tmp$sampleid) <- sample_order
  kstar_res_tmp$Kinase <- factor(kstar_res_tmp$Kinase, levels = kinase_order)

  # Dendrogram rendering disabled by request.
  if(show_significance){
    kstar_ballons <- ggplot2::ggplot(kstar_res_tmp, ggplot2::aes(x = sampleid, y = Kinase, color = Significant, size = log10fpr)) +
      ggplot2::geom_point() +
      cowplot::theme_cowplot() +
      ggplot2::theme(axis.line = element_blank()) +
      ggplot2::theme(axis.text.x = element_text(size = xsize, angle = 90, vjust = 0.5, hjust = 1)) +
      ggplot2::theme(axis.text.y = element_text(size = ysize)) +
      ggplot2::theme(legend.text = element_text(size = legend_size),
            legend.title = element_text(size = legend_size)) +
      ggplot2::ylab("") +
      ggplot2::theme(axis.ticks = element_blank()) +
      ggplot2::labs(x = "Sample", y = "Kinase") +
      ggplot2::scale_color_manual(values = c(
        "FPR>0.05" = "#D3D3D3",
        "FPR<0.05" = "#FFB3BA"
      )) +
      ggplot2::scale_size_continuous(range = c(1, 7), name = "-log10(FPR)")
  }else{
    kstar_ballons <- ggplot2::ggplot(kstar_res_tmp, ggplot2::aes(x = sampleid, y = Kinase, color = score, size = log10fpr)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_viridis_c(name = 'Score') +
      cowplot::theme_cowplot() +
      ggplot2::theme(axis.line = element_blank()) +
      ggplot2::theme(axis.text.x = element_text(size = xsize, angle = 90, vjust = 0.5, hjust = 1)) +
      ggplot2::theme(axis.text.y = element_text(size = ysize)) +
      ggplot2::theme(legend.text = element_text(size = legend_size),
            legend.title = element_text(size = legend_size)) +
      ggplot2::ylab("") +
      ggplot2::theme(axis.ticks = element_blank()) +
      ggplot2::labs(x = "Sample", y = "Kinase") +
      ggplot2::scale_size_continuous(range = c(1, 7), name = "-log10(FPR)")
  }
  if ((by_group)) {
    output <- cowplot::plot_grid(group_ballons, kstar_ballons, ncol = 1, align = "v", rel_heights = c(0.45, 6.55))
  } else{
    output <- kstar_ballons
  }
  return(output)
}

#' make_kstar_barplot
#'
#' @param kstar_res data.frame; KSTAR result table
#' @param residue character vector; choose from "Y" and "ST"
#' @param show_significance logical; if TRUE plot -log10(FPR), else plot score
#' @param nkin integer; number of kinases to include
#' @param color_gradient character vector of length 3; low/mid/high colors
#'
#' @returns ggplot object
#' @export
#'
#' @examples NULL
#' # make_kstar_barplot(kstar_res = kstar_res, residue = "Y", nkin = 30)
make_kstar_barplot <- function(kstar_res,
                               residue = c("Y", "ST"),
                               show_significance = FALSE,
                               nkin = 30,
                               color_gradient = c("darkblue", "whitesmoke", "indianred")) {
  kstar_res$fpr[kstar_res$fpr == 0] <- 0.0001
  kstar_res_tmp <- kstar_res[kstar_res$residue %in% residue, ]

  top_kinases <- kstar_res_tmp %>%
    dplyr::group_by(Kinase) %>%
    dplyr::summarize(best_fpr = min(fpr, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(best_fpr) %>%
    head(nkin) %>%
    dplyr::pull(Kinase)

  kstar_plot_df <- kstar_res_tmp %>%
    dplyr::filter(Kinase %in% top_kinases)

  if (show_significance) {
    p <- ggplot2::ggplot(kstar_plot_df, ggplot2::aes(x = reorder(Kinase, fpr), y = -log10(fpr), fill = -log10(fpr))) +
      ggplot2::geom_col(alpha = 0.9) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~sampleid, scales = "free_x") +
      ggplot2::scale_fill_gradient2(low = color_gradient[1], mid = color_gradient[2], high = color_gradient[3]) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Kinase", y = "-log10(FPR)")
  } else {
    p <- ggplot2::ggplot(kstar_plot_df, ggplot2::aes(x = reorder(Kinase, score), y = score, fill = score)) +
      ggplot2::geom_col(alpha = 0.9) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~sampleid, scales = "free_x") +
      ggplot2::scale_fill_gradient2(low = color_gradient[1], mid = color_gradient[2], high = color_gradient[3]) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Kinase", y = "KSTAR Score")
  }

  return(p)
}

####### Decoupler Kinase Prediction ######

#' load_kinase_network_from_data
#'
#' @param database character; one of "protmapper", "phosphonetworks", "phosphosite", "signor", "kea" 
#' @param db_version character; Omnipath version to load
#' @param db_path character; path to the directory where the database rds files are stored
#'
#' @returns dataframe with columns p_site, enzyme_genesymbol, mor, likelihood, id
#' @export
#'
#' @examples NULL
load_kinase_network_from_data <- function(database = c("protmapper",
                                                       "phosphonetworks",
                                                       "phosphosite",
                                                       "signor",
                                                       "kea"),
                                          db_version = op_version,
                                          db_path = here::here("data", "ptm_databases")) {
  rds_files <- list.files(db_path, pattern = database, full.names = T) %>% .[grep(".rds", .)]
  rds_files <- rds_files[grep(db_version, rds_files)]
  assertthat::assert_that(length(rds_files) == 1)
  omnipath_ptm <- readRDS(rds_files)
  omnipath_ptm$id <- paste(omnipath_ptm$p_site, omnipath_ptm$enzyme_genesymbol, sep =
                             "")
  omnipath_ptm <- omnipath_ptm[!duplicated(omnipath_ptm$id), ]
  omnipath_ptm <- omnipath_ptm[, -5]
  return(omnipath_ptm)
}

### Urgent problem DIA and DDA limma out put requires the right format 
# rename KSN to fit decoupler format
#' dc_kin_activities
#'
#' @param phospho_differential_analysis 
#' @param dc_ptm 
#' @param net_idenitifiers 
#' @param inference_method 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#' # Differential expression result
# tt<-"/mnt/NAS4/user_data/np-dennis/projects/Proteomics/ProteoLab/projects/Sievers_GBMPhopsho_test/results/id_1/DEA/2026-06-09_12-55-32/toptable.xlsx"
# phospho_differential_analysis<-read.xlsx(tt)%>%as.data.frame()
# n_kin<-20
# dc_ptm<-load_kinase_network_from_data(database = "protmapper",db_version = op_version,db_path = here::here("data","ptm_databases"))
# net_idenitifiers<-"genes"
# inference_method<-"viper"
# dc_kin_act<-dc_kin_activities(
#   phospho_differential_analysis = phospho_differential_analysis,
#  dc_ptm = dc_ptm,
#   net_idenitifiers = net_idenitifiers,
#   inference_method = inference_method
# )
#' # Single-sample inference
# test<-dc_kin_activities(
#     phospho_differential_analysis = tt,
#     dc_ptm = dc_ptm,
#     net_idenitifiers = net_idenitifiers,
#     inference_method = inference_method
#   )
dc_kin_activities<-function(
    phospho_differential_analysis,
    dc_ptm,
    net_idenitifiers="genes",
    inference_method=c("viper","wmean","fgsea")
  ){
  
  print(inference_method)
  
  names(dc_ptm)[c(1,2)] <- c("target","tf")
  net<-dc_ptm
  
  if(class(phospho_differential_analysis)=="SummarizedExperiment"){
    dc_mat<-as.matrix(SummarizedExperiment::assay(phospho_differential_analysis))
  }else{
    # Decoupler target = Genename_Residue
    phospho_differential_analysis <- phospho_differential_analysis[!duplicated(phospho_differential_analysis[[net_idenitifiers]]),]
    rownames(phospho_differential_analysis)<-phospho_differential_analysis[[net_idenitifiers]]
    phospho_differential_analysis$targets<-NULL;phospho_differential_analysis$target<-NULL
    contests <- grep("(_t$|^t$)", colnames(phospho_differential_analysis), value = TRUE)
    assertthat::assert_that(length(contests) > 0)
    
    dc_df<-phospho_differential_analysis[, contests, drop = FALSE]
    dc_mat<-as.matrix(dc_df)  
  }
  
  if(inference_method=="viper"){
    message("use ",inference_method)
    contrast_acts <- decoupleR::run_viper(mat=dc_mat,
                                          network = net,
                                          .target = "target",
                                          .source = "tf",
                                          minsize = 5)
    # Filter norm_wmean
    f_contrast_acts <- contrast_acts %>%
      dplyr::filter(statistic == 'viper') %>%
      dplyr::mutate(rnk = NA)
  }
  
  if(inference_method=="wmean"){
    message("use ",inference_method)
    contrast_acts <- decoupleR::run_wmean(mat=dc_mat,
                                          network = net,
                                          .target = "target",
                                          .source = "tf",
                                          times = 1000,
                                          minsize = 5,
                                          seed = 2905)
    
    # Filter norm_wmean
    f_contrast_acts <- contrast_acts %>%
      dplyr::filter(statistic == 'norm_wmean') %>%
      dplyr::mutate(rnk = NA)
  }
  if(inference_method=="fgsea"){
    message("use ",inference_method)
    contrast_acts <- decoupleR::run_fgsea(mat=dc_mat,
                                          network = net,
                                          .target = "target",
                                          .source = "tf",seed = 2905,
                                          minsize = 5)
    contrast_acts$statistic
    # Filter norm_wmean
    f_contrast_acts <- contrast_acts %>%
      dplyr::filter(statistic == 'norm_fgsea') %>%
      dplyr::mutate(rnk = NA)
  }
  
  # Filter top TFs in both signs
  msk <- f_contrast_acts$score > 0
  f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
  f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
  as.data.frame(f_contrast_acts)
}

#' plot_ballon_decoupler
#'
#' @param contrast_acts dataframe output of dc_kin_activities function
#' @param n_kin integer; number of kinases to plot
#' @param color_gradient character vector of length 3; colors for the gradient, low, mid and high
#'
#' @returns ggplot object
#' @export
#'
#' @examples NULL
#test<-dc_kin_activities(
# phospho_differential_analysis = phospho_differential_analysis,
#dc_ptm = dc_ptm,
#net_idenitifiers = net_idenitifiers,
#inference_method = inference_method
#)
# make_dc_dotplot(
#  contrast_acts = test,
#  n_kin = 124,show_significance = F
# )
make_dc_dotplot <- function(contrast_acts,
                                  n_kin = 50,
                                  by_group = (group != "none"),
                                  group = "none",
                                  by_group_color = NULL,
                                  show_significance=T,
                                  cluster_samples = FALSE,
                                  cluster_distance = "euclidean",
                                  cluster_method = "complete",
                                  xsize = 8,
                                  ysize = 10,
                                  legend_size = 10
                                  ) {
  kins_deg <- contrast_acts %>%
    dplyr::group_by(source) %>%
    dplyr::summarize(best_rnk = min(rnk, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(best_rnk) %>%
    head(n_kin) %>%
    dplyr::pull(source)
  
  f_contrast_acts_flt <- contrast_acts %>%
    dplyr::filter(source %in% kins_deg)

  condition_order <- unique(f_contrast_acts_flt$condition)
  kinase_order <- unique(as.character(f_contrast_acts_flt$source))
  if (isTRUE(cluster_samples)) {
    clustered_order <- get_hclust_order_by_score(
      df = f_contrast_acts_flt,
      feature_col = "source",
      sample_col = "condition",
      score_col = "score",
      cluster_distance = cluster_distance,
      cluster_method = cluster_method,
      order_dimension = "sample"
    )
    if (length(clustered_order) > 0) {
      condition_order <- clustered_order
    }

    clustered_kinase_order <- get_hclust_order_by_score(
      df = f_contrast_acts_flt,
      feature_col = "source",
      sample_col = "condition",
      score_col = "score",
      cluster_distance = cluster_distance,
      cluster_method = cluster_method,
      order_dimension = "feature"
    )
    if (length(clustered_kinase_order) > 0) {
      kinase_order <- clustered_kinase_order
    }
  }

  if (by_group) {
    assertthat::assert_that("group" %in% colnames(f_contrast_acts_flt))
    group_tmp <- dplyr::distinct(f_contrast_acts_flt[, c("condition", "group")])
    if (!isTRUE(cluster_samples)) {
      group_tmp <- group_tmp[order(group_tmp$group), ]
    }
    group_tmp$class <- ""
    group_tmp$condition <- factor(group_tmp$condition, levels = condition_order)
    condition_order <- levels(group_tmp$condition)

    group_ballons <- ggplot2::ggplot(group_tmp, ggplot2::aes(x = condition, y = class, fill = group)) +
      ggplot2::geom_tile() +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "right",
        plot.background = ggplot2::element_rect(fill = "white", colour = NA)
      )
    if (!is.null(by_group_color)) {
      group_ballons <- group_ballons + ggplot2::scale_fill_manual(values = by_group_color)
    }
  }

  f_contrast_acts_flt$condition <- factor(f_contrast_acts_flt$condition, levels = condition_order)
  f_contrast_acts_flt$source <- factor(f_contrast_acts_flt$source, levels = kinase_order)
  f_contrast_acts_flt$log10_p_value<-(-log10(f_contrast_acts_flt$p_value))

  # Dendrogram rendering disabled by request.

  if (show_significance) {
    f_contrast_acts_flt$Significant <- ifelse(f_contrast_acts_flt$p_value < 0.05,
                                              "P<0.05",
                                              "P>0.05")
    
    dc_ballons <- ggplot2::ggplot(f_contrast_acts_flt, ggplot2::aes(x = condition, y = source, color = Significant, size = log10_p_value)) +
      ggplot2::geom_point() + 
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0, 0), add = c(0.1, 0.1))) +
      cowplot::theme_cowplot() +
      ggplot2::theme(axis.line  = ggplot2::element_blank()) +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey88", linewidth = 0.35)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = xsize, angle = 90, vjust = 0.5, hjust=1)) +
      ggplot2::theme(axis.text.y =ggplot2::element_text(size = ysize)) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = legend_size),
            legend.title = ggplot2::element_text(size = legend_size)) +
      ggplot2::ylab('') +
      ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
      ggplot2::labs(x = "Condition", y = "Kinase") +
      ggplot2::scale_color_manual(values = c(
        "P>0.05" = "#D3D3D3",
        "P<0.05" = "#FFB3BA"
      )) +
      ggplot2::scale_size_continuous(range = c(1, 9), name = "-log10(P)")
  } else {
    
    dc_ballons <- ggplot2::ggplot(f_contrast_acts_flt,ggplot2::aes(x = condition, y = source, color = score, size = log10_p_value)) + 
      ggplot2::geom_point() + 
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0, 0), add = c(0.1, 0.1))) +
      ggplot2::scale_color_viridis_c(name = 'Score') + 
      cowplot::theme_cowplot() + 
      ggplot2::theme(axis.line  = ggplot2::element_blank()) +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey88", linewidth = 0.35)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = xsize, angle = 90, vjust = 0.5, hjust=1)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = ysize)) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = legend_size),
            legend.title = ggplot2::element_text(size = legend_size)) +
      ggplot2::ylab('') +
      ggplot2::theme(axis.ticks = ggplot2::element_blank()) +
      ggplot2::scale_size_continuous(range = c(1, 9), name = "-log10(P)")
    dc_ballons
  }
  if (by_group) {
    output <- cowplot::plot_grid(group_ballons, dc_ballons, ncol = 1, align = "v", rel_heights = c(0.45, 6.55))
  } else {
    output <- dc_ballons
  }
  output
}

#' plot_barplot_decoupler
#'
#' @param contrast_acts dataframe output of dc_kin_activities function
#' @param n_kin integer; number of kinases to plot
#' @param color_gradient character vector of length 3; colors for the gradient, low, mid and high
#' @param facet_by_condition logical; if TRUE, create faceted plots by condition
#'
#' @returns ggplot object
#' @export
#'
#' @examples NULL
# test<-dc_kin_activities(
#  phospho_differential_analysis = phospho_differential_analysis,
#  dc_ptm = dc_ptm,
#  net_idenitifiers = net_idenitifiers,
#  inference_method = inference_method
# )
# plot_barplot_decoupler(
#  contrast_acts = test,
#  n_kin = 50,
# )
make_dc_barplot <- 
  function(contrast_acts,
                                   n_kin = 50,
                                   color_gradient = c("darkblue", "whitesmoke", "indianred"),
                                   facet_by_condition = TRUE) {
  
  # Calculate log10 adjusted p-value if not present
  if (!"log10_p_value" %in% colnames(contrast_acts)) {
    contrast_acts$log10_p_value <- -log10(contrast_acts$p_value)
  }
  
  # Filter top kinases by rank
  kins_deg <- contrast_acts %>%
    dplyr::arrange(rnk) %>%
    head(n_kin) %>%
    dplyr::pull(source) %>%
    unique(.)
  
  f_contrast_acts_flt <- contrast_acts %>%
    dplyr::filter(source %in% kins_deg) 
  f_contrast_acts_flt$source<-
    factor(f_contrast_acts_flt$source,
           levels = f_contrast_acts_flt$source[order(f_contrast_acts_flt$score,f_contrast_acts_flt$rnk)])
  
  # Create significance column
  f_contrast_acts_flt$significant <- ifelse(f_contrast_acts_flt$p_value < 0.05, 
                                             "p < 0.05", "p >= 0.05")
  # Create barplot
  if (facet_by_condition && dplyr::n_distinct(f_contrast_acts_flt$condition) > 1) {
    barplot <- ggplot2::ggplot(f_contrast_acts_flt, 
                               ggplot2::aes(x = source, y = score, 
                          fill = score, 
                          alpha = log10_p_value)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::facet_wrap(~condition, scales = "free_x") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient2(low = color_gradient[1], high = color_gradient[3],
                           mid = color_gradient[2], midpoint = 0,
                          name = "Kinase activity score") +
      ggplot2::scale_alpha_continuous(range(0,range(f_contrast_acts_flt$log10_p_value)[2]),
                           name = "Log10 adjusted P-value") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        axis.text.y = ggplot2::element_text(hjust = 1, size = 9),
        axis.text.x = ggplot2::element_text(size = 9),
        axis.title.x = ggplot2::element_text(size = 10),
        axis.title.y = ggplot2::element_blank(),
        legend.position = "right") +
      ggplot2::labs(y = "Kinase Activity Score")
    
  } else {
    barplot <- ggplot2::ggplot(f_contrast_acts_flt, 
                               ggplot2::aes(x = source, y = score, 
                          fill = score, 
                          alpha = log10_p_value)) +
      ggplot2::geom_bar(stat = "identity", alpha = 0.8) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient2(low = color_gradient[1], high = color_gradient[3],
                           mid = color_gradient[2], midpoint = 0,
                           name = "Kinase activity score") +
      ggplot2::scale_alpha_continuous(range(0,range(f_contrast_acts_flt$log10_p_value)[2]),
                             name = "Log10 adjusted P-value") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        axis.text.y = ggplot2::element_text(hjust = 1, size = 10),
        axis.text.x = ggplot2::element_text(size = 10),
        axis.title.x = ggplot2::element_text(size = 11),
        axis.title.y = ggplot2::element_blank(),
        legend.position = "right"
      ) +
      ggplot2::labs(y = "Kinase Activity Score")
    barplot
  }
  return(barplot)
}
