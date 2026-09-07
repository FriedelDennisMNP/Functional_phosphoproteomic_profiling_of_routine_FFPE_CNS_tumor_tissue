# prot_import_msdata.R
# Author: Dennis Friedel
# Date: 27.03-2024
# Modification: 29.04-2024
# Description: Utils for Project oriented working
# Detail:
# Turn the output from identification & quantification engines/software Maxquant & DIANN into a summarized experiment
# object
#

"%>%" <- magrittr::"%>%"

# ===================================#
  #### FASTA import ####
# ===================================#

#' read_fasta
#'
#' @author Dennis Friedel
#'
#' @param fasta_path Character; path to a FASTA file (.fasta or .fa).
#'
#' @return Named list where each element is a single character string
#'   containing the concatenated amino-acid sequence for one protein entry.
#'   Names are the full header lines (excluding the leading ">").
#'
#' @details
#' Uses \code{data.table::fread} for fast file reading, then splits the raw
#' lines into entries on every header line (lines starting with ">").
#' Sequence lines between two consecutive headers are concatenated into one
#' string. Empty lines are silently ignored. The list can be subset by any
#' string present in the header, e.g. a UniProt accession.
#'
#' @export
#'
#' @examples
#' # fasta <- read_fasta("path/to/sequences.fasta")
#' # length(fasta)              # number of protein entries
#' # nchar(fasta[["P12345"]])   # length of a specific sequence
read_fasta <- function(fasta_path) {
  if (!is.character(fasta_path) || length(fasta_path) != 1) {
    stop("fasta_path must be a single character string.")
  }
  if (!file.exists(fasta_path)) {
    stop(paste0("File not found: ", fasta_path))
  }
  
  # Read every line as a raw character vector via data.table::fread
  raw <- data.table::fread(
    fasta_path,
    header    = FALSE,
    sep       = "\n",      # treat each line as one field
    colClasses = "character",
    blank.lines.skip = TRUE,
    quote     = ""
  )[[1L]]                  # pull the single column as a character vector
  
  # Locate header lines (start with ">")
  is_header <- startsWith(raw, ">")
  header_idx <- which(is_header)
  
  if (length(header_idx) == 0L) {
    stop("No FASTA header lines ('>') found in the file.")
  }
  
  # For each header, collect all sequence lines until the next header
  n_entries <- length(header_idx)
  sequences <- vector("list", n_entries)
  headers   <- substring(raw[header_idx], 2L)   # strip leading ">"
  
  for (i in seq_len(n_entries)) {
    seq_start <- header_idx[i] + 1L
    seq_end   <- if (i < n_entries) header_idx[i + 1L] - 1L else length(raw)
    sequences[[i]] <- paste(raw[seq_start:seq_end], collapse = "")
  }
  
  names(sequences) <- headers
  sequences
}

# ===================================#
#### Generic import functions  ####
# ===================================#

#' read_ms_file
#'
#' @author Dennis Friedel
#'
#' @param ms_resul_path; Character; path to ms result file.
#' @param clean_column_names logical; TRUE/FALSE if column names in data.frame should follow the same format. here snake_case
#'
#' @return data.frame
#'
#' @description
#' This functions is a generic reading functions which can be applied to all types of tables. It further capplies
#' As input it requires the path of a desired MS-File from MaxQuant, DIA-NN, Spectronaut as input and
#' loads it's content as data.frame into the R environment. Furthermore this function implements the clean column_names function
#' from janitor to provide commonly character format in the data.frame columns.
#'
#' @details
#' fread from data.table  is used to speed up loading of the tabe delimeted file into R.
#'
#' @export
#'
#' @examples NULL
# #mq_txt_path<-"path/to/maxquant/txt/folder"
# #mq_table<-"summary.txt"
# #result_table <-
# # read_ms_file(file.path(mq_txt_path, "/", mq_table))
read_ms_file <- function(mq_file_path, clean_column_names = TRUE, parquet_col_select = NULL) {
  #If file ends with parquet is arrow instead
  if(endsWith(mq_file_path,"parquet")){
    parquet_select <- parquet_col_select
    if (!is.null(parquet_col_select)) {
      parquet_schema_cols <- names(arrow::open_dataset(mq_file_path, format = "parquet"))
      parquet_schema_cols_clean <- janitor::make_clean_names(parquet_schema_cols)
      keep_idx <- which(parquet_schema_cols_clean %in% parquet_col_select)

      if (length(keep_idx) > 0) {
        parquet_select <- parquet_schema_cols[keep_idx]
      } else {
        parquet_select <- NULL
      }
    }

    ms_file <- arrow::read_parquet(mq_file_path, col_select = parquet_select) %>%
      as.data.frame()
    ms_file
  }else{
    ms_file <- data.table::fread(mq_file_path, sep = "\t") %>%
      as.data.frame()
  }
  if (clean_column_names) {
    return(janitor::clean_names(ms_file))
  }
  return(ms_file)
}

#' load_ms_results
#' @author Dennis Friedel
#'
#' @param ms_result_dir character; Directory containing result tables from MaxQuant, DIANN, Spectronaut or any other pre-analysis software.
#' @param file_names character; vector with filenames of the corresponding software.
#'
#' @return List object;
#'
#' @description
#' Reads tables listed in tables  pre-analysis software as named list object into R working space.
#' The user will be infromed about missing files in the directory.
#' @details
#' There is already a pre-selection of files which can be selected via this function.
#' However the user can also customize file names further to load also other files.
#' **MaxQuant files:**
#' "summary.txt","peptides.txt","evidence.txt","proteinGroups.txt","msScans.txt","msmsScans.txt","modificationSpecificPeptides.txt",
#' "Phospho (STY)Sites.txt","parameters.txt","#runningTimes.txt","mqpar.xml".
#'
#' **DIA-NN files:**
#' "result.stats.tsv","result.tsv","result.parquet",report.parquet
#'
#' **Spectronaut files:**
#' "result_spectronaut.tsv"
#'
#' @export
#'
#' @examples NULL
#' # MaxQuant
#' # mq_res<-load_ms_results(ms_result_dir = "/pathtoMaxQuantInput/")
#' # DIA-NN
#' # diann_res<-load_ms_results(ms_result_dir = "/pathtoDIANNInput/")
#' # Spectronaut
#' # spectro_res<-load_ms_results(ms_result_dir = "pathtoSpectronaut")
load_ms_results <-
  function(ms_result_dir,
           file_names = c(
             "summary.txt",
             "peptides.txt",
             "evidence.txt",
             "proteinGroups.txt",
             "Phospho (STY)Sites.txt",
             "result.tsv",
             "report.tsv",
             "result.stats.tsv",
             "report.stats.tsv",
             "result.parquet",
             "report.parquet",
             "result_spectronaut.tsv"
           )) {
    assertthat::assert_that(dir.exists(ms_result_dir))
    assertthat::assert_that(any(file_names%in%list.files(ms_result_dir)))

    available_files <- list.files(ms_result_dir)

    # New DIA-NN runs can contain both files; prefer result.parquet to avoid
    # loading the larger report.parquet twice.
    if ("result.parquet" %in% available_files && "report.parquet" %in% available_files) {
      file_names <- setdiff(file_names, "report.parquet")
    }

    # Columns used by downstream DIA-NN import and SummarizedExperiment creation.
    diann_parquet_cols <- c(
      "run",
      "stripped_sequence",
      "modified_sequence",
      "precursor_normalised",
      "protein_ids",
      "protein_names",
      "genes",
      "proteotypic",
      "q_value",
      "protein_q_value",
      "pg_q_value",
      "gg_q_value",
      "pg_max_lfq",
      "pg_normalised",
      "genes_max_lfq",
      "genes_normalised",
      "genes_max_lfq_unique"
    )
    
    ms_results <- purrr::map(file_names, function(ms_table) {
      if (ms_table %in% available_files) {
        parquet_select <- NULL
        if (ms_table %in% c("result.parquet", "report.parquet")) {
          parquet_select <- diann_parquet_cols
        }

        result_table <-
          read_ms_file(
            file.path(ms_result_dir, "/", ms_table),
            parquet_col_select = parquet_select
          )
        # Summary.txt from MQ comes with an addtional row which we remove here for convenience
        if (ms_table == "summary.txt") {
          result_table <- head(result_table,-1)
          return(result_table)
        } else{
          return(result_table)
        }
      } else {
        warning(paste0(ms_table, " not found, some plots will be missing"))
        print(paste0(ms_table, " not found, some plots will be missing"))
        return(NULL)
      }
    })
    names(ms_results) <- janitor::make_clean_names(file_names)
    ms_results <- plyr::compact(ms_results)
    return(ms_results)
  }

#' check_coldata
#' @author Dennis Friedel
#' @param anno dataframe; dataframe with mandatory columns sample_id, replicate, group (and db_id)
#'
#' @description
#' generic function which controls if the provided annotation dataframe meets the required conditions plus turining
#' characters into common snake_format
#'
#' @return dataframe with snake case column names;
#' @export
#'
#' @examples NULL
check_coldata <- function(anno) {
  anno <- janitor::clean_names(anno)
  if (!("sample_id" %in% colnames(anno))) {
    stop(
      "Sample sheet requires a column name named 'sample_id' which contains the samples in proteinGroups"
    )
  }
  if (!c("replicate" %in% colnames(anno))) {
    stop(
      "Sample sheet requires a column name named 'replicate' indicating if samples are unqiue or replicates"
    )
  }
  if (!c("group" %in% colnames(anno))) {
    stop(
      "Sample sheet requires a column name named 'group' which indicates the primary group you want to investigate"
    )
  }
  if (any(duplicated(colnames(anno)))) {
    stop("anno contains duplicated colnames")
  }
  if (any(duplicated(anno[["sample_id"]]))) {
    stop("sample_id contains duplicated entries")
  }
  anno$sample_id <- janitor::make_clean_names(anno$sample_id)%>%
    gsub("^x","",.)
  rownames(anno)<-anno$sample_id
  anno$group <- janitor::make_clean_names(anno$group, allow_dupes = T)
  anno$replicate <- janitor::make_clean_names(anno$replicate)
  anno
}

#' match_coldata_assay
#'
#' @author Dennis Friedel
#'
#' @param assay_matrix matrix; S3 matrix class object containing sample names in columns and feature names in rows
#' @param anno data.frame; S3 data.frame class dataframe with mandatory columns sample_id, replicate, group (and db_id), while smapl_id contains the sample names
#' present in the colnames of assay_matrix
#'
#' @return list containing coldata and quantmatrix filtered and ordered by samples from both objects.
#'
#' @details
#' 'match_coldata_assay' checks for overlap and order of characters listed in anno column sample_id with the column data. assay_matrix and returns a list containing
#' the quant matrix and coldata matched and ordered according to the samples
#'
#' @export
#'
#' @examples NULL
match_coldata_assay <- function(assay_matrix, anno) {
  anno <- check_coldata(anno)
  samples <- colnames(assay_matrix)
  inter_samples <- intersect(anno$sample_id, samples)
  if (length(inter_samples) == 0) {
    stop(
      "None of the sample_id in anno match the samples in file :",
      "\ncheck column names in 'file' beginning with the intensity as prefix"
    )
  }
  assay_matrix<-assay_matrix[,inter_samples]
  anno_flt <- anno[anno$sample_id %in% colnames(assay_matrix),]
  matched <- match(anno_flt$sample_id, colnames(assay_matrix))
  list("colData"=anno_flt[matched,],"assay"=assay_matrix)
}

#' make_rowdata_ms_se
#'
#' @author Dennis Friedel
#'
#' @param df_ms data.frame; loaded evidence, spectronaut or diann report
#' @param loaded_matrix dataframe; peptide/protein/PTM centric abudnances
#' @param sample_header character; column-name in df_ms containing sample id
#' @param id_header character; column-name in df_ms containing protein ids
#' @param name_header character; column-name in df_ms containing the rownames of loaded matrix
#' @param gene_header character; column-name in df_ms containing the name of the peptide/protein/PTM coding genes
#' @param protein_name_header character; column-name in df_ms containing the full written name of the protein
#'
#' @return NULL
#' @export
#'
#' @examples NULL
make_rowdata_ms_se <- function(df_ms,
                               loaded_matrix,
                               sample_header = "sample_id",
                               id_header = "id",
                               name_header = "id",
                               gene_header = NA,
                               protein_name_header = NA) {
  
  
  if (sum(df_ms[[name_header]] %in% rownames(loaded_matrix)) == 0)
    stop("name header must match rownames of loaded matrix")
  
  df_ms <-
    df_ms[df_ms[[name_header]] %in% rownames(loaded_matrix),] %>%
    data.table::as.data.table()
  
  if (is.null(df_ms[[gene_header]])) {
    genes <- NA
  } else {
    genes <- df_ms[[gene_header]]
  }
  
  if (is.null(df_ms[[protein_name_header]])) {
    protein_name <- NA
  } else {
    protein_name <- df_ms[[protein_name_header]]
  }
  
  default_row_anno <-
    data.frame(
      "id" = gsub("[;].*", "", df_ms[[id_header]]),
      "name" = df_ms[[name_header]],
      "gene" = df_ms[[gene_header]],
      "protein_name" = protein_name
    ) %>% dplyr::distinct()
  #collapse_rows <- default_row_anno[default_row_anno$row_id%in%"AYHEQLSVAEITNACFEPANQMVK",]
  rownames(default_row_anno) <- default_row_anno[["name"]]
  
  #return rowdata with default columns
  default_row_anno[rownames(loaded_matrix),]
}

#' pivot_aggregate
#' @author Dennis Friedel
#' @param df data.table or data.frame; input report table.
#' @param sample_header Character; column name containing sample identifiers.
#' @param id_header Character; column name containing feature identifiers.
#' @param quantity_header Character; column name containing intensity/quantity values.
#' @description
#' modified function from diann R-package
#' 
#' @return NULL
#' @export
#'
#' @examples NULL
pivot_aggregate <-
  function(df,
           sample_header,
           id_header,
           quantity_header) {
    x <-
      data.table::melt.data.table(
        df,
        id.vars = c(sample_header, id_header),
        measure.vars = c(quantity_header)
      )
    x$value[which(x$value == 0)] <- NA
    piv <-
      as.data.frame(
        data.table::dcast.data.table(
          x,
          as.formula(paste0(id_header, "~", sample_header)),
          value.var = "value",
          fun.aggregate = function(x) {
            max(x, na.rm = TRUE)
          }
        )
      )
    rownames(piv) <- piv[[1]]
    piv[[1]] <- NULL
    piv <- piv[order(rownames(piv)),]
    piv <- as.matrix(piv)
    piv[is.infinite(piv)] <- NA
    piv
  }


#' pivot
#' @author Dennis Friedel
#' @param df data.table or data.frame; input report table.
#' @param sample_header Character; column name containing sample identifiers.
#' @param id_header Character; column name containing feature identifiers.
#' @param quantity_header Character; column name containing intensity/quantity values.
#'
#' @return NULL
#' @export
#' @description
#' modified function from diann R-package
#' 
#'
#' @examples NULL
pivot <- function(df,
                  sample_header,
                  id_header,
                  quantity_header) {
  x <-
    data.table::melt.data.table(
      df,
      id.vars = c(sample_header, id_header),
      measure.vars = c(quantity_header)
    )
  x$value[which(x$value == 0)] <- NA
  piv <-
    as.data.frame(data.table::dcast.data.table(x, as.formula(paste0(
      id_header, "~", sample_header
    )), value.var = "value"))
  rownames(piv) <- piv[[1]]
  piv[[1]] <- NULL
  piv <- piv[order(rownames(piv)),]
  as.matrix(piv)
}

# ====================================#
#### Ptm-integration functions #######
# ====================================#
#' get_residue_from_DIANN
#'
#' @param diann_se SummarizedExperimentObject; The S4 Object containing the 
#' intensities of phosphosites in assay 
#' @param fasta_file_path character; path to FASTA file which will be loaded and 
#' used to identified position of modified resiudes 
#'
#' @returns Dataframe with columns showing context between rownames of Diann_se 
#' and Phosphorylated Residues and their positions
#' @export
#'
#' @examples NULL
#' #get_residue_from_DIANN(diann_se = ptm_se_raw,fasta_rds_path = fasta_file_path)
## diann_se<-readRDS("./projects/66_DIA_Phospho_Jurkat_MCF7_120min/rds/ptm_se_id_1_preproc.rds")
## fasta_file_path<-"./data/uniprotkb_Human_AND_reviewed_true_AND_m_2025_01_16.fasta"
## test<-get_residue_from_DIANN(diann_se = diann_se,fasta_file_path = fasta_file_path)
get_residue_from_DIANN<-function(diann_se,fasta_file_path){
  if (grepl("\\.rds$", fasta_file_path, ignore.case = TRUE)) {
    fasta_file <- readRDS(fasta_file_path)
  } else {
    fasta_file <- read_fasta(fasta_file_path)
    gsub("sp[|]","",names(fasta_file))%>%gsub("[|].*","",.)->names(fasta_file)
  }
  
  identifier<-purrr::map(rownames(diann_se),function(pg){
    amino_acid = gsub("[(].*","",pg)%>%stringr::str_sub(., start= -1)
    modified_sequence = gsub("[(].*[)].*","",pg)
    identified_sequence = gsub("[(].*[)]","",pg)
    
    phos_entry = diann_se@elementMetadata[diann_se@elementMetadata$row_id==pg,c("id","gene")]%>%as.data.frame()
    
    if(length(grep(";",phos_entry$id))>0){
      ids<-strsplit(phos_entry$id,";")%>%unlist()
      if(!any(ids%in% names(fasta_file)))
        return(NA)
      which(ids%in% names(fasta_file))%>%ids[.]%>%.[1]->phos_entry$id
    }
    
    ### Return NA if not available
    if(length(grep(phos_entry$id,names(fasta_file)))==0)
      return(NA)
    
    fastaID = fasta_file[names(fasta_file)%in%phos_entry$id][[1]]
    
    ### get start of identified_sequence in FASTA Sequence
    start_position_in_fasta<-nchar(gsub(paste0(identified_sequence,".*"),"",fastaID))
    
    ### get position of modified aminoacid in fasta 
    amino_acids_till_modified_site<-nchar(modified_sequence)
    
    ### add to start position tzo get gobal poisitin in protein sequence & residue
    position_in_protein<-start_position_in_fasta+amino_acids_till_modified_site
    
    residue<-paste0(phos_entry$gene,"_",amino_acid,position_in_protein)
    
    data.frame("residue"=residue,
               "uniprot_id"=phos_entry$id,
               "protein_length"=nchar(fastaID),
               "gene"=phos_entry$gene,
               "modified_sequence"=modified_sequence,
               "start_modified_sequence_in_fasta"=start_position_in_fasta,
               "start_modified_sequence_in_precursor"=nchar(modified_sequence),
               row.names = pg)
  },.progress = TRUE)%>%do.call(rbind,.)
  #%>%plyr::compact()%>%do.call(rbind,.)
  return(identifier)
}

#' extract_highest_prob_sequence
#'
#' @param ptm_probability 
#'
#' @returns NULL
#' @export
#'
#' @examples NULL
#' ptm_probability<-ptm_probability_test<-"ET(0.008)GRQAGVS(0.992)AEMFAMPR"
#' extract_highest_prob_sequence(ptm_probability_test)
extract_highest_prob_sequence <- function(ptm_probability) {
  matches <- gregexpr("[A-Z]\\(([0-9.]+)\\)", ptm_probability, perl = TRUE)
  parts <- regmatches(ptm_probability, matches)[[1]]
  
  letters <- sub("\\(([0-9.]+)\\)", "", parts)
  probs <- as.numeric(sub(".*\\(([0-9.]+)\\)", "\\1", parts))
  
  best<-which.max(probs)
  position_in_nchar<-matches[[1]][best]
  
  chars<-strsplit(ptm_probability,"")%>%unlist()
  chars[position_in_nchar]<-tolower(chars[position_in_nchar])
  ptm_probability_moded<-paste0(chars,collapse = "")
  
  #Clean up
  ptm_probability_cleaned<-gsub("\\(([0-9.]+)\\)", "", ptm_probability_moded)
  
  data.frame("AminoAcid"=gsub("\\(([0-9.]+)\\)", "", parts[best]),
             "peptide"=ptm_probability_cleaned,
             "sequence_probabilities"=ptm_probability)
}

#' make_ptm_anno - Unify PTM Analysis Input for KSTAR,PTM_GSEA and Decoupler
#'
#' @param se_object SummarizedExperimentObject holding Phopphopeptide Intensisites
#' @param data_type character; one of "auto", "DIA-NN", or "MaxQuant". If "auto", the function will attempt to infer the data type based on the row names and rowData of the SummarizedExperiment.
#' @param fasta_path character; path to FASTA file required for DIA-NN conversion. Only used if data_type is "DIA-NN".
#'
#' @returns dataframe with residue information and gene/protein identifiers for PTM analysis.
#' @export
#'
#' @examples NULL
## mq_ptm<-readRDS("/mnt/NAS4/user_data/np-dennis/projects/AG_Sievers/20260717_Phospho_EGFRvariants_Sievers/projects/20260512_RRS_1296_PCF_Phospho/rds/ptm_se_id_1_preproc.rds")
## diann_ptm<-readRDS("./projects/66_DIA_Phospho_Jurkat_MCF7_120min/rds/ptm_se_id_1_preproc.rds")
## fasta_path<-"./data/uniprotkb_Human_AND_reviewed_true_AND_m_2025_01_16.fasta"
## test_diann<-make_ptm_anno(diann_ptm,fasta_path = fasta_path)
make_ptm_anno <- function(se_object,
                          data_type = c("auto", "DIA-NN", "MaxQuant"),
                          object_type = c("auto", "SummarizedExperiment", "Data.Frame"),
                          fasta_path = NULL) {
  
  data_type <- match.arg(data_type)
  object_type <- match.arg(object_type)
  
  if(object_type == "auto"){
    if(methods::is(se_object, "SummarizedExperiment")){
      object_type<-"SummarizedExperiment"
      dt_cols<-colnames(SummarizedExperiment::rowData(se_object))
    }else if(methods::is(se_object, "data.frame")){
      object_type<-"data.frame"
      dt_cols<-colnames((se_object))
    }else{
      stop("se_object must be a SummarizedExperiment or Dataframe.")
    }
  }
  
  if (data_type == "auto") {
    if ("phospho_sty_probabilities" %in% dt_cols) {
      data_type <- "MaxQuant"
    } else {
      data_type <- "DIA-NN"
    }
  }
  
  if (data_type == "MaxQuant") {
    
    if(object_type=="SummarizedExperiment"){
      row_data <- as.data.frame(SummarizedExperiment::rowData(se_object))  
      required_cols <- c("id", "gene", "position", "phospho_sty_probabilities", "sequence_window")
    }else{
      
      row_data <- se_object
      row_data$id<-row_data$proteins
      row_data$gene<-row_data$gene_names
      required_cols <- c("id", "gene", "position", "phospho_sty_probabilities", "sequence_window")
    }
    
    if (!all(required_cols %in% colnames(row_data))) {
      stop('MaxQuant mode requires rowData columns: "id","gene", "position", "phospho_sty_probabilities","sequence_window"')
    }
    ptm_df <- row_data[, required_cols]
    ptm_df$rowid <- rownames(ptm_df)
    ptm_df$id <- gsub(";.*", "", ptm_df$id)
    
    peptide <- purrr::map(ptm_df$phospho_sty_probabilities, function(x) {
      extract_highest_prob_sequence(x)
    }) %>% do.call(rbind, .)
    
    ptm_df$ptm_peptide <- peptide$peptide
    ptm_df$gene_residue <- paste0(gsub(";.*", "", ptm_df$gene), "_", peptide$AminoAcid, ptm_df$position)
    ptm_df$id_residue   <- paste0(ptm_df$id, "_", peptide$AminoAcid, ptm_df$position)
    ptm_df$Name         <- ptm_df$gene_residue
    
  }
  
  if (data_type == "DIA-NN") {
    if (is.null(fasta_path) || !nzchar(fasta_path)) {
      stop("DIA-NN mode requires fasta_path.")
    }
    ptm_df   <- as.data.frame(SummarizedExperiment::rowData(se_object))
    site_info <- get_residue_from_DIANN(diann_se = se_object, fasta_file_path = fasta_path)
    keep     <- intersect(rownames(se_object), rownames(site_info))
    
    ptm_df<-ptm_df[keep,]
    ptm_df$rowid  <- rownames(ptm_df)
    
    ptm_df$ptm_peptdide <- gsub("[(].*[)]", "_", ptm_df$rowid) %>% strsplit(., "_") %>%
      purrr::map(., function(x) {
        mod_seq <- strsplit(x[1], "") %>% unlist()
        mod_seq[length(mod_seq)] <- tolower(mod_seq[length(mod_seq)])
        paste0(paste0(mod_seq, collapse = ""), x[2])
      }) %>% unlist()
    
    # gene_residue: GENE_AMINOACIDPOSITION  e.g. RFX1_S59  (decoupler / PTM-GSEA input)
    ptm_df$gene_residue <- site_info[keep, "residue"]
    
    # id_residue:   UNIPROT_AMINOACIDPOSITION  e.g. P51523_S59  (KSTAR input)
    ptm_df$id_residue   <- paste0(site_info[keep, "uniprot_id"], "_",
                                  gsub(".*_", "", site_info[keep, "residue"]))
    ptm_df$Name         <- ptm_df$gene_residue
    ptm_df
  }
  
  if (!"Name" %in% colnames(ptm_df)) {
    stop("Failed to derive PTM identifiers in format GENE_RESIDUEPOSITION.")
  }
  
  return(ptm_df)
}

#' .se_ptm_aggregate - Aggregate MS Summarized experiment by Column in elementMetadata 
#'
#' @param se_object SummarizedExperiment; ptm_object obtained
#' @param new_id Character; Column name in elementMetadata to aggregate by
#'
#' @returns SummarizedExperiment with aggregated intensities and rowData
#' @export
#'
#' @examples NULL
# se_object<-ptm_se_anno
# new_id<-"Name"
#'#.se_ptm_aggregate(se_object = ptm_se_anno,new_id = "Name")
.se_ptm_aggregate<-function(se_object,new_id){
  
  assertthat::assert_that(methods::is(se_object, "SummarizedExperiment"))
  assertthat::assert_that(is.character(new_id))
  
  new_id_dup<-se_object@elementMetadata[[new_id]][duplicated(se_object@elementMetadata[[new_id]])]
  
  ## Only keep assay intensities in Summarized Experiment for unique gene_residues
  se_object_unique<-se_object[!se_object@elementMetadata[[new_id]]%in%new_id_dup,]
  
  agg_se<-purrr::map(unique(new_id_dup),function(x){
    
    tmp_se<-se_object[se_object@elementMetadata$Name==x,]
    new_assay<-matrixStats::colMaxs(SummarizedExperiment::assay(tmp_se),na.rm=T)%>%as.matrix()%>%t()
    new_assay[is.infinite(new_assay)] <- NA
    rownames(new_assay)<-x
    
    ## Paste0 elementMetadata of first entry of the duplicated gene_resid
    SummarizedExperiment::SummarizedExperiment(
      assays = list("intensity" = new_assay),
      rowData = tmp_se@elementMetadata[1,],
      colData = tmp_se@colData
    )
  },.progress = T)%>%do.call(rbind,.)%>%rbind(se_object_unique,.)
  rownames(agg_se)<-agg_se@elementMetadata[[new_id]]
  return(agg_se)
}

# ====================================#
#### Generic ms stat functions  #####
# ====================================#

#' summarize_mq_identification
#' @author Dennis Friedel
#' @param ms_data data.frame; imported MS result table.
#'
#' @return data.frame with identification numbers and further details
#' @export
#'
#' @examples NULL
summarize_mq_identification <- function(ms_data, mod_data = FALSE) {
  assertthat::assert_that(is.data.frame(ms_data),
                          is.logical(mod_data))
  
  id_table <- ms_data %>%
    dplyr::select(dplyr::contains(
      c(
        "intensity_",
        "reverse",
        "potential",
        "only_identified_by_site"
      )
    ))
  
  # Remove Intensity from name
  colnames(id_table) <- gsub("intensity_", "", colnames(id_table))
  
  table_summary <- data.frame(
    "feature_identified" = nrow(id_table) - colSums(id_table == 0),
    "missing_values" = colSums(id_table == 0),
    "potential_contaminants" = colSums(id_table[id_table$potential_contaminant %in%
                                                  "+",] > 0),
    "reverse" = colSums(id_table[id_table$reverse %in% "+",] >
                          0),
    # "Reverse" = length(which(protein_table$Reverse == '+')),
    "only_identified_by_site" = colSums(id_table[id_table$only_identified_by_site %in%
                                                   "+",] > 0),
    check.names = FALSE
  ) %>%
    .[!(
      row.names(.) %in%
        c(
          "potential_contaminant",
          "reverse",
          "only_identified_by_site"
        )
    ),]
  
  table_summary$experiment <- rownames(table_summary)
  rownames(table_summary) <- NULL
  table_summary <- table_summary[, c(6, 1, 2, 3, 4, 5)]
  
  # check if match between runs was performed and how much it boosts the id rate
  mbr_counts <- ms_data %>%
    dplyr::select(dplyr::contains("identification_type_")) %>%
    stringr::str_count(., "By matching")
  if (sum(mbr_counts) > 0) {
    table_summary[, "by_mbr"] <- mbr_counts
  }
  
  if (length(grep("lfq_", table_summary[, "experiment"])) > 0) {
    lfq_stats <-
      table_summary[grep("lfq_", table_summary[, "experiment"]),]
    colnames(lfq_stats)[-1] <-
      paste0("lfq_", colnames(lfq_stats))[-1]
    stats_out <-
      cbind(table_summary[grep("lfq_", table_summary[, "experiment"],
                               invert = TRUE),], lfq_stats[-1])
    return(stats_out)
  }
  
  if (mod_data) {
    p1_numbers <-
      table_summary[endsWith(x = table_summary[, "experiment"], suffix = "_1"),]
    colnames(p1_numbers) <- paste0("p1_", colnames(p1_numbers))
    p2_numbers <-
      table_summary[endsWith(x = table_summary[, "experiment"], suffix = "_2"),]
    colnames(p2_numbers) <- paste0("p2_", colnames(p2_numbers))
    p3_numbers <-
      table_summary[endsWith(x = table_summary[, "experiment"], suffix = "_3"),]
    colnames(p3_numbers) <- paste0("p3_", colnames(p3_numbers))
    
    remove_experiments <-
      c(p1_numbers$p1_experiment,
        p2_numbers$p2_experiment,
        p3_numbers$p3_experiment)
    table_summary <-
      table_summary[!(table_summary$experiment %in% remove_experiments),]
    table_summary_all <- table_summary[-c(1, 2, 3),]
    stats_out <-
      do.call("cbind",
              list(table_summary_all, p1_numbers[,-1], p2_numbers[,-1], p3_numbers[,-1]))
    return(stats_out)
  }
  table_summary <-
    dplyr::rename(table_summary, c("sample_id" = "experiment"))
  table_summary
}

#' aggreate_ms_states
#' @author Dennis Friedel
#' @param ms_data data.frame; imported MS result table.
#' @param experiment_column Character; column name identifying experiments/runs.
#' @param feature_column Character; column name identifying features (peptides/proteins).
#' @param in_percent Logical; if TRUE return counts as percentages.
#'
#' @return NULL
#' @export
#'
#' @examples NULL
aggregate_ms_stats <-
  function(ms_data,
           experiment_column = "experiment",
           feature_column = c("charge", "precursor_charge", "missed_cleavages"),
           in_percent = TRUE) {
    assertthat::assert_that(
      is.data.frame(ms_data),
      is.character(feature_column),
      is.character(experiment_column),
      experiment_column %in% colnames(ms_data),
      any(
        feature_column %in% c("charge", "precursor_charge", "missed_cleavages", "score")
      ),
      is.logical(in_percent)
    )
    
    
    feature_dcast <-
      ms_data[c(experiment_column, feature_column)] %>%
      reshape2::dcast(., as.formula(paste0(experiment_column, " ~ ", feature_column)), fill = 0)
    
    if (in_percent) {
      feature_dcast <-
        cbind(feature_dcast[, experiment_column],
              feature_dcast[,-1] / rowSums(feature_dcast[,-1]) *
                100)
    }
    
    colnames(feature_dcast) <-
      c(experiment_column, paste0(feature_column, 1:(ncol(feature_dcast) - 1)))
    feature_dcast <-
      dplyr::rename(feature_dcast, c("sample_id" = experiment_column))
    feature_dcast
  }

#' calc_gravy
#' @author Dennis Friedel
#' @param peptides data.frame; peptide-level result table.
#'
#' @return NULL
#' @export
#'
#' @examples NULL
calc_gravy <- function(peptides) {
  assertthat::assert_that(is.data.frame(peptides))
  df <- peptides %>% dplyr::select(dplyr::contains(c(
    "length",
    "count",
    "sequence",
    "experiment"
  )))
  
  # Kyte and Dolittle scale
  # negativ (hydrophil), über Null (neutral), bis positiv (hydrophob)
  df$gravy <- (
    df$a_count * 1.8 +
      df$r_count * -4.5 +
      df$n_count * -3.5 +
      df$d_count * -3.5 +
      df$c_count * 2.5 +
      df$q_count * -3.5 +
      df$e_count * -3.5 +
      df$g_count * -0.4 +
      df$h_count * -3.2 +
      df$i_count * 4.5 +
      df$l_count * 3.8 +
      df$k_count * -3.9 +
      df$m_count * 1.9 +
      df$f_count * 2.8 +
      df$p_count * -1.6 +
      df$s_count * -0.8 +
      df$t_count * -0.7 +
      df$w_count * -0.9 +
      df$y_count * -1.3 +
      df$v_count * 4.2
  ) / df$length
  
  df <-
    df %>% dplyr::select(dplyr::contains(c("gravy", "experiment")))
  df_out <- reshape2::melt(df, id.vars = "gravy")
  
  # Remove value 0,
  df_out$variable <- gsub("experiment_",  "", df_out$variable)
  df_out <- df_out[!is.na(df_out$value),]
  
  # Repeat rows n numbers of times, being n the frequency (value)
  df_expanded <- df_out[rep(rownames(df_out), df_out$value),]
  df_table <- df_expanded %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      Mean = format(round(mean(gravy), 2), nsmall = 1),
      Max = format(round(max(gravy), 2), nsmall = 1),
      Min = format(round(min(gravy), 2), nsmall = 1),
      Median = format(round(median(gravy), 2), nsmall = 1)
    ) %>%
    as.data.frame(df_table)
  colnames(df_table) <-
    c("experiment", paste0("gravy_", colnames(df_table)[-1]))
  
  df_table <- dplyr::rename(df_table, c("sample_id" = "experiment"))
  df_expanded <-
    dplyr::rename(df_expanded, c("sample_id" = "variable"))
  
  list("gravy_summary" = df_table,
       "gravy_scores" = df_expanded)
}

#' gravy_diann
#' @author Dennis Friedel
#' @param diann_res data.frame; DIA-NN result table.
#'
#' @return NULL
#' @export
#'
#' @examples NULL
gravy_diann <- function(diann_res) {
  # split sequence window to obtain aminoacid count
  aa_counts <- unique(diann_res$stripped_sequence) %>%
    strsplit(., "") %>%
    purrr::map(., function(pepseq) {
      a_count <- table(pepseq)
      names(a_count) <- paste0(tolower(names(a_count)), "_count")
      a_count <- as.data.frame.default(a_count) %>%
        t() %>%
        as.data.frame()
    }) %>%
    do.call(plyr::rbind.fill, .)
  
  # substitute NA with 0
  aa_counts[is.na(aa_counts)] <- 0
  # add sequence name and length
  aa_counts$sequence <- unique(diann_res$stripped_sequence)
  aa_counts$length <- nchar(aa_counts$sequence)
  
  # dcast sample and sequence to get a matrix with sequences as rows
  # and samples as columns with the number of occurences of the sequence in each samlpes
  piv <-
    as.data.frame(reshape2::dcast(diann_res, as.formula(
      paste0("stripped_sequence", "~", "run")
    ))) %>%
    tibble::column_to_rownames(., var = "stripped_sequence")
  colnames(piv) <- paste0("experiment_", colnames(piv))
  peptide_diann <- cbind(piv, aa_counts)
  
  # use the calc gravy function to obtain gravy sumamry and expaneded score
  calc_gravy(peptides = peptide_diann)
}



# =========================#
#### MaxQuant  ######
# =========================#

#' load_mq_res_table
#'
#' @param mq_tab ; dataframe from mq_result_list
#' @param remove_low_quality logical; remove low quality features
#' @param is_pg logical; is dataframe proteingroups
#' @param is_ptm logical; is dataframe a STY file
#' @param verbose logical; print messages
#' 
#' @return NULL
#' @export
#'
#' @examples NULL
load_mq_res_table<-function(mq_tab,
                            remove_low_quality = TRUE,
                            is_pg=FALSE,
                            is_ptm =FALSE,
                            verbose =TRUE){
  assertthat::assert_that(
    is.data.frame(mq_tab),
    is.logical(remove_low_quality),
    is.logical(is_pg)
  )
  
  mq_tab <-
    mq_tab %>%
    as.data.frame() %>%
    janitor::clean_names()
  
  # Remove low quality
  if (remove_low_quality == TRUE) {
    if(verbose)
      message("Filter features identified as contaminant/reverse/only_idenitfied_by_site")
    mq_tab <- remove_low_quality_features(ms_table = mq_tab,is_pg = is_pg)
  }
  
  # Add proteotypic column
  if(!is_pg){
    if(verbose)
      message("Add proteotypic column")
    # Adds new column with logicals that are TRUE if the peptide can be assigned
    # to only one protein and FALSE if it can be assigned to multiple
    mq_tab <- mq_tab %>%
      dplyr::mutate(is_proteotypic = stringr::str_detect(
        string = proteins,
        pattern = ";",
        negate = TRUE
      ))
  }
  mq_tab
}

#' remove_low_quality_features
#'
#' @param ms_table dataframe; either maxquant proteingroups table or evidence/sty
#' @param is_pg locigal; is the table proteingrous or peptide based
#' @param is_ptm locigal; is the table a STY file
#' @param thr_localization_probability numeric; value to use for filtering low quality ptm features
#'
#' @return filtered dataframe
#' @export
#'
#' @examples NULL
remove_low_quality_features <- function(ms_table, is_pg = TRUE,is_ptm=FALSE,thr_localization_probability=0) {
  assertthat::assert_that("potential_contaminant" %in% colnames(ms_table) &
                            ("reverse" %in% colnames(ms_table)))
  if (is_pg) {
    ms_table <-
      ms_table[rowSums(ms_table[, c("potential_contaminant",
                                    "reverse",
                                    "only_identified_by_site")] == "+") == 0,]
  } else{
    ms_table <-
      ms_table[rowSums(ms_table[, c("potential_contaminant",
                                    "reverse")] == "+") == 0,]
    if(is_ptm&thr_localization_probability>0){
      ms_table <-
        ms_table[ms_table[["localization_prob"]]>thr_localization_probability,]  
    }
    
  }
  return(ms_table)
}

#' make_ptm_identifier
#' @author Dennis Friedel
#' @param input_dataframe dataframe; Phsophosite table either output from maxquant or DIA-NN
#' @param uniprot character; name of the column in input_data_frame containing the uniprot ids
#' @param gene_name character; name of the column in input_data_frame containing the gene_names ids
#' @param amino_acid character; name of the column in input_data_frame containing the amino_acid modified (STY)
#' @param ptm_position character; name of the column in input_data_frame containing the position of the modified aminoacid wihtin the protein sequence
#' @param sequence character; name of the column in input_data_frame containing the idenitfied sequence
#' @param type character; the required format. Either PTM for generating a unique row_id for downstream analysis or format for PTM_GSEA analysis
#' 
#'
#' @return vector with with the provided ids in the required fromat for downstream analyis or analyiss by ptm_gsea
#' @export
#'
#' @examples NULL
#' # Example call (requires a valid phospho_groups data.frame):
#' # res <- make_ptm_identifier(
#' #   input_dataframe = phospho_groups,
#' #   uniprot = "protein",
#' #   gene_name = "gene_names",
#' #   amino_acid = "amino_acid",
#' #   ptm_position = "position",
#' #   sequence = "sequence_window",
#' #   type = "ptm_gsea"
#' # )
make_ptm_identifier <- function(input_dataframe,
                                uniprot = "protein",
                                gene_name = "gene_names",
                                amino_acid = "amino_acid",
                                ptm_position = "position",
                                sequence = "sequence_window",
                                type = c("ptm", "ptm_gsea")) {
  assertthat::assert_that(
    is.data.frame(input_dataframe),
    type %in% c("ptm", "ptm_gsea"))
  if (type == "ptm") {
    assertthat::assert_that(
      is.character(uniprot)&uniprot%in%colnames(input_dataframe),
      is.character(gene_name)&gene_name%in%colnames(input_dataframe),
      is.character(amino_acid)&amino_acid%in%colnames(input_dataframe),
      is.character(ptm_position)&ptm_position%in%colnames(input_dataframe),
      is.character(sequence)&sequence%in%colnames(input_dataframe)
    )
    id <- paste0(
      input_dataframe[[uniprot]],
      ";",
      input_dataframe[[gene_name]],
      ";",
      input_dataframe[[amino_acid]],
      input_dataframe[[ptm_position]],
      ";",
      input_dataframe[[sequence]]
    )
  }
  if (type == "ptm_gsea") {
    assertthat::assert_that(
      is.character(uniprot)&uniprot%in%colnames(input_dataframe),
      is.character(gene_name)&gene_name%in%colnames(input_dataframe),
      is.character(amino_acid)&amino_acid%in%colnames(input_dataframe)
    )
    id <- paste0(input_dataframe[[uniprot]],
                 ";",
                 input_dataframe[[ptm_position]],
                 input_dataframe[[amino_acid]],
                 "-p")
  }
  return(id)
}

#' create_maxquant_se_list
#' @author Dennis Friedel
#' @param mq_results List; List object containing maxquant tables.
#' @param anno dataframe; Sample annotation. mandatory columns are sample_id, group, replicate
#' @param remove_low_quality logical;
#' Remove feautres which are marekd with + in columns reverse, contaminant, only identified by site
#' @param ptm_probability_flt numeric; Threshold for localization probabilities of ptm features.
#' @param verbose logical; print progress messages
#' @param db_pg character; Name of the database to use for mapping protein ids to gene symbols. Default is "org.Hs.eg.db" for human.
#' 
#' @return List with summarized experiment objects. Conatining peptide,protein or PTM quantities.
#' @description
#' Wrapper to load intensitiy/lfq intensities with feature information
#' into summarized experiment format. Prot_se for proteingroups,
#' pep_se contains information from evidence.txt, phos_se contains information from phospho-peptides.
#' @export
#'
#' @examples NULL
## test<-create_maxquant_se_list(
##   mq_results = proteoLab::load_ms_results(
##     ms_result_dir = "/mnt/add50/PATHO-PROTEOMICS/bioinformatics/results/DDA/20260709_RRS_1296PCF/MaxQuant2.4.2.0/",
##     file_names = c(
##       "summary.txt",
##       "peptides.txt",
##       "evidence.txt",
##       "proteinGroups.txt",
##       "Phospho (STY)Sites.txt"
##     )
##   ),
##   remove_low_quality = TRUE,
##   ptm_probability_flt = 0.75,
##   verbose = TRUE,
##   db_pg = "org.Hs.eg.db"
## )
### Phospho
## test<-create_maxquant_se_list(
##   mq_results = proteoLab::load_ms_results(
##     ms_result_dir = "/mnt/add50/PATHO-PROTEOMICS/bioinformatics/results/DDA/20251010_phospho_EGFR_EGF_woEGF/MaxQuant2.4.2.0/",
##     file_names = c(
##      "summary.txt",
##       "peptides.txt",
##       "evidence.txt",
##       "proteinGroups.txt",
##       "Phospho (STY)Sites.txt"
##     )
##   ),
##   remove_low_quality = TRUE,
##   ptm_probability_flt = 0.75,
##   verbose = TRUE,
##   db_pg = "org.Hs.eg.db"
## )
##
create_maxquant_se_list <- function(mq_results,
                                    remove_low_quality = TRUE,
                                    ptm_probability_flt = 0.75,
                                    verbose=TRUE,
                                    db_pg="org.Hs.eg.db") {
  assertthat::assert_that(is.list(mq_results),
                          #is.data.frame(anno),
                          any(names((mq_results))%in%c("summary_txt")),
                          any(
                            names((mq_results)) %in% c(
                              "protein_groups_txt",
                              "evidence_txt",
                              "phospho_sty_sites_txt"
                            )
                          ))
  mq_se_list <- list()
  
  # =========#
  ##### Anno
  # =========#
  anno<-mq_results[["summary_txt"]]
  anno$sample_id<-anno$group<-anno$replicate<-anno$experiment
  
  # =========#
  ##### Protein
  # =========#
  if ("protein_groups_txt" %in% names((mq_results))) {
    prot_se <-
      make_se_mq_proteins(
        protein_groups_txt = mq_results$protein_groups_txt,
        anno,
        gene_names = "gene_names",
        protein_id = "protein_i_ds",
        remove_low_quality,
        db = db_pg,
        delim = ";"
      )
    
    mq_se_list[[length(mq_se_list) + 1]] <- prot_se[["pg_se"]]
    names(mq_se_list)[[length(mq_se_list)]] <- "pg_se"
    
    mq_se_list[[length(mq_se_list) + 1]] <- prot_se[["gg_se"]]
    names(mq_se_list)[[length(mq_se_list)]] <- "gg_se"
    
  } else {
    warning("protein_groups_txt not available, prot_se will be missing")
    print("protein_groups_txt not available, prot_se will be missing")
  }
  
  # =========#
  ##### Peptide
  # =========#
  if ("evidence_txt" %in% names(mq_results)) {
    pep_se <-
      make_se_mq_peptides(evidence_txt = mq_results$evidence_txt,
                          anno,
                          remove_low_quality)
    mq_se_list[[length(mq_se_list) + 1]] <- pep_se
    names(mq_se_list)[[length(mq_se_list)]] <- "pep_se"
  } else {
    warning("evidence-txt not found, pep_se will be missing")
    print("evidence_txt not found, pep_se will be missing")
  }
  
  # =========#
  ##### PTM
  # =========#
  if ("phospho_sty_sites_txt" %in% names(mq_results)) {
    phos_se <-
      make_se_mq_phos(
        phospho_sty_sites_txt = mq_results$phospho_sty_sites_txt,
        anno,
        remove_low_quality = remove_low_quality,
        localization_prob_threshold = ptm_probability_flt
      )
    mq_se_list[[length(mq_se_list) + 1]] <- phos_se
    names(mq_se_list)[[length(mq_se_list)]] <- "ptm_se"
  } else {
    warning("phospho_sty_sites_txt not found, ptm_se will be missing")
    print("phospho_sty_sites_txt not found, ptm_se will be missing")
  }
  return(mq_se_list)
}

#' make_se_mq_proteins
#' @author Dennis Friedel
#' @param protein_groups_txt dataframe; proteingroups.txt from maxquant txt folder
#' @param anno dataframe; Annotation for samples in proteinGroups needs columns named SampleID,group and replicate
#' @param protein_id character; Name of the column which contains the protein ID
#' @param gene_names Character vector; gene symbols to use for mapping.
#' @param remove_low_quality logical;
#' Remove feautres which are marekd with + in columns reverse, contaminant, only identified by site
#' @param delim character; Delimter of name and protein id column
#'
#' @return A object of class SummarizedExperiment
#' @export
#'
#' @examples NULL
#' # Data in this example is derived from 61_Jurkat_MCF7_120min
#' # Get SummarizedExperiment with protein intensities
## protein_groups_txt = mq_results$protein_groups_txt
## anno
## gene_names = "gene_names"
## protein_id = "protein_i_ds"
## remove_low_quality
## db = db_pg
## delim = ";"
## intensity_only=T
make_se_mq_proteins <- function(protein_groups_txt,
                                anno,
                                remove_low_quality = TRUE,
                                gene_names = "gene_names",
                                protein_id = "protein_i_ds",
                                delim = ";",
                                intensity_only=TRUE,
                                db="org.Hs.eg.db",
                                verbose = TRUE) {
  
  # =========================================#
  # Load and preproc res table for SummarizedExperimentConversion
  # =========================================#
  assertthat::assert_that(
    is.data.frame(protein_groups_txt),
    is.data.frame(anno),
    length(gene_names) == 1 & is.character(gene_names),
    length(protein_id) == 1 & is.character(protein_id),
    length(delim) == 1
  )
  anno <- check_coldata(anno)
  proteins <-
    load_mq_res_table(mq_tab = protein_groups_txt,
                      remove_low_quality = TRUE,
                      is_pg = TRUE)
  
  if (!gene_names %in% colnames(proteins)) {
    stop("gene_names does not match in column names of file")
  }
  if (!protein_id %in% colnames(proteins)) {
    stop("protein_id does not match in column names of file")
  }
  
  # =========================================#
  # Preprocess
  # =========================================#
  # make row_id from gene_names/proteinid  
  if(verbose)
    message("Making unique names")
  
  double_na <-
    apply(proteins[, c(gene_names, protein_id)], 1, function(x) {
      all(is.na(x))
    })
  if (any(double_na)) {
    stop("NAs in both the 'names' and 'ids' columns")
  }
  
  proteins$id <-
    gsub(paste0(delim, ".*"), "", proteins[, protein_id])
  proteins$row_id <-
    gsub(paste0(delim, ".*"), "", proteins[, gene_names])
  
  proteins$row_id <-
    make.unique(ifelse(
      proteins$row_id == "" |
        is.na(proteins$row_id),
      proteins$id,
      proteins$row_id
    ))
  rownames(proteins)<-proteins$row_id
  
  # =========================================#
  # Assay data 
  # =========================================#
  # Get Assays from PG Table
  sample_name<-anno$sample_id
  sample_quantities <- proteins %>%
    dplyr::select(dplyr::contains(c(sample_name)))
  
  quant_cols<-colnames(sample_quantities)
  for(i in 1:length(sample_name)){
    quant_cols<-gsub(paste0("_", sample_name[i],".*"), "", quant_cols) %>% unique()  
  }
  assay_list<-purrr:::map(quant_cols,function(x){
    res<-as.data.frame(sample_quantities[,grep(paste0("^",x),colnames(sample_quantities))])
    colnames(res)<-gsub(paste0(x,"_"),"",colnames(res)) %>% 
      gsub(paste0("_percent"),"",.)
    if(x%in%c("intensity", "lfq_intensity", "ibaq_intensity")){
      res <- log2(1 + res)
      res[res == 0] <- NA
    }
    res
  })
  names(assay_list)<-quant_cols
  if(intensity_only){
    assay_list<-assay_list[grep("intensity",names(assay_list))]
  }
  
  #========================================#
  # Experiment Coldata - Column Annotation #
  #========================================#
  if(length(assay_list)<1)
    stop("No overlap found between samples in anno and proteingroups")
  experiment_coldata <- match_coldata_assay(assay_matrix = assay_list[[1]],anno = anno)[[1]]
  rownames(experiment_coldata)<-experiment_coldata$sample_id
  
  #========================================#
  # Experiment Rowdata - Feature Annotation#
  #========================================#
  pre_row_data <-
    proteins[,-which(colnames(proteins) %in% colnames(sample_quantities))]
  
  # add ensembl ids
  # message("Add ENTREZID")
  # ensembl_ids<-clusterProfiler::bitr(pre_row_data$protein_i_ds, fromType="UNIPROT", toType="ENTREZID", OrgDb=db)
  # mapped_ensembl_ids<-plyr::mapvalues(gsub(";.*","",pre_row_data$protein_i_ds),ensembl_ids[["UNIPROT"]],ensembl_ids[["ENTREZID"]],warn_missing = FALSE)
  # mapped_ensembl_ids[mapped_ensembl_ids%in%pre_row_data$protein_i_ds]<-NA
  
  row_data <- data.frame(
    "row_id" = pre_row_data$row_id,
    "uniprot" = pre_row_data$id,
    "gene_symbol" = pre_row_data$gene_names,
    #"ensembl_id" = mapped_ensembl_ids,
    pre_row_data
  ) %>% dplyr::distinct()
  row_data <-
    row_data[, c("row_id", "uniprot", "gene_symbol", sort(setdiff(
      colnames(row_data),c("row_id", "uniprot", "gene_symbol")
    )))]
  
  #========================================#
  # Experiment Summarization #
  #========================================#
  if(verbose)
    message("Convert proteins into SummarizedExperiment object")
  
  assertthat::assert_that(identical(rownames(assay_list[[1]]),rownames(row_data)))
  se_gene<-se_prot <-
    SummarizedExperiment::SummarizedExperiment(assays = assay_list,
                                               rowData = row_data,
                                               colData = experiment_coldata)
  rownames(se_prot)<-rownames(se_prot@elementMetadata)<-se_prot@elementMetadata$uniprot
  
  return(list("pg_se"=se_prot,"gg_se"=se_gene))
}

#' make_se_mq_peptides
#' @author Dennis Friedel
#' @param evidence_txt dataframe; evidence.txt from maxquant txt folder
#' @param anno dataframe; Sample annotation. mandatory columns are sample_id, group, replicate
#' @param remove_low_quality logical;
#' Remove feautres which are marked with + in columns reverse, contaminant, only identified by site
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#' #make_se_mq_peptides(file = "./data/evidence_mq.txt", read.xlsx("./data/anno_mq.xlsx"))
# #mq_results <- load_ms_results(ms_result_dir = mq_txt_dir)
# #evidence_txt<-mq_results$evidence_txt
make_se_mq_peptides <- function(evidence_txt = "",
                                anno,
                                remove_low_quality = TRUE,
                                verbose =TRUE) {
  
  # =========================================#
  # Load and preproc res table for SummarizedExperimentConversion
  # =========================================#
  if(verbose)
    message("Load and control data")
  assertthat::assert_that(
    is.data.frame(evidence_txt),
    is.data.frame(anno),
    is.logical(remove_low_quality)
  )
  anno <- check_coldata(anno = anno)
  
  evidence<-load_mq_res_table(mq_tab = evidence_txt,
                              remove_low_quality = TRUE,
                              is_pg = FALSE,
                              is_ptm = FALSE,
                              verbose = verbose)
  evidence<-dplyr::rename(evidence,c("sample_id" = "experiment"))
  # =========================================#
  # Assay data 
  # =========================================#
  if(verbose)
    message("Get peptide intensities")
  
  ev_quanities<-dplyr::select_if(evidence,is.numeric) %>%
    colnames()
  
  # Return Peptide abundances
  peptide_abundance <- aggregate(
    evidence$intensity,
    by = list(
      "sample_id" = evidence$sample_id,
      "sequence" = evidence$sequence
    ),
    FUN = "sum"
  ) %>%
    tidyr::spread(.,
                  key = sample_id,
                  value = x) %>%
    tibble::column_to_rownames(., var = "sequence")
  peptide_abundance <- log2(1 + peptide_abundance)
  colnames(peptide_abundance) <-
    janitor::make_clean_names(colnames(peptide_abundance)) %>% 
    gsub("^x","",.)
  
  #========================================#
  # Experiment ColData - Sample #
  #========================================#
  experiment_coldata <-
    match_coldata_assay(assay_matrix = peptide_abundance, anno = anno)
  
  #========================================#
  # Experiment Rowdata - Feature Annotation#
  #========================================#
  # keep relevant columns
  experiment_row_data <- data.frame(
    "row_id" = evidence$sequence,
    "id" = evidence$proteins,
    "gene" = evidence$gene_names,
    "protein_names" = evidence$protein_names,
    evidence[, c(
      "sequence",
      "length",
      "gene_names",
      "proteins",
      "is_proteotypic",
      "potential_contaminant",
      "reverse"
    )]
  ) %>% dplyr::distinct()
  
  # Convert to summarized experiment
  se_peptide <-
    SummarizedExperiment::SummarizedExperiment(
      assays = list("intensity" = experiment_coldata$assay),
      rowData = experiment_row_data,
      colData = experiment_coldata$colData
    )
  se_peptide
}

#' make_se_mq_phos
#' @author Dennis Friedel
#' @param phospho_sty_sites_txt dataframe; Phospho(STY)Sites.txt from maxquant txt folder
#' @param localization_prob_threshold numeric;
#' @param remove_low_quality Logical; if TRUE remove low-quality precursors.
#' @param anno data.frame; Annotation for samples in proteinGroups needs columns named SampleID,group and replicate
#'
#' @return SummarizedExperiment object.
#' @export
#'
#' @examples NULL
# phospho_sty_sites_txt = ms_res$phospho_sty_sites_txt
# anno = ms_res$summary_txt
# remove_low_quality = TRUE
# ptm_probability_flt = 0.75
# #make_se_mq_phos(phospho_sty_sites_txt = mq_results$phospho_sty_sites_txt
#  # anno = anno
#  # remove_low_quality = TRUE
#  # ptm_probability_flt = 0.75
# #)
make_se_mq_phos <- function(phospho_sty_sites_txt,
                            anno,
                            localization_prob_threshold = 0.75,
                            remove_low_quality = TRUE,
                            verbose=TRUE) {
  # =========================================#
  # Load and preproc res table for SummarizedExperimentConversion
  # =========================================#
  if(verbose)
    message("Load and control data")
  assertthat::assert_that(
    is.data.frame(phospho_sty_sites_txt),
    is.data.frame(anno),
    is.logical(remove_low_quality)
  )
  if(verbose)
    message("control anno")
  anno <- check_coldata(anno = anno)
  
  if(verbose)
    message("Load STY-File")
  phospho_groups <- phospho_sty_sites_txt %>%
    as.data.frame() %>%
    janitor::clean_names()
  
  if(verbose)
    message("Create Residue Annotation")
  
  # Make a unique Identifier
  rownames(phospho_groups) <-
    make_ptm_identifier(input_dataframe = phospho_groups, type = "ptm")
  
  # Make a full annotation with uniprot_residue, gene_residue 
  ptm_df<-make_ptm_anno(se_object = phospho_groups)
  
  keep<-intersect(rownames(ptm_df),rownames(phospho_groups))
  phospho_groups_annotated<-cbind(phospho_groups[keep, ], 
        ptm_df[keep, !colnames(ptm_df) %in% colnames(phospho_groups)])%>%
    as.data.frame()
  
  # Filter features by localization probability
  if (remove_low_quality == TRUE) {
    phospho_groups_annotated <-
      remove_low_quality_features(
        ms_table = phospho_groups_annotated,
        is_pg = FALSE,
        is_ptm = TRUE,
        thr_localization_probability = localization_prob_threshold
      )
  }
  # =========================================#
  # Assay data 
  # =========================================#
  
  # Get peptide quantities
  prefix <- "intensity_"
  sample_id_match <-
    unique(grep(
      paste(paste0("intensity_", anno$sample_id), collapse = "|"),
      colnames(phospho_groups),
      value = FALSE
    ))
  
  # Filter features by localization probability
  if(!any(duplicated(phospho_groups_annotated$gene_residue))){
    message("Duplicated gene Residues, fall back to row_id")
    rownames(phospho_groups_annotated)<-phospho_groups_annotated$gene_residue
  }
  
  phospho_groups_abundances <- phospho_groups_annotated[, sample_id_match]
  colnames(phospho_groups_abundances) <-
    gsub(paste0("intensity_"),
         "",
         colnames(phospho_groups_abundances))
  if (sum(anno$sample_id %in% colnames(phospho_groups_abundances)) < 2) {
    stop(
      "specified quantities prefix ('",
      prefix,
      "') does not indicate >1 columns",
      "\nRun make_se_mq_proteins() with the appropriate prefix argument",
      call. = FALSE
    )
  }
  if (any(!apply(phospho_groups_abundances, 2, is.numeric))) {
    stop(
      "specified quantities prefix ('",
      prefix,
      "') does not indicate numeric columns",
      "\nRun import_MaxQuant() with the appropriate quantities argument",
      call. = FALSE
    )
  }
  
  message("Set up assays based on phosphorylation status all, 1-fold, 2-fold, 3-fold")
  phos_intensity_list <-
    purrr::map(c("", "_1", "_2", "_3"), function(suffix) {
      phos_intensity <-
        phospho_groups_abundances[, paste0(anno$sample_id, suffix)]
      phos_intensity <- log2(1 + (phos_intensity))
      # Duno why this was zero but its now solved
      phos_intensity[phos_intensity==0]<-NA
      colnames(phos_intensity) <- anno$sample_id
      phos_intensity
    })
  names(phos_intensity_list) <-
    paste0("intensity", c("", "_1", "_2", "_3"))
  
  #========================================#
  # Experiment Rowdata - Feature Annotation#
  #========================================#
  pre_row_data <- phospho_groups_annotated[,-sample_id_match]
  
  row_data <- data.frame(
    "id" = pre_row_data$proteins,
    pre_row_data
  ) %>% dplyr::distinct()
  
  #========================================#
  # Experiment ColData - Sample #
  #========================================#
  experiment_coldata <-
    match_coldata_assay(assay_matrix = phos_intensity_list$intensity, anno = anno)$colData
  
  # Find position macth between SampleID and sample in protein
  message("Convert phos_peptides into SummarizedExperiment object")
  se_phos <-
    SummarizedExperiment::SummarizedExperiment(assays = phos_intensity_list,
                                               rowData = row_data,
                                               colData = experiment_coldata)
  se_phos
}

#' mq_stats_report
#' @author Dennis Friedel
#' @param mq_result List; List object containing maxquant tables.
#'
#' @return List; List object with summary of different MQ metrics as list elements
#' @description
#' Evaluates and returns additional metrics from mq_results and returns a list with MQ metrics
#' reporting about overall identification, charge distribution, hydrophobicity of peptides.
#'
#' @export
#'
#' @examples NULL
mq_stats_report <- function(mq_result) {
  assertthat::assert_that(is.list(mq_result))
  stat_tables <- list()
  
  # ============================================#
  ##### mq parameters & summary
  # ============================================#
  if ("summary_txt" %in% names((mq_result))) {
    mq_summary <- mq_result$summary_txt %>%
      dplyr::select(dplyr::contains(
        c(
          "file",
          "experiment",
          "fraction",
          "enzyme",
          "modifications",
          "ms_",
          "peaks",
          "lc_ms_run_type"
        )
      ))
    stat_tables[[length(stat_tables) + 1]] <- mq_summary
    names(stat_tables)[[length(stat_tables)]] <- "summary"
  } else {
    warning("summary_txt not available, some plots will be missing")
    print("summary_txt not available, some plots will be missing")
  }
  
  if ("parameters_txt" %in% names((mq_result))) {
    stat_tables[[length(stat_tables) + 1]] <- mq_result$parameters_txt
    names(stat_tables)[[length(stat_tables)]] <- "paramaters"
  } else {
    warning("parameters_txt not available, some plots will be missing")
    print("parameters_txt not available, some plots will be missing")
  }
  
  # ============================================#
  ##### protein/peptide/PTM identification
  # ============================================#
  if ("protein_groups_txt" %in% names((mq_result))) {
    protein_identification <-
      summarize_mq_identification(ms_data = mq_result$protein_groups_txt)
    
    stat_tables[[length(stat_tables) + 1]] <- protein_identification
    names(stat_tables)[[length(stat_tables)]] <-
      "protein_identification"
  } else {
    warning("protein_groups_txt not available, protein_identification will be missing")
    print("protein_groups_txt not available, protein_identification will be missing")
  }
  if ("peptides_txt" %in% names((mq_result))) {
    # overall identification
    peptide_identification <-
      summarize_mq_identification(ms_data = mq_result$peptides_txt)
    stat_tables[[length(stat_tables) + 1]] <- peptide_identification
    names(stat_tables)[[length(stat_tables)]] <-
      "peptide_identification"
    # hydrophobcity
    gravy_summary <-
      calc_gravy(peptides = mq_result$peptides_txt)
    stat_tables[[length(stat_tables) + 1]] <- gravy_summary
    names(stat_tables)[[length(stat_tables)]] <- "gravy_summary"
  } else {
    warning(
      "peptides_txt not available, peptide_identification & gravy summary will be missing"
    )
    print(
      "peptides_txt not available, peptide_identification & gravy summary will be missing"
    )
  }
  
  if ("phospho_sty_sites_txt" %in% names((mq_result))) {
    ptm_identification <-
      summarize_mq_identification(ms_data = mq_result$phospho_sty_sites_txt,
                                  mod_data = TRUE)
    stat_tables[[length(stat_tables) + 1]] <- ptm_identification
    names(stat_tables)[[length(stat_tables)]] <-
      "ptm_identification"
  } else {
    warning("phospho_sty_sites_txt not available, ptm_identification will be missing")
    print("phospho_sty_sites_txt not available, ptm_identification will be missing")
  }
  
  # ============================================#
  ##### peptide/protein traits
  # ============================================#
  if ("evidence_txt" %in% names((mq_result))) {
    # charge state of precursor ions
    charge <-
      aggregate_ms_stats(
        ms_data = mq_result$evidence_txt,
        experiment_column = "experiment",
        feature_column = "charge",
        in_percent = TRUE
      )
    stat_tables[[length(stat_tables) + 1]] <- charge
    names(stat_tables)[[length(stat_tables)]] <- "charge"
    
    # protease specifity
    protease_specfity <-
      aggregate_ms_stats(
        ms_data = mq_result$evidence_txt,
        experiment_column = "experiment",
        feature_column = "missed_cleavages",
        in_percent = TRUE
      )
    stat_tables[[length(stat_tables) + 1]] <- protease_specfity
    names(stat_tables)[[length(stat_tables)]] <- "protease_specfity"
  } else {
    warning("evidence_txt not available, charge & protease specifity will be missing")
    print("evidence_txt not available, charge & protease specifity will be missing")
  }
  
  return(stat_tables)
}


# ============================#
#### DIA-NN #####
# ============================#

#' create_diann_se_list
#'
#' @param result_tsv 
#' @param anno 
#' @param remove_low_quality 
#' @param threshold_precursor 
#' @param protein_theshold 
#' @param pq_q_threshold 
#' @param gg_q_threshold 
#' @param fasta_file 
#' @param quantities 
#' @param verbose 
#' @param number_cores 
#'
#' @return NULL
#' @export
#'
#' @examples NULL
##create_diann_se_list(
# # result_tsv = test_tmp$report_parquet,
# # stats_tsv= test_tmp$report_stats_tsv,
# # remove_low_quality = TRUE,
# # use_only_proteotypic = TRUE,
# # threshold_precursor = 0.01,
# # protein_theshold = 0.01,
# # pq_q_threshold = 0.01,
# # gg_q_threshold = 0.01,
# #pg_quantity=c("pg_max_lfq"),#"pg_normalised"
# #gg_quantity=c("genes_max_lfq"),#"genes_normalised"
# #quantities = c("peptide","ptm","protein","gene")
# #verbose=TRUE
# #)
#### PTM Analysis
##ms_results <- load_ms_results(ms_result_dir = "./data/testDIANN181_phospho_jurkatcells120min/")
##  diann_se_list<-create_diann_se_list(
##  result_tsv = ms_results$result_tsv,
##  stats_tsv = ms_results$result_stats_tsv,
##  fasta_path = "./data/uniprotkb_Human_AND_reviewed_true_AND_m_2025_01_16.fasta",
##  verbose = TRUE
##)
create_diann_se_list<-function(result_tsv,
                               stats_tsv,
                               remove_low_quality = TRUE,
                               use_only_proteotypic = TRUE,
                               threshold_precursor = 0.01,
                               protein_theshold = 0.01,
                               pq_q_threshold = 0.01,
                               gg_q_threshold = 0.01,
                               pg_quantity=c("pg_max_lfq"),#"pg_normalised"
                               gg_quantity=c("genes_max_lfq"),#"genes_normalised"
                               quantities = c(
                                 "peptide",
                                 "ptm",
                                 "protein",
                                 "gene"
                               ),
                               fasta_path = NULL,
                               verbose=TRUE){
  stats_tsv<-as.data.frame(stats_tsv)
  ## Assert right input
  assertthat::assert_that(is.data.frame(result_tsv),
                          is.data.frame(stats_tsv),
                          any(quantities %in% c(
                            "peptide",
                            "ptm",
                            "protein",
                            "gene"
                          ),
                          is.numeric(threshold_precursor),
                          is.numeric(pq_q_threshold),
                          is.logical(verbose)
                          ))
  diann_se_list <- list()
  
  # =========#
  ##### Peptide
  # =========#
  if ("peptide" %in% quantities) {
    pep_se<-make_se_diann(result_tsv = result_tsv,
                          anno = stats_tsv,
                          row_header = "stripped_sequence",
                          q_id_header = "precursor_normalised",
                          remove_low_quality = remove_low_quality,
                          only_proteotypic =use_only_proteotypic,
                          threshold_precursor = threshold_precursor,
                          protein_theshold = 1,
                          pq_q_threshold = 1,
                          gg_q_threshold = 1)
    diann_se_list[[length(diann_se_list) + 1]] <- pep_se
    names(diann_se_list)[[length(diann_se_list)]] <- "pep_se"
  }
  #=================#
  #### PTM
  #=================#
  result_ptm<-result_tsv[grep("UniMod[:]21", result_tsv[, "modified_sequence"]),]
  if (("ptm" %in% quantities) & (nrow(result_ptm)!=0)) {
      ptm_se<-make_se_diann(result_tsv = result_ptm,
                            anno = stats_tsv,
                            row_header = "modified_sequence",
                            q_id_header = "precursor_normalised",
                            remove_low_quality = remove_low_quality,
                            only_proteotypic =use_only_proteotypic,
                            threshold_precursor = threshold_precursor,
                            protein_theshold = 1,
                            pq_q_threshold = 1,
                            gg_q_threshold = 1)
      if(!is.null(fasta_path)&file.exists(fasta_path)){
        message("Convert PTM Se to residue level using fasta file: ",fasta_path)
        ptm_se<-convert_ptmse_residue(ptm_se = ptm_se,fasta_path = fasta_path)
      }
      diann_se_list[[length(diann_se_list) + 1]] <- ptm_se
      names(diann_se_list)[[length(diann_se_list)]] <- "ptm_se"  
  }
  
  # =========#
  ##### Protein
  # =========#
  if ("protein" %in% quantities) {
    assertthat::assert_that(pg_quantity%in%c("pg_max_lfq","pg_normalised"))
    protein_se<-make_se_diann(result_tsv = result_tsv,
                              anno = stats_tsv,
                              row_header = "protein_ids",
                              q_id_header = pg_quantity,
                              remove_low_quality = remove_low_quality,
                              only_proteotypic =use_only_proteotypic,
                              threshold_precursor = threshold_precursor,
                              protein_theshold = protein_theshold,
                              pq_q_threshold = pq_q_threshold,
                              gg_q_threshold = gg_q_threshold)
    diann_se_list[[length(diann_se_list) + 1]] <- protein_se
    names(diann_se_list)[[length(diann_se_list)]] <- "pg_se"
    
  }
  # =========#
  ##### Gene summary
  # =========#
  if ("gene" %in% quantities) {
    assertthat::assert_that(gg_quantity%in%c("genes_normalised","genes_max_lfq"))
    assertthat::assert_that(gg_quantity%in%colnames(result_tsv))
    gene_se<-make_se_diann(result_tsv = result_tsv,
                           anno = stats_tsv,
                           row_header = "genes",
                           q_id_header = gg_quantity,
                           remove_low_quality = remove_low_quality,
                           only_proteotypic =use_only_proteotypic,
                           threshold_precursor = threshold_precursor,
                           protein_theshold = protein_theshold,
                           pq_q_threshold = pq_q_threshold,
                           gg_q_threshold = gg_q_threshold)
    gene_se
    diann_se_list[[length(diann_se_list) + 1]] <- gene_se
    names(diann_se_list)[[length(diann_se_list)]] <- "gg_se"
  }
  # =========#
  ##### Gene summary unique
  # =========#
  if ("gene" %in% quantities) {
    assertthat::assert_that("genes_max_lfq_unique"%in%colnames(result_tsv))
    gene_se<-make_se_diann(result_tsv = result_tsv,
                           anno = stats_tsv,
                           row_header = "genes",
                           q_id_header = "genes_max_lfq_unique",
                           remove_low_quality = remove_low_quality,
                           only_proteotypic =use_only_proteotypic,
                           threshold_precursor = threshold_precursor,
                           protein_theshold = protein_theshold,
                           pq_q_threshold = pq_q_threshold,
                           gg_q_threshold = gg_q_threshold)
    gene_se
    diann_se_list[[length(diann_se_list) + 1]] <- gene_se
    names(diann_se_list)[[length(diann_se_list)]] <- "ggu_se"
  }
  return(diann_se_list)
}  

#' get_diann_quantity
#'
#' @param result_table Dataframe; DIA-NN result table
#' @param row_id_header Character; Name of the column which contains the row ID
#' @param quantity_id_header Character; Name of the column which contains the quantity ID
#' @param is_ptm Boolean; if TRUE, the function will filter for PTM features
#' @param filter_low_quality Boolean; if TRUE, the function will filter out low-quality features based on q-values
#' @param only_proteotypic Boolean; if TRUE, the function will filter for proteotypic features
#' @param precursor_q numeric; threshold for precursor q-value between 0 and 1 . Default is 0.01
#' @param protein_q numeric; threshold for protein q-value between 0 and 1 Default is 1
#' @param protein_group_q numeric; threshold for protein group q-value between 0 and 1 Default is 1
#' @param gene_group_q numeric; threshold for gene group q-value between 0 and 1 Default is 1
#'
#' @return Dataframe; A data frame containing the summarized quantity information, with rows representing features and columns representing samples.
#' @export
#'
#' @examples NULL
#' # Aquire Peptide Intensity
#' se_assay <- get_diann_quantity(
    # #result_table = result_tsv,
    # #row_id_header = "stripped_sequence",
    # #quantity_id_header = "precursor_normalised",
    # #filter_low_quality = TRUE,
    # #only_proteotypic = TRUE,
    # #precursor_q = 0.01,
    # #protein_q = 1,
    # #protein_group_q = 1,
    # #gene_group_q = 1
    # ) %>%
    #   log2(.)
get_diann_quantity <- function(result_table,
                               row_id_header = "stripped_sequence",
                               quantity_id_header = "precursor_normalised",
                               filter_low_quality = TRUE,
                               only_proteotypic = TRUE,
                               precursor_q = 0.01,
                               protein_q = 1,
                               protein_group_q = 1,
                               gene_group_q = 1) {
  assertthat::assert_that(
    is.data.frame(result_table),
    is.character(row_id_header),
    is.character(quantity_id_header),
    is.logical(only_proteotypic),
    is.numeric(precursor_q),
    is.numeric(protein_group_q),
    is.numeric(gene_group_q)
  )
  result_table <- janitor::clean_names(result_table)
  df <- data.table::as.data.table(result_table)
  if (only_proteotypic) {
    df <- df[which(df[["proteotypic"]] != 0),]
  }
  if (filter_low_quality) {
    df <- unique(df[which(
      df[[row_id_header]] != "" &
        df[[quantity_id_header]] > 0 &
        df[["q_value"]] <= precursor_q &
        df[["protein_q_value"]] <= protein_q &
        df[["pg_q_value"]] <= protein_group_q &
        df[["gg_q_value"]] <= gene_group_q
    ),
    c("run", row_id_header, quantity_id_header), with = FALSE])
  }
  dups <-
    any(duplicated(paste0(df[["run"]], ":", df[[quantity_id_header]])))
  if (dups) {
    warning("Multiple quantities per id: the maximum of these will be calculated")
    qtable <-
      pivot_aggregate(df, "run", row_id_header, quantity_id_header)
  }else {
    qtable <- pivot(df, "run", row_id_header, quantity_id_header)
  }
  colnames(qtable) <- janitor::make_clean_names(colnames(qtable))
  return(qtable)
}

#' make_se_diann - create SummarizedExperiment object from DIA-NN result table
#'
#' @param result_tsv Dataframe; DIA-NN result table
#' @param anno Dataframe; Sample annotation. mandatory columns are sample_id, group, replicate
#' @param row_header chracter; Name of the column which contains the row ID
#' @param q_id_header character; Name of the column which contains the quantity ID
#' @param remove_low_quality Boolean; if TRUE, the function will filter out low-quality features based on q-values
#' @param only_proteotypic Boolean; if TRUE, the function will filter for proteotypic features
#' @param threshold_precursor Numeric; threshold for precursor q-value between 0 and 1 . Default is 0.01
#' @param protein_theshold Numeric; threshold for protein q-value between 0 and 1 Default is 1
#' @param pq_q_threshold Numeric; threshold for protein group q-value between 0 and 1 Default is 1
#' @param gg_q_threshold Numeric; threshold for gene group q-value between 0 and 1 Default is 1
#'
#' @return SummarizedExperiment; A SummarizedExperiment object containing the summarized quantity information, with rows representing features and columns representing samples.
#' @export
#' @description
#' This fucntion summarizes the quantitative sample and feature information into
#' SummarizedExperiment Format required for further downstream analysis
#' NOTE: in this Version it is avoided to use Annotation to further 
#' to enable later annotation of the experiment which comes with less constraints 
#' 
#' @examples NULL
# anno=stats_tsv
# row_header = "stripped_sequence"
# q_id_header = "precursor_normalised"
# remove_low_quality = TRUE
# only_proteotypic =TRUE
# threshold_precursor = 0.01
# protein_theshold = 1
# pq_q_threshold = 1
# gg_q_threshold = 1
make_se_diann <-
  function(result_tsv,
           anno = stats_tsv,
           row_header = "stripped_sequence",
           q_id_header = "precursor_normalised",
           remove_low_quality = TRUE,
           only_proteotypic = TRUE,
           threshold_precursor = 0.01,
           protein_theshold = 1,
           pq_q_threshold = 1,
           gg_q_threshold = 1) {
    
    ## Assert right input
    assertthat::assert_that(
      is.data.frame(result_tsv),
      is.data.frame(anno),
      is.numeric(threshold_precursor),
      is.numeric(pq_q_threshold),
      is.numeric(protein_theshold),
      is.numeric(pq_q_threshold),
      is.numeric(gg_q_threshold)
    )
    #========================================#
    # Experiment Assay #
    #========================================#
    se_assay <- get_diann_quantity(
      result_table = result_tsv,
      row_id_header = row_header,
      quantity_id_header = q_id_header,
      filter_low_quality = remove_low_quality,
      only_proteotypic = TRUE,
      precursor_q = threshold_precursor,
      protein_q = protein_theshold,
      protein_group_q = pq_q_threshold,
      gene_group_q = gg_q_threshold
    ) %>%
      log2(.)
    
    #========================================#
    # Experiment ColData - Sample #
    #========================================#
    experiment_coldata <-
      data.frame("sample_id" = colnames(se_assay),
                 row.names = colnames(se_assay))
    rownames(anno) <- gsub(".*/", "", anno$file_name) %>% gsub("[.]d", "", .) %>%
      janitor::make_clean_names()
    
    # Combine with anno 
    assertthat::assert_that(!any(duplicated(rownames(anno))))
    if(!all(rownames(anno)%in%rownames(experiment_coldata))){
      anno<-anno[rownames(anno)%in%rownames(experiment_coldata),]
    }
    assertthat::assert_that(all(rownames(anno)%in%rownames(experiment_coldata))&
                              all(rownames(experiment_coldata)%in%rownames(anno)))
    
    # Evalaute Anno matching columns of DIANN
    experiment_coldata <- cbind(experiment_coldata, anno[rownames(experiment_coldata), ])
    
    #========================================#
    # Experiment Rowdata - Feature Annotation#
    #========================================#
    # keep relevant columns
    experiment_row_data <- data.frame(
      "row_id" = rownames(se_assay),
      "id" = plyr::mapvalues(rownames(se_assay), result_tsv[[row_header]], result_tsv[["protein_ids"]], warn_missing = FALSE),
      "gene" = plyr::mapvalues(rownames(se_assay), result_tsv[[row_header]], result_tsv[["genes"]], warn_missing = FALSE),
      "protein_names" = plyr::mapvalues(rownames(se_assay), result_tsv[[row_header]], result_tsv[["protein_names"]], warn_missing = FALSE)
    ) %>% dplyr::distinct()
    
    # ## add ensembl ids
    # ensembl_ids <-
    #   clusterProfiler::bitr(
    #     experiment_row_data[["id"]],
    #     fromType = "UNIPROT",
    #     toType = "ENTREZID",
    #     OrgDb = "org.Hs.eg.db"
    #   )
    # experiment_row_data$ensembl_ids <-
    #   plyr::mapvalues(experiment_row_data$id,
    #                   ensembl_ids[["UNIPROT"]],
    #                   ensembl_ids[["ENTREZID"]],
    #                   warn_missing = FALSE)
    
    #========================================#
    # Convert to summarized experiment #
    #========================================#
    se <-
      SummarizedExperiment::SummarizedExperiment(
        assays = list("intensity" = se_assay[experiment_row_data$row_id, experiment_coldata$sample_id]),
        rowData = experiment_row_data,
        colData = experiment_coldata
      )
    se
  }

#' convert_PTMse_residue
#'
#' @param ptm_se SummarizedExperiment; ptm_object obtained  
#' @param fasta_path character; Path to fasta file 
#'
#' @returns SummarizedExperiment with gene residue as rownames and Intensities 
#' @export
#'
#' @examples NULL
##ptm_se<-test_diann$ptm_se
##fasta_path<-"./data/uniprotkb_Human_AND_reviewed_true_AND_m_2025_01_16.fasta"
##test<-convert_ptmse_residue(ptm_se = ptm_se,fasta_path = fasta_path)
convert_ptmse_residue<-function(ptm_se,fasta_path){
  
  assertthat::assert_that(methods::is(ptm_se, "SummarizedExperiment"))
  assertthat::assert_that(file.exists(fasta_path))
  
  ### Get residue-information
  ptm_anno<-
    make_ptm_anno(se_object = ptm_se,
                  fasta_path = fasta_path)
  
  plyr::join(as.data.frame(ptm_se@elementMetadata),
             ptm_anno,
             by = "row_id") %>%
    .[,grep("[.]",colnames(.),invert = T)]%>%
    S4Vectors::DataFrame(.) ->ptm_se@elementMetadata
  
  ptm_se_anno<-ptm_se[!is.na(ptm_se@elementMetadata$Name),]
  
  #### Residue centric Intensities 
  ptm_residue_se<-.se_ptm_aggregate(se_object = ptm_se_anno,new_id = "Name")
  return(ptm_residue_se)
}

# ============================#
####  Spectronaut  ######
# ============================#

#' import_spectronaut
#' @author Dennis Friedel
#' @param file character; file to spectronaut report tsv
#' @param anno data.frame; Annotation for samples in result.tsv needs columns named label,condition and replicate
#' @param pg_q_threshold numeric;
#' @param eg_q_threshold numeric;
#' @param quantities character; vector containing the quanitites to generate summarized experiements from
#' @param fasta_file character; path to the fasta file
#' @param number_cores integer; number of cores to use to find position of modified amino acid in fasta sequence
#' @param save_data boolean; save results
#' @param output_dir character; directory in which the data will be saved
#' @param output_name character; name that will be used for the object
#'
#' @return NULL
#' @export
#'
#' @examples NULL
#' 
# #test <- import_spectronaut(
#  # file = "./tests/testthat/data/ms_data/result_spectronaut.tsv",
#   #anno = openxlsx::read.xlsx("./tests/testthat/data/ms_data/anno_spectronaut.xlsx"),
#   #quantities = c("proteins", "phospho_peptides"),
#   #fasta_file = "./tests/testthat/data/ms_data/UP000005640_9606.fasta",
#   #number_cores = 2)
import_spectronaut <- function(file,
                               anno,
                               pg_q_threshold = 0.01,
                               eg_q_threshold = 0.01,
                               quantities = c("peptides",
                                              "phospho_peptides",
                                              "proteins"),
                               fasta_file = NULL,
                               number_cores = 2,
                               save_data = FALSE,
                               output_dir = "./output",
                               output_name = "unproc_spectronaut_se") {
  # =========================================#
  # Control validity of paramters
  # =========================================#
  assertthat::assert_that(
    file.exists(file),
    is.data.frame(anno),
    is.numeric(pg_q_threshold),
    is.character(quantities),
    is.logical(save_data),
    is.character(output_name)
  )
  
  # =========================================#
  # Control if required columns for SummarizedExperiment are in anno
  # =========================================#
  print("Load control annotation-file")
  anno <- control_anno_df(anno)
  
  # =========================================#
  # Load and control data
  # =========================================#
  print("Load and control data")
  evidence <-
    data.table::fread(file, sep = "\t", header = TRUE) %>%
    as.data.frame() %>%
    janitor::clean_names()
  
  evidence <-
    evidence %>% dplyr::rename(sample_id = r_file_name,
                               group = r_condition,
                               replicate = r_replicate)
  evidence$sample_id <- plyr::mapvalues(
    evidence$sample_id,
    unique(evidence$sample_id),
    janitor::make_clean_names(unique(evidence$sample_id), allow_dupes = TRUE),
    warn_missing = FALSE
  )
  evidence$id <- gsub(";.*", "", evidence$pg_protein_accessions)
  
  if (all(!anno$sample_id %in% evidence$sample_id)) {
    stop("No entry in anno column sample_id is matching run in Spectronaut result")
  }
  
  # ===============================================#
  # Get quantities and create summarized experiments
  # ===============================================#
  se_list <- purrr::map(quantities, function(x) {
    list()
  })
  names(se_list) <- quantities
  
  #================================#
  # peptides
  #================================#
  if ("peptides" %in% quantities) {
    message("import peptides")
    peptides <- load_spectronaut_matrix(
      x = evidence,
      id_header = "pep_stripped_sequence",
      quantity_header = "pep_quantity",
      eg_q = eg_q_threshold,
      pg_q = pg_q_threshold
    )
    
    row_data <- make_rowdata_ms_se(
      df_ms = evidence,
      loaded_matrix = peptides,
      sample_header = "sample_id",
      id_header = "id",
      name_header = "pep_stripped_sequence",
      gene_header = NA,
      protein_name_header = NA
    )
    
    se_pep <- SummarizedExperiment::SummarizedExperiment(
      assays = peptides,
      rowData = row_data,
      colData = make_coldata_se(assay_matrix = peptides, anno = anno)
    )
    se_list$peptides <- se_pep
  }
  
  #================================#
  # PTM
  #================================#
  if ("phospho_peptides" %in% quantities) {
    message("import phospho-peptides")
    if (is.null(fasta_file))
      stop("FASTA File cannot be NULL if phospho-peptides are wanted")
    if (!file.exists(fasta_file))
      stop("FASTA File not found")
    
    # Filter only for phospho entries
    evidence$eg_phos_sequence <-
      gsub("[[]Oxi.*[]]", "", evidence$eg_modified_sequence) %>%
      gsub("[[]Ace.*[]]", "", .) %>%
      gsub("_", "", .)
    evidence_phos <-
      evidence[grep("Phospho", evidence$eg_phos_sequence),]
    
    phos_abundance <- load_spectronaut_matrix(
      x = evidence_phos,
      id_header = "eg_phos_sequence",
      quantity_header = "pep_quantity",
      eg_q = eg_q_threshold,
      pg_q = pg_q_threshold
    )
    
    row_data <- make_rowdata_ms_se(
      df_ms = evidence_phos,
      loaded_matrix = phos_abundance,
      sample_header = "sample_id",
      id_header = "id",
      name_header = "eg_phos_sequence",
      gene_header = NA,
      protein_name_header = NA
    )
    
    row_data_detail <-
      add_phos_detail_rd(
        row_data = row_data,
        fasta_file = fasta_file,
        fasta_cores = 2,
        delim = "[]"
      )
    
    se_phos <- SummarizedExperiment::SummarizedExperiment(
      assays = phos_abundance,
      rowData = row_data_detail,
      colData = make_coldata_se(assay_matrix = phos_abundance, anno = anno)
    )
    se_list$phospho_peptides <- se_phos
  }
  
  #===============================#
  # Proteins
  #===============================#
  if ("proteins" %in% quantities) {
    message("import proteins")
    proteins <- load_spectronaut_matrix(
      x = evidence,
      id_header = "id",
      quantity_header = "pg_quantity",
      eg_q = eg_q_threshold,
      pg_q = pg_q_threshold
    )
    
    row_data <- make_rowdata_ms_se(
      df_ms = evidence,
      loaded_matrix = proteins,
      sample_header = "sample_id",
      id_header = "id",
      name_header = "id",
      gene_header = NA,
      protein_name_header = NA
    )
    
    se_prot <- SummarizedExperiment::SummarizedExperiment(
      assays = proteins,
      rowData = row_data,
      colData = make_coldata_se(assay_matrix = proteins, anno = anno)
    )
    se_list$proteins <- se_prot
  }
  
  if (save_data) {
    out_cache_dir <- paste0(output_dir, "/import_spectronaut/cache/")
    dir.create(out_cache_dir, recursive = TRUE)
    output_file <- paste0(out_cache_dir, "/", output_name, ".rds")
    saveRDS(se_list, file = output_file)
  }
  se_list
}



#===============#
# assay_data
#===============#

#' load_spectronaut_matrix
#' @author Dennis Friedel
#' @param x data.table; input data for pivoting.
#' @param id_header Character; column name containing feature identifiers.
#' @param quantity_header Character; column name containing intensity/quantity values.
#' @param proteotypic.only Logical; retain only proteotypic peptides.
#' @param eg_q Numeric; EG.Q-value threshold.
#' @param pq_q Numeric; PG.Q-value threshold.
#'
#' @return NULL
#' @export
#'
#' @examples NULL
load_spectronaut_matrix <-
  function(x,
           id_header = "id",
           quantity_header = "pg_quantity",
           eg_q = 0.01,
           pg_q = 0.01) {
    df <- data.table::as.data.table(x)
    df <-
      unique(df[which(df[[id_header]] != "" &
                        df[[quantity_header]] >
                        0 &
                        df[["eg_qvalue"]] <= eg_q &
                        df[["pg_qvalue"]] <= pg_q &
                        df[["f_excluded_from_quantification"]] == FALSE &
                        df[["f_frg_loss_type"]] == "noloss"),
                c("sample_id", id_header, quantity_header), with = FALSE])
    
    is_duplicated <- any(duplicated(paste0(df[["sample_id"]],
                                           ":", df[[id_header]])))
    if (is_duplicated) {
      warning("Multiple quantities per id: the maximum of these will be calculated")
      pivot_aggregate(df, "sample_id", id_header, quantity_header)
    } else {
      pivot(df, "sample_id", id_header, quantity_header)
    }
  }

#' add_phos_detail_rd
#' @author Dennis Friedel
#' @param row_data dataframe; required columns are "id" "name","gene","protein_name"
#' @param fasta_file character; path to fasta file
#' @param fasta_cores integer; number of cores to use to find position of modified amino acid in fasta sequence
#'
#'
#' @return Data.frame with required default columns plus addtional columns required for phosphorproteomic analysis
#' @export
#'
#' @examples NULL
add_phos_detail_rd <- function(row_data,
                               fasta_file = "./tests/testthat/data/ms_data/UP000005640_9606.fasta",
                               fasta_cores = 2,
                               delim = c("()", "[]")) {
  assertthat::assert_that(file.exists(fasta_file),
                          is.character(delim),
                          is.data.frame(row_data))
  delim1 <- substr(delim, start = 1, stop = 1)
  delim2 <- substr(delim, start = 2, stop = 2)
  # Additional annotation for modified seqences
  row_data$sequence_window <-
    gsub(paste0("[", delim1, "].*[", delim2, "]"), "", row_data$row_id)
  
  row_data$amino_acid <-
    strsplit(row_data$row_id, paste0("[", delim1, "]")) %>%
    purrr::map(., function(x) {
      substr(x[1], start = nchar(x[1]), stop = nchar(x[1]))
    }) %>%
    unlist()
  
  #Add downstream relevant columns to rowdata and change rowname of se including position of PTM in protein
  fasta<-read_fasta("/mnt/NAS4/user_data/np-dennis/projects/Proteomics/ProteoLab/data/uniprotkb_Human_AND_reviewed_true_AND_m_2025_01_16.fasta")
  names(fasta) <-
    gsub("sp[|]", "", names(fasta_new)) %>% gsub("[|].*", "", .)
  
  row_data$positions_within_proteins <-
    furrr::future_map(seq(from = 1, to = nrow(row_data)), function(pg) {
      phos_entry <- row_data[pg, ]
      if (length(grep(phos_entry$id, names(fasta))) == 0)
        return(NA)
      
      fasta_seq <-
        fasta[names(fasta) %in% phos_entry$id][[1]]
      position_in_aa <-
        gregexpr(phos_entry$sequence_window, fasta_seq, perl = TRUE)[[1]][1]
      
      position_mod_in_seq <-
        gregexpr(paste0("[", delim1, "]"), phos_entry$row_id)[[1]][1] - 1 # minus one because +1 offset
      position_mod_in_protein <-
        position_in_aa + position_mod_in_seq
      position_mod_in_protein
    }, .progress = TRUE, .options = furrr::furrr_options(globals = "fasta_cores")) %>% unlist()
  
  
  # general identified
  row_data$idenitfier <- make_ptm_identifier(
    input_dataframe = row_data,
    type = "ptm",
    uniprot = "id",
    gene_name = "gene",
    amino_acid = "amino_acid",
    ptm_position = "positions_within_proteins",
    sequence = "sequence_window"
  )
  
  row_data$idenitfier_ptmgsea <- make_ptm_identifier(
    input_dataframe = row_data,
    type = "ptm_gsea",
    uniprot = "id",
    gene_name = "gene",
    amino_acid = "amino_acid",
    ptm_position = "positions_within_proteins",
    sequence = "sequence_window"
  )
  row_data
}