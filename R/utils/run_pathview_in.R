#' run_pathview_in
#'
#' @param out_dir 
#' @param ... 
#'
#' @returns
#' @export
#'
#' @examples
#' run_pathview_in("path/to/your/folder",
  #' gene.data = kinase_scores,
  #' pathway.id = "04630",
  #' species = "hsa",
  #' gene.idtype = "SYMBOL",
  #' limit = list(gene = max(abs(kinase_scores))))
run_pathview_in <- function(out_dir, ...) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(out_dir)
  pathview::pathview(...)
}

