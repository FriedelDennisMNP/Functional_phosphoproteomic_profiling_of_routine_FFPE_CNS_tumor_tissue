# rip_functions.R
# Author: Dennis Friedel
# Date: 11.04-2024
# Modification: 11.04-2024
# Description: Utils for Project oriented working
#
# Detail:
#' - function to start a project with rip
#'  - building subdirectroy structure for input and output
#' - functions to add new dataset
#' - function to start new analysis


#' create_rip_envir
#'
#' @param rip_dir character; relative or absolute path in
#' which the folder structure should be generated
#'
#' @return
#' @export
#' @description
#' Creates directories and sub-directories which
#' can be used as backbone for project oriented working
#'
#' @examples
#' create_rip_envir(rip_dir = "./")
create_rip_envir <- function(rip_dir = "") {
  # Create a set of directories with some example subdirectories
  scripts_subdir <- paste0(
    paste0(rip_dir, "/data/"),
    c(
      "/sample_sheets/",
      "/omics/panel/",
      "/omics/exome/",
      "/omics/mrna/",
      "/omics/ms/",
      "/omics/methylation/"
    )
  )
  lapply(scripts_subdir, function(x) {
    dir.create(x, recursive = TRUE, showWarnings = FALSE)
  })

  # Create a set of directories in which R-Scripts/Python-scripts can be stored
  scripts_subdir <- paste0(paste0(rip_dir, "/R/"),
                           c("analysis/",
                             "utils/",
                             "diagnosis/"))
  lapply(scripts_subdir, function(x) {
    dir.create(x, recursive = TRUE, showWarnings = FALSE)
  })

  # Create a set of directories in which objects,
  # tables and plots resulting from scripts can be stored
  output_subdir <- paste0(paste0(rip_dir, "/output/"))
  lapply(output_subdir, function(x) {
    dir.create(x, recursive = TRUE, showWarnings = FALSE)
  })

  # Create a set of directories in
  # which genesets and databases can be stored can be any datatype
  library_subdir <-
    paste0(paste0(paste0(rip_dir, "/library/")),
           c("genesets",
             "tables"))
  lapply(library_subdir, function(x) {
    dir.create(x, recursive = TRUE, showWarnings = FALSE)
  })

  # Create a subdirectory for reporting
  # results should be ignored by git
  results_subdir <- paste0(rip_dir, ("/results"))
  lapply(results_subdir, function(x) {
    dir.create(x, recursive = TRUE, showWarnings = FALSE)
  })

  # Create a sub-directory for archive to keep old stuff..for reasons?
  archive_subdir <- paste0(rip_dir, ("/archive"))
  lapply(archive_subdir, function(x) {
    dir.create(x, recursive = TRUE, showWarnings = FALSE)
  })
  message("Created rip oriented enviroment")
  return()
}


#' add_new_data
#'
#' @param new_file character; path to new file which should
#'  be used as data in rip envir
#' @param data_dir character; path in project dir to which
#' new file should be copied. Default is data in the current wd
#'
#' @return
#' @export
#'
#' @examples
#'
add_new_data <-
  function(new_file,
           data_dir = "./data/") {
    assertthat::assert_that(file.exists(new_file))
    file.copy(new_file, to = data_dir)
    return()
  }

#' add_new_data_directory
#'
#' @param new_directory character; path to new directory
#' which should be used as data in rip envir
#' @param data_dir character; path in project dir to
#' which new dir should be copied. Default is data in the current wd
#'
#' @return
#' @export
#'
#' @examples
#'
add_new_data_directory <-
  function(new_directory,
           data_dir = "./data/") {
    assertthat::assert_that(dir.exists(new_directory))
    R.utils::copyDirectory(new_directory, to = data_dir)
    return()
  }

# Add analysis
#' add_analysis
#'
#' @param script_dir character;
#' name of the analysis directory and where it  should be created
#' @param analysis_name character;
#' name of the analysis
#' @param author_name character;
#' name of the author how is planning to di the analysis
#'
#' @return
#' @export
#'
#' @examples
#' add_analysis(analysis_name = "test",author_name = "Dennis Friedel")
add_analysis <-
  function(script_dir = "./R", analysis_name, author_name) {
    assertthat::assert_that(is.character(analysis_name),
                            is.character(author_name),
                            dir.exists(script_dir))

    # create subdir in script dir
    # with readme where the users can describe their plans
    analysis_dir <- paste0(script_dir, "/", analysis_name)
    lapply(analysis_dir, function(x) {
      dir.create(x, recursive = TRUE, showWarnings = FALSE)
    })
    # create RMD file with smart goals
    # Create a templacte file for project oriented working
    rmd_analyis_summary <-
      paste0(analysis_dir, "/", analysis_name, "_summary.rmd")

    file.create(rmd_analyis_summary, showWarnings = FALSE)
    writeLines(
      text = c(
        paste0("# ", analysis_name),
        paste0("Author: ", author_name, "; "),
        paste0("Date: ", Sys.Date()),
        "\n",
        paste0("## Description: "),
        "\n",
        "\n",
        "Aim of the analysis:\n",
        "Condition for completeness:\n",
        "Which datasets will be used:\n",
        "Contributors:\n",
        "Aprroximated Timewindow:\n",
        ""
      ),
      sep = "\n",
      con = rmd_analyis_summary
    )
  }

# Create a templacte file for project oriented working
#' Title
#'
#' @param analysis_dir
#' @param analysis_name
#' @param script_name
#' @param author_name
#' @param dataname
#' @param short_description
#'
#' @return
#' @export
#'
#' @examples
#' add_analysis_script(analysis_dir = "./R/test/",
#'  analysis_name = "test",
#'  script_name = "test_script",
#'  author_name = "Dennis",
#'  dataname = "test_data",
#'  short_description = "This script is just a
#'  test script to show how add analayis script works")
add_analysis_script <-
  function(analysis_dir = "./R",
           analysis_name,
           script_name,
           author_name,
           dataname,
           short_description) {
    assertthat::assert_that(
      is.character(analysis_dir),
      dir.exists(analysis_dir),
      is.character(script_name),
      is.character(author_name),
      is.character(dataname),
      is.character(short_description)
    )
    # Create a templacte file for project oriented working
    template_file <- paste0(analysis_dir,
                            "/",
                            script_name,
                            ".R")

    file.create(template_file, showWarnings = FALSE)
    writeLines(
      text = c(
        paste0(
          "####=============== Rscript: ",
          script_name,
          "===============####"
        ),
        paste0("# Author:", author_name),
        paste0("# Date: ", Sys.Date()),
        "# Modification:",
        paste0("# Description: ", short_description),
        "# Detail:",
        "#'",
        "####===================================####",
        "",
        paste0(
          "dataset  = '",
          dataname,
          "' # - Name of the dataset that is going to be analysed. "
        ),
        paste0(
          "analysis = '",
          analysis_name,
          "' # - Name of the analysis e.g Marker Identification"
        ),
        ""
      ),
      sep = "\n",
      con = template_file
    )
  }

#' create_output_path
#'
#' @param project_dir character;
#' absolute/realtive path to current working dir
#' @param dataset_name character;
#' name of the dataset
#' @param analysis_name character;
#' descritpion of the analysis proteomics/metabolimics etc
#'
#' @return
#' @export
#'
#' @examples
create_output_path <-
  function(output_dir = "./",
           dataset_name = get(x = "dataset"),
           analysis_name = get(x = "analysis")) {
    assertthat::assert_that(
      is.character(output_dir),
      dir.exists(output_dir),
      is.character(dataset_name),
      is.character(analysis_name)
    )
    outdir <- paste0(output_dir,
                    "/output/",
                    dataset_name,
                    "/",
                    analysis_name,
                    "/")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    outdir
  }

#' save_here
#'
#' @param ...
#' @param output_dir
#' @param dataset_name
#' @param analysis_name
#' @param object_name
#'
#' @return returns path to which the data will be saved
#' @export
#' @details
#' Uses here and information provided by the user to save datatypes into their corresponding subfolders in output
#' eg. rds object are directed into cache while figures to plots and excel tables to tables
#' If no datatype or name is provided, then only a directory will be created
#' @examples
#'
save_here = function(output_dir = "./",
                     dataset_name = get(x = "dataset"),
                     analysis_name = get(x = "analysis"),
                     object_name = NULL,
                     ...) {

  assertthat::assert_that(
    is.character(output_dir),
    dir.exists(output_dir),
    is.character(dataset_name),
    is.character(analysis_name),
    is.character(object_name)
  )

  # Create ouptut path
  outpath <-
    create_output_path(
      output_dir = output_dir,
      dataset_name = dataset_name,
      analysis_name = analysis_name
    )

  path_destiny = paste0(outpath,"/",object_name)
  message("new file path:", path_destiny)
  path_destiny
}

# ===========================================================================#
#### Color Functions
# ===========================================================================#

# # Create list of colors
# brewer_colors <- c("Accent", "Dark2", "Pastel1", "Pastel2")
# selected_color_brewer <-
#   (RColorBrewer::brewer.pal.info)[brewer_colors, ]
# selected_color_brewer <-
#   purrr::map(rownames(selected_color_brewer), function(x) {
#     RColorBrewer::brewer.pal(n = selected_color_brewer[x, ]$maxcolor, name = x)
#   })
# names(selected_color_brewer) <- brewer_colors

# all_colors <- colors()
# decent <- all_colors[grep("4", all_colors)]
# nphd_categorial_colors <-
#   list("decent1" = decent[seq(1, (length(decent) / 2))],
#        "decent2" = decent[seq(length(decent) / 2, length(decent))])
# nphd_categorial_colors <-
#   c(selected_color_brewer, nphd_categorial_colors)

#' seq_color_scale
#'
#' @param col_limit
#' @param palette_set
#'
#' @return
#' @export
#' @description
#' Set color gardient in accordance to provided set from brewer pal
#'
#'
#' @examples
seq_color_scale <- function(col_limit = 5,
                            palette_set = c("PRGn", "RdGy")) {
  circlize::colorRamp2(seq(-col_limit, col_limit, (col_limit / 5)),
                       rev(RColorBrewer::brewer.pal(11, palette_set)))
}

#' col_fun
#'
#' @param min integer
#' @param mid integer
#' @param max integer
#' @param mincolorm character; color to use for negative value
#' @param midcolor character; color to use for values close to 0
#' @param maxcolor character; color to use for high values
#'
#' @return
#' @export
#'
#' @examples
col_fun <-
  function(min = 0,
           mid = 2,
           max = 3,
           mincolorm = "purple4",
           midcolor = "yellow4",
           maxcolor = "green4") {
    circlize::colorRamp2(c(min, mid, max), c(mincolorm, midcolor, maxcolor))
  }


#' cat_color_code
#'
#' @param ggsci_colorpalette Color palettes from ggsci
#' @param categorial_vector character; vector with group
#' @param use_ggsci
#'
#' @return
#' @export
#'
#' @examples
cat_color_code <-
  function(ggsci_colorpalette,
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
