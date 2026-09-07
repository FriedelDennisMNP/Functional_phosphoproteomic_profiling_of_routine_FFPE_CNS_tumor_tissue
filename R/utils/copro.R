#### CoPro Function ####

##### Load Genome positions
anno_genome_size<-"./data/hs_genome/anno.genome.size.xlsx" %>%
  openxlsx::read.xlsx(., colNames = F)
anno_genome_pq<-"./data/hs_genome/anno.genome.pq.xlsx" %>%
  openxlsx::read.xlsx(., colNames = F)
anno_genome_chr<-"./data/hs_genome/anno.genome.chr.xlsx" %>%
  openxlsx::read.xlsx(., colNames = F)
chr_line_pos <-
  cumsum(as.numeric(head(anno_genome_size$X1, -1)))
names(chr_line_pos) <- c((2:22), "X", "Y")
pq_line_pos <-
  cumsum(c(0, head(anno_genome_size$X1, -1))) + anno_genome_pq$X1
anno_xtick_pos <-
  (cumsum(c(0, head(
    anno_genome_size$X1, -1
  ))) + anno_genome_size$X1 / 2)


#' roll_align_all
#'
#' @param numdata
#' @param window_size
#' @param fun
#'
#' @return
#' @export
#'
#' @examples
#'
roll_align_all<-function(numdata,window_size,fun){
  assertthat::assert_that(is.numeric(numdata))
  roll_center <-
    data.table::frollapply(numdata, align = "center", n = window_size,FUN = fun)
  rolls_right <-
    data.table::frollapply(numdata, align = "right", n = window_size,FUN = fun)
  rolls_left <-
    data.table::frollapply(numdata, align = "left", n = window_size,FUN = fun)
  keep_rolling<-data.frame(roll_center,rolls_right,rolls_left) %>% rowMeans(.,na.rm=T)
  keep_rolling
}

#' smooth_roller
#'
#' @param vals numeric; values to smooth, should be ordered for example by genomic position
#' @param iteration integer; number of iteration in which the rolling mean will be applied
#' @param align character; parameter passed down to align in frollmean().
#' define if rolling window covers preceding rows ("right"), following rows ("left") or centered ("center"). Defaults to "right"
#' @param window_size integer; size of the rolling window.
#'
#' @return
#' @export
#'
#' @examples
smooth_roller <-
  function(vals,
           iteration,
           window_size = 3,
           fun = c("mean", "median")) {
    # Repetitive smoothing
    itrolls <- vals
    for (i in 1:iteration) {
      itrolls <-
        roll_align_all(numdata = itrolls,
                       window_size = window_size,
                       fun = fun)
    }
    itrolls
  }

#' plot_chromprot
#'
#' @param long_df
#' @param intensity_value
#' @param plot_title
#' @param col_high
#' @param col_mid
#' @param col_low
#' @param breaks
#' @param break_length
#' @param cumu_chr_pos
#' @param cumu_pq_pos
#' @param xtick_pos
#' @param genome_size
#' @param genome_anno
#' @param dot_alpha
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_chromprot<-function(long_df,
                         intensity_value,
                         plot_title,
                         col_high = "orange3",
                         col_mid = "grey70",
                         col_low = "dodgerblue4",
                         dot_alpha=1,
                         breaks = 3,
                         break_length = 13,
                         cumu_chr_pos = get("chr_line_pos"),
                         cumu_pq_pos = get("pq_line_pos"),
                         xtick_pos = get("anno_xtick_pos"),
                         genome_size = get("anno_genome_size"),
                         genome_anno = get("anno_genome_chr")
){
  
  assertthat::assert_that(
    is.data.frame(long_df),
    intensity_value %in% colnames(long_df),
    "start_position" %in% colnames(long_df),
    all(long_df[["start_position"]][long_df[["chromosome_name"]] == 2] >
          long_df[["start_position"]][long_df[["chromosome_name"]] == 1])
  )
  assertthat::assert_that(is.numeric(cumu_chr_pos),
                          assertthat::are_equal(length(cumu_chr_pos), 23),
                          all(names(cumu_chr_pos) %in% c(2:22, "X", "Y")))
  
  chrom_plot <-
    ggplot(long_df,
           aes(x = start_position,
               y = long_df[[intensity_value]],
               color = long_df[[intensity_value]],
               alpha = dot_alpha)) +
    geom_point() +
    scale_color_gradient2(
      high = col_high,
      mid = col_mid,
      low = col_low,
      name = intensity_value,
      breaks = round(seq(
        from = -breaks,
        to = breaks,
        length.out = break_length
      )),
      na.value = NA,
      midpoint = 0,
      guide = "colourbar"
    ) +
    geom_hline(yintercept = 0,
               color = "black",
               size = .5) +
    geom_vline(xintercept = cumu_chr_pos,
               color = "black",
               size = .3) +
    geom_vline(
      xintercept = cumu_pq_pos,
      color = "black",
      size = .3,
      linetype = "dotted"
    ) +
    geom_vline(xintercept = 0,
               color = "black",
               size = 1) +
    geom_vline(
      xintercept = sum(as.numeric(genome_size$X1)),
      color = "black",
      size = 1
    ) +
    ylim(-breaks, breaks) +
    xlab(label = "chromosomes") +
    ylab(label = intensity_value) +
    scale_x_continuous(breaks=xtick_pos,labels = genome_anno$X1) +
    theme(axis.text.x= element_text(size=9,angle = 45),
          panel.background = element_rect(fill = "white", colour = "white")) +
    ggtitle(plot_title)
}

