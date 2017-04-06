#' @title Computes flanking base composition for a mutation type in an MFF file.
#'
#' @description This function computes the composition of nucleotides at the
#' left and right flanking bases for a particular mutation type (say "C->T").
#'
#' @param file an MFF file.
#' @param flanking_bases The number of flanking bases in the signatures
#' @param dist_from_ends The distance from ends of the reads taken to filter the mutations.
#'                       A mutation is not considered if it occurs at a position farther than
#'                       this value from either ends of the read.
#' @param pattern The type of mutation. The default is C->T.
#' @param plot if TRUE, then the function plots the nucleotide composition in a logo plot.
#'             Defaults to TRUE.
#' @param logo.control The control parameters for the logo plot representation. Defaults to
#'                     \code{list()}
#'
#' @return Returns composition of nucleotides at the left flanking base and the right flanking
#'         base. If plot is TRUE, then it plots a logo plot representation of the nucleotides.
#'
#' @keywords flanking_base_composition
#' @export
#'


flanking_base_composition <- function(file,
                                      flanking_bases = 1,
                                      dist_from_ends = 20,
                                      pattern = "C->T",
                                      plot = TRUE,
                                      logo.control = list()
){

  data <- read.csv(file, header=FALSE)
  data1 <- data[which(data[,6] == "+"),]
  data2 <- data1[union(which(data1[,2] < 20 ), which(data1[,3] < 20)),]
  sig_split <- do.call(rbind, lapply(data2[,1], function(x) return(strsplit(as.character(x), "")[[1]])))
  subs <- apply(sig_split, 1, function(x) return(paste0(x[(flanking_bases+1):(flanking_bases+4)], collapse = "")))
  left_flank <- apply(sig_split, 1, function(x) return(paste0(x[flanking_bases], collapse = "")))
  right_flank <- apply(sig_split, 1, function(x) return(paste0(x[(flanking_bases+5)], collapse = "")))

  data3 <- dplyr::tbl_df(data2) %>% dplyr::mutate(subs = subs) %>%
    dplyr::mutate(left_flank = left_flank) %>%
    dplyr::mutate(right_flank = right_flank) %>% filter(subs == pattern) %>%
    select(V7, left_flank) %>% group_by(left_flank)  %>% as.data.frame()

  tab <- tapply(data3[,1], data3[,2], sum)
  tab_filtered = tab[match(c("A", "C", "G", "T"), names(tab))]
  tab1prop <- tab_filtered/sum(tab_filtered)

  data4 <- dplyr::tbl_df(data2) %>% dplyr::mutate(subs = subs) %>%
    dplyr::mutate(left_flank = left_flank) %>%
    dplyr::mutate(right_flank = right_flank) %>% filter(subs == pattern) %>%
    select(V7, right_flank) %>% group_by(right_flank)  %>% as.data.frame()

  tab <- tapply(data4[,1], data4[,2], sum)
  tab_filtered = tab[match(c("A", "C", "G", "T"), names(tab))]
  tab2prop <- tab_filtered/sum(tab_filtered)

  if(plot){
    left_strand_tab <- matrix(0, 4, 1)
    left_strand_tab[,1] <- tab1prop

    right_strand_tab <- matrix(0, 4, 1)
    right_strand_tab[,1] <- tab2prop

    tab <- cbind(left_strand_tab, right_strand_tab)

    colnames(tab) <- c("left flank", "right flank")
    rownames(tab) <- c("A", "C", "G", "T")

    color_profile <- list("type" = "per_row",
                          "col" = RColorBrewer::brewer.pal(dim(tab)[1],name ="Spectral"))


    logo.control.default <- list(ic = NULL, hist = TRUE, frame_width =1,
                                 ic.scale = FALSE,
                                 xaxis_fontsize = 10, xlab_fontsize = 15,
                                 y_fontsize = 15, main_fontsize = 16, start = 0.001,
                                 yscale_change = TRUE, pop_name = paste0("logo plot flanking base:",pattern), xlab = "X",
                                 ylab = "composition", col_line_split = "grey80", scale0 = 0.01,
                                 scale1 = 0.99, newpage = TRUE)

    logo.control <- modifyList(logo.control.default, logo.control)

    do.call(Logolas::logomaker, c(list(table=tab, color_profile = color_profile),
                       logo.control))
  }

  ll <- list("left-flank" = tab1prop,
             "right-flank" = tab2prop)
  return(ll)
}
