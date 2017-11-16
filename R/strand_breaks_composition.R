#' @title Computes strand break composition for mutation signatures for an MFF file
#'
#' @description This function computes the composition of nucleotides at the
#' strand breaks positions corresponding to mutation signatures (either all signatures
#' or the ones containing C->T only).
#'
#' @param file The file name for which the strand break composition is sought
#' @param flanking_bases The number of flanking bases in the signatures
#' @param dist_from_ends The distance from ends of the reads taken to filter the mutations.
#'                       A mutation is not considered if it occurs at a position farther than
#'                       this value from either ends of the read.
#' @param CtoT if TRUE, it focuses on reads with C->T mutations only, else looks at
#'             reads with all types of mutations in determining strand break composition.
#'
#' @return Returns a list with two elements, each being a vector of counts of nucleotides
#'         mapped to the left and right strand breaks
#' @keywords strand_breaks_composition
#' @export
#'

strand_breaks_composition <- function(file,
                                      flanking_bases=1,
                                      dist_from_ends = 10,
                                      CtoT = TRUE){

  data <- read.csv(file, header=FALSE)
 # data1 <- data[ which(data[,2] < dist_from_ends | data[,3] < dist_from_ends), ]
  data1 <- data[ which(data[,2] < dist_from_ends), ]

  sig_split <- sapply(data1[,1],
         function(x) return(paste0(strsplit(as.character(x), "")[[1]][(flanking_bases+1):(flanking_bases+4)], collapse = "")))
  plus_indices <-  which(data1[,6] == "+")

  if (CtoT){
    CtoT_indices <- which(sig_split == "C->T")
    total_indices <- intersect(CtoT_indices, plus_indices)
  }else{
    total_indices <- plus_indices
  }

  filtered_plus_data <- data1[total_indices,]
  tab1 = table(filtered_plus_data[,4])
  tab1_filtered = tab1[match(c("A", "C", "G", "T"), names(tab1))]
  tab1prop <- tab1_filtered/sum(tab1_filtered)
  left_strand_tab <- rbind(c(tab1prop[c(1,3)], 0, 0),
                           c(0,0, tab1prop[c(2,4)]))
  rownames(left_strand_tab) <- c("purine", "pyrimidine")
  colnames(left_strand_tab) <- c("A", "G", "C", "T")

  left_strand_tab <- t(left_strand_tab)

  color_profile <- list("type" = "per_column",
                        "col" = c("blue", "red"))


  logomaker(left_strand_tab,
            color_profile = color_profile,
            hist=TRUE,
            frame_width = 1,
            ic.scale = TRUE,
            yscale_change = TRUE,
            pop_name = "left strand break",
            xlab = "",
            ylab = "",
            yaxis=FALSE,
            main_fontsize = 15,
            xaxis_fontsize = 15,
            col_line_split="black",
            newpage=TRUE)


  data2 <- data[ which(data[,3] < dist_from_ends), ]
  minus_indices <-  which(data2[,6] == "-")

  if (CtoT){
    GtoA_indices <- which(sig_split == "G->A")
    total_indices <- intersect(GtoA_indices, minus_indices)
  }else{
    total_indices <- minus_indices
  }

  filtered_minus_data <- data1[total_indices,]
  tab2 = table(filtered_minus_data[,5])
  tab2_filtered = tab2[match(c("A", "C", "G", "T"), names(tab2))]
  tab2prop <- tab2_filtered/sum(tab2_filtered)
  right_strand_tab <- rbind(c(tab2prop[c(1,3)], 0, 0),
                           c(0,0, tab2prop[c(2,4)]))
  rownames(right_strand_tab) <- c("purine", "pyrimidine")
  colnames(right_strand_tab) <- c("A", "G", "C", "T")
  right_strand_tab <- t(right_strand_tab)

  color_profile <- list("type" = "per_column",
                        "col" = c("blue", "red"))


  logomaker(right_strand_tab,
            color_profile = color_profile,
            hist=TRUE,
            frame_width = 1,
            ic.scale = TRUE,
            yscale_change = TRUE,
            pop_name = "right strand break",
            xlab = "",
            ylab = "",
            yaxis=FALSE,
            main_fontsize = 15,
            xaxis_fontsize = 15,
            col_line_split="black",
            newpage=TRUE)


  ll <- list("left-strand-break" = tab1prop,
             "right-strand-break" = tab2prop)
  return(ll)
}
