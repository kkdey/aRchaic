#' @title Read length distribution plots for an MFF file
#'
#' @description For an MFF file, we plot the read length distribution under three settings -
#'  all reads, reads containing C to T mutations and reads containing C to T towards the
#'  ends of the fragments.
#'
#' @param file an MFF file
#' @param end_break The position on the read from the end of the fragment such that if
#' C->T change occurs beyond that end break, it is considered in the third class of
#' reads (potentially damaged reads)
#' @param cols The color profiles of the three plots. Defaults to "red", "green" and "blue".
#' @param cex_legend The cex of the legend. Defaults to 0.5.
#' @param cex_main The cex of the title of the figure. Defaults to 1.
#'
#' @return A 3-fold density plot of an MFF file under three scenarios - all reads,
#' reads containing C to T mutations and reads containing C to T towards the ends
#' of the fragments.
#'
#' @keywords read_length
#'
#' @export

read_length_distribution <- function(file,
                                     end_break = 5,
                                     cols = c("red", "green", "blue"),
                                     cex_legend = 0.5,
                                     cex.main = 1){

  ###  read in all the files in the directory

    tab1 <- read.csv(paste0(file), header=FALSE)
    read_length <- tab1$V2 + tab1$V3
    indices <- grep("C->T", tab1$V1)
    indices <- c(indices, grep("G->A", tab1$V1))
    read_length_CtoT <- tab1[indices, ]$V2 + tab1[indices, ]$V3 ## read lengths for all C to T
    indices2 <- which(tab1$V2 < end_break | tab1$V3 < end_break)
    indices_matched <- intersect(indices, indices2)
    read_length_3 <- (tab1[indices_matched, ]$V2 + tab1[indices_matched, ]$V3)
    plot(table(read_length)/sum(table(read_length)), type="o", col=cols[1],
         main=paste0(ancient_names[l]), cex.main = cex.main, xlab="read pos", ylab="prop of occur")
    points(table(read_length_CtoT)/sum(table(read_length_CtoT)), type="o", col=cols[2])
    points(table(read_length_3)/ sum(table(read_length_3)), type="o", col=cols[3])
    legend("topright", fill=c("red", "green", "blue"),
           legend = c("all", "CtoT", paste0("CtoT < ", end_break)), cex=cex_legend)
    cat("we are at sample", l, "\n")
}
