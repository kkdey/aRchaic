#' @title Decompose the theta matrix output from model_archaic into probabilities
#' of different mismatch features
#'
#' @description Decomposes the theta matrix outout from \code{model_archaic}
#'              into component probability distributions for the mismatch type,
#'              flanking bases, strand break base and position of strand
#'              composition. These are the probabilities that get represented
#'              in the logo plot in \code{plot_archaic}.
#' @param theta_out A theta matrix from a GoM model fit, with columns
#'                  representing clusters and rows representing the mutational
#'                  signatures.
#' @param max_pos The maximum distance from the end of the read that is used for
#'                filtering.
#' @return Returns a list with number of items equal to number of clusters.
#'         Each item of this list is another comprising of the probability
#'         distribution of types of mismatch, flanking base,
#'         distance from end of read and strand break composition
#'         separately.
#' @keywords decompose_theta
#' @import magrittr
#' @import gridBase
#' @import grid
#' @import ggplot2
#' @export



decompose_theta = function (theta_out, max_pos = 20) {
  flanking_bases <- 1
  signature_set <- rownames(theta_out)
  signature_patterns <- substring(signature_set, 1, 4+2*flanking_bases)

  breakbase <- substring(signature_set, 6+2*flanking_bases,  6+2*flanking_bases)
  library(dplyr)
  theta_break <- dplyr::tbl_df(data.frame(theta_out)) %>%
                    dplyr::mutate(sig = breakbase) %>%
                        dplyr::group_by(sig) %>%
                            dplyr::summarise_all(funs(sum)) %>%
                                    as.data.frame()
  rownames(theta_break) <- theta_break[,1]
  theta_break <- theta_break[,-1]

  theta_break <- theta_break[match(c("A", "C", "G", "T"), rownames(theta_break)),]
  breaks_theta <- theta_break

  if(flanking_bases%%1 != 0){
    stop("flanking bases not evenly distributed")
  }

  theta <- dplyr::tbl_df(data.frame(theta_out)) %>%
    dplyr::mutate(sig = signature_set) %>%
    dplyr::group_by(sig) %>%
    dplyr::summarise_all(funs(sum)) %>%
    as.data.frame()

  rownames(theta) <-  theta[,1]
  theta <- theta[,-1]

  sig_names <- rownames(theta)

  prob_mutation <- filter_signatures_only_location(t(theta_out), max_pos = max_pos, flanking_bases = flanking_bases)
  prob_mutation <- t(apply(prob_mutation, 1, function(x) {
    y <- x[!is.na(x)];
    return(y/sum(y))
  }))

  sig_split <- do.call(rbind,
                       lapply(sig_names,
                              function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))

  ncol_sig <- (4+2*flanking_bases)

  sub_pattern <- sapply(1:dim(sig_split)[1],
                        function(x) paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse=""))

  new_sig_split <- cbind(sig_split[,1:flanking_bases], sub_pattern, sig_split[,((ncol_sig - flanking_bases +1):ncol_sig)])
  colnames(new_sig_split) = NULL


  prop_patterns_list <- list()

  for(l in 1:dim(theta)[2]){
    prop_patterns_list[[l]] <- numeric();
    for(j in 1:ncol(new_sig_split)){
      temp <- tapply(theta[,l], factor(new_sig_split[,j], levels=c("A", "C", "G", "T",
                                                                   "C->T", "C->A", "C->G",
                                                                   "T->A", "T->C", "T->G")), sum)

      temp[is.na(temp)]=0
      prop_patterns_list[[l]] <- cbind(prop_patterns_list[[l]], temp)
    }
  }


  ll <- list()

  for(l in 1:dim(theta)[2]){
    ll[[l]] <- list()
    colnames(prop_patterns_list[[l]]) <- c(paste0("flank:", -1), "mismatch", paste0("flank:", 1))
    ll[[l]][["mismatch-flank"]] <- prop_patterns_list[[l]]
    ll[[l]][["strand-break-composition"]] <- breaks_theta[, l, drop=FALSE]
    ll[[l]][["mismatch-trend-read"]] <- prob_mutation[l,]
  }

  return(ll)
}
