#' @title Breaking down the theta matrix obtained from GoM model by different mutational features.
#'
#' @description Breaking down a theta matrix into component probabilities for strand composition (
#'         probability that a mutation comes from positive or negative strand for each cluster
#'         represented by the theta matrix), strand break composition (probability it is A, C, G or T),
#'         the distance from the end of the read and the mutation and flanking base.
#' @param theta_pool A theta matrix from a GoM model fit, with columns representing clusters
#'                   and rows representing the mutational signatures.
#' @param flanking_bases The numbe rof flanking bases to the mutation considered. Defaults to 1.
#' @param max_pos The maximum distance from the end of the read that is used for filtering.
#' @return Returns a list with number of items equal to number of clusters. Each item of this
#'         list is another list which comprises of the individual probability distribution of types of
#'         mutation, flanking base, strand break, distance from end of read and strand break composition
#'         separately.
#' @keywords theta_breakdown
#' @import Logolas, gridBase, grid, ggplot2
#' @export



theta_breakdown = function (theta_pool, flanking_bases = 1, max_pos = 20) {
  signature_set <- rownames(theta_pool)
  signature_patterns <- substring(signature_set, 1, 4+2*flanking_bases)

  indices_minus <- grep("_-_", signature_set)
  strand_theta <- data.frame("minus" = colSums(theta_pool[indices_minus,]),
                             "plus" = colSums(theta_pool[-indices_minus,]))

  breakbase <- substring(signature_set, 8+2*flanking_bases,  8+2*flanking_bases)

  theta_break <- dplyr::tbl_df(data.frame(theta_pool)) %>%
                    dplyr::mutate(sig = breakbase) %>%
                        dplyr::group_by(sig) %>%
                            dplyr::summarise_each(funs(sum)) %>%
                                    as.data.frame()
  rownames(theta_break) <- theta_break[,1]
  theta_break <- theta_break[,-1]

  theta_break <- theta_break[match(c("A", "C", "G", "T"), rownames(theta_break)),]
  breaks_theta <- theta_break


  sub_pattern <- sapply(1:dim(sig_split)[1],
                        function(x) paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse=""))

  new_sig_split <- cbind(sig_split[,1:flanking_bases], sub_pattern, sig_split[,((ncol_sig - flanking_bases +1):ncol_sig)])
  colnames(new_sig_split) = NULL

  if(flanking_bases%%1 != 0){
    stop("flanking bases not evenly distributed")
  }

  theta <- dplyr::tbl_df(data.frame(theta_pool)) %>% dplyr::mutate(sig = signature_set) %>% dplyr::group_by(sig) %>% dplyr::summarise_each(funs(sum)) %>% as.data.frame()
  rownames(theta) <-  theta[,1]
  theta <- theta[,-1]

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

  if(is.null(sig_names))
    sig_names <- rownames(theta)

  #prob_mutation <- filter_by_pos(t(theta_pool), max_pos = max_pos)
  prob_mutation <- filter_signatures_only_location(t(theta_pool), max_pos = max_pos, flanking_bases = flanking_bases)
  prob_mutation <- t(apply(prob_mutation, 1, function(x) {
    y <- x[!is.na(x)];
    return(y/sum(y))
  }))

  ll <- list()

  for(l in 1:dim(theta)[2]){
    ll[[l]] <- list()
    colnames(prop_patterns_list[[l]]) <- c(paste0("flank:", -1), "mismatch", paste0("flank:", 1))
    ll[[l]][["mismatch-flank"]] <- prop_patterns_list[[l]]
    ll[[l]][["strand-orientation"]] <- strand_theta[l,]
    ll[[l]][["strand-break-composition"]] <- breaks_theta[, l, drop=FALSE]
    ll[[l]][["mismatch-trend-read"]] <- prob_mutation[l,]
  }

  return(ll)
}
