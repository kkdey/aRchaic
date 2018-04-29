#' @title Function to extract signature counts for position along the read
#'
#' @description This function can be used by the user to extract the counts data for positions of the mutation signatures
#' along the read length.
#'
#' @param mat The matrix of the signature counts obtained by \code{aggregate_bin_counts} or
#'  \code{club_signature_counts} functions.
#' @param max_pos The highest position along the read from the start for which the counts are to be extracted.
#' @param flanking_bases The number of flanking bases in the signatures
#'
#' @return Returns a filtered matrix of signature counts for positions along the read length.
#' @keywords filter-signatures
#' @export
#'


filter_signatures_only_location <-  function(mat, max_pos = 20, flanking_bases=1){
  input_pos <- 0:max_pos;
  pos <- as.numeric(sapply(as.character(colnames(mat)), function(x) return (tail(strsplit(x, "_")[[1]], 1))))
  reduced_dat <- mat[,which(!is.na(match(pos, input_pos)))]
  pos2 <- pos[which(!is.na(match(pos, input_pos)))];
  pos2fac <- factor(pos2, levels=input_pos)

  pos_mat <- as.numeric();
  for(m in 1:dim(mat)[1]){
    pos_mat <- rbind(pos_mat, tapply(reduced_dat[m,], pos2fac, sum));
  }
  return(pos_mat)
}
