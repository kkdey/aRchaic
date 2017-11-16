#' @title Filters only mutations from mutational signatures in a matrix of signature counts.
#'
#' @description This function filters the mutations from mutational signatures and generates
#'              signature counts for the mutations.
#' @param mat  The matrix of counts of signatures in samples.
#' @param flanking_bases The number of flanking bases for each mutation. Defaults to 1.
#'
#' @return Returns a new matrix of filtered signatures, where the filtered signatures represent
#'         just the mutation types (6 in this case C->T, C->A, C->G, T->A, T->C and T->G). The
#'         other types of mutations are clubbed with above 6. (see \code{club_signature_counts} for more
#'         details on clubbing)
#' @keywords filter_by_mutation
#' @export
#'


filter_by_mutation <-  function(mat, flanking_bases=1){
  mutations <- as.character(sapply(as.character(colnames(mat)), function(x) return (paste0(strsplit(x, "")[[1]][(flanking_bases + 1): (flanking_bases + 4)], collapse=""))))

  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat[l,], mutations, sum))
  }
  rownames(mat_filtered) <- rownames(mat)

  return(mat_filtered)
}
