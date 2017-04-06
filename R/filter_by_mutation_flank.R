#' @title Filters mutations and flanking bases from mutational signatures in a
#' matrix of signature counts.
#'
#' @description This function filters the mutations and flanking basesfrom mutational signatures and
#'              generates signature counts for the filtered signatures in each sample.
#' @param mat  The matrix of counts of signatures in samples.
#'
#' @return Returns a new matrix of filtered signatures, where the filtered signatures represent
#'         mutation and flanking bases.
#' @keywords filter_by_mutation_flank
#' @export
#'


filter_by_mutation_flank <-  function(mat){
  mutation_flank <- as.character(sapply(as.character(colnames(mat)), function(x) return (strsplit(x,"_")[[1]][1])))

  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat[l,], mutation_flank, sum))
  }
  rownames(mat_filtered) <- rownames(mat)

  return(mat_filtered)
}
