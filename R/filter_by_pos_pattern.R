#' @title Filters position of a specific type of mutation from ends of read in a matrix of signature counts.
#'
#' @description This function filters position of a specific type of mismatch
#'              given by \code{pattern} (for example C->T) from ends of the read,
#'              and generates signature counts for the filtered signatures in each sample.
#' @param mat  The matrix of counts of signatures in samples.
#' @param max_pos The maximum distance from ends of the read that is used for filtering.
#'                Defaults to 20.
#' @param pattern The pattern of the mutation type. Defaults to C->T.
#' @param flanking_bases the number of flanking bases for the mutation.
#'
#' @return Returns a new matrix of filtered signatures, where the filtered signatures
#'         represent position of a specific type of mutational type from ends of the read.
#' @keywords filter_by_pos_pattern
#' @export
#'


filter_by_pos_pattern <-  function(mat, max_pos = 20, pattern = "C->T", flanking_bases=1){

  input_pos <- 1:max_pos
  mutations <- as.character(sapply(as.character(colnames(mat)), function(x) return (paste0(strsplit(x, "")[[1]][(flanking_bases + 1): (flanking_bases + 4)], collapse=""))))
  CtoTindices <- which(mutations == pattern)

  pos <- as.numeric(sapply(as.character(colnames(mat)), function(x) return (tail(strsplit(x,"_")[[1]],1))))
  pos_indices <- which(!is.na(match(pos, input_pos)))

  matched_indices <- intersect(CtoTindices, pos_indices)

  mutations_pos <- paste0(mutations, "_", pos)
  mutations_pos_reduced <- mutations_pos[matched_indices]
  mat_reduced <- mat[,matched_indices]

  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat_reduced[l,], mutations_pos_reduced, sum))
  }
  rownames(mat_filtered) <- rownames(mat)
  return(mat_filtered)
}
