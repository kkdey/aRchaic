#' @title Filters mutation and position from ends of read from
#' mutational signatures in a matrix of signature counts.
#'
#' @description This function filters the mutations and position
#'              of mismatch from ends of the read, from mutational signatures and
#'              generates signature counts for the filtered signatures in each sample.
#' @param mat  The matrix of counts of signatures in samples.
#' @param max_pos The maximum distance from ends of the read that is used for filtering.
#'                Defaults to 20.
#' @param flanking_bases The number of flanking bases to the mutation.Defaults to 1
#'
#' @return Returns a new matrix of filtered signatures, where the filtered signatures
#'         represent mutation and position of mutation from ends of the read.
#' @keywords filter_by_mutation_pos
#' @export
#'

filter_by_mutation_pos <-  function(mat, max_pos = 10, flanking_bases=1){
  input_pos = 1:max_pos
  mutations <- as.character(sapply(as.character(colnames(mat)), function(x) return (paste0(strsplit(x, "")[[1]][(flanking_bases + 1): (flanking_bases + 4)], collapse=""))))
  pos <- as.numeric(sapply(as.character(colnames(mat)), function(x) return (strsplit(x,"_")[[1]][4])))
  indices <- which(!is.na(match(pos, input_pos)))

  mutations_pos <- paste0(mutations, "_", pos)
  mutations_pos_reduced <- mutations_pos[indices]
  mat_reduced <- mat[,indices]

  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat_reduced[l,], mutations_pos_reduced,  sum))
  }
  rownames(mat_filtered) <- rownames(mat)

  return(mat_filtered)
}
