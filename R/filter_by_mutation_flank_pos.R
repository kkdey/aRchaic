#' @title Filters mutations, flanking bases and position from ends of read from
#' mutational signatures in a matrix of signature counts.
#'
#' @description This function filters the mutations, flanking bases, and position
#'              of mismatch from ends of the read, from mutational signatures and
#'              generates signature counts for the filtered signatures in each sample.
#' @param mat  The matrix of counts of signatures in samples.
#' @param max_pos The maximum distance from ends of the read that is used for filtering.
#'                Defaults to 20.
#' @param flanking_bases The number of flanking bases to the mutation.Defaults to 1
#'
#' @return Returns a new matrix of filtered signatures, where the filtered signatures
#'         represent mutation, flanking bases and position of mutation from ends of the read.
#' @keywords filter_by_mutation_flank_pos
#' @export
#'

filter_by_mutation_flank_pos <-  function(mat, max_pos = 20, flanking_bases=1){

  input_pos <- 1:max_pos
  mutation_flank <- as.character(sapply(as.character(colnames(mat)), function(x) return (strsplit(x,"_")[[1]][1])))
  pos <- as.numeric(sapply(as.character(colnames(mat)), function(x) return (tail(strsplit(x,"_")[[1]],1))))
  pos_indices <- which(!is.na(match(pos, input_pos)))
  mutation_flank_pos <- paste0(mutation_flank, "_", pos)
  mutation_flank_pos_reduced <- mutation_flank_pos[pos_indices]
  mat_reduced <- mat[,pos_indices]

  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat_reduced[l,], mutation_flank_pos_reduced,  sum))
  }
  rownames(mat_filtered) <- rownames(mat)

  return(mat_filtered)

}


