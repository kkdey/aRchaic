#' @title Filters position of mutation from ends of read from
#' mutational signatures in a matrix of signature counts.
#'
#' @description This function filters position of mismatch from ends of the read,
#'              from mutational signatures and generates signature counts for the filtered
#'              signatures in each sample.
#' @param mat  The matrix of counts of signatures in samples.
#' @param max_pos The maximum distance from ends of the read that is used for filtering.
#'                Defaults to 20.
#'
#' @return Returns a new matrix of filtered signatures, where the filtered signatures
#'         represent position of mutation from ends of the read.
#' @keywords filter_by_pos
#' @export
#'

filter_by_pos <- function(mat, max_pos = 20){
  input_pos <- 1:max_pos
  pos <- as.numeric(sapply(as.character(colnames(mat)), function(x) return (tail(strsplit(x,"_")[[1]],1))))
  pos_indices <- which(!is.na(match(pos, input_pos)))
  pos_reduced <- pos[pos_indices]

  mat_reduced <- mat[,pos_indices]
  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat_reduced[l,], pos_reduced, sum))
  }
  rownames(mat_filtered) <- rownames(mat)

  return(mat_filtered)
}
