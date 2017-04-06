#' @title Filters out strand break information from signatures in a matrix of signature counts.
#'
#' @description This function filters out strand break information from signatures
#'              and generates signature counts for the filtered signatures in each sample.
#' @param mat  The matrix of counts of signatures in samples.
#' @return Returns a new matrix of filtered signatures, where the filtered signatures
#'         do not contain the strand break information.
#' @keywords filter_out_strand_break
#' @export
#'

filter_out_strand_break <- function(mat){
  strand_break_out <- as.character(sapply(as.character(colnames(mat)), function(x) return (paste0(strsplit(x, "_")[[1]][-3], collapse="_"))))
  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat[l,], strand_break_out, sum))
  }
  rownames(mat_filtered) <- rownames(mat)

  return(mat_filtered)
}
