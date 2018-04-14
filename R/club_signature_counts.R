#' @title Club signatures in an MFF file to account for strand effect
#'
#' @description Function for clubbing signatures from the two strands to remove
#' strand bias. An example is G->A and C->T are clubbed into C->T.
#'
#' @param signature_counts The matrix of counts of all signatures as produced by
#'                         \code{aggregate_bin_counts()}.
#' @param flanking_bases The number of flanking bases to the mutation. Defaults to 1.
#'
#' @return Returns a matrix of clubbed signatures. The default choice of the mutation
#' signatures are C->T, C->A, C->G, T->A, T->C and T->G. Other mutations are clubbed to
#' above signatures (for example G->A is clubbed with C->T).
#'
#' @keywords club_signature_counts
#'
#' @export
#'


club_signature_counts <- function(signature_counts, flanking_bases=1){

  # leftflank <- grep("left", colnames(signature_pos_str_counts))
  # rightflank <- grep("right", colnames(signature_pos_str_counts))
  # signature_counts_flank <- signature_pos_str_counts[,c(leftflank, rightflank)]
  # signature_counts <- signature_pos_str_counts[, -c(leftflank, rightflank)]

  signature_set <- colnames(signature_counts)
  signature_set_2 <- signatureclub2(signature_set, flanking_bases = flanking_bases)

  # signature_set_3 <- signature_set_2
  #
  # temp <- sapply(signature_set_2[indices], function(x) {
  #                                               bases <- strsplit(as.character(x), split="")[[1]]
  #                                               left_bases <- bases[1:(flanking_bases)]
  #                                               right_bases <- bases[(flanking_bases+4+1):(4+2*flanking_bases)]
  #                                               other_bases <- bases[(4+2*flanking_bases+1):nchar(x)]
  #                                               if(other_bases[2] == "-") {other_bases[2] = "+"} else {other_bases[2] = "-"}
  #                                               newsig <- paste0(c(rev(right_bases), bases[(flanking_bases+1):(flanking_bases+4)], rev(left_bases), other_bases), collapse="")
  #                                               return(newsig)
  #                                          })
  # signature_set_3[indices] <-  temp
  #
  #
  signature_counts_pooled <- do.call(rbind, lapply(1:dim(signature_counts)[1], function(x) tapply(signature_counts[x,], signature_set_2, sum)))
  rownames(signature_counts_pooled) <- rownames(signature_counts)
  temp_split <- do.call(rbind, lapply(colnames(signature_counts_pooled), function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))

  if(length(which(temp_split[,(flanking_bases+1)]=="G")) !=0 || length(which(temp_split[,(flanking_bases+1)]=="A"))!=0){
    stop("G->A conversion did not fully work; aborting")
  }
  return(signature_counts_pooled)
}

gsub2 <- function(pattern, replacement, x, ...) {
  for(i in 1:length(pattern))
    x <- chartr(pattern[i], replacement[i], x, ...)
  x
}

signatureclub2 <- function(signature_set, flanking_bases){
  from <- c('A','T','G','C')
  to <- c('t','a','c','g');
  signature_set_mod <- array(0, length(signature_set));
  for(m in 1:length(signature_set)){
    if(substring(signature_set[m], (1+flanking_bases), (1+flanking_bases)) == "A" | substring(signature_set[m], (1+flanking_bases), (1+flanking_bases)) == "G"){
      temp_split <- strsplit(as.character(signature_set[m]), split="")[[1]]
      bases_flanked <- toupper(gsub2(from, to, substring(signature_set[m], 1, (4+2*flanking_bases))))
      temp_split[1:(4+2*flanking_bases)] <- strsplit(as.character(bases_flanked), split="")[[1]]
      side1 <- temp_split[1:flanking_bases]
      side2 <- temp_split[(5+flanking_bases):(4+2*flanking_bases)]
      temp_split[(5+flanking_bases):(4+2*flanking_bases)] <- rev(side1)
      temp_split[1:flanking_bases] <- rev(side2)
      sign <- temp_split[6+2*flanking_bases]
    }else{
      temp_split <- strsplit(as.character(signature_set[m]), split="")[[1]]
    }
    temp_new <- paste0(temp_split, collapse = "")
    signature_set_mod[m] <- temp_new
  }
  return(signature_set_mod)
}

