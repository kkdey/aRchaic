#' @title Aggregating signature counts from MFF files in one or more folders.
#'
#' @description Aggregate the signatures into bins for each MFF file in each of the folders
#' listed in the vector \code{folders} and then create a matrix of counts of signatures, where
#' each row corresponds to one sample MFF file and each column to one signature (binned from
#' \code{aggregate_bin_counts ()}).
#'
#' @param folders a vector of folder names hosting the MFF files. Each folder may represent
#'                MFF files from same source sequenced and analysed together. (see vignette for
#'                examples)
#' @param pattern If the user wants to incorporate specific files in the matrix and not all in
#'                files in the folder, she can put the names as a vector in \code{pattern}.
#' @param breaks The breaks used for binning the distance from ends of the reads of the mutations, when
#'               processing the MFF files.
#' @param flanking_bases The numbe rof flanking bases to the mutation considered. Defaults to 1.
#'
#' @return The function creates a matrix of counts of MFF samples and the mutational signatures.
#'         Also if the processed files in a folder are not present, it creates them as .rda
#'         files in the folders.
#' @keywords aRchaic_pool
#' @import e1071, classtpx
#' @export



aRchaic_pool = function(folders,
                        pattern = NULL,
                        breaks = c(-1, seq(1,20,1)),
                        flanking_bases = 1){

  message("Checking if the folders exist")

  if(!is.null(pattern)){

  }

  for(i in 1:length(folders)){
    if(!file.exists(folders[i]))
      stop("A folder in the folder list does not exist:  aborting")
  }

  datalist <- vector("list", length(folders))

  for(i in 1:length(folders)){
    if(!file.exists(paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda"))){
      message (paste0("Processing the MutationFeatureFormat files in the directory", folders[i]))
      out <- aggregate_signature_counts(dir = paste0(folders[i]),
                                          pattern = NULL,
                                          breaks = breaks,
                                          flanking_bases = flanking_bases)
        clubbed_data <- club_signature_counts(out, flanking_bases = flanking_bases)
        save(clubbed_data, file = paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda"))
        datalist[[i]] <- clubbed_data
    }else{
        datalist[[i]] <- get(load(paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda")))
    }
  }

  sig_names <- colnames(datalist[[1]])
  row_names_pool <- rownames(datalist[[1]])
  if(length(datalist) >= 2){
    for(num in 2:length(datalist)){
      sig_names <- union(sig_names, colnames(datalist[[num]]))
      row_names_pool <- c(row_names_pool, rownames(datalist[[num]]))
    }
  }

  pooled_data <- matrix(0, length(row_names_pool), length(sig_names))
  rownames(pooled_data) <- row_names_pool
  colnames(pooled_data) <- sig_names

  for(num in 1:length(datalist)){
    pooled_data[match(rownames(datalist[[num]]), rownames(pooled_data)), match(colnames(datalist[[num]]), sig_names)] <- datalist[[num]]
  }

  if(is.null(pattern)){
    return(pooled_data)
  }

  if(!is.null(pattern)){
    temp <- as.numeric()
    for(l in 1:length(pattern)){
      temp <- rbind(temp, pooled_data[grep(pattern = pattern[l], paste0(rownames(pooled_data), ".csv")),, drop=FALSE])
    }
    return(temp)
  }

}
