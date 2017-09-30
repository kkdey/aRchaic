#' @title Bin signatures and aggregate counts of signatures from an MFF file.
#'
#' @description For each file in directory \code{dir}, aggregates the data from MFF files
#' into bins determined by distance from ends of the read based on \code{breaks} provided.
#' This function is an Internal data processing function for \code{aRchaic_cluster} and related
#' fucntions.
#'
#' @param dir The directory containing the MFF files.
#' @param breaks The breaks used for binning the distance from ends of the reads of the mutations.
#' @param flanking_bases The number of flanking bases. Defaults to 1.
#' @param output_rda If non-NULL, the output matrix is stored in the path provided (should end with .rda).
#'
#' @return Returns a matrix with rows being the samples (each MFF file), columns representing
#' the binned mutational signatures, where bins are determined by \code{breaks} in user input.
#' The cells contain counts of the number of mutational signatures observed in that MFF file.
#'
#' @keywords aggregate_counts
#' @export


aggregate_signature_counts <- function(dir,
                                       pattern = NULL,
                                       breaks = NULL,
                                       flanking_bases = 1,
                                       output_rda = NULL){
  if(is.null(pattern)){
    ancient_files <- setdiff(list.files(dir, pattern = ".csv"), list.files(dir, pattern = ".csv#"))
  }else{
    ancient_files <- list.files(dir, pattern = pattern)
  }

  signature_ancient <- vector(mode="list")
  signature_counts_ancient <- vector(mode="list")


  for(num in 1:length(ancient_files)){
    tmp_dat <- damage_build_bin_counts(paste0(dir, ancient_files[num]),
                                       breaks=breaks,
                                       type=2)
    signature_counts_ancient[[num]] <- tmp_dat[,2];
    signature_ancient[[num]] <- as.character(tmp_dat[,1]);
    cat("Reading file ", paste0(dir, ancient_files[num]), "\n")
  }

  merged_signature_ancient <- signature_ancient[[1]]

  if(length(ancient_files) >= 2){
    for(num in 2:length(ancient_files)){
      merged_signature_ancient <- union(merged_signature_ancient, signature_ancient[[num]])
    }
  }

  ancient_counts <- matrix(0, length(ancient_files), length(merged_signature_ancient))

  for(num in 1:length(ancient_files)){
    ancient_counts[num, match(signature_ancient[[num]], merged_signature_ancient)] <- signature_counts_ancient[[num]]
  }

  ancient_files_filt <- as.character(sapply(ancient_files, function(x) strsplit(x, ".csv")[[1]][1]))

  rownames(ancient_counts) <- ancient_files_filt

  signature_split <- do.call(rbind, lapply(merged_signature_ancient, function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases+6)]))

  indices1 <- which(signature_split[,(flanking_bases+1)]==signature_split[,(flanking_bases+4)])

  wrong_letters <- c("B", "D", "E", "F", "H", "I", "J", "K", "L", "M", "N", "O",
                     "P", "Q", "R", "S", "U", "V", "W", "X", "Y", "Z")
  temp <- list()
  for(l in 1:length(wrong_letters)){
    temp[[l]] <- grep(paste0(wrong_letters[l]), merged_signature_ancient)
  }

  indices2 <- Reduce(union, temp)

  # temp1 <- grep("N", sign)
  # indices2 <- numeric()
  # for(m in 1:(4+2*flanking_bases)){
  #   temp1 <- gre
  #   indices2 <- c(indices2, which(signature_split[,m]=="N" | signature_split[,m]=="R"));
  # }
  #
  # indices3 <- numeric()
  # indices3 <- c(indices3, which(signature_split[,(4+2*flanking_bases + 4)]=="N" | signature_split[,(4+2*flanking_bases + 4)]=="R"))
  # indices3 <- c(indices3, which(signature_split[,(4+2*flanking_bases + 6)]=="N" | signature_split[,(4+2*flanking_bases + 4)]=="R"))
  #

  # indices <- union(indices1, union(indices2, indices3))

  indices <- union(indices1, indices2)

  if(length(indices) > 0) {
    ancient_counts_filtered <- ancient_counts[,-indices]
  }else{
    ancient_counts_filtered <- ancient_counts
  }

  ancient_counts_filtered <- matrix(ancient_counts_filtered, nrow = nrow(ancient_counts))

  rownames(ancient_counts_filtered) <- ancient_files_filt
  if(length(indices) > 0){
    colnames(ancient_counts_filtered) <- merged_signature_ancient[-indices]
  }else{
    colnames(ancient_counts_filtered) <- merged_signature_ancient
  }

  if(is.null(output_rda)){
    return(ancient_counts_filtered)
  }else{
    save(ancient_counts_filtered, output_rda)
  }
}

damage_build_bin_counts =  function(file,
                                    breaks=NULL,
                                    type=1)
{
  file <-  read.csv(file=file, header=TRUE);
  file <- file[,-7] ## remove the read name

  colnames(file) <- c("mut", "leftpos", "rightpos", "lsb", "rsb", "strand")
  library(dplyr)
  tab <- dplyr::tbl_df(file) %>% dplyr::group_by(mut, leftpos, rightpos, lsb, rsb, strand) %>% dplyr::summarize(n= n())
  tab2 <- as.data.frame(tab)

  if(length(which(tab2[,3]==-1)) !=0){
    tab2[which(tab2[,3]==-1), 3] <- 0
  }
  if(length(which(tab2[,2]==-1)) !=0){
    tab2[which(tab2[,2]==-1), 2] <- 0
  }

  min_dist_from_end <- as.numeric(apply(tab2[,2:3], 1, function(x) return(min(x))))


  if(is.null(breaks)){
    bins <- c(0, 5, 10, 15, 20, max(min_dist_from_end))
  }else{
   # message("breaks values provided : adjusting")
    bins <- c(intersect(breaks, (-1):max(min_dist_from_end)),max(min_dist_from_end))
  }

  bin_values <- .bincode(min_dist_from_end, bins, TRUE)

  modified_file <- cbind.data.frame(tab2[,1], tab2[,4], tab2[,5], tab2[,6], tab2[,7], bin_values)
  colnames(modified_file) <- c("pattern", "base.lsb", "base.rsb", "strand", "counts", "bin_values")

  library(plyr)
  df1 <- plyr::ddply(modified_file, .(pattern, strand, base.lsb, base.rsb, strand, bin_values), plyr::summarise, newvar = sum(counts))
  colnames(df1) <- c("pattern", "strand", "base.lsb", "base.rsb", "bin", "counts")

  if(type==2){
    df2 <- cbind.data.frame(paste0(df1[,1], "_", df1[,2], "_", df1[,3], "_", df1[,4], "_", df1[,5]), df1[,6])
    colnames(df2) <- c("pattern-strand-breaks-bin", "counts")
    out <- df2
  }else{
    out <- df1
  }
}

