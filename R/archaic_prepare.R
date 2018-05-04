#' @title Prepare and aggregate the mismatch feature counts from one or multiple
#' MFF files
#'
#' @description For each file in the vector of directories \code{dirs},
#' aggregates the mismatch feature counts data from the MFF files
#'
#' @param dirs The directory/directories containing the MFF files.
#' @param max_pos The maximum position from the ends of the reads for which
#'  mismatches are considered in aRchaic.
#' @param output_rda If non-NULL, the processed data for each
#'  directory in \code{dirs} is saved as a .Rdata file.
#'
#' @return Returns a matrix with rows being the samples (each MFF file),
#' columns representing the mismatch signatures (comprising of features like
#' mismatch type, flanking bases and strand break information).
#' The cells contain counts of the number of mutational signatures observed in
#' that MFF file.
#'
#' @keywords aggregate_counts
#' @export

archaic_prepare <- function(dirs,
                            max_pos = 20,
                            from_scratch = FALSE,
                            delete = FALSE,
                            output_rda = TRUE){
  if(class(dirs) != "character") stop("the dirs input must be names of the folders - hence of class character")
    for(numdir in 1:length(dirs)){
      if(regmatches(dirs[numdir],regexpr(".$", dirs[numdir])) != "/") dirs[numdir] <- paste0(dirs[numdir], "/")
    }
    dat_list  <- list()
    for(numdir in 1:length(dirs)){
      model_filename <- paste0(tail(strsplit(dirs[numdir], "/")[[1]],1), ".rda")
      csvnames <- list.files(dirs[numdir], pattern = ".csv")

####################   what if no CSV -MFF files are present in directories? ######################################

      if(length(csvnames) == 0 & file.exists(paste0(dirs[numdir], model_filename))){
          cat("no csv files present in directory, ", dirs[numdir], " to process,
              but the .RData file is present, loading and
              returning the .RData file \n")
          dat_list[[numdir]] <- get(load(paste0(dirs[numdir], model_filename)))
      }

      if(length(csvnames) == 0 & !file.exists(paste0(dirs[numdir], model_filename))){
          cat("no csv or .RData files present in the directory ", dirs[numdir], " : returning NULL, \n")
        dat_list[[numdir]] <- NULL
      }


###################  what if no .RData file is present and CSV files are present  #########################

      if(length(csvnames) > 0 & !file.exists(paste0(dirs[numdir], model_filename)) | length(csvnames) > 0 & file.exists(paste0(dirs[numdir], model_filename)) & from_scratch){
         files_to_process <- csvnames
         csvnames2 <- substr(csvnames, 1, nchar(csvnames) - 4)
         signature_counts_from_file <- vector(mode="list")
         signatures_file <- vector(mode="list")


         for(numfile in 1:length(files_to_process)){
           tmp_dat <- damage_build_bin_counts(file = paste0(dirs[numdir], files_to_process[numfile]),
                                              max_pos = max_pos,
                                              type=2)
           signature_counts_from_file[[numfile]] <- tmp_dat[,2]
           signatures_file[[numfile]] <- tmp_dat[,1]
           cat("Reading file ", paste0(dirs[numdir], files_to_process[numfile]), "\n")
         }

         merged_signatures <- signatures_file[[1]]

         if(length(files_to_process) >= 2){
           for(num in 2:length(files_to_process)){
             merged_signatures <- union(merged_signatures, signatures_file[[num]])
           }
         }

         ancient_counts <- matrix(0, length(files_to_process),
                                  length(merged_signatures))
         for(num in 1:length(files_to_process)){
           ancient_counts[num, match(signatures_file[[num]],
                                     merged_signatures)] <- signature_counts_from_file[[num]]
         }

         flanking_bases = 1
         signature_split <- do.call(rbind, lapply(merged_signatures,
                                                  function(x) strsplit(as.character(x),
                                                                       split="")[[1]][1:(4+2*flanking_bases+6)]))
         indices1 <- which(signature_split[,(flanking_bases+1)]==signature_split[,(flanking_bases+4)])
         wrong_letters <- c("B", "D", "E", "F", "H", "I", "J", "K", "L", "M", "N", "O",
                            "P", "Q", "R", "S", "U", "V", "W", "X", "Y", "Z")
         temp <- list()
         for(l in 1:length(wrong_letters)){
           temp[[l]] <- grep(paste0(wrong_letters[l]), merged_signatures)
         }

         indices2 <- Reduce(union, temp)
         indices <- union(indices1, indices2)
         ancient_counts_filtered <- matrix(ancient_counts[, -indices],
                                           nrow = nrow(ancient_counts))

         rownames(ancient_counts_filtered) <- csvnames2
         if(length(indices) > 0){
           colnames(ancient_counts_filtered) <- merged_signatures[-indices]
         }else{
           colnames(ancient_counts_filtered) <- merged_signatures
         }

         outdat <- club_signature_counts(ancient_counts_filtered)

         dat_list[[numdir]] <- outdat
         if(output_rda){
           save(outdat, file = paste0(dirs[numdir],
                          tail(strsplit(dirs[numdir], "/")[[1]],1), ".rda"))
         }
         next
      }


##################  if CSV files and .RData file both are present, and not scratch   #################################

      if(length(csvnames) > 0 & file.exists(paste0(dirs[numdir], model_filename)) & !from_scratch){
          model <- get(load(paste0(dirs[numdir], model_filename)))
          csvnames2 <- substr(csvnames, 1, nchar(csvnames) - 4)
          deleted_samples <- setdiff(rownames(model), csvnames2)
          if(delete){
            if(length(deleted_samples) > 0){
              row_ids_to_delete <- match(deleted_samples, rownames(model))
              cat("Deleting files : ", deleted_samples, " from the .RData file \n")
              model <- model[-row_ids_to_delete,]
            }
          }else{
            if(length(deleted_samples) > 0){
              message("some of the row names in model file do not appear as
                      csv files: may be they are from older files , if you
                      want to remove them, try delete = TRUE")
            }}
          matched_ids <- match(csvnames2, rownames(model))
          absent_ids <- which(is.na(matched_ids))

          if(length(absent_ids) == 0 & length(deleted_samples) > 0){
            dat_list[[numdir]] <- model
            if(output_rda){
              save(model, file = paste0(dirs[numdir],
                                            tail(strsplit(dirs[numdir], "/")[[1]],1), ".rda"))
            }
          }

          if(length(absent_ids) == 0 & length(deleted_samples) == 0){
            dat_list[[numdir]] <- model
            if(output_rda){
              save(model, file = paste0(dirs[numdir],
                                            tail(strsplit(dirs[numdir], "/")[[1]],1), ".rda"))
            }
          }

          if(length(absent_ids) > 0){
          files_to_process <- csvnames[absent_ids]
          signature_counts_from_file <- vector(mode="list")
          signatures_file <- vector(mode="list")


          for(numfile in 1:length(files_to_process)){
            tmp_dat <- damage_build_bin_counts(file = paste0(dirs[numdir], files_to_process[numfile]),
                                               max_pos = max_pos,
                                               type=2)
            signature_counts_from_file[[numfile]] <- tmp_dat[,2]
            signatures_file[[numfile]] <- tmp_dat[,1]
            cat("Reading file ", paste0(dirs[numdir], files_to_process[numfile]), "\n")
          }

          merged_signatures <- signatures_file[[1]]

          if(length(files_to_process) >= 2){
            for(num in 2:length(files_to_process)){
              merged_signatures <- union(merged_signatures, signatures_file[[num]])
            }
          }

          ancient_counts <- matrix(0, length(files_to_process),
                                   length(merged_signatures))
          for(num in 1:length(files_to_process)){
            ancient_counts[num, match(signatures_file[[num]],
                                      merged_signatures)] <- signature_counts_from_file[[num]]
          }

          flanking_bases = 1
          signature_split <- do.call(rbind, lapply(merged_signatures,
                                                   function(x) strsplit(as.character(x),
                                                                        split="")[[1]][1:(4+2*flanking_bases+6)]))
          indices1 <- which(signature_split[,(flanking_bases+1)]==signature_split[,(flanking_bases+4)])
          wrong_letters <- c("B", "D", "E", "F", "H", "I", "J", "K", "L", "M", "N", "O",
                             "P", "Q", "R", "S", "U", "V", "W", "X", "Y", "Z")
          temp <- list()
          for(l in 1:length(wrong_letters)){
            temp[[l]] <- grep(paste0(wrong_letters[l]), merged_signatures)
          }

          indices2 <- Reduce(union, temp)
          indices <- union(indices1, indices2)
          ancient_counts_filtered <- matrix(ancient_counts[, -indices],
                                            nrow = nrow(ancient_counts))

          rownames(ancient_counts_filtered) <- csvnames2[absent_ids]
          if(length(indices) > 0){
            colnames(ancient_counts_filtered) <- merged_signatures[-indices]
          }else{
            colnames(ancient_counts_filtered) <- merged_signatures
          }

          outdat <- club_signature_counts(ancient_counts_filtered)

          total_sigs <- union(colnames(model), colnames(outdat))
          final_dat <- matrix(0, nrow(model) +  nrow(outdat), length(total_sigs))
          final_dat[1:nrow(model), match(colnames(model), total_sigs)] = model
          final_dat[(nrow(model)+1):(nrow(model)+nrow(outdat)),
                    match(colnames(outdat), total_sigs)] = outdat
          rownames(final_dat) = c(rownames(model), rownames(ancient_counts_filtered))
          colnames(final_dat) = total_sigs

          dat_list[[numdir]] <- final_dat

          if(output_rda){
            save(final_dat, file = paste0(dirs[numdir],
                             tail(strsplit(dirs[numdir], "/")[[1]],1), ".rda"))
          }

        }
      }
    }
    names(dat_list) <- sapply(dirs, function(x) return(tail(strsplit(x, "/")[[1]],1)))
    return(dat_list)
}


damage_build_bin_counts =  function(filename,
                                    max_pos = 20,
                                    type=2)
{
  breaks = c(-1, seq(1,max_pos,1))
  file <- read.csv(filename, header = FALSE, stringsAsFactors = FALSE)
  file[,7] <- rep(1, dim(file)[1])
  #file <- file[,-7] ## remove the read name

  colnames(file) <- c("mut", "leftpos", "rightpos", "lsb", "rsb", "strand", "counts")

  library(dplyr)
  tab <- dplyr::tbl_df(file) %>% dplyr::group_by(mut, leftpos, rightpos, lsb, rsb, strand) %>% dplyr::summarize(n= n())
  tab2 <- as.data.frame(tab)
  min_dist_from_end_pre <- as.numeric(apply(tab2[,2:3], 1, function(x) return(min(x))))
  idx <- which(min_dist_from_end_pre <= tail(breaks, 1))
  tab2 <- tab2[idx,]
  min_dist_from_end <- min_dist_from_end_pre[idx]

  if(length(which(tab2[,3]==-1)) !=0){
    tab2[which(tab2[,3]==-1), 3] <- 0
  }
  if(length(which(tab2[,2]==-1)) !=0){
    tab2[which(tab2[,2]==-1), 2] <- 0
  }

  bases <- apply(tab2, 1, function(x) {
    idx <- which.min(c(x[2],x[3]))
    if(idx == 1)
      return (x[4])
    else if (idx == 2)
      return (gsub2("ACGT", "TGCA", x[5]))
  })

  breakbases <- bases
  modified_file <- cbind.data.frame(tab2[,1], breakbases, tab2[,6], tab2[,7], min_dist_from_end)
  colnames(modified_file) <- c("pattern", "breakbase", "strand", "counts", "bin_values")

  library(plyr)
  df1 <- plyr::ddply(modified_file, .(pattern, strand, breakbases, strand, bin_values), plyr::summarise, newvar = sum(as.numeric(counts)))
  colnames(df1) <- c("pattern", "strand", "breakbase", "bin", "counts")

  if(type==2){
    df2 <- cbind.data.frame(paste0(df1[,1], "_", df1[,2], "_", df1[,3], "_", df1[,4]), df1[,5])
    colnames(df2) <- c("pattern-strand-breaks-bin", "counts")
    out <- df2
  }else{
    out <- df1
  }
}

club_signature_counts <- function(signature_counts, flanking_bases=1){

  signature_set <- colnames(signature_counts)
  signature_set_2 <- signatureclub2(signature_set,
                                    flanking_bases = flanking_bases)
  signature_counts_pooled <- do.call(rbind, lapply(1:dim(signature_counts)[1],
              function(x) tapply(signature_counts[x,], signature_set_2, sum)))
  rownames(signature_counts_pooled) <- rownames(signature_counts)
  temp_split <- do.call(rbind, lapply(colnames(signature_counts_pooled),
  function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))
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
    }else{
      temp_split <- strsplit(as.character(signature_set[m]), split="")[[1]]
    }
    temp_new <- paste0(temp_split, collapse = "")
    signature_set_mod[m] <- temp_new
  }
  return(signature_set_mod)
}


