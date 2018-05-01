#' @title Grade of Membership (GoM) model clustering of aDNA samples using DNA damage patterns
#'
#' @description Performs GoM model clustering of aDNA samples using DNA damage 
#' patterns- mutation, flanking base, distance from read end and strand 
#' break information. The default implementation of this model follows the 
#' modeling framework of Shiraishi et al (2015). 
#' 
#' 
#' @param K the number of clusters to fit to the model.
#' @param tol The tolerance level of convergence of the GoM model fit
#' @param labs The factor of labels used to group the samples in visualization. 
#'             May be used to distinguish
#'             samples from different labs, or different library prep.
#' @param gom_method The GoM method type. Defaults to \code{independent} model proposed by
#'                   Y. Shiraichi and M. Stephens. The other option is to use the \code{full}
#'                   model which is uses the \code{maptpx} package by Matt Taddy.
#' @param gom.control Control parameters for the GoM model fit.
#' @param output_dir The output directory where the model is saved. 
#'                   If NULL, it picks the current working directory.      
#'                   
#' @return Fits the GoM model on the aggregated data from \code{prepare_archaic}
#'   and then saves the model output as a .RData file. Also outputs model
#'   assessment scores like the BIC. 
#'   
#' @importFrom CountClust compGoM
#' @import maptpx 
#' @export             


model_archaic <- function(dat,
                          K,
                          tol=0.1,
                          labs = NULL,
                          gom_method = "independent",
                          gom.control = list(),
                          output_dir = NULL){
  
  
  gom.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                                 nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                                 light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)
  
  gom.control <- modifyList(gom.control.default, gom.control)
  if(class(dat) == "list"){
    if(is.null(labs)) labs <- c()
    cat("The data is read as a list of matrices - processed by prepare_archaic() \n")
    datalist <- dat
    for(numdir in 1:length(datalist)){
      labs <- c(labs, rep(names(datalist)[numdir], dim(datalist[[numdir]])[1]))
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
      pooled_data[match(rownames(datalist[[num]]), rownames(pooled_data)), 
                  match(colnames(datalist[[num]]), sig_names)] <- as.matrix(datalist[[num]])
    }
  }else if(class(dat) == "character"){
    if(is.null(labs)) labs <- c()
    message("The data is read as names of folders")
    folders <- dat
    datalist <- list()
    for(numdir in 1:length(folders)){
      if(file.exists(paste0(folders[numdir], tail(strsplit(folders[numdir], "/")[[1]],1), ".rda"))){
        datalist[[numdir]] <- get(load(paste0(folders[numdir], tail(strsplit(folders[numdir], "/")[[1]],1), ".rda")))
        cat("Successfully read .RData file from the folder, ", folders[numdir], "CHECK : \n")
      }else{
        message(".RData file not found in folder", folders[numdir], "running 
             prepare_archaic on the MFF files in the folder")
        proc_out <- prepare_aRchaic(folders[numdir])
        datalist[[numdir]] <- proc_out
      }
      labs <- c(labs, rep(tail(strsplit(folders[numdir], "/")[[1]],1), dim(datalist[[numdir]])[1]))
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
      pooled_data[match(rownames(datalist[[num]]), rownames(pooled_data)), 
                  match(colnames(datalist[[num]]), sig_names)] <- as.matrix(datalist[[num]])
    }
  }else if (class(dat) == "matrix" | class(dat) == "data.frame"){
    pooled_data <- dat
    if(is.null(labs)) labs <- 1:nrow(pooled_data)
    if(is.null(levels)) levels <- labs
    if (min( pooled_data %% 1) !=0 ) stop ("dat must be a counts data matrix ")
    if (min(pooled_data) < 0 ) stop ("dat must contain only positive integers ")
  }else{
    stop("class of the dat input must be a character or a vector of characters
          indicating dorectory names, or a 
           matrix of data frame of counts")
  }
  
  if(length(labs) != dim(pooled_data)[1]) stop(paste0("the length of labs is :", length(labs),
          "which does not match the dimension of the data matrix compiled: ", dim(pooled_data)[1]))
  
  zero_sum_rows <- which(rowSums(pooled_data) == 0)
  if(length(zero_sum_rows) > 0){
    pooled_data <- pooled_data[-zero_sum_rows, ]
    labs <- labs[-zero_sum_rows]
  }
  
  pooled_data_2 <- filter_out_strand(pooled_data)
  
  
  if(gom_method == "full"){
    message("Fitting the Grade of Membership Model - full version - due to Matt Taddy")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = pooled_data_2, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
  }else if(gom_method == "independent"){
    message("Fitting the Grade of Membership Model - independent version - due to Y. Shiraichi and M. Stephens")
    
    signature_set <- colnames(pooled_data_2)
    sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
    new_sig_split <- matrix(0, dim(sig_split)[1], 3);
    new_sig_split[,1] <- sig_split[,flanking_bases]
    new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse="")))
    new_sig_split[,3] <- sig_split[,(flanking_bases+5)]
    
    levels(new_sig_split[,1]) <- c("0", "1", "2", "3", "4")
    
    pos <- t(sapply(1:length(signature_set), function(x)
    {
      y = strsplit(signature_set[x], "")[[1]]
      return(paste(y[10:length(y)], collapse=""))
    }))
    
    
    
    mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
    for(k in 1:dim(new_sig_split)[2]){
      temp <- as.factor(new_sig_split[,k])
      mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp),
                                          to = 0:(length(levels(temp))-1))))
    }
    
    pos <- as.numeric(pos)
    pos <- pos - min(pos)
    pos <- factor(pos, levels = 0:21)
    
    signatures <- mat;
    signature_pos <- cbind.data.frame(signatures, pos)
    
    suppressWarnings(topic_clus <- do.call(maptpx::topics, 
        append(list(counts = pooled_data_2, K=K, tol=tol, 
              model = "full", signatures = signatures), gom.control)))
  }
  
  model_assessment <- CountClust::compGoM(pooled_data_2, topic_clus)
  ll <- list("omega" = topic_clus$omega,
             "theta" = topic_clus$theta,
             "assessment" = model_assessment,
             "labs" = labs)
  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}
  save(ll, file = paste0(output_dir, "model.rda"))
  return(ll)
}


filter_out_strand <- function(mat){
  strand_out <- as.character(sapply(as.character(colnames(mat)), function(x) return (paste0(strsplit(x, "_")[[1]][-2], collapse="_"))))
  mat_filtered <- as.numeric()
  for(l in 1:dim(mat)[1]){
    mat_filtered <- rbind(mat_filtered, tapply(mat[l,], strand_out, sum))
  }
  rownames(mat_filtered) <- rownames(mat)
  
  return(mat_filtered)
}

