#' @title Classification (model-based and non-model based) of aDNA samples using
#' filtered DNA damage patterns
#'
#' @description Performs classification of aDNA samples using filtered DNA damage patterns,
#' which comprises of some subset of the 5 features (mutation,
#' flanking base, distance from read end, strand and strand break) but not all of them.
#' We provide two classification options - SVM and classification GoM (classtpx). We fit
#' soft classification for both cases, thereby obtaining probabilities of belonging to the
#' classes for each test sample.
#'
#' @param folders a vector of folder names hosting the MFF files. Each folder may represent
#'                MFF files from same source sequenced and analysed together. (see vignette for
#'                examples)
#' @param type The type of filtering mechanism used. The possible options are
#'                   - \code{mutation} : just the mutations (C->T, C->G, C->A, T->A, T->C and T->G)
#'                   - \code{mutation-flank}: the mutation and flanking base features
#'                   - \code{mutation-pos}: the mutation and position of the mutation from the nearest end of read.
#'                   - \code{mutation-flank-pos}: mutation + flanking base + position from nearest end of read
#'                   - \code{specific-mutation-pos}: specific mutation (say C->T) + position of that mutation from nearest end of read
#'                   - \code{wo-strand}: all features except strand feature
#'                   - \code{wo-strand-break}: all features except the strand break feature
#' @param class_labs The class labels for the samples. Takes values from 1 to K (where K is the
#'                   number of classes) and NA. The NA values correspond to test samples whereas
#'                   the others are training samples assigned to a class between 1 and K.
#' @param class_method The method used for performing supervised learning. The options are \code{SVM}
#'                     and \code{classtpx}.
#' @param run_from Can take one of two options - \code{class} and \code{start}.  For \code{start},
#'                 the function will perform all steps - processing the MFF files and then
#'                 classification. for \code{class}, the function will first check
#'                 if the processed files are present already from previous run or otherwise.
#'                  If present, it will use the processed files to run the classification model.
#' @param run_index The index vector of files to be included for each folder in the vector \code{folders}.
#'                  Defaults to using all files in the folder.
#' @param breaks The breaks used for binning the distance from ends of the reads of the mutations, when
#'               processing the MFF files. Defaults to assigning each base as one bin from
#'               end of the read to 20 bases.
#' @param flanking_bases The numbe rof flanking bases to the mutation considered. Defaults to 1.
#' @param pattern type of mutation (say "C->T"). Comes into play when \code{type} is \code{specific-mutation-pos}.
#'                Defaults to NULL.
#' @param max_pos The maximum distance from ends of read that a mismatch would be considered.
#'                The rest will be filtered out.
#' @param normalize If SVM option is chosen by the user, he can run it on the original data (for \code{normalize = FALSE})
#'                  or on log CPM normalized data (for \code{normalize = TRUE})
#' @param svm.control Control parameters for the SVM function in e1071.
#' @param classtpx.control Control parameters for the supervised GoM model application.
#'
#' @return For the \code{start} option of \code{run_from}, the function processes the MFF files
#' in each folder, aggregates them into a matrix and saves it in the folder as a .rda file.
#' Then it filters the signatures based on the \code{type} selected. Subsequently,
#' it fits SVM or classtpx to the aggregated data from one or more folders as determined by
#' the vector \code{folders} and generates the probabilities of belonging to each class for the test
#' samples (which correspond to those samples that have NA in \code{class_labs}).
#'
#' @keywords aRchaic_class_beta
#' @export


aRchaic_class_beta = function(folders,
                              type = c("mutation", "mutation-flank", "mutation-pos","mutation-flank-pos",
                                  "specific-mutation-pos", "wo-strand", "wo-strand-break"),
                              class_labs = NULL,
                              class_method = c("SVM", "classtpx"),
                              run_from = "class",
                              run_index = NULL,
                              breaks = c(-1, seq(1,20,1)),
                              flanking_bases =1,
                              pattern = NULL,
                              max_pos = 20,
                              normalize = TRUE,
                              svm.control = list(),
                              classtpx.control = list()){

  classtpx.control.default <- list(method="theta.fix", shrink=FALSE,
                                   shrink.method = 1, tol=0.001,
                                   ord=FALSE)
  svm.control.default <- list(scale = TRUE, type = NULL, kernel ="radial",
                              degree = 3,
                              coef0 = 0, cost = 1, nu = 0.5,
                              class.weights = NULL, cachesize = 40, tolerance = 0.001, epsilon = 0.1,
                              shrinking = TRUE, cross = 0, fitted = TRUE)


  classtpx.control <- modifyList(classtpx.control.default, classtpx.control)
  svm.control <- modifyList(svm.control.default, svm.control)

  if(is.null(class_labs)){
    stop("Classification labels not provided: please provide a vector of labels in class_labs")
  }

  if(run_from == "start"){
    if(is.null(run_index)){
      run_index <- 1:length(folders)
    }
    if(sum(run_index - 1:length(folders))^2 == 0){
      for(i in 1:length(folders)){
        file.remove(paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda"))
      }
    }else{
      folders1 <- folders[run_index]
      for(i in 1:length(folders1)){
        file.remove(paste0(folders1[i], tail(strsplit(folders1[i], "/")[[1]],1), ".rda"))
      }
    }
  }

  datalist <- vector("list", length(folders))

  for(i in 1:length(folders)){
    if(!file.exists(paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda"))){
      message (paste0("Processing the MutationFeatureFormat files in the directory", folders[i]))
      out <- aggregate_signature_counts(dir = paste0(folders[i]),
                                        pattern = NULL,
                                        breaks = breaks,
                                        flanking_bases = flanking_bases)
      clubbed_data <- club_signature_counts(out, flanking_bases = 1)
      save(clubbed_data, file = paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda"))
      datalist[[i]] <- clubbed_data
    }else{
      datalist[[i]] <- get(load(paste0(folders[i], tail(strsplit(folders[i], "/")[[1]],1), ".rda")))
    }
  }

  ######################  Pooling data from different folders if present  ############################################

  if(run_from == "start"){
    message("Pooling the data from multiple sources")
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

  mat <- pooled_data

  if(type == "mutation"){
    mat_reduced <- filter_by_mutation(mat)
  }else if (type == "mutation-pos"){
    mat_reduced <- filter_by_mutation_pos(mat, max_pos, flanking_bases)
  }else if (type == "mutation-flank"){
    mat_reduced <- filter_by_mutation_flank(mat)
  }else if (type == "mutation-flank-pos"){
    mat_reduced <- filter_by_mutation_flank_pos(mat, max_pos, flanking_bases)
  } else if (type == "specific-mutation-pos"){
    if(is.null(pattern)){
      stop("for this option, user needs to provide a mutation type: C->T, C->A C->G, T->A, T->C or T->G")
    }
    mat_reduced <- filter_by_pos_pattern(mat, max_pos, pattern = pattern, flanking_bases)
  }else if (type == "wo-strand"){
    mat_reduced <- filter_out_strand(mat)
  }else if (type == "wo-strand-break"){
    mat_reduced <- filter_out_strand_break(mat)
  }else {
    stop("the type of filtering does not match with the possible options: see documentation")
  }

  rownames(mat_reduced) <- rownames(mat)

  pooled_data <- mat_reduced

  if(length(class_labs) != dim(pooled_data)[1]){
    stop("number of samples not same as number of class labels, aborting")
  }

  test_indices <- which(is.na(class_labs))
  if(length(test_indices)==0){
    stop("no test sample provided, aborting")
  }

  train_indices <- which(!is.na(class_labs))
  if(length(train_indices)==0){
    stop("no training sample provided, aborting")
  }


  if (class_method == "SVM"){
    if(normalize){
      voom_pooled_data <- t(limma::voom(t(pooled_data))$E);
      trainX <- voom_pooled_data[train_indices,]
      testX <- voom_pooled_data[test_indices,]
      y = factor(class_labs[train_indices])
      data <- cbind.data.frame(y, trainX);
      library(e1071)
      model_SVM <- do.call(e1071::svm, append(list(formula = y ~ .,
                                                   data=data,
                                                   probability=TRUE), svm.control))
      prob  = predict(model_SVM, testX, probability=TRUE)
    }else{
      trainX <- pooled_data[train_indices,]
      testX <- pooled_data[test_indices,]
      y = factor(class_labs[train_indices])
      data <- cbind.data.frame(y, trainX);
      library(e1071)
      model_SVM <- do.call(e1071::svm, append(list(formula = y ~ .,
                                                   data=data,
                                                   probability=TRUE), svm.control))

      prob  = predict(model_SVM, testX, probability=TRUE)
    }
    return(list("model" = model_SVM, "test_class_prob" = prob))
  }

  if(class_method == "classtpx"){
    known_samples <- train_indices
    class_labs_train <- class_labs[train_indices]
    model_classtpx <- do.call(classtpx::class_topics, append(list(counts = as.matrix(pooled_data),
                                                                  K = length(unique(class_labs_train)),
                                                                  known_samples = known_samples,
                                                                  class_labs = class_labs_train), classtpx.control))

    omega_test <- model_classtpx$omega[test_indices,]
    return(list("model" = model_classtpx, "test_class_prob" = omega_test))
  }

}
