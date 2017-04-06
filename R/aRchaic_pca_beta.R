#' @title Performs Principal Component Analysis (PCA) of MFF samples using filtered aDNA damage
#' patterns
#'
#' @description Performs principal component analysis on the aDNA samples using
#' filtered DNA damage patterns,which comprises of some subset of the 5 features (mutation,
#' flanking base, distance from read end, strand and strand break) but not all of them.
#' The function fits the PCA and saves the PCA output and the plot in the output directory.
#'
#' @param folders a vector of folder names hosting the MFF files. Each folder may represent
#'                MFF files from same source sequenced and analysed together. (see vignette for
#'                examples)
#' @param labs The labels used to group the samples in visualization. May be used to distinguish
#'             samples from different labs, or different library prep.
#' @param type The type of filtering mechanism used. The possible options are
#'                   - \code{mutation} : just the mutations (C->T, C->G, C->A, T->A, T->C and T->G)
#'                   - \code{mutation-flank}: the mutation and flanking base features
#'                   - \code{mutation-pos}: the mutation and position of the mutation from the nearest end of read.
#'                   - \code{mutation-flank-pos}: mutation + flanking base + position from nearest end of read
#'                   - \code{specific-mutation-pos}: specific mutation (say C->T) + position of that mutation from nearest end of read
#'                   - \code{wo-strand}: all features except strand feature
#'                   - \code{wo-strand-break}: all features except the strand break feature
#' @param run_from Can take one of two terms - \code{start} and \code{pca}. For \code{start},
#'                 the function will perform all steps - processing the MFF files, PCA fitting
#'                 and visualization from scratch. For \code{pca}, the function will first check
#'                 if the processed files are present already from previous run or from aRchaic_pool().
#'                 If present, it will use the processed files to run PCA, save the model
#'                 as \code{pca.rda} file in the \code{output_dir} and then do subsequent
#'                 visualization.
#' @param run_index The index vector of files to be included for each folder in the vector \code{folders}.
#'                  Defaults to using all files in the folder.
#' @param breaks The breaks used for binning the distance from ends of the reads of the mutations, when
#'               processing the MFF files. Defaults to assigning each base as one bin from
#'               end of the read to 20 bases.
#' @param flanking_bases The numbe rof flanking bases to the mutation considered. Defaults to 1.
#' @param pattern type of mutation (say "C->T"). Comes into play when \code{type} is \code{specific-mutation-pos}.
#'                Defaults to NULL.
#' @param max_pos The maximum distance from ends of read that a mismatch would be considered.
#'                The rest will be filtered out. Defaults to 20.
#' @param normalize If \code{normalize = FALSE}), original data is used for fitting PCA, else
#'                  log CPM normalized data is used. Defaults to TRUE.
#' @param cols colors attributed to each label in the \code{labs} vector.
#' @param pcs_to_plot a vector of the indicators of which PCs to plot. For example
#'                    \code{pcs_to_plot = c("PC1", "PC2", "PC3")} will plot three scatter
#'                    plots of PC1 vs PC2, PC1 vs PC3 and PC2 vs PC3.
#' @param lay A layout matrix. Defaults to to the upper diagonal layout. For other
#'            options, check \code{grid.arrange} page
#'            (https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html).
#' @param plot_width The width of the PCA plot figure.
#' @param plot_height The height of the PCA plot figure.
#' @param output_dir  The output directory where the model and the PCA plot figure are saved.
#'                    If NULL, it picks the current working directory.
#' @param output_name The name of the PCA figure.
#'
#' @return For the \code{start} option of \code{run_from}, the function processes the MFF files
#' in each folder, aggregates them into a matrix, filters the signatures by the \code{type} of
#' the mutational signature and then fits PCA and saves the PCA plot.
#' @keywords aRchaic_pca
#' @import ggplot2, gridBase, grid
#' @export



aRchaic_pca_beta =  function(folders,
                             labs = NULL,
                             type = c("mutation", "mutation-flank", "mutation-pos","mutation-flank-pos",
                                      "specific-mutation-pos", "wo-strand", "wo-strand-break"),
                             run_from = c("start", "pca"),
                             run_index = NULL,
                             breaks = c(-1, seq(1,20,1)),
                             flanking_bases = 1,
                             pattern = NULL,
                             max_pos = 20,
                             normalize=TRUE,
                             cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                 "hotpink","burlywood","darkkhaki","yellow","darkgray","deepskyblue",
                                 "brown4","darkorchid","magenta", "azure1","azure4"),
                             pcs_to_plot = c("PC1", "PC2", "PC3"),
                             lay = NULL,
                             plot_width = 10,
                             plot_height = 7,
                             output_dir = NULL,
                             output_name = "pca"){

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}

  if(is.null(labs)){
    labs <- c()
    for(i in 1:length(folders)){
      temp <- setdiff(list.files(folders[i], pattern = ".csv"), list.files(folders[i], pattern = ".csv#"))
      labs <- c(labs, rep(tail(strsplit(folders[i], "/")[[1]],1), length(temp)))
    }
  }

  ##########################  Check if the folder names actually exist #######################################

  message("Checking if the folders exist")

  for(i in 1:length(folders)){
    if(!file.exists(folders[i]))
      stop("A folder in the folder list does not exist:  aborting")
  }

  datalist <- vector("list", length(folders))

  #########################  If the user wants to run from scratch  ##########################################

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

  #########################  Run aggregation functions on MutationFeatureFormat #############################################


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
  if(normalize){
    voom_mat_reduced <- t(limma::voom(t(mat_reduced))$E);
    pr <- prcomp(voom_mat_reduced)
  }else{
    pr <- prcomp(mat_reduced)
  }

 if(is.null(labs)){
   labs <- rep(1, dim(voom_mat_reduced)[1])
 }

  pc_data_frame <- cbind.data.frame(pr$x, labs)

  graphList <- vector(mode="list")

  indices <- array(0, length(pcs_to_plot))

  pc_data_frame <- cbind.data.frame(pr$x, labs)

  for(num in 1:length(pcs_to_plot)){
    indices[num] <- match(paste0(pcs_to_plot[num]), colnames(pc_data_frame))
  }

  pc_data_frame_filtered <- pc_data_frame[,indices]

  total_comb <- length(pcs_to_plot)*(length(pcs_to_plot) - 1)/2

  if(is.null(lay)){
    lay <- matrix(0, length(pcs_to_plot) - 1, length(pcs_to_plot) - 1)
    l <- 1
    for(m in 1:nrow(lay)){
      for(n in 1:ncol(lay)){
        if(m > n){
          lay[m,n] <- NA
        }else{
          lay[m,n] = l
          l <- l+1
        }
      }
    }
  }

  lnum <- 1
  for(m in 1:(length(pcs_to_plot)-1)){
    for(n in (m+1):(length(pcs_to_plot))){
      a <- ggplot2::ggplot(pc_data_frame_filtered, ggplot2::aes_string(x = paste0("PC", indices[m]), y = paste0("PC", indices[n]))) +
        geom_point(aes(colour = factor(labs))) + ggplot2::scale_colour_manual(values = cols, guide = ggplot2::guide_legend(title = "Groups")) +
        ggplot2::xlab(paste0(pcs_to_plot[m])) + ggplot2::ylab(paste0(pcs_to_plot[n]))
      graphList[[lnum]] <- a
      lnum = lnum +1
    }
  }

  plot.new()
  grid.newpage()
  a <- do.call("grid.arrange",
               args = list(grobs=graphList,
                           layout_matrix = lay))

  ggplot2::ggsave(paste0(output_dir, output_name, ".png"), a, width = plot_width,
                  height = plot_height)

  return(pr)
}
