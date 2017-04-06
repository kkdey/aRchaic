#' @title Performs non-linear dimension reduction using t-SNE of MFF samples using aDNA damage
#' patterns
#'
#' @description Performs t-SNE on the aDNA samples using the mutational
#' signature patterns- mutation, flanking base, distance from read end, strand and
#' strand break information.
#'
#' @param folders a vector of folder names hosting the MFF files. Each folder may represent
#'                MFF files from same source sequenced and analysed together. (see vignette for
#'                examples)
#' @param labs The labels used to group the samples in visualization. May be used to distinguish
#'             samples from different labs, or different library prep.
#' @param run_from Can take one of three terms - \code{start}, \code{tsne} and \code{plot}. For \code{start},
#'                 the function will perform all steps - processing the MFF files, t-SNE fitting
#'                 and visualization from scratch. For \code{tsne}, the function will first check
#'                 if the processed files are present already from previous run or from aRchaic_pool().
#'                 If present, it will use the processed files to run t-SNE, save the model
#'                 as \code{tsne.rda} file in the \code{output_dir} and then do subsequent
#'                 visualization. For \code{plot}, it will first check if there is a
#'                 \code{tsne.rda} in \code{output_dir}. If present, it will use it to do the
#'                 visualization. If not, will switch to \code{tsne}. \code{plot} and \code{tsne} are
#'                 designed to save time if user is repeating the run with the same data with
#'                 updates to the plot or to model respectively.
#' @param run_index The index vector of files to be included for each folder in the vector \code{folders}.
#'                  Defaults to using all files in the folder.
#' @param breaks The breaks used for binning the distance from ends of the reads of the mutations, when
#'               processing the MFF files.
#' @param flanking_bases The numbe rof flanking bases to the mutation considered. Defaults to 1.
#' @param normalize If SVM option is chosen by the user, he can run it on the original data (for \code{normalize = FALSE})
#'                  or on log CPM normalized data (for \code{normalize = TRUE})
#' @param cols colors attributed to each label in the \code{labs} vector.
#' @param dims The dimensions of the manifold on which the data is projected.
#' @param lay A layout matrix. Defaults to to the upper diagonal layout. For other
#'            options, check \code{grid.arrange} page
#'            (https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html).
#' @param filter_indices An option to select a subset of features to run tSNE. t-SNE is expensive
#'                       on around 15,000 features, so this can be a easy way to test the model with
#'                       fewer features. Defaults to NULL, meaning all features are taken.
#' @param plot_width The width of the t-SNE plot figure. Defaults to 10.
#' @param plot_height The height of the t-SNE plot figure. Defaults to 7.
#' @param output_dir  The output directory where the model and the t-SNE plot figure are saved.
#'                    If NULL, it picks the current working directory.
#'
#' @return For the \code{start} option of \code{run_from}, the function processes the MFF files
#' in each folder, aggregates them into a matrix and saves it in the folder as a .rda file. Then
#' it fits t-SNE to the aggregated data from one or more folders as determined by
#' the vector \code{folders} and saves the result as \code{tsne.rda} file in the
#' \code{output_dir}. Then the t-SNE plot corresponding to the number of dimensions
#' to fit \code{dims} is saved as \code{tsne.png} file in the \code{output_dir}.For option \code{tsne}, the function skips
#' the data aggregation and processing step if the .rda files in the folders are present. For
#' \code{plot}, the functions skips the t-SNE fitting if the \code{tsne.rda} file is already
#' present in \code{output_dir}. The \code{tsne} and \code{plot} options are meant to reduce
#' time on repeated jobs.
#'
#' @keywords aRchaic_tsne
#' @import ggplot2
#' @import gridBase
#' @import grid
#' @export


aRchaic_tsne =  function(folders,
                        labs = NULL,
                        run_from = c("start", "tsne", "plot"),
                        run_index = NULL,
                        breaks = c(-1, seq(1,20,1)),
                        flanking_bases = 1,
                        normalize=TRUE,
                        cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                 "hotpink","burlywood","darkkhaki","yellow","darkgray","deepskyblue",
                                 "brown4","darkorchid","magenta", "azure1","azure4"),
                        dims = 3,
                        lay = NULL,
                        filter_indices = NULL,
                        plot_width = 10,
                        plot_height = 7,
                        output_dir = NULL){

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

  if(is.null(filter_indices)){
    filter_indices <- 1:dim(pooled_data)[2]
  }
  pooled_data <- pooled_data[,filter_indices]

  if(run_from == "tsne" | run_from == "start"){
    message("Fitting t-SNE and saving it to file")
    if(normalize){
      voom_pooled_data <- t(limma::voom(t(pooled_data))$E);
      tsne_out <- tsne::tsne(voom_pooled_data, k=dims)
    }else{
      tsne_out <- tsne::tsne(pooled_data, k=dims)
    }
    save(tsne_out, file = paste0(output_dir, "tsne.rda"))
  }else if(run_from == "plot"){
    if(!file.exists(paste0(output_dir, "tsne.rda"))){
      message("Fitting t-SNE and saving it to file")
      if(normalize){
        voom_pooled_data <- t(limma::voom(t(pooled_data))$E);
        tsne_out <- tsne::tsne(voom_pooled_data, k=dims)
      }else{
        tsne_out <- tsne::tsne(pooled_data, k = dims)
      }
      save(tsne_out, file = paste0(output_dir, "tsne.rda"))
    }else{
      message("Reading t-SNE fit from previously saved file")
      tsne_out <- get(load(paste0(output_dir, "tsne.rda")))
    }
  }else{
    stop("run from must be either from start (which clears everything out and restarts) or from t-SNE (which does t-SNE and visualization)
         or from plot (which just does the t-SNE visualization on fitted model if present)")
  }

  colnames(tsne_out) <- paste0("tSNE", 1:dims)
  dims_to_plot <- colnames(tsne_out)
  tsne_data_frame <- cbind.data.frame(tsne_out, labs)

  graphList <- vector(mode="list")

  indices <- array(0, length(dims_to_plot))

  tsne_data_frame <- cbind.data.frame(tsne_out, labs)

  total_comb <- length(dims_to_plot)*(length(dims_to_plot) - 1)/2

  if(is.null(lay)){
    lay <- matrix(0, length(dims_to_plot) - 1, length(dims_to_plot) - 1)
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
  for(m in 1:(length(dims_to_plot)-1)){
    for(n in (m+1):(length(dims_to_plot))){
      a <- ggplot2::ggplot(tsne_data_frame, ggplot2::aes_string(x = paste0("tSNE",m), y = paste0("tSNE",n))) +
        geom_point(aes(colour = factor(labs))) + ggplot2::scale_colour_manual(values = cols, guide = ggplot2::guide_legend(title = "Groups")) +
        ggplot2::xlab(paste0(dims_to_plot[m])) + ggplot2::ylab(paste0(dims_to_plot[n]))
      graphList[[lnum]] <- a
      lnum = lnum +1
    }
  }

  plot.new()
  grid.newpage()
  a <- do.call("grid.arrange",
               args = list(grobs=graphList,
                           layout_matrix = lay))

  ggplot2::ggsave(paste0(output_dir, "tsne.png"), a, width = plot_width,
                  height = plot_height)

  }








