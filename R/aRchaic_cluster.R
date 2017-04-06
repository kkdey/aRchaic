#' @title Grade of Membership (GoM) model clustering of aDNA samples using DNA damage patterns
#'
#' @description Performs GoM model clustering of aDNA samples using DNA damage patterns- mutation,
#' flanking base, distance from read end, strand and strand break information. Upon performing the
#' model fit, it performs model visualizations and saves the figures in a user defined directory.
#'
#' @param folders a vector of folder names hosting the MFF files. Each folder may represent
#'                MFF files from same source sequenced and analysed together. (see vignette for
#'                examples)
#' @param K the number of clusters to fit to the model.
#' @param tol The tolerance level of convergence of the GoM model fit
#' @param labs The labels used to group the samples in visualization. May be used to distinguish
#'             samples from different labs, or different library prep.
#' @param levels The levels of the \code{labs} vector. Defaults to all unique values of \code{labs} vector.
#' @param run_from Can take one of three terms - \code{start}, \code{gom} and \code{plot}. For \code{start},
#'                 the function will perform all steps - processing the MFF files, clustering
#'                 and visualization from scratch. For \code{gom}, the function will first check
#'                 if the processed files are present already from previous run or from aRchaic_pool().
#'                 If present, it will use the processed files to run clustering and then do subsequent
#'                 visualization. For \code{plot}, it will first check if there is a clustering output
#'                 already present from previous clustering fit, if present, it will use it to do the
#'                 visualization. If not, will switch to \code{gom}. \code{plot} and \code{gom} are designed
#'                 to save time if user is repeating the run with the same data with updates to the plot or to
#'                 model respectively.
#' @param run_index The index vector of files to be included for each folder in the vector \code{folders}.
#'                  Defaults to using all files in the folder.
#' @param breaks The breaks used for binning the distance from ends of the reads of the mutations, when
#'               processing the MFF files. Defaults to assigning each base as one bin from
#'               end of the read to 20 bases.
#' @param flanking_bases The numbe rof flanking bases to the mutation considered. Defaults to 1.
#' @param gom_method The GoM method type. Defaults to \code{independent} model proposed by
#'                   Y. Shiraichi and M. Stephens. The other option is to use the \code{full}
#'                   model which is uses the \code{maptpx} package by Matt Taddy.
#' @param topic_cols colors attributed to each cluster used for the Structure plot visualization
#'                   of the grades of membership.
#' @param structure.control Control parameters for the Structure plot representation in
#'                          \code{StructureGGplot()} fucntion.
#' @param logo.control Control parameters for the logo plot representation \code{damageLogo5()}
#'                     function.
#' @param topics.control Control parameters for the maptpx GoM model fit.
#' @param output_dir The output directory where the model, Structure plot and the logo plots
#'                   of the clusters will be saved. If NULL, it picks the current working directory.
#' @param structure_width The width of the image of the Structure plot representation.
#' @param structure_height The height of the image of the Structure plot representation.
#'
#' @return For the \code{start} option of \code{run_from}, the function processes the MFF files
#' in each folder, aggregates them into a matrix and saves it in the folder as a .rda file. Then
#' upon GoM model fitting, it saves model output in a model.rda file, the Structure plot in a
#' structure.png file and the logo plots in \code{logo_clus_k.png} files, with \code{k}
#' representing the cluster k -  all in the \code{output_dir} provided. For \code{gom} and
#' \code{plot} options of the function, each or both rda generation steps may be skipped
#' respectively, if they are already present.
#'
#' @keywords aRchaic_cluster
#' @import gridBase
#' @import ggplot2
#' @import Logolas
#' @import CountClust
#' @export



aRchaic_cluster = function(folders,
                           K,
                           tol=100,
                           labs = NULL,
                           levels = NULL,
                           run_from = c("start", "gom", "plot"),
                           run_index = NULL,
                           breaks = c(-1, seq(1,20,1)),
                           flanking_bases = 1,
                           gom_method = "independent",
                           topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                      "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                      "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                           structure.control = list(),
                           logo.control = list(),
                           topics.control = list(),
                           output_dir = NULL,
                          # save_plot = TRUE,
                           structure_width = 5,
                           structure_height = 8){

  structure.control.default <- list(yaxis_label = "aRchaic pops",
                                 order_sample = FALSE,
                                 figure_title = paste0("  StructurePlot: K=", K,""),
                                 axis_tick = list(axis_ticks_length = .1,
                                                  axis_ticks_lwd_y = .1,
                                                  axis_ticks_lwd_x = .1,
                                                  axis_label_size = 10,
                                                  axis_label_face = "bold"),
                                 legend_title_size = 10,
                                 legend_key_size = 0.7,
                                 legend_text_size = 8)

  logo.control.default <- list(sig_names = NULL, ic.scale=TRUE,
                               max_pos = 20, flanking_bases=1,
                               yscale_change = TRUE, xaxis=TRUE,
                               yaxis=TRUE, xlab = " ", xaxis_fontsize=30,
                               xlab_fontsize=10, title_aligner = 11,
                               y_fontsize=27, title_fontsize = 35,
                               mut_width=2, start=0.0001,
                               renyi_alpha = 5, inflation_factor = c(3,1,3),
                               pop_names = paste0("Cluster : ", 1:K),
                               logoport_x = 0.25, logoport_y= 0.50, logoport_width= 0.25, logoport_height= 0.50,
                               lineport_x = 0.9, lineport_y=0.40, lineport_width=0.32, lineport_height=0.28,
                               breaklogoport_x = 0.94, breaklogoport_y = 0.40, breaklogoport_width=0.30, breaklogoport_height=0.45,
                               barport_x = 0.60, barport_y=0.60, barport_width=0.25, barport_height=0.35,
                               output_width = 1200, output_height = 700)

  topics.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                        nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                        light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)

  structure.control <- modifyList(structure.control.default, structure.control)
  logo.control <- modifyList(logo.control.default, logo.control)
  topics.control <- modifyList(topics.control.default, topics.control)



  if(is.null(labs)){
    labs <- c()
    for(i in 1:length(folders)){
      temp <- setdiff(list.files(folders[i], pattern = ".csv"), list.files(folders[i], pattern = ".csv#"))
      labs <- c(labs, rep(tail(strsplit(folders[i], "/")[[1]],1), length(temp)))
    }
  }

  if(is.null(levels)){
    levels <- unique(labs)
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

  zero_sum_rows <- which(rowSums(pooled_data) == 0)
  if(length(zero_sum_rows) > 0){
    pooled_data <- pooled_data[-zero_sum_rows, ]
    labs <- labs[-zero_sum_rows]
  }


  ########################  Grade of Membership Model  ######################################################

  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/")
  }

  if(run_from == "gom" | run_from == "start"){
    if(gom_method == "full"){
        message("Fitting the Grade of Membership Model - full version - due to Matt Taddy")
        suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = pooled_data, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
        save(topic_clus, file = paste0(output_dir, "model.rda"))
    }else if(gom_method == "independent"){
      message("Fitting the Grade of Membership Model - independent version - due to Y. Shiraichi and M. Stephens")

      signature_set <- colnames(pooled_data)
      sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
      new_sig_split <- matrix(0, dim(sig_split)[1], 3);
      new_sig_split[,1] <- sig_split[,flanking_bases]
      new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse="")))
      new_sig_split[,3] <- sig_split[,(flanking_bases+5)]

      levels(new_sig_split[,1]) <- c("0", "1", "2", "3", "4")

      pos <- t(sapply(1:length(signature_set), function(x)
      {
        y = strsplit(signature_set[x], "")[[1]]
        return(paste(y[12:length(y)], collapse=""))
      }))



      mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
      for(k in 1:dim(new_sig_split)[2]){
        temp <- as.factor(new_sig_split[,k])
        mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp), to = 0:(length(levels(temp))-1))))
      }

      pos <- as.numeric(pos)
      pos <- pos - min(pos)
      pos <- factor(pos, levels = 0:21)

      signatures <- mat;
      signature_pos <- cbind.data.frame(signatures, pos)

      suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = pooled_data, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
      save(topic_clus, file = paste0(output_dir, "model.rda"))
    }
  }else if(run_from == "plot"){
    if(gom_method == "full"){
      if(!file.exists(paste0(output_dir, "model.rda"))){
       # message("Fitting the Grade of Membership Model - full version - due to Matt Taddy")
        suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = pooled_data, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
        save(topic_clus, file = paste0(output_dir, "model.rda"))
      }else{
        topic_clus <- get(load(paste0(output_dir, "model.rda")))
      }
    }else if(gom_method == "independent"){
     # message("Fitting the Grade of Membership Model - independent version - due to Y. Shiraichi and M. Stephens")
      signature_set <- colnames(pooled_data)
      sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:8])))
      new_sig_split <- matrix(0, dim(sig_split)[1], 3);
      new_sig_split[,1] <- sig_split[,1]
      new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,2:5], collapse="")))
      new_sig_split[,3] <- sig_split[,6]

      levels(new_sig_split[,1]) <- c("0", "1", "2", "3", "4")

      pos <- t(sapply(1:length(signature_set), function(x)
      {
        y = strsplit(signature_set[x], "")[[1]]
        return(paste(y[12:length(y)], collapse=""))
      }))



      mat <- matrix(0, dim(new_sig_split)[1], dim(new_sig_split)[2])
      for(k in 1:dim(new_sig_split)[2]){
        temp <- as.factor(new_sig_split[,k])
        mat[,k] <- as.numeric(as.matrix(plyr::mapvalues(temp, from = levels(temp), to = 0:(length(levels(temp))-1))))
      }

      pos <- as.numeric(pos)
      pos <- pos - min(pos)
      pos <- factor(pos, levels = 0:21)

      signatures <- mat;
      signature_pos <- cbind.data.frame(signatures, pos)

      if(!file.exists(paste0(output_dir, "model.rda"))){
        suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = pooled_data, K=K, tol=tol, model = "independent", signatures = signature_pos), topics.control)))
        save(topic_clus, file = paste0(output_dir, "model.rda"))
      }else{
        topic_clus <- get(load(paste0(output_dir, "model.rda")))
      }
    }

  }else {
    stop("run from must be either from start (which clears everything out and restarts) or from gom (which does clustering and follow up)
         or from plot (which just does the plots)")
  }


  message ("Structure plot and Logo plot representations of clusters")

  omega <- topic_clus$omega
  annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(labs, levels = levels)
  )

    if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}
    plot.new()
    grid::grid.newpage()
    do.call(CountClust::StructureGGplot, append(list(omega= omega,
                                         annotation = annotation,
                                         palette = topic_cols),
                                                structure.control))
    ggplot2::ggsave(paste0(output_dir, "structure.png"), width = structure_width,
                    height = structure_height)

  # if(save_plot){
    if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}
    plot.new()
    do.call(damageLogo_five, append(list(theta_pool = topic_clus$theta,
                                output_dir = output_dir),
            logo.control))
    graphics.off()
  # }else if(!save_plot){
  #   plot.new()
  #   if(is.null(output_dir)){ output_dir <- paste0(getwd(), "/")}
  #   do.call(StructureGGplot, append(list(omega= omega, annotation = annotation, palette = topic_cols), structure.control))
  #   do.call(damageLogo_five, append(list(theta_pool = topic_clus$theta, output_dir = output_dir),
  #           logo.control))
  # }

  message ("Finished")
}
