#' @title Grade of Membership (GoM) model clustering of aDNA samples using filtered DNA damage patterns
#'
#' @description Performs GoM model clustering of aDNA samples using filtered DNA damage
#' patterns, which comprises of some subset of the 5 features (mutation,
#' flanking base, distance from read end, strand and strand break). An example is mutation-flank,
#' which would mean we just used mutation and flanking base features to do the clustering.
#' Upon performing the model fit, it performs model visualizations and saves the figures
#' in a user defined directory.
#'
#'
#' @param folders a vector of folder names hosting the MFF files. Each folder may represent
#'                MFF files from same source sequenced and analysed together. (see vignette for
#'                examples)
#' @param K the number of clusters to fit to the model.
#' @param type The type of filtering mechanism used. The possible options are
#'                   - \code{mutation} : just the mutations (C->T, C->G, C->A, T->A, T->C and T->G)
#'                   - \code{mutation-flank}: the mutation and flanking base features
#'                   - \code{mutation-pos}: the mutation and position of the mutation from the nearest end of read.
#'                   - \code{mutation-flank-pos}: mutation + flanking base + position from nearest end of read
#'                   - \code{specific-mutation-pos}: specific mutation (say C->T) + position of that mutation from nearest end of read
#'                   - \code{wo-strand}: all features except strand feature
#'                   - \code{wo-strand-break}: all features except the strand break feature
#' @param tol The tolerance level of convergence of the GoM model fit
#' @param pattern type of mutation (say "C->T"). Comes into play when \code{type} is \code{specific-mutation-pos}.
#'                Defaults to NULL.
#' @param max_pos The maximum distance from ends of read that a mismatch would be considered.
#'                The rest will be filtered out.
#' @param labs The labels used to group the samples in visualization. May be used to distinguish
#'             samples from different labs, or different library prep.
#' @param levels The levels of the \code{labs} vector. Defaults to all unique values of \code{labs} vector.
#' @param flanking_bases The numbe rof flanking bases to the mutation considered. Defaults to 1.
#' @param gom_method The GoM method type. Defaults to \code{independent} model proposed by
#'                   Y. Shiraichi and M. Stephens. The other option is to use the \code{full}
#'                   model which is uses the \code{maptpx} package by Matt Taddy.
#' @param topic_cols colors attributed to each cluster used for the Structure plot visualization
#'                   of the grades of membership.
#' @param structure.control Control parameters for the Structure plot representation in
#'                          \code{StructureGGplot()} fucntion. Defaults to \code{list()}
#' @param graph.control Control parameters for the graphical representation of the clusters.
#'                      Defaults to \code{list()}.
#' @param topics.control Control parameters for the maptpx GoM model fit.
#' @param output_dir The output directory where the model, Structure plot and the logo plots
#'                   of the clusters will be saved. If NULL, it picks the current working directory.
#' @param structure_width The width of the image of the Structure plot representation.
#' @param structure_height The height of the image of the Structure plot representation.
#' @param inflation The inflation of flanking base logo plot if needed. Defaults to c(2,1,2),
#'                   which means the flanking base heights will be inflated 2 times with respect
#'                   to mutation in the middle in the logo plot.
#' @param output_width the width of the output plots of the graphs for cluster representation.
#' @param output_height the height of the output plots of the graphs for cluster representation.
#'
#' @return The function first filters the data based on the filtering \code{type} provided. Then
#' it performs GoM model fitting on the filtered data and saves model output in a model.rda file
#' Then it performs visualization of the grades of memebrships of the clusters and graphs for
#' cluster representation and save them as  Structure plot in a
#' structure.png file and the logo plots in \code{logo_clus_k.png} files, with \code{k}
#' representing the cluster k respectively-  all in the \code{output_dir} provided.
#'
#' @keywords aRchaic_cluster_beta
#' @import gridBase
#' @import ggplot2
#' @import Logolas
#' @import CountClust
#' @export



aRchaic_cluster_beta = function(folders,
                                K,
                                type = c("mutation", "mutation-flank", "mutation-pos", "mutation-flank-pos",
                                         "specific-mutation-pos", "wo-strand", "wo-strand-break"),
                                tol=0.01,
                                run_from = c("start", "gom"),
                                pattern = NULL,
                                max_pos = 20,
                                labs = NULL,
                                levels = NULL,
                                breaks = c(-1, seq(1,20,1)),
                                flanking_bases = 1,
                                gom_method = "independent",
                                topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                               "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                                structure.control = list(),
                                graph.control = list(),
                                topics.control = list(),
                                output_dir = NULL,
                                structure_width = 5,
                                structure_height = 8,
                                inflation = rep(1,1,1),
                                output_width = 1200,
                                output_height = 700){


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


 mat <- pooled_data


if(type == "mutation"){
  aRchaic_cluster_beta_mutation(mat=mat, K=K, tol=tol,
                                labs = labs, levels = levels, flanking_bases = flanking_bases,
                                gom_method = gom_method, topic_cols = topic_cols,
                                structure.control = structure.control,
                                logo.control = graph.control, topics.control = topics.control,
                                output_dir = output_dir, structure_width = structure_width,
                                structure_height = structure_height,
                                output_width = output_width, output_height = output_height)
}
else if(type == "mutation-flank"){
  aRchaic_cluster_beta_mutation_flank(mat=mat, K=K, tol=tol,
                                labs = labs, levels = levels, flanking_bases = flanking_bases,
                                gom_method = gom_method, topic_cols = topic_cols,
                                structure.control = structure.control,
                                logo.control = graph.control, topics.control = topics.control,
                                output_dir = output_dir, structure_width = structure_width,
                                structure_height = structure_height, inflation = inflation,
                                output_width = output_width, output_height = output_height)
}
else if(type == "mutation-pos"){
  aRchaic_cluster_beta_mutation_pos(mat=mat, K=K, tol=tol, max_pos = max_pos,
          labs = labs, levels = levels, flanking_bases = flanking_bases,
          gom_method = gom_method, topic_cols = topic_cols, structure.control = structure.control,
          logo.control = graph.control, topics.control = topics.control,
          output_dir = output_dir, structure_width = structure_width,
          structure_height = structure_height,
           output_width = output_width, output_height = output_height)
}
else if(type == "mutation-flank-pos"){
  aRchaic_cluster_beta_mutation_flank_pos(mat=mat, K=K, tol=tol, max_pos = max_pos,
                                    labs = labs, levels = levels, flanking_bases = flanking_bases,
                                    gom_method = gom_method, topic_cols = topic_cols, structure.control = structure.control,
                                    logo.control = graph.control, topics.control = topics.control,
                                    output_dir = output_dir, structure_width = structure_width,
                                    structure_height = structure_height, inflation = inflation,
                                    output_width = output_width, output_height = output_height)
}
else if(type == "specific-mutation-pos"){
  if(is.null(pattern)){
    stop("for this option, user needs to provide a mutation type: C->T, C->A C->G, T->A, T->C or T->G")
  }
  aRchaic_cluster_beta_pos(mat=mat, pattern = pattern, K=K, tol=tol, max_pos = max_pos,
                           labs = labs, levels = levels, flanking_bases = flanking_bases,
                           gom_method = gom_method, topic_cols = topic_cols,
                           structure.control = structure.control,
                           graph.control = graph.control,
                           topics.control = topics.control,
                           output_dir = output_dir, structure_width = structure_width,
                           structure_height = structure_height,
                           output_width = output_width, output_height = output_height)
}
else if(type == "wo-strand"){
  aRchaic_cluster_beta_wo_strand(mat=mat, K=K, tol=tol, max_pos = max_pos,
                                labs = labs, levels = levels, flanking_bases = flanking_bases,
                                gom_method = gom_method, topic_cols = topic_cols, structure.control = structure.control,
                                logo.control = graph.control, topics.control = topics.control,
                                output_dir = output_dir, structure_width = structure_width,
                                structure_height = structure_height, inflation = inflation,
                                output_width = output_width, output_height = output_height)
}
else if(type == "wo-strand-break"){
    aRchaic_cluster_beta_wo_strand_break(mat=mat, K=K, tol=tol, max_pos = max_pos,
                                        labs = labs, levels = levels, flanking_bases = flanking_bases,
                                        gom_method = gom_method, topic_cols = topic_cols, structure.control = structure.control,
                                        logo.control = graph.control, topics.control = topics.control,
                                        output_dir = output_dir, structure_width = structure_width,
                                        tructure_height = structure_height, inflation = inflation,
                                        output_width = output_width, output_height = output_height)
}else{
  stop("the type of filtering does not match with the possible options: see documentation")
  }
}


######################  aRchaic_cluster_beta : Type 1  ####################################

aRchaic_cluster_beta_mutation = function(mat,
                                         K,
                                         tol=0.01,
                                         labs = NULL,
                                         levels = NULL,
                                         flanking_bases = 1,
                                         gom_method = "independent",
                                         topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                                        "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                                        "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                                         structure.control = list(),
                                         logo.control = list(),
                                         topics.control = list(),
                                         output_dir = NULL,
                                         structure_width = 5,
                                         structure_height = 8,
                                         output_width = 1200,
                                         output_height = 700){

  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/")
  }

  topics.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                                 nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                                 light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)
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

  logo.control.default <- list(ic = NULL, hist = TRUE, frame_width =1,
                               ic.scale = FALSE,
                               xaxis_fontsize = 10, xlab_fontsize = 15,
                               y_fontsize = 15, main_fontsize = 16, start = 0.001,
                               yscale_change = TRUE, pop_name = paste0("logo plot mutation"), xlab = "",
                               ylab = "composition", col_line_split = "grey80", scale0 = 0.01,
                               scale1 = 0.99, newpage = TRUE)

  structure.control <- modifyList(structure.control.default, structure.control)
  logo.control <- modifyList(logo.control.default, logo.control)
  topics.control <- modifyList(topics.control.default, topics.control)



  mat_reduced <- filter_by_mutation(mat)
  signature_set <- colnames(mat_reduced)



  if(gom_method == "full"){
    message("Fitting the Grade of Membership Model - full version - due to Matt Taddy")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(gom_method == "independent"){
    message("Fitting the Grade of Membership Model - full version - due to Y. Shiraichi and M. Stephens")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "independent", signatures = signature_set), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(is.null(labs)){
    labs <- rownames(mat_reduced)
  }

  if(is.null(levels)){
    levels <- unique(labs)
  }

  message ("Structure plot and Logo plot representations of clusters")


  omega <- topic_clus$omega
  annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(labs, levels = levels)
  )

  rownames(omega) <- annotation$sample_id

  plot.new()
  grid.newpage()
  do.call(CountClust::StructureGGplot, append(list(omega= omega,
                                       annotation = annotation,
                                       palette = topic_cols),
                                  structure.control))
  ggplot2::ggsave(paste0(output_dir, "structure.png"), width = structure_width,
                  height = structure_height)

  theta_pool <- topic_clus$theta
  signature_set <- rownames(theta_pool)
  theta <- dplyr::tbl_df(data.frame(theta_pool)) %>% dplyr::mutate(sig = signature_set) %>% dplyr::group_by(sig) %>% dplyr::summarise_each(funs(sum)) %>% as.data.frame()
  rownames(theta) <-  theta[,1]
  theta <- theta[,-1]

  mutations <- as.character(sapply(rownames(theta), function(x) return(paste0(strsplit(x,"")[[1]][-2], collapse = ""))))
  rownames(theta) <- mutations

  for(m in 1:dim(theta)[2]){
    tab <- matrix(theta[,m], nrow=length(theta[,1]))
    rownames(tab) <- rownames(theta)
    colnames(tab) <- "mutation"

    cols = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, cols$maxcolors, rownames(cols)))

    total_chars = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O",
                    "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "zero", "one", "two",
                    "three", "four", "five", "six", "seven", "eight", "nine", "dot", "comma",
                    "dash", "colon", "semicolon", "leftarrow", "rightarrow")

    cols = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, cols$maxcolors, rownames(cols)))
    color_code = sample(col_vector, length(total_chars), replace=FALSE)

    color_code[c(1,3,7,20, 43)] <- c("green", "blue", "orange", "red", "gray")

    color_profile <- list("type" = "per_symbol",
                          "col" = color_code)

    plot.new()
    png(paste0(output_dir, "logo_clus_", m, ".png"), width=output_width, height = output_height)
    do.call(Logolas::logomaker, c(list(table=tab, color_profile = color_profile),
                                  logo.control))
    dev.off()
  }

  graphics.off()
  message("Finished")
}

#######################  aRchaic_cluster_beta type 2  #####################################

aRchaic_cluster_beta_mutation_flank = function(mat,
                                               K,
                                               tol=0.01,
                                               labs = NULL,
                                               levels = NULL,
                                               flanking_bases = 1,
                                               gom_method = "independent",
                                               topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                                              "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                                              "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                                               structure.control = list(),
                                               logo.control = list(),
                                               topics.control = list(),
                                               output_dir = NULL,
                                               structure_width = 5,
                                               structure_height = 8,
                                               inflation = c(1,1,1),
                                               output_width = 1200,
                                               output_height = 700){
  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/")
  }

  topics.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                                 nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                                 light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)
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

  logo.control.default <- list(hist = FALSE, frame_width =1,
                               ic.scale = TRUE,
                               xaxis_fontsize = 10, xlab_fontsize = 15,
                               y_fontsize = 15, main_fontsize = 16, start = 0.001,
                               yscale_change = TRUE, pop_name = paste0("logo plot mutation"), xlab = "",
                               ylab = "composition", col_line_split = "grey80", scale0 = 0.01,
                               scale1 = 0.99, newpage = TRUE)

  structure.control <- modifyList(structure.control.default, structure.control)
  logo.control <- modifyList(logo.control.default, logo.control)
  topics.control <- modifyList(topics.control.default, topics.control)

  mat_reduced <- filter_by_mutation_flank(mat)
  signature_set <- colnames(mat_reduced)
  sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:(flanking_bases+5)])))
  new_sig_split <- matrix(0, dim(sig_split)[1], 3);
  new_sig_split[,1] <- sig_split[,1]
  new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse="")))
  new_sig_split[,3] <- sig_split[,(flanking_bases+5)]

  if(gom_method == "full"){
    message("Fitting the Grade of Membership Model - full version - due to Matt Taddy")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(gom_method == "independent"){
    message("Fitting the Grade of Membership Model - full version - due to Y. Shiraichi and M. Stephens")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "independent", signatures = signature_set), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(is.null(labs)){
    labs <- rownames(mat_reduced)
  }

  if(is.null(levels)){
    levels <- unique(labs)
  }

  message ("Structure plot and Logo plot representations of clusters")

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}

  omega <- topic_clus$omega
  annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(labs, levels = levels)
  )

  rownames(omega) <- annotation$sample_id

  plot.new()
  grid.newpage()
  do.call(CountClust::StructureGGplot, append(list(omega= omega,
                                       annotation = annotation,
                                       palette = topic_cols),
                                  structure.control))
  ggplot2::ggsave(paste0(output_dir, "structure.png"), width = structure_width,
                  height = structure_height)

  theta_pool <- topic_clus$theta
  signature_set <- rownames(theta_pool)
  sig_split <- do.call(rbind,
                       lapply(signature_set,
                              function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))
  ncol_sig <- (4+2*flanking_bases)

  if(flanking_bases%%1 != 0){
    stop("flanking bases not evenly distributed")
  }


  sub_pattern <- sapply(1:dim(sig_split)[1],
                        function(x) paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse=""))
  new_sig_split <- cbind(sig_split[,1:flanking_bases], sub_pattern, sig_split[,((ncol_sig - flanking_bases +1):ncol_sig)])
  colnames(new_sig_split) = NULL

  prop_patterns_list <- list()

  theta <- dplyr::tbl_df(data.frame(theta_pool)) %>% dplyr::mutate(sig = signature_set) %>% dplyr::group_by(sig) %>% dplyr::summarise_each(funs(sum)) %>% as.data.frame()
  rownames(theta) <-  theta[,1]
  theta <- theta[,-1]

  for(l in 1:dim(theta)[2]){
    prop_patterns_list[[l]] <- numeric();
    for(j in 1:ncol(new_sig_split)){
      temp <- tapply(theta[,l], factor(new_sig_split[,j], levels=c("A", "C", "G", "T",
                                                                   "C->T", "C->A", "C->G",
                                                                   "T->A", "T->C", "T->G")), sum)

      temp[is.na(temp)]=0
      prop_patterns_list[[l]] <- cbind(prop_patterns_list[[l]], temp)
      rownames(prop_patterns_list[[l]]) <- gsub("->", ">", rownames(prop_patterns_list[[l]]))

    }
  }

  ic <- damage.ic(prop_patterns_list, alpha=1, inflation_factor = c(1,1,1))


  for(l in 1:dim(theta)[2]){
    tab <- prop_patterns_list[[l]]
    colnames(tab) <- c("left \n flank", "mutation", "right \n flank")

    cols = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, cols$maxcolors, rownames(cols)))

    total_chars = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O",
                    "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "zero", "one", "two",
                    "three", "four", "five", "six", "seven", "eight", "nine", "dot", "comma",
                    "dash", "colon", "semicolon", "leftarrow", "rightarrow")

    cols = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, cols$maxcolors, rownames(cols)))
    color_code = sample(col_vector, length(total_chars), replace=FALSE)

    color_code[c(1,3,7,20, 43)] <- c("green", "blue", "orange", "red", "gray")

    color_profile <- list("type" = "per_symbol",
                          "col" = color_code)

    png(paste0(output_dir, "logo_clus_", l, ".png"), width=output_width, height = output_height)
    do.call(Logolas::logomaker, c(list(table=tab, ic = ic[,l], color_profile = color_profile),
                                  logo.control))
    dev.off()
  }

  graphics.off()
  message("Finished")
}


######################  aRchaic_cluster_beta Type 3  ###################################

aRchaic_cluster_beta_mutation_pos = function(mat,
                                             K,
                                             tol=0.01,
                                             max_pos = 20,
                                             labs = NULL,
                                             levels = NULL,
                                             flanking_bases = 1,
                                             gom_method = "independent",
                                             topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                                            "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                                            "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                                             structure.control = list(),
                                             logo.control = list(),
                                             topics.control = list(),
                                             output_dir = NULL,
                                             structure_width = 5,
                                             structure_height = 8,
                                             output_width = 1200,
                                             output_height = 700){

  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/")
  }

  topics.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                                 nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                                 light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)
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

  logo.control.default <- list(sig_names = NULL, ic.scale=FALSE,
                               max_pos = max_pos, flanking_bases=flanking_bases,
                               yscale_change = TRUE, xaxis=TRUE,
                               yaxis=TRUE, xlab = " ", xaxis_fontsize=27,
                               xlab_fontsize=10, title_aligner = 11,
                               y_fontsize=27, title_fontsize = 35,
                               mut_width=2, start=0.0001,
                               renyi_alpha = 5, inflation_factor = c(2,1,2),
                               pop_names = paste0("Cluster : ", 1:K),
                               logoport_x = 0.24, logoport_y= 0.50, logoport_width= 0.15,
                               logoport_height= 0.50,
                               lineport_x = 0.95, lineport_y=0.73, lineport_width=0.35,
                               lineport_height=0.48,
                               output_width = output_width, output_height = output_height)

  structure.control <- modifyList(structure.control.default, structure.control)
  logo.control <- modifyList(logo.control.default, logo.control)
  topics.control <- modifyList(topics.control.default, topics.control)

  mat_reduced <- filter_by_mutation_flank_pos(mat, max_pos = 20)

  signature_set <- colnames(mat_reduced)
  sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:(flanking_bases+5)])))
  new_sig_split <- matrix(0, dim(sig_split)[1], 3);
  new_sig_split[,1] <- sig_split[,1]
  new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse="")))
  new_sig_split[,3] <- sig_split[,(flanking_bases+5)]

  pos <- sapply(signature_set, function(x) return(strsplit(x, "_")[[1]][2]))

  pos <- as.numeric(pos)
  pos <- pos - min(pos)
  pos <- factor(pos, levels = 0:21)


  signatures <- new_sig_split;
  signature_pos <- cbind.data.frame(signatures, pos)

  if(gom_method == "full"){
    message("Fitting the Grade of Membership Model - full version - due to Matt Taddy")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(gom_method == "independent"){
    message("Fitting the Grade of Membership Model - full version - due to Y. Shiraichi and M. Stephens")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "independent", signatures = signature_pos), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(is.null(labs)){
    labs <- rownames(mat_reduced)
  }

  if(is.null(levels)){
    levels <- unique(labs)
  }

  message ("Structure plot and Logo plot representations of clusters")

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}

  omega <- topic_clus$omega
  annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(labs, levels = levels)
  )

  rownames(omega) <- annotation$sample_id

  plot.new()
  grid.newpage()
  do.call(CountClust::StructureGGplot, append(list(omega= omega,
                                       annotation = annotation,
                                       palette = topic_cols),
                                  structure.control))
  ggplot2::ggsave(paste0(output_dir, "structure.png"), width = structure_width,
                  height = structure_height)

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}
  plot.new()
  do.call(damageLogo_four, append(list(theta_pool = topic_clus$theta,
                                       output_dir = output_dir),
                                  logo.control))
  graphics.off()

  message("Finished")
}

##################### aRchaic_cluster_beta type 4  ##################################

aRchaic_cluster_beta_mutation_flank_pos = function(mat,
                                                   K,
                                                   tol=0.01,
                                                   max_pos = 20,
                                                   labs = NULL,
                                                   levels = NULL,
                                                   flanking_bases = 1,
                                                   gom_method = "independent",
                                                   topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                                                  "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                                                  "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                                                   structure.control = list(),
                                                   logo.control = list(),
                                                   topics.control = list(),
                                                   output_dir = NULL,
                                                   structure_width = 5,
                                                   structure_height = 8,
                                                   inflation = rep(2,1,2),
                                                   output_width = 1200,
                                                   output_height = 700){

  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/")
  }

  topics.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                                 nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                                 light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)
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
                               max_pos = max_pos, flanking_bases=flanking_bases,
                               yscale_change = TRUE, xaxis=TRUE,
                               yaxis=TRUE, xlab = " ", xaxis_fontsize=20,
                               xlab_fontsize=10, title_aligner = 11,
                               y_fontsize=20, title_fontsize = 35,
                               mut_width=2, start=0.0001,
                               renyi_alpha = 5, inflation_factor = c(2,1,2),
                               pop_names = paste0("Cluster : ", 1:K),
                               logoport_x = 0.25, logoport_y= 0.50, logoport_width= 0.25,
                               logoport_height= 0.50,
                               lineport_x = 0.9, lineport_y=0.73, lineport_width=0.35,
                               lineport_height=0.48,
                               output_width = output_width, output_height = output_height)

  structure.control <- modifyList(structure.control.default, structure.control)
  logo.control <- modifyList(logo.control.default, logo.control)
  topics.control <- modifyList(topics.control.default, topics.control)

  mat_reduced <- filter_by_mutation_flank_pos(mat, max_pos = 20)

  signature_set <- colnames(mat_reduced)
  sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:(flanking_bases+5)])))
  new_sig_split <- matrix(0, dim(sig_split)[1], 3);
  new_sig_split[,1] <- sig_split[,1]
  new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse="")))
  new_sig_split[,3] <- sig_split[,(flanking_bases+5)]

  pos <- sapply(signature_set, function(x) return(strsplit(x, "_")[[1]][2]))

  pos <- as.numeric(pos)
  pos <- pos - min(pos)
  pos <- factor(pos, levels = 0:21)


  signatures <- new_sig_split;
  signature_pos <- cbind.data.frame(signatures, pos)

  if(gom_method == "full"){
    message("Fitting the Grade of Membership Model - full version - due to Matt Taddy")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(gom_method == "independent"){
    message("Fitting the Grade of Membership Model - full version - due to Y. Shiraichi and M. Stephens")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "independent", signatures = signature_pos), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(is.null(labs)){
    labs <- rownames(mat_reduced)
  }

  if(is.null(levels)){
    levels <- unique(labs)
  }

  message ("Structure plot and Logo plot representations of clusters")

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}

  omega <- topic_clus$omega
  annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(labs, levels = levels)
  )

  rownames(omega) <- annotation$sample_id

  plot.new()
  grid.newpage()
  do.call(CountClust::StructureGGplot, append(list(omega= omega,
                                       annotation = annotation,
                                       palette = topic_cols),
                                  structure.control))
  ggplot2::ggsave(paste0(output_dir, "structure.png"), width = structure_width,
                  height = structure_height)

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}
  plot.new()
  do.call(damageLogo_one, append(list(theta_pool = topic_clus$theta,
                                      output_dir = output_dir),
                                 logo.control))
  graphics.off()

  message("Finished")
}


##################  aRchaic_clus_beta  Type 5  ########################################

aRchaic_cluster_beta_pos = function(mat,
                                    pattern = "C->T",
                                    K,
                                    tol=0.01,
                                    labs = NULL,
                                    levels = NULL,
                                    max_pos = 20,
                                    flanking_bases = 1,
                                    gom_method = "independent",
                                    topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                                   "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                                   "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                                    structure.control = list(),
                                    graph.control = list(),
                                    topics.control = list(),
                                    output_dir = NULL,
                                    structure_width = 5,
                                    structure_height = 8,
                                    output_width = 500,
                                    output_height = 700){

  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/")
  }

  topics.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                                 nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                                 light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)
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

  graph.control.default <- list(col="red",
                                cex=unit(1, "npc"), pch=unit(16,"npc"),
                                xlab="position", ylab="prob. of mutation",
                                main=paste0("mutation trend:", pattern),
                                cex.axis=unit(1, "npc"),
                                cex.main=unit(1, "npc"))

  structure.control <- modifyList(structure.control.default, structure.control)
  graph.control <- modifyList(graph.control.default, graph.control)
  topics.control <- modifyList(topics.control.default, topics.control)

  mat_reduced <- filter_by_pos_pattern(mat, max_pos = max_pos,
                                       pattern = pattern)

  signature_set <- colnames(mat_reduced)

  if(gom_method == "full"){
    message("Fitting the Grade of Membership Model - full version - due to Matt Taddy")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(gom_method == "independent"){
    message("Fitting the Grade of Membership Model - full version - due to Y. Shiraichi and M. Stephens")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "independent", signatures = signature_set), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(is.null(labs)){
    labs <- rownames(mat_reduced)
  }

  if(is.null(levels)){
    levels <- unique(labs)
  }

  message ("Structure plot and Logo plot representations of clusters")

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}

  omega <- topic_clus$omega
  annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(labs, levels = levels)
  )

  rownames(omega) <- annotation$sample_id

  plot.new()
  grid.newpage()
  do.call(CountClust::StructureGGplot, append(list(omega= omega,
                                       annotation = annotation,
                                       palette = topic_cols),
                                  structure.control))
  ggplot2::ggsave(paste0(output_dir, "structure.png"), width = structure_width,
                  height = structure_height)

  theta_pool <- topic_clus$theta
  pos <- sapply(rownames(theta_pool), function(x) return(strsplit(x, "_")[[1]][2]))
  rownames(theta_pool) <- pos
  max_prob <- max(theta_pool)

  theta_pool_mod <- theta_pool[match(1:max_pos, as.numeric(rownames(theta_pool))),]

  plot(1, type="n", axes=F, xlab="", ylab="")

  for(l in 1:dim(theta_pool_mod)[2]){
    png(paste0(output_dir, "plot_clus_", l, ".png"), width=output_width, height = output_height)
    do.call(plot_graph_2, c(list(probs = theta_pool_mod[,l],
                                 max_pos = max_pos,
                                 max_prob = max_prob),
                            graph.control))
    dev.off()
  }

  graphics.off()
  message("Finished")
}



plot_graph_2 <- function(probs, max_pos, max_prob, col="red",
                         cex=unit(1, "npc"), pch=unit(16,"npc"),
                         xlab="position", ylab="prob. of mutation",
                         main="",
                         cex.axis=unit(1, "npc"),
                         cex.main=unit(1, "npc")){
  # if (length(probs) != max_pos){
  #   stop(cat('probability vector must be of length ', max_pos))
  # }
  par(font.axis = 2)
  plot(as.numeric(names(probs)), probs/max_prob, xlim = c(0, max_pos), ylim=c(0,1),
       type = "b", xaxt = "n", yaxt = "n", cex = cex, pch=pch, col=col, main=main,
       cex.main=cex.main, ylab="", xlab="")
  axis(side = 1, at = floor(seq(1, max_pos, length.out=5)), cex.axis = cex.axis, lwd.ticks = 1, tck=-0.05,
       cex.lab=2, mgp=c(2.5, 0.5, 0))
  title(xlab = xlab, mgp=c(2.5,1.5,0), cex.lab=1.8)
  ylimit <- c(0.0, 0.5, 1.0)*max_prob
  axis(side = 2, at = c(0.0, 0.5, 1.0), labels = round(ylimit,2), cex.axis = cex.axis, lwd.ticks=1, tck=-0.05,
       cex.lab=2, mgp=c(2.5, 0.5, 0))
  title(ylab = ylab, mgp=c(2.5,1,0), cex.lab=1.8)
}

#########################  aRchaic_cluster_beta  type 6  #############################

aRchaic_cluster_beta_wo_strand = function(mat,
                                          K,
                                          tol=10,
                                          max_pos = 20,
                                          labs = NULL,
                                          levels = NULL,
                                          flanking_bases = 1,
                                          gom_method = "independent",
                                          topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                                         "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                                         "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                                          structure.control = list(),
                                          logo.control = list(),
                                          topics.control = list(),
                                          output_dir = NULL,
                                          structure_width = 5,
                                          structure_height = 8,
                                          inflation = rep(2,1,2),
                                          output_width = 1200,
                                          output_height = 700){
  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/")
  }

  topics.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                                 nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                                 light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)
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
                               renyi_alpha = 5, inflation_factor = c(2,1,2),
                               pop_names = paste0("Cluster : ", 1:K),
                               logoport_x = 0.26,
                               logoport_y= 0.50,
                               logoport_width= 0.28,
                               logoport_height= 0.40,
                               lineport_x = 0.95,
                               lineport_y=0.40,
                               lineport_width=0.32,
                               lineport_height=0.28,
                               breaklogoport_x = 0.95,
                               breaklogoport_y = 0.40,
                               breaklogoport_width=0.35,
                               breaklogoport_height=0.40,
                               output_width = output_width,
                               output_height = output_height)

  structure.control <- modifyList(structure.control.default, structure.control)
  logo.control <- modifyList(logo.control.default, logo.control)
  topics.control <- modifyList(topics.control.default, topics.control)

  mat_reduced <- filter_out_strand(mat)

  signature_set <- colnames(mat_reduced)
  sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:(flanking_bases+5)])))
  new_sig_split <- matrix(0, dim(sig_split)[1], 3);
  new_sig_split[,1] <- sig_split[,1]
  new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse="")))
  new_sig_split[,3] <- sig_split[,(flanking_bases+5)]

  pos <- sapply(signature_set, function(x) return(strsplit(x, "_")[[1]][3]))

  pos <- as.numeric(pos)
  pos <- pos - min(pos)
  pos <- factor(pos, levels = 0:21)

  strand_break <- sapply(signature_set, function(x) return(strsplit(x, "_")[[1]][2]))
  strand_break <- factor(strand_break)
  signature <- new_sig_split

  signature_pos_strand_break <- cbind.data.frame(new_sig_split, strand_break, pos)

  if(gom_method == "full"){
    message("Fitting the Grade of Membership Model - full version - due to Matt Taddy")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(gom_method == "independent"){
    message("Fitting the Grade of Membership Model - full version - due to Y. Shiraichi and M. Stephens")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "independent", signatures = signature_pos_strand_break), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(is.null(labs)){
    labs <- rownames(mat_reduced)
  }

  if(is.null(levels)){
    levels <- unique(labs)
  }

  message ("Structure plot and Logo plot representations of clusters")

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}

  omega <- topic_clus$omega
  annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(labs, levels = levels)
  )

  rownames(omega) <- annotation$sample_id

  plot.new()
  grid.newpage()
  do.call(CountClust::StructureGGplot, append(list(omega= omega,
                                       annotation = annotation,
                                       palette = topic_cols),
                                  structure.control))
  ggplot2::ggsave(paste0(output_dir, "structure.png"), width = structure_width,
                  height = structure_height)

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}
  plot.new()
  do.call(damageLogo_three, append(list(theta_pool = topic_clus$theta,
                                        output_dir = output_dir),
                                   logo.control))
  graphics.off()

  message("Finished")

}

########################  aRchaic_cluster_beta  Type 7  #############################

aRchaic_cluster_beta_wo_strand_break = function(mat,
                                                K,
                                                tol=10,
                                                max_pos = 20,
                                                labs = NULL,
                                                levels = NULL,
                                                flanking_bases = 1,
                                                gom_method = "independent",
                                                topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                                               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                                               "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                                                structure.control = list(),
                                                logo.control = list(),
                                                topics.control = list(),
                                                output_dir = NULL,
                                                structure_width = 5,
                                                structure_height = 8,
                                                inflation = rep(2,1,2),
                                                output_width = 1200,
                                                output_height = 700){

  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/")
  }

  topics.control.default <- list(bf = FALSE, kill = 2, ord = TRUE, verb = 1, admix = TRUE,
                                 nbundles = 1, use_squarem = FALSE, init.adapt = FALSE, type = "full",
                                 light = 1, method_admix = 1, sample_init = TRUE, tmax = 10000)
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
                               max_pos = 20, flanking_bases=flanking_bases,
                               yscale_change = TRUE, xaxis=TRUE,
                               yaxis=TRUE, xlab = " ", xaxis_fontsize=20,
                               xlab_fontsize=10, title_aligner = 11,
                               y_fontsize=20, title_fontsize = 35,
                               mut_width=2, start=0.0001,
                               renyi_alpha = 5, inflation_factor = c(2,1,2),
                               pop_names = paste0("Cluster : ", 1:K),
                               logoport_x = 0.25,
                               logoport_y= 0.50,
                               logoport_width= 0.28,
                               logoport_height= 0.40,
                               lineport_x = 0.9,
                               lineport_y=0.40,
                               lineport_width=0.32,
                               lineport_height=0.28,
                               barport_x = 0.72,
                               barport_y=0.60,
                               barport_width=0.35,
                               barport_height=0.25,
                               output_width = output_width,
                               output_height = output_height)

  structure.control <- modifyList(structure.control.default, structure.control)
  logo.control <- modifyList(logo.control.default, logo.control)
  topics.control <- modifyList(topics.control.default, topics.control)

  mat_reduced <- filter_out_strand_break(mat)

  signature_set <- colnames(mat_reduced)
  sig_split <- t(sapply(1:length(signature_set), function(x) return(strsplit(signature_set[x], "")[[1]][1:(flanking_bases+5)])))
  new_sig_split <- matrix(0, dim(sig_split)[1], 3);
  new_sig_split[,1] <- sig_split[,1]
  new_sig_split[,2] <- sapply(1:length(signature_set), function(x) return(paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse="")))
  new_sig_split[,3] <- sig_split[,(flanking_bases+5)]

  pos <- sapply(signature_set, function(x) return(strsplit(x, "_")[[1]][3]))

  pos <- as.numeric(pos)
  pos <- pos - min(pos)
  pos <- factor(pos, levels = 0:21)

  strand <- sapply(signature_set, function(x) return(strsplit(x, "_")[[1]][2]))
  strand <- factor(strand)
  signature <- new_sig_split

  signature_pos_strand <- cbind.data.frame(new_sig_split, strand, pos)


  if(gom_method == "full"){
    message("Fitting the Grade of Membership Model - full version - due to Matt Taddy")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "full", signatures = NULL), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(gom_method == "independent"){
    message("Fitting the Grade of Membership Model - full version - due to Y. Shiraichi and M. Stephens")
    suppressWarnings(topic_clus <- do.call(maptpx::topics, append(list(counts = mat_reduced, K=K, tol=tol, model = "independent", signatures = signature_pos_strand), topics.control)))
    save(topic_clus, file = paste0(output_dir, "model.rda"))
  }

  if(is.null(labs)){
    labs <- rownames(mat_reduced)
  }

  if(is.null(levels)){
    levels <- unique(labs)
  }

  message ("Structure plot and Logo plot representations of clusters")

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}

  omega <- topic_clus$omega
  annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(labs, levels = levels)
  )

  rownames(omega) <- annotation$sample_id

  plot.new()
  grid.newpage()
  do.call(CountClust::StructureGGplot, append(list(omega= omega,
                                       annotation = annotation,
                                       palette = topic_cols),
                                  structure.control))
  ggplot2::ggsave(paste0(output_dir, "structure.png"), width = structure_width,
                  height = structure_height)

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}
  plot.new()
  do.call(damageLogo_two, append(list(theta_pool = topic_clus$theta,
                                      output_dir = output_dir),
                                 logo.control))
  graphics.off()

  message("Finished")

}

damage.ic<-function(pwm, alpha=1, inflation_factor = c(1,1,1)) {
  if(length(inflation_factor) != ncol(pwm[[1]])){
    stop("inflation factor vector size
         must equal to the number of sites - flanking bases + mismatch")
  }
  npos<-ncol(pwm[[1]])
  ic<- matrix(0, npos, length(pwm))

  for(i in 1:npos){
    mat <- numeric()
    for(j in 1:length(pwm)){
      mat <- cbind(mat, pwm[[j]][,i])
    }
    mat_clean <- mat[rowSums(mat) != 0,]
    ic[i,] <- inflation_factor[i]*ic_computer_2(mat_clean, alpha)
  }

  return(ic)
  }
