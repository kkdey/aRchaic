


#############   aRchaic_cluster_beta_1 ()   ############################

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

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}

  omega <- topic_clus$omega
  annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(labs, levels = levels)
  )

  rownames(omega) <- annotation$sample_id

  plot.new()
  grid.newpage()
  do.call(StructureGGplot, append(list(omega= omega,
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
