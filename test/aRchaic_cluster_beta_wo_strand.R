

##############  aRchaic cluster beta without strand  #####################

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
                               yaxis=TRUE, xlab = " ", xaxis_fontsize=20,
                               xlab_fontsize=10, title_aligner = 11,
                               y_fontsize=20, title_fontsize = 35,
                               mut_width=2, start=0.0001,
                               renyi_alpha = 5, inflation_factor = c(3,1,3),
                               pop_names = paste0("Cluster : ", 1:K),
                               logoport_x = 0.25,
                               logoport_y= 0.50,
                               logoport_width= 0.28,
                               logoport_height= 0.40,
                               lineport_x = 0.9,
                               lineport_y=0.40,
                               lineport_width=0.32,
                               lineport_height=0.28,
                               breaklogoport_x = 0.90,
                               breaklogoport_y = 0.40,
                               breaklogoport_width=0.35,
                               breaklogoport_height=0.50,
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
  do.call(StructureGGplot, append(list(omega= omega,
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
