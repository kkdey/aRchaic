

#############     aRchaic_cluster_beta_pos ()   ############################


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
  do.call(StructureGGplot, append(list(omega= omega,
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
