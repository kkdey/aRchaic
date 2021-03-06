#' @title STRUCTURE plot and logo plot representation of the clusterring
#' model from archaic_fit
#'
#' @description Takes the clustering model fit  from \code{archaic_fit}
#' as an input and plots the clusters using as EDLogo plots (Dey et al 2018)
#' and the proportional mixing of the clusters for each sample using a
#' STRUCTURE plot (Pritchard et al 2000, Dey2017) representation.
#'
#' @param model Fitted model from \code{archaic_fit}.
#' @param topic_cols A vector of color assignment to the clusters/topics used
#'                   for cluster representation in the STRUCTURE (Rosenberg2002)
#'                   representation in \code{archaic_plot()}.
#' @param background if equals "modern", as in the default, compares enrichment
#'                   of mismatch features against a modern background - else
#'                   uses a background with equal probability of all mismatch
#'                   features.
#' @param structure.control The control or tuning parameters for the STRUCTURE
#'                          plot representation \code{structure.pdf} output
#'                          of \code{archaic_plot} (Dey2017)
#' @param logo.control The control or tuning parameters for the EDLogo plots
#'                          representation \code{logo_clus_*.pdf} output
#'                          of \code{archaic_plot} (Dey2018)
#' @param output_dir The path/directory where to save the output plots
#'
#' @references
#'     Rosenberg2002.
#'     Rosenberg, N.A., Pritchard, J.K., Weber, J.L., Cann, H.M., Kidd, K.K.,
#'     Zhivotovsky, L.A. and Feldman, M.W., 2002. Genetic structure of
#'     human populations. science, 298(5602), pp.2381-2385.
#'
#'     Dey2017.
#'     Dey, K.K., Hsiao, C.J. and Stephens, M., 2017. Visualizing the structure
#'     of RNA-seq expression data using grade of membership models.
#'     PLoS genetics, 13(3), p.e1006599.
#'
#'     Dey2018.
#'     Dey, K.K., Xie, D. and Stephens, M., 2017. A new sequence logo plot to
#'     highlight enrichment and depletion. bioRxiv, p.226597.
#'
#'     Pritchard2002.
#'     Pritchard, J.K., Stephens, M. and Donnelly, P., 2000.
#'     Inference of population structure using multilocus genotype data.
#'     Genetics, 155(2), pp.945-959.
#'
#' @return Returns a \code{structure.pdf} and as many logo plots of the form
#'         \code{logo_clus_k.pdf} for each cluster k, in the output path
#'         provided \code{output_dir}.
#'
#' @keywords structure EDLogo
#' @importFrom CountClust StructureGGplot
#' @importFrom Logolas get_viewport_logo nlogomaker
#' @import ggplot2
#' @export

archaic_plot <- function(model,
                         topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                        "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                        "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                         background = "modern",
                         structure.control = list(),
                         logo.control = list(),
                         output_dir = NULL){

  labs <- model$labs
  K <- dim(model$omega)[2]
  if(is.null(levels(labs))) {
    levels <- unique(labs)
  }else{
    levels <- levels(labs)
  }
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
                                    legend_text_size = 8,
                                    structure_width = 5,
                                    structure_height = 8)

  logo.control.default <- list(max_pos = 20, flanking_bases=1,
                               base_probs_list = NULL,
                               clip = 0,
                               mut_ranges = c(0,0),
                               break_ranges = c(0,0),
                               logoport_x = 0.7,
                               logoport_y= 0.50,
                               logoport_width= 1, logoport_height= 1.3,
                               breaklogoport_x = 0.60,
                               breaklogoport_y = 0.467, breaklogoport_width=0.76,
                               breaklogoport_height=1.25, lineport_x = 0.65, lineport_y=0.53,
                               lineport_width=0.8, lineport_height=1.4, panelname_x = 0.75,
                               panelname_y= 0.6, panelname_width= 0.3, panelname_height= 0.3,
                               mutlogo.control = list(main_fontsize = 25,
                                                      control = list(npc_units_main = 0.985,
                                                                     lines_units_main = 1)),
                               breaklogo.control = list(main_fontsize = 25,
                                                        control = list(npc_units_main = 0.98,
                                                                       lines_units_main = 1)),
                               output_width = 20, output_height = 7)

  if(background == "null"){
    logo.control.default$base_probs_list = NULL
  }else{
    data("base_probs_moderns")
    logo.control.default$base_probs_list = base_probs_moderns
    logo.control.default$mut_ranges = c(1, 1)
    logo.control.default$break_ranges = c(1, 1)
  }


  structure.control <- modifyList(structure.control.default, structure.control)
  logo.control <- modifyList(logo.control.default, logo.control)


  structure.control.two <- structure.control
  structure.control.two$structure_height = NULL
  structure.control.two$structure_width = NULL


  message ("Structure plot and Logo plot representations of clusters")

  omega <- model$omega
  annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(labs, levels = levels)
  )

  if(is.null(output_dir)){ output_dir <- paste0(getwd(),"/")}else{
    if(regmatches(output_dir,regexpr(".$", output_dir)) != "/"){output_dir <- paste0(output_dir, "/")}
  }

  plot.new()
  grid::grid.newpage()
  do.call(CountClust::StructureGGplot, append(list(omega= omega,
                                                   annotation = annotation,
                                                   palette = topic_cols),
                                              structure.control.two))
  ggplot2::ggsave(paste0(output_dir, "structure.pdf"),
                  width =  structure.control$structure_width,
                  height = structure.control$structure_height)

  ###################   Logo plot representation  #########################


  plot.new()
  do.call(Logo_aRchaic_cluster, append(list(theta_pool = model$theta,
                                            output_dir = output_dir,
                                            topic_cols = topic_cols),
                                       logo.control))
  graphics.off()

  message ("Finished")
}


# @title Builds damage logo plots based on full mutation signatures (mutation, flanking base, position, strand and strand break)
#
# @description Damage Logo plots for each cluster from the GoM model fit. It showcases the
# different mutational features - for example, mutation, flanking base, position on read,
# strand and strand break information for each cluster.
#
# @param theta_pool The theta matrix obtained from running the grade of membership model that stores for each cluster, the
# probability distribution over all the mutational signatures.
# @param max_pos The maximum distance from the end of the read upto which mutations are considered.
# @param flanking_bases The number of flanking bases of the mutational signature.
# @param mutlogo.control The control parameters for the mismatch and flanking bases logo.
# @param breaklogo.control The control parameters for the logo for strand break.
# @param logoport_x the X-axis position of the plot window for the logo plot
# @param logoport_y the Y-axis position of the plot window for the logo plot
# @param logoport_width the width of the plot window for the logo plot
# @param logoport_height the width of the plot window for the logo plot
# @param lineport_x the X-axis position of the plot window for the mutational profile line plot.
# @param lineport_y the Y-axis position of the plot window for the mutational profile line plot.
# @param lineport_width the width of the plot window for the mutational profile line plot.
# @param lineport_height the width of the plot window for the mutational profile line plot.
# @return Returns logo plots for each cluster
# @param breaklogoport_x the X-axis position of the plot window for strand break logo plot.
# @param breaklogoport_y the Y-axis position of the plot window for the strand break logo plot.
# @param breaklogoport_width the width of the plot window for the strand break logo plot.
# @param breaklogoport_height the width of the plot window for the strand break logo plot.
# @param output_dir The directory where the logo plot will be saved.
# @param output_width The width of the logo plot figure.
# @param output_height the height of the logo plot figure.
#
# @return Returns logo plot for each cluster
#
# @import grid
# @import gridBase
#
# @export


Logo_aRchaic_cluster <- function(theta_pool,
                                 max_pos = 20,
                                 flanking_bases=1,
                                 mutlogo.control = list(),
                                 breaklogo.control = list(),
                                 base_probs_list = NULL,
                                 clip = 0,
                                 mut_ranges  = c(0, 0),
                                 break_ranges = c(0, 0),
                                 logoport_x = 0.7,
                                 logoport_y= 0.5,
                                 logoport_width= 1.2,
                                 logoport_height= 1.1,
                                 breaklogoport_x = 0.5,
                                 breaklogoport_y = 0.4,
                                 breaklogoport_width=0.7,
                                 breaklogoport_height=1,
                                 lineport_x = 0.4,
                                 lineport_y=0.5,
                                 lineport_width=1,
                                 lineport_height=1,
                                 panelname_x = 0.8,
                                 panelname_y= 0.6,
                                 panelname_width= 0.3,
                                 panelname_height= 0.3,
                                 topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                                "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                                "brown4","darkorchid","magenta","yellow", "azure1","azure4"),
                                 output_dir = NULL,
                                 filename = NULL,
                                 output_width = 18,
                                 output_height = 7){

  library(grid)
  library(gridBase)
  library(ggplot2)


  if(is.null(output_dir)){output_dir <- getwd();}
  flag <- 0
  if(dim(theta_pool)[2] == 1){
    flag = 1
    theta_pool <- cbind(theta_pool, theta_pool)
    colnames(theta_pool) <- c("sample1", "sample2")
  }
  signature_set <- rownames(theta_pool)
  signature_patterns <- substring(signature_set, 1, 4+2*flanking_bases)
  library(dplyr)
  theta2 <- dplyr::tbl_df(data.frame(theta_pool)) %>%
    dplyr::mutate(sig = signature_patterns) %>%
    dplyr::group_by(sig) %>%
    dplyr::summarise_all(funs(sum)) %>%
    as.data.frame()
  rownames(theta2) <-  theta2[,1]
  theta2 <- theta2[,-1, drop=FALSE]

  breakbase <- substring(signature_set, 6+2*flanking_bases,  6+2*flanking_bases)

  theta_break <- dplyr::tbl_df(data.frame(theta_pool)) %>%
    dplyr::mutate(sig = breakbase) %>%
    dplyr::group_by(sig) %>%
    dplyr::summarise_all(funs(sum)) %>%
    as.data.frame()
  rownames(theta_break) <- theta_break[,1]
  theta_break <- theta_break[,-1]

  theta_break <- theta_break[match(c("A", "C", "G", "T"), rownames(theta_break)),]
  breaks_theta <- theta_break

  sig_names <- rownames(theta_pool)

  # prob_mutation <- filter_by_pos(t(theta_pool), max_pos = max_pos)
  prob_mutation <- filter_signatures_only_location(t(theta_pool),
                                                   max_pos = max_pos, flanking_bases = flanking_bases)

  prob_mutation <- t(apply(prob_mutation, 1, function(x) {
    y <- x[!is.na(x)];
    return(y/sum(y))
  }))
  clipped_bases <- setdiff(0:20, as.numeric(colnames(prob_mutation)))

  max_prob <- max(prob_mutation);
 # clipped_bases <- setdiff(0:20, as.numeric(colnames(prob_mutation)))

  if(is.null(base_probs_list)){
    prob_limits = c(round(min(prob_mutation, na.rm=TRUE), 2)-0.01, round(max(prob_mutation, na.rm=TRUE), 2) + 0.01)
    prob_breaks = c(0, round(min(prob_mutation),2)-0.01,
                    round(0.5*(min(prob_mutation, na.rm=TRUE)+max(prob_mutation, na.rm=TRUE)), 2),
                    round(max(prob_mutation), 2)+0.01)
  }else{
    if(length(clipped_bases) > 0){
      prob1_mutation <- prob_mutation - t(replicate(dim(prob_mutation)[1], as.numeric(base_probs_list[[(2 * flanking_bases + 3)]])[-(clipped_bases+1)]))
    }else{
      prob1_mutation <- prob_mutation - t(replicate(dim(prob_mutation)[1], as.numeric(base_probs_list[[(2 * flanking_bases + 3)]])))
    }
    colnames(prob1_mutation) <- colnames(prob_mutation)
    prob_mutation_after_clipping <- prob1_mutation[,(clip+1):(dim(prob1_mutation)[2])]
    prob_limits = c(round(min(prob_mutation_after_clipping), 2)-0.01,
                    round(max(prob_mutation_after_clipping), 2) + 0.01)
    prob_breaks = c(0, round(min(prob_mutation_after_clipping),2)-0.01,
                    round(0.5*(min(prob_mutation_after_clipping)+max(prob_mutation_after_clipping)), 2),
                    round(max(prob_mutation_after_clipping), 2)+0.01)
  }

  sig_split <- do.call(rbind,
                       lapply(sig_names,
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

  for(l in 1:dim(theta_pool)[2]){
    prop_patterns_list[[l]] <- numeric();
    for(j in 1:ncol(new_sig_split)){
      temp2 <- tapply(theta_pool[,l], factor(new_sig_split[,j], levels=c("A", "C", "G", "T",
                                                                         "C->T", "C->A", "C->G",
                                                                         "T->A", "T->C", "T->G")), sum)

      prop_patterns_list[[l]] <- cbind(prop_patterns_list[[l]], temp2)
    }
  }

  grob_list <- list()
  for(l in 1:length(prop_patterns_list)){
    pdf(paste0(output_dir, "logo_clus_", l, ".pdf"), width=output_width, height = output_height)
    damageLogo_six.skeleton(pwm = prop_patterns_list[[l]],
                            probs = prob_mutation[l,],
                            breaks_theta_vec = breaks_theta[,l, drop=FALSE],
                            prob_limits = prob_limits,
                            prob_breaks = prob_breaks,
                            mutlogo.control = mutlogo.control,
                            breaklogo.control = breaklogo.control,
                            background = base_probs_list,
                            clip = clip,
                            mut_ranges = mut_ranges,
                            break_ranges = break_ranges,
                            flanking_bases = flanking_bases,
                            logoport_x = logoport_x,
                            logoport_y= logoport_y,
                            logoport_width= logoport_width,
                            logoport_height= logoport_height,
                            breaklogoport_x = breaklogoport_x,
                            breaklogoport_y = breaklogoport_y,
                            breaklogoport_width=breaklogoport_width,
                            breaklogoport_height=breaklogoport_height,
                            lineport_x = lineport_x,
                            lineport_y= lineport_y,
                            lineport_width=lineport_width,
                            lineport_height=lineport_height,
                            panelname_x = panelname_x,
                            panelname_y= panelname_y,
                            panelname_width= panelname_width,
                            panelname_height= panelname_height,
                            panelname_title = paste0("cluster ", l),
                            panelname_color = topic_cols[l])
    dev.off()
  }
}

damageLogo_six.skeleton <- function(pwm,
                                    probs,
                                    breaks_theta_vec,
                                    prob_limits,
                                    prob_breaks,
                                    mutlogo.control = list(),
                                    breaklogo.control = list(),
                                    background = NULL,
                                    clip = 0,
                                    mut_ranges = c(0, 0),
                                    break_ranges = c(0, 0),
                                    flanking_bases = 1,
                                    logoport_x = 0.7,
                                    logoport_y= 0.5,
                                    logoport_width= 1.2,
                                    logoport_height= 1.1,
                                    breaklogoport_x = 0.5,
                                    breaklogoport_y = 0.4,
                                    breaklogoport_width=0.7,
                                    breaklogoport_height=1,
                                    lineport_x = 0.4,
                                    lineport_y=0.5,
                                    lineport_width=1,
                                    lineport_height=1,
                                    panelname_x = 0.8,
                                    panelname_y= 0.6,
                                    panelname_width= 0.3,
                                    panelname_height= 0.3,
                                    panelname_title = "",
                                    panelname_color = "black"){

  mut_lowrange = mut_ranges[1]
  mut_uprange =  mut_ranges[2]

  break_lowrange = break_ranges[1]
  break_uprange =  break_ranges[2]

  mutlogo.control.default <- list(ic = FALSE,
                                  score = "log",
                                  total_chars = c("A", "B", "C", "D", "E", "F", "G",
                                                  "H", "I", "J", "K", "L", "M", "N", "O",
                                                  "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y",
                                                  "Z", "zero", "one", "two",
                                                  "three", "four", "five", "six", "seven", "eight",
                                                  "nine", "dot", "comma",
                                                  "dash", "colon", "semicolon", "leftarrow", "rightarrow"),
                                  frame_width=c(1,2,1), yscale_change=TRUE,
                                  pop_name = "Mismatch and \n flanking base composition",
                                  addlogos = NULL, addlogos_text = NULL, newpage = FALSE,
                                  yrange = NULL, xaxis=TRUE, yaxis=TRUE, xaxis_fontsize=23,
                                  y_fontsize=22, main_fontsize = 25,
                                  xlab_fontsize=22,
                                  start=0.001, xlab = "", ylab = "Enrichment Score",
                                  col_line_split="grey80", control = list(epsilon=0.25,gap_ylab=3.5, gap_xlab = 4,
                                                                          round_off = 1, posbins = 3,
                                                                          negbins = 3,
                                                                          lowrange = mut_lowrange,
                                                                          uprange = mut_uprange,
                                                                          size_port = 1, symm = FALSE,
                                                                          npc_units_main = 1.5,
                                                                          lines_units_main = 1))

  breaklogo.control.default <- list( ic = FALSE,
                                     score = "log",
                                     total_chars = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O",
                                                     "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "zero", "one", "two",
                                                     "three", "four", "five", "six", "seven", "eight", "nine", "dot", "comma",
                                                     "dash", "colon", "semicolon", "leftarrow", "rightarrow"),
                                     frame_width=NULL, yscale_change=TRUE,
                                     pop_name = "Base composition \n 5' of  strand break \n",
                                     addlogos = NULL, addlogos_text = NULL, newpage = FALSE,
                                     yrange = NULL, xaxis=FALSE, yaxis=TRUE, xaxis_fontsize=10,
                                     xlab_fontsize=22, y_fontsize=22, main_fontsize=25,
                                     start=0.001, xlab = "", ylab = "Enrichment Score",
                                     col_line_split="white", control = list(gap_ylab=3.5, epsilon = 0.01,
                                                                            round_off = 1, symm = TRUE,
                                                                            npc_units_main = 1.5,
                                                                            lines_units_main = 1,
                                                                            lowrange = break_lowrange,
                                                                            uprange = break_uprange))
  mutlogo.control <- modifyList(mutlogo.control.default, mutlogo.control)
  breaklogo.control <- modifyList(breaklogo.control.default, breaklogo.control)


  cols = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category ==
                                         'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, cols$maxcolors, rownames(cols)))
  total_chars = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O",
                  "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "zero", "one", "two",
                  "three", "four", "five", "six", "seven", "eight", "nine", "dot", "comma",
                  "dash", "colon", "semicolon", "leftarrow", "rightarrow")

  set.seed(20)
  cols2 <- sample(col_vector, length(total_chars), replace=FALSE)
  cols2[match(c("A", "C", "G", "T"), total_chars)] <- c(RColorBrewer::brewer.pal(4,name ="Spectral"))
  color_profile_1 <- list("type" = "per_symbol",
                          "col" = cols2)

  if(!is.null(background)){
    base_probs_list <- background
    base_probs_mat <- matrix(NA, 10, 3)
    rownames(base_probs_mat) <- c("A", "C", "G", "T", "C->A", "C->G", "C->T", "T->A", "T->C", "T->G")
    for(l in 1:(2*flanking_bases+1)){
      base_probs_mat[match(names(base_probs_list[[l]]), rownames(base_probs_mat)),l] <- as.numeric(base_probs_list[[l]])
    }
    base_probs_mat <- base_probs_mat[match(rownames(pwm), rownames(base_probs_mat)),]
  }

  pwm1 <- pwm
  rownames(pwm1)[match(c("C->A", "C->G", "C->T",
                         "T->A", "T->C", "T->G"), rownames(pwm1))] <- c("C>A", "C>G", "C>T",
                                                                        "T>A", "T>C", "T>G")
  colnames(pwm1) <- c("5' \n flank", "mismatch", "3' \n flank")


  if(!is.null(background)){
    rownames(base_probs_mat)[match(c("C->A", "C->G", "C->T",
                                     "T->A", "T->C", "T->G"), rownames(base_probs_mat))] <- c("C>A", "C>G", "C>T",
                                                                                              "T>A", "T>C", "T>G")
    colnames(base_probs_mat) <- colnames(pwm1)
  }


  color_profile_2 = list("type" = "per_row",
                         "col" = RColorBrewer::brewer.pal(4,name ="Spectral"))



  Logolas::get_viewport_logo(1, 4, widths.val = c(3,5,5,5))

  seekViewport(paste0("plotlogo", 1))
  vp1 <- viewport(x=panelname_x, y=panelname_y, width=panelname_width, height=panelname_height)
  df <- data.frame(x1 = -0.5, x2 = 0.5, y1 = -1, y2 = 1)
  p <- ggplot() +
    geom_rect(data=df, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),  fill=panelname_color, color = panelname_color, alpha=1) +
    ggtitle(panelname_title) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          plot.title = element_text(size = 30, face = "bold", hjust = 0.5))
  print(p, vp = vp1)
  upViewport(0)

  seekViewport(paste0("plotlogo", 2))
  vp2 <- viewport(x=logoport_x, y=logoport_y, width=logoport_width, height=logoport_height)
  pushViewport(vp2)
  if(!is.null(background)){
    do.call(Logolas::nlogomaker, append(list(table = pwm1,
                                    color_profile = color_profile_1,
                                    bg = base_probs_mat),
                               mutlogo.control))
  }else{
    do.call(Logolas::nlogomaker, append(list(table = pwm1,
                                    color_profile = color_profile_1,
                                    bg = NULL),
                               mutlogo.control))
  }
  upViewport(0)


  if(is.null(background)){

    pos_data <- data.frame(position = as.numeric(names(probs)),
                           val = as.numeric(probs))

    seekViewport(paste0("plotlogo", 3))
    vp3 = viewport(x = lineport_x, y = lineport_y, width=lineport_width, height=lineport_height)
    p <- ggplot(data=pos_data, aes(x=position,y=val)) +
      geom_point(size = 3, aes(colour = "red")) +
      geom_line(aes(colour = "red"))+
      ggtitle("Location of \n mismatch in read" ) +
      theme(plot.title = element_text(lineheight=1.2, margin=margin(0,0,20,0),
                                      hjust = 0.5, size = 25),
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 22), axis.title.y = element_text(size = 22),
            axis.text.x = element_text(colour="black", hjust=0.8, size = 18),
            axis.text.y = element_text(size = 22, hjust=0.8, colour = "black"),
            legend.position="none",
            axis.ticks.length=unit(0.3,"cm"))+
      labs(x="position in read",y="probability of mismatch") +
      scale_x_continuous(limits = c(0, 20))  +
      scale_y_continuous(limits = prob_limits,
                         breaks = prob_breaks)
     # theme_bw() + theme() +
     # theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
     # theme(axis.text.x = element_text(colour="black", hjust=0.8, size = 18),
     #       axis.text.y = element_text(size = 18, hjust=0.8, colour = "black")) +
     # theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
     # theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
     # theme(plot.title = element_text(size = 25)) +
     # theme(plot.title = element_text(margin=margin(b = 30, unit = "pt"))) +
     # theme(plot.title = element_text(lineheight=2, hjust = 0.5))+
     # theme(legend.position="none") +
     # theme(axis.ticks.length=unit(0.3,"cm"))
    # geom_hline(yintercept=0, linetype="dashed")
    print(p, vp = vp3)

  }else{
    bg_pos_vec <- base_probs_list[[(2*flanking_bases+3)]]
    probs1 <- (probs+1e-10) - (bg_pos_vec[match(names(probs), names(bg_pos_vec))]+1e-10)
    # probs1 <- probs1 - median(probs1)
    num_pos <- length(as.numeric(probs1))
    pos_data <- data.frame(position = as.numeric(names(probs)[(clip+1):num_pos]),
                           val = as.numeric(probs1)[(clip+1):num_pos])
    seekViewport(paste0("plotlogo", 3))
    vp3 = viewport(x = lineport_x, y = lineport_y, width=lineport_width, height=lineport_height)
    p <- ggplot(data=pos_data, aes(x=position,y=val)) +
      geom_point(size = 3, aes(colour = "red")) +
      geom_line(aes(colour = "red"))+
      ggtitle("Location of \n mismatch in read" ) +
      theme(plot.title = element_text(lineheight=1.2, margin=margin(0,0,20,0),
                                      hjust = 0.5, size = 25),
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 22), axis.title.y = element_text(size = 22),
            axis.text.x = element_text(colour="black", hjust=0.8, size = 18),
            axis.text.y = element_text(size = 22, hjust=0.8, colour = "black"),
            legend.position="none",
            axis.ticks.length=unit(0.3,"cm"))+
      labs(x="position in read",y="Enrichment in probability") +
      scale_x_continuous(limits = c(0, 20))  +
      scale_y_continuous(limits = prob_limits,
                         breaks = prob_breaks,
                         expand = c(0,0)) +
   #   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
   #   theme(axis.title.x = element_text(size = 22), axis.title.y = element_text(size = 22)) +
  #    theme(axis.text.x = element_text(colour="black", hjust=0.8, size = 22),
  #          axis.text.y = element_text(size = 22, hjust=0.8, colour = "black")) +
  #    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  #    theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  #    theme(plot.title = element_text(size = 25)) +
  #    theme(plot.title = element_text(margin=margin(b = 40, unit = "pt"))) +
  #    theme(plot.title = element_text(hjust = 0.5))+
  #    theme(legend.position="none") +
  #    theme(axis.ticks.length=unit(0.3,"cm")) +
      geom_hline(yintercept=0, linetype="dashed")
    print(p, vp = vp3)

  }


  seekViewport(paste0("plotlogo", 4))
  vp4 <- viewport(x=breaklogoport_x, y=breaklogoport_y, width=breaklogoport_width, height=breaklogoport_height)
  pushViewport(vp4)
  if(!is.null(background)){
    bg_breaks_theta_vec <- matrix(base_probs_list[[(2*flanking_bases+2)]], ncol=1)
    rownames(bg_breaks_theta_vec) <- names(base_probs_list[[(2*flanking_bases+2)]])
    colnames(bg_breaks_theta_vec) <- colnames(breaks_theta_vec)
    do.call(Logolas::nlogomaker, append(list(table = breaks_theta_vec,
                                             color_profile = color_profile_2,
                                             bg = bg_breaks_theta_vec),
                                        breaklogo.control))
  }else{
    do.call(Logolas::nlogomaker, append(list(table = breaks_theta_vec,
                                             color_profile = color_profile_2,
                                             bg = NULL),
                                        breaklogo.control))
  }
  upViewport(0)
}
