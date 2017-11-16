#' @title Builds damage logo plots based on full mutation signatures (mutation, flanking base, position, strand and strand break)
#'
#' @description Damage Logo plots for each cluster from the GoM model fit. It showcases the
#' different mutational features - for example, mutation, flanking base, position on read,
#' strand and strand break information for each cluster.
#'
#' @param theta_pool The theta matrix obtained from running the grade of membership model that stores for each cluster, the
#' probability distribution over all the mutational signatures.
#' @param sig_names The mutational signature names. Defaults to the rownames of the theta matrix above.
#' @param ic.scale A binary variable indicating whether the height of the bars for substitution and flanking bases should be
#'        adjusted by the information criterion.
#' @param max_pos The maximum distance from the end of the read upto which mutations are considered.
#' @param flanking_bases The number of flanking bases of the mutational signature.
#' @param yscale_change A binary variable indicating whether the Y axis scale should be adjusted based on the size of the
#'        logos, defaults to TRUE.
#' @param xlab The labels for X axis.
#' @param xaxis A binary indicating whether the X axis of the logo plot should be shown
#' @param yaxis A binary indicating whether the Y axis of the logo plot should be shown
#' @param xaxis_fontsize The fontsize of the X axis ticks.
#' @param xlab_fontsize The fontsize of the X axis labels.
#' @param y_fontsize The fontsize of the Y axis ticks.
#' @param mut_width Thw width of the bar for the mutation at the center.
#' @param start The starting point of the stacking of logos on the Y axis. Should be close to 0, defau;ts to 0.0001.
#' @param renyi_alpha The entropy scale for the Renyi entropy on the flanking bases and mutations.
#' @param inflation The inflation scale of flanking bases with respect to mutation.
#'                  Will be a 3 length vector. Defaults to c(1,1,1) implying no inflation.
#'                  c(2,1,2) will mean the flanking bases are 2 times inflated compared to mutation.
#' @param pop_names The title of the plot. Defaults to the cluster labels.
#' @param logoport_x the X-axis position of the plot window for the logo plot
#' @param logoport_y the Y-axis position of the plot window for the logo plot
#' @param logoport_width the width of the plot window for the logo plot
#' @param logoport_height the width of the plot window for the logo plot
#' @param lineport_x the X-axis position of the plot window for the mutational profile line plot.
#' @param lineport_y the Y-axis position of the plot window for the mutational profile line plot.
#' @param lineport_width the width of the plot window for the mutational profile line plot.
#' @param lineport_height the width of the plot window for the mutational profile line plot.
#' @return Returns logo plots for each cluster
#' @param breaklogoport_x the X-axis position of the plot window for strand break logo plot.
#' @param breaklogoport_y the Y-axis position of the plot window for the strand break logo plot.
#' @param breaklogoport_width the width of the plot window for the strand break logo plot.
#' @param breaklogoport_height the width of the plot window for the strand break logo plot.
#' @param output_dir The directory where the logo plot will be saved.
#' @param output_width The width of the logo plot figure.
#' @param output_height the height of the logo plot figure.
#'
#' @return Returns logo plot for each cluster
#'
#' @import grid
#' @import gridBase
#'
#' @export


damageLogo_six <- function(theta_pool,
                            sig_names = NULL,
                            max_pos = 20,
                            flanking_bases=1,
                            mutlogo.control = list(),
                            breaklogo.control = list(),
                            base_probs_list = NULL,
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
                            output_dir = NULL,
                            output_width = 1200,
                            output_height = 700){
  library(grid)
  library(gridBase)

  if(length(inflation_factor)==1){
    inflation_factor <- rep(inflation_factor, dim(theta_pool)[2])
  }
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

  indices_minus <- grep("_-_", signature_set)
  strand_theta <- data.frame("minus" = colSums(theta_pool[indices_minus,]),
                             "plus" = colSums(theta_pool[-indices_minus,]))

  if(flag == 1){
    strand_theta <- data.frame("minus" = colSums(matrix(theta_pool[indices_minus,])),
                               "plus" = colSums(matrix(theta_pool[-indices_minus,])))
    strand_theta <- strand_theta/2;
  }
  breakbase <- substring(signature_set, 8+2*flanking_bases,  8+2*flanking_bases)

  theta_break <- dplyr::tbl_df(data.frame(theta_pool)) %>%
    dplyr::mutate(sig = breakbase) %>%
    dplyr::group_by(sig) %>%
    dplyr::summarise_all(funs(sum)) %>%
    as.data.frame()
  rownames(theta_break) <- theta_break[,1]
  theta_break <- theta_break[,-1]

  theta_break <- theta_break[match(c("A", "C", "G", "T"), rownames(theta_break)),]
  breaks_theta <- theta_break


  if(is.null(sig_names))
    sig_names <- rownames(theta)

  # prob_mutation <- filter_by_pos(t(theta_pool), max_pos = max_pos)
  prob_mutation <- filter_signatures_only_location(t(theta_pool),
                                                   max_pos = max_pos, flanking_bases = flanking_bases)
  prob_mutation <- t(apply(prob_mutation, 1, function(x) {
    y <- x[!is.na(x)];
    return(y/sum(y))
  }))
  max_prob <- max(prob_mutation);

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

  for(l in 1:dim(theta)[2]){
    prop_patterns_list[[l]] <- numeric();
    for(j in 1:ncol(new_sig_split)){
      temp2 <- tapply(theta_pool[,l], factor(new_sig_split[,j], levels=c("A", "C", "G", "T",
                                                                         "C->T", "C->A", "C->G",
                                                                         "T->A", "T->C", "T->G")), sum)

      prop_patterns_list[[l]] <- cbind(prop_patterns_list[[l]], temp2)
    }
  }


  grob_list <- list()
  if(flag == 1){
    l = 1
    png(paste0(output_dir, "logo_", pop_names[l], ".png"), width=output_width, height = output_height)
    damageLogo_six.skeleton(pwm = prop_patterns_list[[l]],
                            probs = prob_mutation[l,],
                            breaks_theta_vec = breaks_theta[,l, drop=FALSE],
                            mutlogo.control = mutlogo.control,
                            breaklogo.control = breaklogo.control,
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
                            lineport_height=lineport_height)
    dev.off()
  }else{
    for(l in 1:length(prop_patterns_list)){
      png(paste0(output_dir, "logo_clus_", l, ".png"), width=output_width, height = output_height)
      damageLogo_six.skeleton(pwm = prop_patterns_list[[l]],
                              probs = prob_mutation[l,],
                              breaks_theta_vec = breaks_theta[,l, drop=FALSE],
                              mutlogo.control = mutlogo.control,
                              breaklogo.control = breaklogo.control,
                              bg = base_probs_mat,
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
                              lineport_height=lineport_height)
      dev.off()
    }
  }
}


damagelogo_six.skeleton <- function(pwm,
                                    probs,
                                    breaks_theta_vec,
                                    mutlogo.control = list(),
                                    breaklogo.control = list(),
                                    bg = NULL,
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
                                    lineport_height=1){

  cols = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category ==
                                         'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, cols$maxcolors, rownames(cols)))
  total_chars = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O",
                  "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "zero", "one", "two",
                  "three", "four", "five", "six", "seven", "eight", "nine", "dot", "comma",
                  "dash", "colon", "semicolon", "leftarrow", "rightarrow")

  set.seed(20)
  color_profile_1 <- list("type" = "per_symbol",
                          "col" = sample(col_vector, length(total_chars), replace=FALSE))

  rownames(tab)[match(c("C->A", "C->G", "C->T",
                        "T->A", "T->C", "T->G"), rownames(tab))] <- c("C>A", "C>G", "C>T",
                                                                      "T>A", "T>C", "T>G")

  color_profile_2 = list("type" = "per_row",
                       "col" = RColorBrewer::brewer.pal(4,name ="Spectral"))

  pos_data <- data.frame(position = as.numeric(names(probs)),
                         val = as.numeric(probs))


  Logolas::get_viewport_logo(1, 3)

  seekViewport(paste0("plotlogo", 1))
  vp1 <- viewport(x=logoport.x, y=logoport.y, width=logoport.width, height=logoport.height)
  pushViewport(vp1)
  do.call(nlogomaker, append(list(table = pwm,
                             color_profile = color_profile_1,
                             bg = base_probs_list,
                             newpage = FALSE),
                             mutlogo.control))
  upViewport(0)

  seekViewport(paste0("plotlogo", 2))
  vp2 <- viewport(x=breaklogoport_x, y=breaklogoport_y, width=breaklogoport_width, height=breaklogoport_height)
  pushViewport(vp2)
  do.call(nlogomaker, append(list(table = breaks_theta_vec,
                                  color_profile = color_profile_2,
                                  col_line_split = "white",
                                  newpage = FALSE),
                             breaklogo.control))
  upViewport(0)

  seekViewport(paste0("plotlogo", 3))

  vp3 = viewport(x = lineport_x, y = lineport_y, width=lineport_width, height=lineport_height)
  p <- ggplot(data=pos_data, aes(x=position,y=val)) + geom_point() +
    ggtitle("Probability of mismatch along the read") +
    labs(x="position from end of read",y="probability of mismatch") +
    scale_x_continuous(limits = c(0, 20))
  print(p, vp = vp3)

}
