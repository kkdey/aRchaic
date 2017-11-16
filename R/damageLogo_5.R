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
#' @param barport_x the X-axis position of the plot window for strand bar plot.
#' @param barport_y the Y-axis position of the plot window for the strand bar plot.
#' @param barport_width the width of the plot window for the strand bar plot.
#' @param barport_height the width of the plot window for the strand bar plot.
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



damageLogo_five <- function(theta_pool,
                            sig_names = NULL,
                            ic.scale=TRUE,
                            max_pos = 20,
                            flanking_bases=1,
                            yscale_change = TRUE,
                            xaxis=TRUE,
                            yaxis=TRUE,
                            xlab = " ",
                            xaxis_fontsize=10,
                            xlab_fontsize=20,
                            title_aligner = 15,
                            y_fontsize=20,
                            title_fontsize = 20,
                            mut_width=2,
                            start=0.0001,
                            renyi_alpha = 1,
                            inflation_factor = c(2,1,2),
                            base_probs_list = NULL,
                            pop_names=paste0("Cluster ",1:dim(theta_pool)[2]),
                            logoport_x = 0.25,
                            logoport_y= 0.50,
                            logoport_width= 0.28,
                            logoport_height= 0.40,
                            lineport_x = 0.9,
                            lineport_y=0.40,
                            lineport_width=0.32,
                            lineport_height=0.28,
                            breaklogoport_x = 1.00,
                            breaklogoport_y = 0.40,
                            breaklogoport_width=0.40,
                            breaklogoport_height=0.50,
                            barport_x = 0.58,
                            barport_y=0.60,
                            barport_width=0.25,
                            barport_height=0.40,
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

  ic <- damage.ic(prop_patterns_list, alpha=renyi_alpha,
                  inflation_factor = inflation_factor, base_probs_list = base_probs_list)

  grob_list <- list()
  if(flag == 1){
    l = 1
    png(paste0(output_dir, "logo_", pop_names[l], ".png"), width=output_width, height = output_height)
    damageLogo.pos.str.skeleton(pwm = prop_patterns_list[[l]],
                                probs = prob_mutation[l,],
                                breaks_theta_vec = breaks_theta[,l, drop=FALSE],
                                strand_theta_vec = strand_theta[1,],
                                ic = ic[,l],
                                max_pos = max_pos,
                                max_prob = max_prob,
                                ic.scale = ic.scale,
                                yscale_change = yscale_change,
                                xlab = xlab,
                                xaxis=xaxis,
                                yaxis=yaxis,
                                xaxis_fontsize=xaxis_fontsize,
                                xlab_fontsize=xlab_fontsize,
                                title_aligner = title_aligner,
                                y_fontsize=y_fontsize,
                                title_fontsize = title_fontsize,
                                mut_width=mut_width,
                                start=start,
                                pop_name = pop_names[l],
                                logoport_x = logoport_x,
                                logoport_y= logoport_y,
                                logoport_width= logoport_width,
                                logoport_height= logoport_height,
                                lineport_x = lineport_x,
                                lineport_y= lineport_y,
                                lineport_width=lineport_width,
                                lineport_height=lineport_height,
                                breaklogoport_x = breaklogoport_x,
                                breaklogoport_y = breaklogoport_y,
                                breaklogoport_width=breaklogoport_width,
                                breaklogoport_height=breaklogoport_height,
                                barport_x = barport_x,
                                barport_y = barport_y,
                                barport_width = barport_width,
                                barport_height = barport_height)
    dev.off()
    # }else{
    #   l=1
    #   par(new=TRUE)
    #   grid.newpage()
    #   damageLogo.pos.str.skeleton(pwm = prop_patterns_list[[l]],
    #                               probs = prob_mutation[l,],
    #                               breaks_theta_vec = breaks_theta[,l, drop=FALSE],
    #                               strand_theta_vec = strand_theta[l,],
    #                               ic = ic[,l],
    #                               max_pos = max_pos,
    #                               max_prob = max_prob,
    #                               ic.scale = ic.scale,
    #                               yscale_change = yscale_change,
    #                               xlab = xlab,
    #                               xaxis=xaxis,
    #                               yaxis=yaxis,
    #                               xaxis_fontsize=xaxis_fontsize,
    #                               xlab_fontsize=xlab_fontsize,
    #                               title_aligner = title_aligner,
    #                               y_fontsize=y_fontsize,
    #                               title_fontsize = title_fontsize,
    #                               mut_width=mut_width,
    #                               start=start,
    #                               pop_name = pop_names[l],
    #                               logoport_x = logoport_x,
    #                               logoport_y= logoport_y,
    #                               logoport_width= logoport_width,
    #                               logoport_height= logoport_height,
    #                               lineport_x = lineport_x,
    #                               lineport_y= lineport_y,
    #                               lineport_width=lineport_width,
    #                               lineport_height=lineport_height,
    #                               breaklogoport_x = breaklogoport_x,
    #                               breaklogoport_y = breaklogoport_y,
    #                               breaklogoport_width=breaklogoport_width,
    #                               breaklogoport_height=breaklogoport_height,
    #                               barport_x = barport_x,
    #                               barport_y = barport_y,
    #                               barport_width = barport_width,
    #                               barport_height = barport_height)
    # }
  } else {
    for(l in 1:length(prop_patterns_list)){
        png(paste0(output_dir, "logo_clus_", l, ".png"), width=output_width, height = output_height)
        damageLogo.pos.str.skeleton(pwm = prop_patterns_list[[l]],
                                    probs = prob_mutation[l,],
                                    breaks_theta_vec = breaks_theta[,l, drop=FALSE],
                                    strand_theta_vec = strand_theta[1,],
                                    ic = ic[,l],
                                    max_pos = max_pos,
                                    max_prob = max_prob,
                                    ic.scale = ic.scale,
                                    yscale_change = yscale_change,
                                    xlab = xlab,
                                    xaxis=xaxis,
                                    yaxis=yaxis,
                                    xaxis_fontsize=xaxis_fontsize,
                                    xlab_fontsize=xlab_fontsize,
                                    title_aligner = title_aligner,
                                    y_fontsize=y_fontsize,
                                    title_fontsize = title_fontsize,
                                    mut_width=mut_width,
                                    start=start,
                                    pop_name = pop_names[l],
                                    logoport_x = logoport_x,
                                    logoport_y= logoport_y,
                                    logoport_width= logoport_width,
                                    logoport_height= logoport_height,
                                    lineport_x = lineport_x,
                                    lineport_y= lineport_y,
                                    lineport_width=lineport_width,
                                    lineport_height=lineport_height,
                                    breaklogoport_x = breaklogoport_x,
                                    breaklogoport_y = breaklogoport_y,
                                    breaklogoport_width=breaklogoport_width,
                                    breaklogoport_height=breaklogoport_height,
                                    barport_x = barport_x,
                                    barport_y = barport_y,
                                    barport_width = barport_width,
                                    barport_height = barport_height)
        dev.off()
      # }else{
      #   damageLogo.pos.str.skeleton(pwm = prop_patterns_list[[l]],
      #                               probs = prob_mutation[l,],
      #                               breaks_theta_vec = breaks_theta[,l, drop=FALSE],
      #                               strand_theta_vec = strand_theta[l,],
      #                               ic = ic[,l],
      #                               max_pos = max_pos,
      #                               max_prob = max_prob,
      #                               ic.scale = ic.scale,
      #                               yscale_change = yscale_change,
      #                               xlab = xlab,
      #                               xaxis=xaxis,
      #                               yaxis=yaxis,
      #                               xaxis_fontsize=xaxis_fontsize,
      #                               xlab_fontsize=xlab_fontsize,
      #                               title_aligner = title_aligner,
      #                               y_fontsize=y_fontsize,
      #                               title_fontsize = title_fontsize,
      #                               mut_width=mut_width,
      #                               start=start,
      #                               pop_name = pop_names[l],
      #                               logoport_x = logoport_x,
      #                               logoport_y= logoport_y,
      #                               logoport_width= logoport_width,
      #                               logoport_height= logoport_height,
      #                               lineport_x = lineport_x,
      #                               lineport_y= lineport_y,
      #                               lineport_width=lineport_width,
      #                               lineport_height=lineport_height,
      #                               breaklogoport_x = breaklogoport_x,
      #                               breaklogoport_y = breaklogoport_y,
      #                               breaklogoport_width=breaklogoport_width,
      #                               breaklogoport_height=breaklogoport_height,
      #                               barport_x = barport_x,
      #                               barport_y = barport_y,
      #                               barport_width = barport_width,
      #                               barport_height = barport_height)
      # }
    }
  }
}

damageLogo.pos.str.skeleton <- function(pwm,
                                        probs,
                                        breaks_theta_vec,
                                        strand_theta_vec,
                                        ic,
                                        max_pos,
                                        max_prob,
                                        ic.scale=TRUE,
                                        xlab = "",
                                        xaxis=TRUE,
                                        yaxis=TRUE,
                                        xaxis_fontsize=10,
                                        xlab_fontsize=15,
                                        title_aligner = 8,
                                        y_fontsize=15,
                                        title_fontsize = 20,
                                        mut_width=2,
                                        start=0.0001,
                                        yscale_change=TRUE,
                                        pop_name = NULL,
                                        logoport_x=0.25,
                                        logoport_y= 0.5,
                                        logoport_width= 0.3,
                                        logoport_height= 0.9,
                                        lineport_x = 0.6,
                                        lineport_y=0.25,
                                        lineport_width=0.30,
                                        lineport_height=0.30,
                                        breaklogoport_x = 1,
                                        breaklogoport_y = 0.50,
                                        breaklogoport_width=0.40,
                                        breaklogoport_height=0.35,
                                        barport_x = 0.5,
                                        barport_y=0.80,
                                        barport_width=0.25,
                                        barport_height=0.25){

  if (class(pwm) == "data.frame"){
    pwm <- as.matrix(pwm)
  }else if (class(pwm) != "matrix"){
    stop("pwm must be of class matrix or data.frame")
  }

  if (any(abs(1 - apply(pwm,2,sum)) > 0.01))
    stop("Columns of PWM must add up to 1.0")

  chars <- c("A", "C", "G", "T", "X",
             "C->T", "C->A", "C->G",
             "T->A", "T->C", "T->G")

  letters <- list(x=NULL,y=NULL,id=NULL,fill=NULL)
  npos <- ncol(pwm)

  if (ic.scale){
    if(yscale_change){
      if(max(ic)<1){ylim <- 1
      #facs <- ic + 1 - max(ic)
      facs <- ic/max(ic)
      }
      if(max(ic)>1){ylim <- 2
      facs <- ic}
    }else{
      ylim <- ceiling(max(ic))
      facs <- ic
    }
    ylab <- "Information content"
  }else{
    ylim <- 1
    ylab <- "Probability"
    facs <- rep(1, npos)
  }

  wt <- c(rep(1, floor(npos/2)),mut_width,rep(1, floor(npos/2)))
  flanked_coord <- c((-floor(npos/2)):(-1), 0, 1:floor(npos/2))
  x.pos <- 0
  for (j in 1:npos){

    column <- pwm[,j]
    hts <- as.numeric(0.99*column*facs[j])
    letterOrder <- order(hts)

    y.pos <- 0
    for (i in 1:length(chars)){
      letter <- chars[letterOrder[i]]
      ht <- hts[letterOrder[i]]
      if (ht>0) letters <- addLetter2(letters,letter,x.pos,y.pos,ht,wt[j])
      y.pos <- y.pos + ht + start
    }
    x.pos <- x.pos + wt[j]
  }

  xlim <- cumsum(wt) - wt/2;
  # xlim <- c(wt[1]/2, wt[1] + wt[2]/2, wt[1]+wt[2]+wt[3]/2, wt[1]+wt[2]+wt[3], 5.5)
  low_xlim <- c(xlim - 0.5*wt, xlim[length(xlim)]+0.5*wt[length(xlim)])
  ylim_scale <- seq(0, ylim, length.out=6);
  if(ic.scale){
    ic_lim_scale <- seq(0, max(ic), length.out=6)
  }else{
    ic_lim_scale <- ylim_scale
  }
  if(ic.scale){
    if(yscale_change){
      if(ylim  > 1){
        letters$y <- letters$y*(ylim/max(ic));
      }
    }}

  plot.new()
  grid::grid.newpage()
  vp <- viewport(x=logoport_x, y=logoport_y, width=logoport_width, height=logoport_height)
  pushViewport(vp)
  #pushViewport(plotViewport(c(bottomMargin,leftMargin,max(ylim)+1.5,max(ylim)+ 1.5)))
  # pushViewport(viewport(layout = grid.layout(2, 2),
  #              x = bottomMargin,
  #              y = leftMargin,
  #              width = max(xlim/2)+0.5,
  #              height = max(ylim/2)+0.5))
  pushViewport(dataViewport(0:ncol(pwm),0:ylim,name="vp1"))
  grid.polygon(x=unit(letters$x,"native"), y=unit(letters$y,"native"),
               id=letters$id, gp=gpar(fill=letters$fill,col="transparent"))
  grid.polygon(x=unit(letters$x,"native"), y=unit(letters$y,"native"),
               id=letters$id,
               gp=gpar(fill=letters$fill,col="transparent"))

  xlim <- cumsum(wt) - wt/2;
  # xlim <- c(wt[1]/2, wt[1] + wt[2]/2, wt[1]+wt[2]+wt[3]/2, wt[1]+wt[2]+wt[3], 5.5)
  low_xlim <- c(xlim - 0.5*wt, xlim[length(xlim)]+0.5*wt[length(xlim)])
  ylim_scale <- seq(0, ylim, length.out=6);
  ic_lim_scale <- seq(0, max(ic), length.out=6)

  for(n in 1:length(xlim)){
    grid.lines(x = unit(low_xlim[n], "native"),
               y = unit(c(0, ylim), "native"),
               gp=gpar(col="grey80"))
  }

  if(is.null(pop_name)){
    grid.text("Logo plot", x = unit(1.5, "npc"), y = unit(1, "npc") + unit(title_aligner, "lines"),
              gp = gpar(fontsize = title_fontsize))
  }else{
    grid.text(paste0(pop_name), x = unit(1.5, "npc"), y = unit(title_aligner, "lines"),
              gp = gpar(fontsize = title_fontsize, col="black"))
  }

  if (xaxis){
    grid.xaxis(at=xlim,
               label=c(paste0("\n left \n flank \n"),
                        "\n mismatch",
                        paste0("\n right \n flank \n")
               ),
               gp=gpar(fontsize=xaxis_fontsize))
    grid.text(paste0(xlab),y=unit(-3,"lines"),
              gp=gpar(fontsize=xlab_fontsize))
  }
  if (yaxis){
    if(yscale_change==TRUE){
      grid.yaxis(at = ylim_scale,
                 label = round(ic_lim_scale,2),
                 gp=gpar(fontsize=y_fontsize))
    }else{
      grid.yaxis(gp=gpar(fontsize=y_fontsize))
    }
    grid.text(ylab,x=unit(-4,"lines"),rot=90,
              gp=gpar(fontsize=y_fontsize))
  }

  upViewport(0)
  vp2 <- viewport(x=lineport_x,y=lineport_y,width=lineport_width, height=lineport_height, just=c("right","top"))
  pushViewport(vp2)
  par(plt = gridPLT(), new=TRUE)
  plot_graph(probs, max_pos=max_pos, max_prob=max_prob)
  upViewport(0)
  #pushViewport(vps$figure)
  # vp3 <- viewport(x=breaklogoport_x, y=breaklogoport_y, width=breaklogoport_width, height=breaklogoport_height, just=c("right","bottom"))
  # pushViewport(vp3)
  # par(plt = gridPLT(), new=TRUE)

  vp4 <- viewport(width = barport_width, height = barport_height, x = barport_x, y = barport_y)
  p = plot_bar(strand_theta_vec_2 = strand_theta_vec)
  print(p, vp = vp4)

  vp3 <- viewport(width = breaklogoport_width, height = breaklogoport_height, x = breaklogoport_x,
                  y = breaklogoport_y, just = c("right","bottom"))
  pushViewport(vp3)
  plot_logo(breaks_theta_vec = breaks_theta_vec)
  upViewport(0)
  #print(p, vp = vp3)
  #plot(1:4, col="red", pch=20)
  popViewport(0)
  par(ask=FALSE)
}

plot_bar <- function(strand_theta_vec_2){
  df <- data.frame("value" = as.numeric(t(strand_theta_vec_2)), "strand"=plyr::revalue(factor(names(strand_theta_vec_2)), c("plus" = "+",  "minus" = "\U2012")))
  ggplot2::ggplot(df, aes(x=strand, y=value, fill = c("lightpink", "green")))  +
    geom_bar(stat='identity', width=0.8, position = position_dodge(width=1.5),
             colour = 'black') + labs(title="   strand \n orientation ") +
    xlab("") + ylim(0,1) + ylab("") + theme(legend.position="none",
                                            axis.ticks = element_blank(),
                                            axis.text = element_blank(),
                                            panel.background = element_rect(fill = "white"),
                                            panel.border = element_rect(colour = "white"),
                                            title = element_text(size = 27)) +
    #  geom_text(aes(x=strand, y=value+0.1, label=value)) +
    geom_text(aes(x=strand, y=value*0.5, label=as.character(strand)), colour="white", size=25) + coord_flip()
    # theme(axis.text = element_blank(),
    #       axis.ticks = element_blank(),
    #       panel.grid  = element_blank(),
    #       panel.border = element_rect(colour = "white")) +
}

plot_logo <- function(breaks_theta_vec){
  if(any(rownames(breaks_theta_vec) != c("A", "C", "G", "T"))){
    stop("the rownames must be A, C, G, T in that sequence")
  }
  mat <- matrix(0, nrow=4, ncol=2)
  mat[c(1,3),1] <- breaks_theta_vec[c(1,3),1]
  mat[c(2,4),2] <- breaks_theta_vec[c(2,4),1]
  rownames(mat) <- c("A", "C", "G", "T")
  colnames(mat) <- c("purine", "pyrimidine")

  cols = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, cols$maxcolors, rownames(cols)))

  color_profile <- list("type" = "per_column",
                        "col" = c("blue", "red"))


  Logolas::logomaker(mat,
            color_profile = color_profile,
            frame_width = 1,
            ic.scale = TRUE,
            yscale_change = TRUE,
            pop_name = "strand-break composition",
            xlab = "",
            ylab = "",
            yaxis=FALSE,
            main_fontsize = 30,
            xaxis_fontsize = 27,
            col_line_split="black",
            control = list(hist=TRUE),
            newpage=FALSE)
}

pwm2ic<-function(pwm) {
  npos<-ncol(pwm)
  ic<-numeric(length=npos)
  for (i in 1:npos) {
    ic[i]<- log(length(which(pwm[,i]!=0.00)), base=2) + sum(sapply(pwm[, i], function(x) {
      if (x > 0) { x*log2(x) } else { 0 }
    }))
  }
  ic
}

ic_computer_2 <-function(mat, alpha, base_probs = NULL) {

  mat <- apply(mat, 2, function(x) return(x/sum(x)))
  npos<-ncol(mat)
  ic <-numeric(length=npos)
  for (i in 1:npos) {

    if(is.null(base_probs)){
      ll <- length(which(mat[,i]!=0.00))
      probs <- rep(1/ll, ll)
    } else{
      probs <- base_probs
    }

    if(alpha == 1){
      ic[i] <- sum(sapply(1:length(mat[,i]), function(x) {
        if (x > 0) { mat[x, i]*log2(mat[x,i]) } else { 0 }
      })) - sum(sapply(1:length(mat[,i]), function(x) {
        if (x > 0) { mat[x, i]*log2(probs[x]) } else { 0 }
      }))
    }
    else if(alpha == Inf){
      ic[i] <- log(length(which(mat[,i]!=0.00)), base=2) + log(max(mat[,i]))
    }
    else if(alpha <= 0){
      stop("alpha value must be greater than 0")
    }
    else{
      ic[i] <-  abs((log(length(which(mat[,i] !=0.00)), base=2) - (1/(1-alpha))* log2(sum(mat[,i]^{alpha}))) -
                      (log(length(which(mat[,i] !=0.00)), base=2) - (1/(1-alpha))* log2(sum(probs^{alpha}))))
    }
  }
  return(ic)
}


damage.ic<-function(pwm, alpha=1, inflation_factor = c(1,1,1), base_probs_list = NULL) {
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

    if(!is.null(base_probs_list)){
      base_probs <- base_probs_list[[i]]
      idx <- match(names(base_probs), rownames(mat))
      mat <- mat[idx[!is.na(idx)],]
      mat_clean <- mat
    }else{
      base_probs = NULL
      mat_clean <- mat[rowSums(mat) != 0,]
    }

    ic[i,] <- inflation_factor[i]*ic_computer_2(mat_clean, alpha, base_probs)
  }

  return(ic)
}


###########   skeleton damage logo    ###################



###################################################################
#####################  Damage Logos  ##############################
###################################################################



###################  letter  A   #################################


letterA <- function(x.pos,y.pos,ht,wt,id=NULL){

  x <- c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6)
  y <- c(0,10,10,0,0,3,3,0,0,4,7.5,4,4)
  x <- 0.1*x
  y <- 0.1*y

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- c(rep(1,9),rep(2,4))
  }else{
    id <- c(rep(id,9),rep(id+1,4))
  }

  fill <- c("green","white")

  list(x=x,y=y,id=id,fill=fill)
}

#grid.newpage()

#grid.polygon(x=a_let$x, y=a_let$y, gp=gpar(fill=a_let$fill,
#             col="transparent"))

################   letter  T   ##################################

letterT <- function(x.pos,y.pos,ht,wt,id=NULL){

  x <- c(0,10,10,6,6,4,4,0)
  y <- c(10,10,9,9,0,0,9,9)
  x <- 0.1*x
  y <- 0.1*y

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- rep(1,8)
  }else{
    id <- rep(id,8)
  }

  fill <- "red"

  list(x=x,y=y,id=id,fill=fill)
}

#################### letter  C #####################################

letterC <- function(x.pos,y.pos,ht,wt,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- rep(1,length(x))
  }else{
    id <- rep(id,length(x))
  }

  fill <- "blue"

  list(x=x,y=y,id=id,fill=fill)
}


##################  letter G  ######################################

letterG <- function(x.pos,y.pos,ht,wt,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))

  h1 <- max(y.l1)
  r1 <- max(x.l1)

  h1 <- 0.4
  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)



  if (is.null(id)){
    id <- c(rep(1,length(x)),rep(2,length(x.add)))
  }else{
    id <- c(rep(id,length(x)),rep(id+1,length(x.add)))
  }

  x <- c(rev(x),x.add)
  y <- c(rev(y),y.add)

  x <- x.pos + wt*x
  y <- y.pos + ht*y


  fill <- c("orange","orange")

  list(x=x,y=y,id=id,fill=fill)

}

letterX <- function(x.pos,y.pos,ht,wt,id=NULL){

  x <- c( 0, 0.4, 0, 0.2, 0.5, 0.8, 1, 0.6, 1, 0.8, 0.5, 0.2)
  x <- 0.05 + 0.90*x
  y <- c( 0, 0.5, 1, 1, 0.6, 1, 1, 0.5, 0, 0, 0.4, 0)

  if (is.null(id)){
    id <- rep(1,length(x))
  }else{
    id <- rep(id,length(x))
  }

  fill <- "pink"

  x <- x.pos + wt*x
  y <- y.pos + ht*y


  ll <- list("x"= x,
             "y"= y,
             "id" = id,
             "fill" = fill)
  return(ll)
}





################  letter C to T  ###############################

letter_C_to_T <- function(x.pos,y.pos,ht,wt,id=NULL){

  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  x1 <- 0.4*x
  y1 <- 1*y


  x <- c( 0.4, 0.4, 0, 0, 1, 1, 0.6, 0.6)
  y <- c( 0, 0.85, 0.85, 1, 1, 0.85, 0.85, 0)
  x2 <- 0.6 + 0.4*x
  y2 <- 1*y


  x3 <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3 <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)

  xpool <- c(x1, x2, x3)
  ypool <- c(y1, y2, y3)

  if(!is.null(id)){
    id_pool <- id + c(rep(1,length(x1)), rep(3, length(x2)),
                      rep(4, length(x3)))
  }else{
    id_pool <- c(rep(1,length(x1)), rep(3, length(x2)),
                 rep(4, length(x3)))
  }

  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool

  fill=c("blue","red","grey80")

  list(x=x,y=y,id=id_pool,fill=fill)

}

###############  letter C to G  ###################################

letter_C_to_G <- function(x.pos,y.pos,ht,wt,id=NULL){

  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  x1_pool <- 0.4*x
  y1_pool <- 1*y

  id1_pool <- rep(1,length(x1_pool))



  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))

  h1 <- max(y.l1)
  r1 <- max(x.l1)

  h1 <- 0.4
  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)



  id2_pool <- c(rep(2,length(x)),rep(3,length(x.add)))


  x2_pool <- 0.6 + 0.4*c(rev(x),x.add)
  y2_pool <- c(rev(y),y.add)


  x3_pool <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3_pool <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)

  xpool <- c(x1_pool, x2_pool, x3_pool)
  ypool <- c(y1_pool, y2_pool, y3_pool)

  if(!is.null(id)){
    id_pool <- id + c(id1_pool, id2_pool, rep(4, length(x3_pool)))
  }else{
    id_pool <- c(id1_pool, id2_pool, rep(4, length(x3_pool)))
  }

  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool

  fill <- c("blue","orange", "orange", "grey80")

  list(x=x,y=y,id=id_pool,fill=fill)
}

##############  letter C to A  ####################################

letter_C_to_A <- function(x.pos,y.pos,ht,wt,id=NULL){

  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  x1 <- 0.4*x
  y1 <- 1*y


  x <- 0.1* (c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6))
  y <- 0.1*(c(0,10,10,0,0,3,3,0,0,4,7.5,4,4))
  x2 <- 0.6 + 0.4*x
  y2 <- 1*y

  x3 <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3 <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)

  xpool <- c(x1, x2, x3)
  ypool <- c(y1, y2, y3)

  if(!is.null(id)){
    id_pool <- id +  c(rep(1,length(x1)), c(rep(2,9),rep(3,4)),
                       rep(4, length(x3)))
  }else{
    id_pool <- c(rep(1,length(x1)), c(rep(2,9),rep(3,4)),
                 rep(4, length(x3)))
  }

  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool

  fill=c("blue","green", "white", "grey80")

  list(x=x,y=y,id=id_pool,fill=fill)
}

################  letter T to G  ################################


letter_T_to_G <- function(x.pos,y.pos,ht,wt,id=NULL){

  x <- c( 0.4, 0.4, 0, 0, 1, 1, 0.6, 0.6)
  y <- c( 0, 0.85, 0.85, 1, 1, 0.85, 0.85, 0)
  x1_pool <- 0.4*x
  y1_pool <- 1*y
  id1_pool <- rep(1,length(x1_pool))


  x3_pool <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3_pool <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)


  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))

  h1 <- max(y.l1)
  r1 <- max(x.l1)

  h1 <- 0.4
  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)



  id2_pool <- c(rep(2,length(x)),rep(3,length(x.add)))


  x2_pool <- 0.6 + 0.4*c(rev(x),x.add)
  y2_pool <- c(rev(y),y.add)

  if(!is.null(id)){
    id_pool <- id + c(id1_pool, id2_pool, rep(4, length(x3_pool)))
  }else{
    id_pool <- c(id1_pool, id2_pool, rep(4, length(x3_pool)))
  }

  xpool <- c(x1_pool, x2_pool, x3_pool)
  ypool <- c(y1_pool, y2_pool, y3_pool)

  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool

  fill=c("red","orange","orange", "grey80")

  list(x=x,y=y,id=id_pool,fill=fill)
}

################ letter  T to C  ################################

letter_T_to_C <- function(x.pos,y.pos,ht,wt,id=NULL){

  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  x2 <- 0.6 + 0.4*x
  y2 <- 1*y


  x <- c( 0.4, 0.4, 0, 0, 1, 1, 0.6, 0.6)
  y <- c( 0, 0.85, 0.85, 1, 1, 0.85, 0.85, 0)
  x1 <- 0.4*x
  y1 <- 1*y


  x3 <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3 <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)

  xpool <- c(x1, x2, x3)
  ypool <- c(y1, y2, y3)

  if(!is.null(id)){
    id_pool <- id +c(rep(1,length(x1)), rep(3, length(x2)),
                     rep(4, length(x3)))
  }else{
    id_pool <- c(rep(1,length(x1)), rep(3, length(x2)),
                 rep(4, length(x3)))
  }

  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool

  fill=c("red","blue","grey80")

  list(x=x,y=y,id=id_pool,fill=fill)

}

#################  letter T to A  ####################################

letter_T_to_A <- function(x.pos,y.pos,ht,wt,id=NULL){
  x <- c( 0.4, 0.4, 0, 0, 1, 1, 0.6, 0.6)
  y <- c( 0, 0.85, 0.85, 1, 1, 0.85, 0.85, 0)
  x1 <- 0.4*x
  y1 <- 1*y


  x <- 0.1* (c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6))
  y <- 0.1*(c(0,10,10,0,0,3,3,0,0,4,7.5,4,4))
  x2 <- 0.6 + 0.4*x
  y2 <- 1*y


  x3 <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3 <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)

  xpool <- c(x1, x2, x3)
  ypool <- c(y1, y2, y3)

  if(!is.null(id)){
    id_pool <- id + c( rep(3, length(x1)), c(rep(1,9),rep(2,4)),
                       rep(4, length(x3)))
  }else{
    id_pool <- c( rep(3, length(x1)), c(rep(1,9),rep(2,4)),
                  rep(4, length(x3)))
  }

  fill=c("green","white", "red", "grey80")

  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool

  list(x=x,y=y,id=id_pool,fill=fill)
}

################# add letters to logo plot #########################


addLetter2 <- function(letters,which,x.pos,y.pos,ht,wt){

  if (which == "A"){
    letter <- letterA(x.pos,y.pos,ht,wt)
  }else if (which == "C"){
    letter <- letterC(x.pos,y.pos,ht,wt)
  }else if (which == "G"){
    letter <- letterG(x.pos,y.pos,ht,wt)
  }else if (which == "X"){
    letter <- letterX(x.pos,y.pos,ht,wt)
  }else if (which == "T"){
    letter <- letterT(x.pos,y.pos,ht,wt)
  }else if (which == "C->T"){
    letter <- letter_C_to_T(x.pos,y.pos,ht,wt)
  }else if (which == "C->A"){
    letter <- letter_C_to_A(x.pos,y.pos,ht,wt)
  }else if (which == "C->G"){
    letter <- letter_C_to_G(x.pos,y.pos,ht,wt)
  }else if (which == "T->A"){
    letter <- letter_T_to_A(x.pos,y.pos,ht,wt)
  }else if (which == "T->G"){
    letter <- letter_T_to_G(x.pos,y.pos,ht,wt)
  }else if (which == "T->C"){
    letter <- letter_T_to_C(x.pos,y.pos,ht,wt)
  }else{
    stop("which must be one of A,C,G,T,X, C->T,
         C->G, C->A, T->A, T->C, T->G")
  }

  letters$x <- c(letters$x,letter$x)
  letters$y <- c(letters$y,letter$y)

  lastID <- ifelse(is.null(letters$id),0,max(letters$id))
  letters$id <- c(letters$id,lastID+letter$id)
  letters$fill <- c(letters$fill,letter$fill)
  letters
  }


plot_graph <- function(probs, max_pos, max_prob, col="red",
                       cex=unit(1.3, "npc"), pch=unit(16,"npc"),
                       xlab="position", ylab="prob. of mismatch",
                       main="",
                       cex.axis=unit(1.5, "npc"),
                       cex.main=unit(1.5, "npc")){
  # if (length(probs) != max_pos){
  #   stop(cat('probability vector must be of length ', max_pos))
  # }
  par(font.axis = 2)
  plot(as.numeric(names(probs)), probs/max_prob, xlim = c(0, max_pos), ylim=c(0,1),
       type = "b", xaxt = "n", yaxt = "n", cex = cex, pch=pch, col=col, main=main,
       cex.main=cex.main, ylab="", xlab="")
  axis(side = 1, at = floor(seq(1, max_pos, length.out=5)), cex.axis = cex.axis, lwd.ticks = 1, tck=-0.05,
       cex.lab=2.8, mgp=c(3.5, 1, 0))
  title(xlab = xlab, mgp=c(3,1.5,0), cex.lab=2.3)
  ylimit <- c(0.0, 0.5, 1.0)*max_prob
  axis(side = 2, at = c(0.0, 0.5, 1.0), labels = round(ylimit,2),
       cex.axis = cex.axis, lwd.ticks=1, tck=-0.05,
       cex.lab=2, mgp=c(3.5, 1, 0))
  title(ylab = ylab, mgp=c(3,1,0), cex.lab=2.3)
}


