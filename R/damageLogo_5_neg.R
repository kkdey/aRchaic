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
#' @param use_log TRUE/FALSE depending on whether to use log transform for calculating
#'                enrichment scores.
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



damageLogo_five_neg <- function(theta_pool,
                            sig_names = NULL,
                            ic.scale=TRUE,
                            use_log = FALSE,
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
  theta <- dplyr::tbl_df(data.frame(theta_pool)) %>% dplyr::mutate(sig = signature_patterns) %>% dplyr::group_by(sig) %>% dplyr::summarise_each(funs(sum)) %>% as.data.frame()
  rownames(theta) <-  theta[,1]
  theta <- theta[,-1, drop=FALSE]

  indices_minus <- grep("_-_", signature_set)
  strand_theta <- data.frame("minus" = colSums(theta_pool[indices_minus,]),
                             "plus" = colSums(theta_pool[-indices_minus,]))

  if(flag == 1){
    strand_theta <- data.frame("minus" = colSums(matrix(theta_pool[indices_minus,])),
                               "plus" = colSums(matrix(theta_pool[-indices_minus,])))
    strand_theta <- strand_theta/2;
  }
  breakbase <- substring(signature_set, 8+2*flanking_bases,  8+2*flanking_bases)

  theta_break <- dplyr::tbl_df(data.frame(theta_pool)) %>% dplyr::mutate(sig = breakbase) %>% dplyr::group_by(sig) %>% dplyr::summarise_each(funs(sum)) %>% as.data.frame()
  rownames(theta_break) <- theta_break[,1]
  theta_break <- theta_break[,-1]

  theta_break <- theta_break[match(c("A", "C", "G", "T"), rownames(theta_break)),]
  breaks_theta <- theta_break


  if(is.null(sig_names))
    sig_names <- rownames(theta)

  # prob_mutation <- filter_by_pos(t(theta_pool), max_pos = max_pos)
  prob_mutation <- filter_signatures_only_location(t(theta_pool), max_pos = max_pos, flanking_bases = flanking_bases)
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
      temp2 <- tapply(theta[,l], factor(new_sig_split[,j], levels=c("A", "C", "G", "T", "X",
                                                                    "C->T", "C->A", "C->G",
                                                                    "T->A", "T->C", "T->G")), sum)

      temp2[is.na(temp2)]=0
      prop_patterns_list[[l]] <- cbind(prop_patterns_list[[l]], temp2)
    }
  }

  ic <- damage.ic(prop_patterns_list, alpha=renyi_alpha,
                  inflation_factor = inflation_factor, base_probs_list = base_probs_list)

  grob_list <- list()
  if(flag == 1){
    l = 1
    png(paste0(output_dir, "logo_", pop_names[l], ".png"), width=output_width, height = output_height)
    damageLogo.pos.str.skeleton.neg(pwm = prop_patterns_list[[l]],
                                probs = prob_mutation[l,],
                                breaks_theta_vec = breaks_theta[,l, drop=FALSE],
                                strand_theta_vec = strand_theta[1,],
                                ic = ic[,l],
                                max_pos = max_pos,
                                max_prob = max_prob,
                                ic.scale = ic.scale,
                                use_log = use_log,
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
      damageLogo.pos.str.skeleton.neg(pwm = prop_patterns_list[[l]],
                                  probs = prob_mutation[l,],
                                  breaks_theta_vec = breaks_theta[,l, drop=FALSE],
                                  strand_theta_vec = strand_theta[1,],
                                  ic = ic[,l],
                                  max_pos = max_pos,
                                  max_prob = max_prob,
                                  ic.scale = ic.scale,
                                  use_log = use_log,
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

damageLogo.pos.str.skeleton.neg <- function(pwm,
                                        probs,
                                        breaks_theta_vec,
                                        strand_theta_vec,
                                        ic,
                                        max_pos,
                                        max_prob,
                                        ic.scale=TRUE,
                                        use_log = FALSE,
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
                                        barport_height=0.25,
                                        col_line_split="gray",
                                        main_fontsize = 20){

  pwm[-(1:4), c(1,3)] <- NA
  pwm[1:5, 2] <- NA

  if (class(pwm) == "data.frame"){
    pwm <- as.matrix(pwm)
  }else if (class(pwm) != "matrix"){
    stop("pwm must be of class matrix or data.frame")
  }

  if (any(abs(1 - apply(pwm,2,function(x) return(sum(x, na.rm=TRUE)))) > 0.01))
    stop("Columns of PWM must add up to 1.0")

  chars <- c("A", "C", "G", "T", "X",
             "C->T", "C->A", "C->G",
             "T->A", "T->C", "T->G")



  pwm_adj <- apply(pwm, 2, function(x)
  {
    indices <- which(is.na(x))
    if(length(indices) == 0){
      y = x
      if(use_log){
        y[y < 0.01] = 0.01
        z <- log(y) - median(log(y))
      }else{
        z <- y - median(y)
      }
      return(z)
    }else{
      y <- x[!is.na(x)]
      if(use_log){
        y[y < 0.01] = 0.01
        z <- log(y) - median(log(y))
      }else{
        z <- y - median(y)
      }
      zext <- array(0, length(x))
      zext[indices] <- 0
      zext[-indices] <- z
      return(zext)
    }
  })

  pwm_mat_pos <- pwm_adj
  pwm_mat_pos[pwm_mat_pos<= 0] = 0
  pwm_mat_pos_norm  <- apply(pwm_mat_pos, 2, function(x) return(x/sum(x)))
  pwm_mat_pos_norm[pwm_mat_pos_norm == "NaN"] = 0


  pwm_mat_neg <- pwm_adj
  pwm_mat_neg[pwm_mat_neg >= 0] = 0
  pwm_mat_neg_norm  <- apply(pwm_mat_neg, 2, function(x) return(x/sum(x)))
  pwm_mat_neg_norm[pwm_mat_neg_norm == "NaN"] = 0

  tab_neg <- apply(pwm_adj, 2, function(x) {
    y = x[x < 0]
    if(length(y) == 0){
      return(0)
    }else{
      return(abs(sum(y)))
    }
  })

  tab_pos <- apply(pwm_adj, 2, function(x) {
    y = x[x > 0]
    y[y < 0.01] = 0.01
    if(length(y) == 0){
      return(0)
    }else{
      return(abs(sum(y)))
    }
  })

  pos_neg_scaling <- apply(rbind(tab_pos, tab_neg), 2, function(x) return(x/sum(x)))
  pos_ic <- pos_neg_scaling[1, ] * ic
  neg_ic <- pos_neg_scaling[2, ] * ic


  letters <- list(x=NULL,y=NULL,id=NULL,fill=NULL)
  npos <- ncol(pwm)
  facs <- pos_ic

  ylim <- ceiling(max(pos_ic))
  x.pos <- 0
  slash_inds <- grep("/", chars)

  # if (ic.scale){
  #   if(yscale_change){
  #     if(max(ic)<1){ylim <- 1
  #     #facs <- ic + 1 - max(ic)
  #     facs <- ic/max(ic)
  #     }
  #     if(max(ic)>1){ylim <- 2
  #     facs <- ic}
  #   }else{
  #     ylim <- ceiling(max(ic))
  #     facs <- ic
  #   }
  #   ylab <- "Information content"
  # }else{
  #   ylim <- 1
  #   ylab <- "Probability"
  #   facs <- rep(1, npos)
  # }

  wt <- c(rep(1, floor(npos/2)),mut_width,rep(1, floor(npos/2)))
  flanked_coord <- c((-floor(npos/2)):(-1), 0, 1:floor(npos/2))
  x.pos <- 0

  for (j in 1:npos){

    column <- pwm_mat_pos_norm[,j]
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
  low_xlim <- c(xlim - 0.5*wt, xlim[length(xlim)]+0.5*wt[length(xlim)])
  letters$y <- letters$y + max(abs(neg_ic))

  y1 <- min(letters$y)
  max1 <- max(letters$y)
  ylim_scale <- seq(0, ylim, length.out=6);

  negbins <- ceiling((y1/max1)*6)
  posbins <- 6 - negbins
  ic_lim_scale <- c(seq(0, y1, length.out = negbins),
                    seq(y1, ylim, length.out = posbins))

  letters$y <- letters$y/ylim

  markers <- round(ic_lim_scale,2) - round(y1,2)

  # ylim_scale <- seq(0, ylim, length.out=6);
  # if(ic.scale){
  #   ic_lim_scale <- seq(0, max(ic), length.out=6)
  # }else{
  #   ic_lim_scale <- ylim_scale
  # }
  # if(ic.scale){
  #   if(yscale_change){
  #     if(ylim  > 1){
  #       letters$y <- letters$y*(ylim/max(ic));
  #     }
  #   }}

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
  pushViewport(dataViewport(0:ncol(pwm),0:1,name="vp1"))
  grid.polygon(x=grid::unit(letters$x,"native"), y=grid::unit(letters$y,"native"),
                     id=letters$id, gp=grid::gpar(fill=letters$fill,col="transparent"))
  grid.polygon(x=grid::unit(letters$x,"native"), y=grid::unit(letters$y,"native"),
                     id=letters$id,
                     gp=grid::gpar(fill=letters$fill,col="transparent"))


  for(n in 2:length(xlim)){
    grid::grid.lines(x = grid::unit(low_xlim[n], "native"),
                     y = grid::unit(c(0, max(markers)/ylim), "native"),
                     gp=grid::gpar(col="gray"))
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
    ylab <- "Enrichment Score"
    # if(yscale_change==TRUE){
    grid::grid.yaxis(at = ic_lim_scale/ylim,
                     label = round(ic_lim_scale,2) - round(y1,2),
                     gp=grid::gpar(fontsize=y_fontsize))
    # }else{
    #   grid::grid.yaxis(gp=grid::gpar(fontsize=y_fontsize))
    # }
    grid::grid.text(ylab,x=grid::unit(-3.5,"lines"),rot=90,
                    gp=grid::gpar(fontsize=y_fontsize))
  }


  x.pos <- 0
  letters <- list(x=NULL,y=NULL,id=NULL,fill=NULL)
  npos <- ncol(pwm)

  letters <- list(x=NULL,y=NULL,id=NULL,fill=NULL)
  facs <- neg_ic

  slash_inds <- grep("/", chars)

  for (j in seq_len(npos)){

    column <- pwm_mat_neg_norm[,j]
    hts <- as.numeric(0.99*column*facs[j])
    letterOrder <- rev(order(hts))

    y.pos <- - neg_ic[j]
    for (i in seq_along(chars)){
      letter <- chars[letterOrder[i]]
      ht <- hts[letterOrder[i]]
      if (ht>0) letters <- addLetter2(letters,letter,x.pos,y.pos,ht,wt[j])
      y.pos <- y.pos + ht + start
    }
    x.pos <- x.pos + wt[j]
  }

  letters$y <- letters$y + max(abs(neg_ic))
  #  letters$y <- 0.8*letters$y*(ylim/max(max(pos_ic), max(neg_ic)))
  letters$y <- letters$y/ylim

  xlim <- cumsum(wt) - wt/2;
  # xlim <- c(wt[1]/2, wt[1] + wt[2]/2, wt[1]+wt[2]+wt[3]/2, wt[1]+wt[2]+wt[3], 5.5)
  low_xlim <- c(xlim - 0.5*wt, xlim[length(xlim)]+0.5*wt[length(xlim)])
  ylim_scale <- seq(0, ylim, length.out=6);
  ic_lim_scale <- seq(-max(neg_ic), 0, length.out=6)


  grid::grid.polygon(x=grid::unit(letters$x,"native"), y=grid::unit(letters$y,"native"),
                     id=letters$id, gp=grid::gpar(fill=letters$fill,col="transparent"))
  grid::grid.polygon(x=grid::unit(letters$x,"native"), y=grid::unit(letters$y,"native"),
                     id=letters$id,
                     gp=grid::gpar(fill=letters$fill,col="transparent"))

  grid::grid.lines(x = grid::unit(c(0, (xlim+0.5*wt)), "native"),
                   y = grid::unit(y1/ylim, "native"),
                   gp=grid::gpar(col="black"))


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
                     hist=TRUE,
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
                     newpage=FALSE)
}

