#' @title Mutation trends along the read at the  5' and 3' ends.
#'
#' @description Plots mutation trends along the read at the 5' or at the 3' end
#'
#' @param file a MFF file.
#' @param pattern a vector of mutation types for which the trends are to be plotted.
#' @param plot_type if "left", the trend at the 5' end is calculated, else at the 3' end.
#' @param sample_name The name of the sample corresponding to the file. Defaults to filename
#' @param strand can take either +, - or both, representing whether the function consider positive,
#'               negative or both strand mutations to calculate the trends.
#' @param numBases The number of bases from the end upto which the trend is to be plotted.
#'                 Defaults to 20.
#' @param cols Colors of the mutation types in the plot.
#' @param legend_cex The size of the legend. Defaults to 0.5
#' @param main_cex The size of the title of the figure. Defaults to 0.7.
#' @param lab_cex The size of the labels of the figure. Defaults to 0.5.
#' @param layout_left The layout of the legend if the plot is at 5' end. Defaults to "topright".
#' @param layout_right The layout of the legend if the plot is at 3' end. Defaults to "topleft".
#'
#' @return The function creates a logo plot representation of the mutational features for a single
#'         MFF file.
#'
#' @keywords aRchaic_mutation_trend
#' @import  gridBase
#' @import  grid
#' @import  ggplot2
#' @export



aRchaic_mutation_trend <- function(file,
                                   pattern = c("C->T", "C->A", "C->G",
                                               "T->A", "T->C", "T->G",
                                               "G->A", "G->C", "G->T",
                                               "A->C", "A->G", "A->C"
                                   ),
                                   plot_type=c("left", "right"),
                                   sample_name = NULL,
                                   strand = "both",
                                   numBases = 20,
                                   cols= c("black", "yellow", "green",
                                           "blue", "orange", "purple",
                                           "red", "pink", "brown",
                                           "gray", "darkseagreen", "magenta"),
                                   legend_cex = 0.5,
                                   main_cex = 0.7,
                                   lab_cex = 0.7,
                                   layout_left = "topright",
                                   layout_right = "topleft"){
  data <- read.csv(file, header=FALSE)
  totsum <- sum(data[,7]);

  ######## putting G->A and C->T t the beginning ######################

  if(plot_type == "left"){
    index2 <- which(pattern == "G->A")
    if(length(index2) > 0){
      cols <- c(cols[index2], cols[-index2])
      pattern <- c(pattern[index2], pattern[-index2])
    }

    index1 <- which(pattern == "C->T")
    if(length(index1) > 0){
      cols <- c(cols[index1], cols[-index1])
      pattern <- c(pattern[index1], pattern[-index1])
    }
  }

  if(plot_type == "right"){
    index2 <- which(pattern == "C->T")
    if(length(index2) > 0){
      cols <- c(cols[index2], cols[-index2])
      pattern <- c(pattern[index2], pattern[-index2])
    }

    index1 <- which(pattern == "G->A")
    if(length(index1) > 0){
      cols <- c(cols[index1], cols[-index1])
      pattern <- c(pattern[index1], pattern[-index1])
    }
  }


  ###########   if strand information, focus on only one strand ###########


  if(strand == "+"){
    data <- data[data[,6] == "+",]
  }
  if(strand =="_"){
    data <- data[data[,6] == "-",]
  }

  flag <- as.numeric();

  for(i in 1:length(pattern)){
    if(length(grep(pattern[i], data[,1])) > 0){
      flag <- c(flag, i);
      pattern_data <- data[grep(pattern[i], data[,1]),]
      pattern_data[,7] <- pattern_data[,7]/totsum;

      if(length(which(pattern_data[,3]==-1)) !=0){
        pattern_data[which(pattern_data[,3]==-1), 3] <- 0
      }

      if(min(pattern_data[,3]) == 1){ pattern_data[,3] = pattern_data[,3] - 1}

      pattern_data <- pattern_data[pattern_data[,2] <= numBases | pattern_data[,3] <= numBases, ]

      tab_pattern_left <- tapply(pattern_data[pattern_data[,2] <= numBases,7], pattern_data[pattern_data[,2] <= numBases,2], sum);

      lo_left <- loess(as.numeric(tab_pattern_left)~as.numeric(names(tab_pattern_left)))


      tab_pattern_right <- tapply(pattern_data[pattern_data[,3] <= numBases,7], pattern_data[pattern_data[,3] <= numBases,3], sum);

      lo_right <- loess(as.numeric(tab_pattern_right)~as.numeric(names(tab_pattern_right)))

      if(i==1){

        if(plot_type=="left"){
          if(is.null(sample_name)){
            temp2 <- predict(lo_left)
            values <- as.numeric(names(tab_pattern_left))
            temp <- array(0, numBases)
            temp[match(values, 0:numBases)] <- temp2
            temp[temp < 0] =0
            if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
            plot(temp, col=cols[i], ylim= c(0,max(temp)+0.0005),
                 xlim = c(1, numBases), lwd=2, type="l",
                 ylab=paste0("predicted no. of damages: "),
                 xlab="read positions (from left)",
                 xaxs="i",yaxs="i",
                 main=paste0("Damage plot (5' end) : "),
                 cex.main = main_cex,
                 cex.lab = lab_cex
            )
          }
          if(!is.null(sample_name)){
            temp2 <- predict(lo_left)
            values <- as.numeric(names(tab_pattern_left))
            temp <- array(0, numBases)
            temp[match(values, 0:numBases)] <- temp2
            temp[temp < 0] =0
            if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
            plot(temp, col=cols[i], ylim= c(0,max(temp)+0.0005),
                 xlim = c(1,numBases), lwd=2, type="l",
                 ylab=paste0("predicted no. of damages: "),
                 xlab="read positions (from left)",
                 xaxs="i",yaxs="i",
                 main=paste0("Damage plot (5' end): ", sample_name),
                 cex.main = main_cex,
                 cex.lab = lab_cex
            )}
        }

        if(plot_type=="right"){
          if(is.null(sample_name)){
            temp2 <- predict(lo_right)
            values <- as.numeric(names(tab_pattern_right))
            temp <- array(0, numBases)
            temp[match(values, 0:numBases)] <- temp2
            temp[temp < 0] =0
            if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
            plot(temp, col=cols[i], lwd=2, ylim= c(0,max(temp)+0.0005),
                 type="l",
                 ylab=paste0("predicted no. of damages: "),
                 xlab="read positions (from right)",
                 main=paste0("Damage plot (3' end): "),
                 cex.main = main_cex,
                 cex.lab = lab_cex,
                 xaxs="i",yaxs="i",
                 #xlim=rev(range(as.numeric(names(tab_pattern_right))))
                 xlim = c(numBases, 1))
          }
          if(!is.null(sample_name)){
            temp2 <- predict(lo_right)
            values <- as.numeric(names(tab_pattern_right))
            temp <- array(0, numBases)
            temp[match(values, 0:numBases)] <- temp2
            temp[temp < 0] =0
            if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
            plot(temp, col=cols[i], lwd=2, ylim= c(0,max(temp)+0.0005), type="l",
                 ylab=paste0("predicted no. of damages: "),
                 xlab="read positions (from right)",
                 main=paste0("Damage plot (3' end): ", sample_name),
                 cex.main = main_cex,
                 cex.lab = lab_cex,
                 xaxs="i",yaxs="i",
                 # xlim=rev(range(as.numeric(names(tab_pattern_right))))
                 xlim = c(numBases, 1))
          }
        }

      }else{

        if(plot_type=="left"){
          if(is.null(sample_name)){
            temp2 <- predict(lo_left)
            values <- as.numeric(names(tab_pattern_left))
            temp <- array(0, numBases)
            temp[match(values, 0:numBases)] <- temp2
            temp[temp < 0] =0
            if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
            lines(temp, col=cols[i], lwd=2)
          }
          if(!is.null(sample_name)){
            temp2 <- predict(lo_left)
            values <- as.numeric(names(tab_pattern_left))
            temp <- array(0, numBases)
            temp[match(values, 0:numBases)] <- temp2
            temp[temp < 0] =0
            if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
            lines(temp, col=cols[i], lwd=2)
          }
        }

        if(plot_type=="right"){
          if(is.null(sample_name)){
            temp2 <- predict(lo_right)
            values <- as.numeric(names(tab_pattern_right))
            temp <- array(0, numBases)
            temp[match(values, 0:numBases)] <- temp2
            temp[temp < 0] =0
            if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
            lines(temp, col=cols[i], lwd=2)
          }
          if(!is.null(sample_name)){
            temp2 <- predict(lo_right)
            values <- as.numeric(names(tab_pattern_right))
            temp <- array(0, numBases)
            temp[match(values, 0:numBases)] <- temp2
            temp[temp < 0] =0
            if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
            lines(temp, col=cols[i], lwd=2)
          }
        }

      }
    }
  }
  if(plot_type == "left"){
    legend(layout_left, legend=pattern[flag], fill=cols[flag], cex=legend_cex);
  }else if (plot_type == "right"){
    legend(layout_right, legend=pattern[flag], fill=cols[flag], cex=legend_cex);
  }
}
