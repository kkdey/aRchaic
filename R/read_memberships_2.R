#' @title Read level cluster membership assignment from fitted GoM model
#'
#' @description Performs a read level cluster membership assignment using the
#' fitted model from \code{aRchiac_cluster()} for each read with mismatch in a
#' BAM file.
#'
#' @param fit a Grade of Membership model saved under \code{model.rda} produced
#' as output by \code{aRchaic_cluster()} function.
#'
#' @param reads_data The new counts matrix of signature patterns for each read
#' in the BAM file.
#'
#' @param method The method used for fitting cluster memberships to the read.
#' There are four options - \code{independent}, \code{independent-nostrand},
#' \code{lik} and \code{map}.
#' The method \code{independent} assumes the topic features (flanking bases,
#' mutations, strand, strand breaks etc) to be independent and uses the
#' probability weighting of each feature to determine read memberships.
#' The \code{independent-nostrand} is same as the \code{independent} method
#' but does not take account of strand and strand break information.
#' The method \code{lik} assumes a multinomial likelihood over the mismatch
#' patterns to determine read memberships.
#' The method \code{map} uses \code{predict.topics} function from maptpx
#' package for the prediction. The default choice of method is \code{independent}.
#'
#'
#' @param add_sig A Boolean stating whether or not to add the mismatch
#' signature patterns to the table as an additional column. Defaults to TRUE.
#'
#' @param round_off A number stating how many places to round off after decimal
#' for the cluster membership proportions
#'
#' @return Produces a matrix of the cluster assignments per read, with the reads
#' marked by read identifiers along the rows and the columns representing the
#' different clusters. The row sums should sum up to 1.
#'
#' @keywords read_memberships_2
#' @importFrom maptpx predict.topics
#' @export

read_memberships_2 = function (fit, reads_data, method = "independent",
                             add_sig=TRUE, round_off = 5){
  if(class(fit) != "topics"){
    stop("The class of the fit object must be topics - generated from aRchaic_cluster()")
  }
  omega <- fit$omega;
  theta <- fit$theta;

  if(is.null(rownames(theta))){
    stop("rownames of the ")
  }

  rownames <- rownames(theta)

  if(method == "map"){
    table <- maptpx::predict.topics(fit, reads_data)
    if(add_sig){
      sig <- apply(reads_data, 1, function(x)
        {
         tmp <- paste0(rownames[which(x > 0)], collapse = " ; ")
         return(tmp)
      })
      table <- cbind.data.frame(table, sig)
    }
  }else if(method == "lik"){
    table <- t(apply(reads_data, 1, function(x){
      tmp1 <- array(0, dim(theta)[2])
      for(k in 1:length(tmp1)){
        tmp1[k] <-  prod(theta[which(x>0),k])
      }
      tmp2 <- tmp1/sum(tmp1)
      tmp2 <- round(tmp2, round_off)
      tmp2 <- tmp2/sum(tmp2)
      if(add_sig){
       tmp3 <- c(tmp2, paste0(rownames[which(x > 0)], collapse = " ; "))
      }else{
        tmp3 <- tmp2
      }
      return(tmp3)
    }))
    table <- as.data.frame(table)
  }else if (method == "independent"){
    message("This method is only applicable when the mutation signatures are provided as row names to fit$theta")

    topic_break <- theta_breakdown(theta)
    table <- t(apply(reads_data, 1, function(x){

      symbol_vec <- rownames(theta)[which(x == 1)]

      if(length(symbol_vec) == 0){
        sym_prob <- rep(1/dim(theta)[2], dim(theta)[2])
      }else{
        sym_prob <- array(1, length(topic_break))
        for(v in 1:length(symbol_vec)){
          sym_prob <- sym_prob * get_prob_symbol(topic_break, as.character(symbol_vec[v]))
        }
        sym_prob <- sym_prob/sum(sym_prob)
        sym_prob <- round(sym_prob, round_off)
        sym_prob <- sym_prob/sum(sym_prob)
      }

      if(add_sig){
        sym_prob_2 <- c(sym_prob, paste0(rownames[which(x == 1)], collapse = " ; "))
      }else{
        sym_prob_2 <- sym_prob
      }
      return(sym_prob_2)
    }))
    table <- as.data.frame(table)
  }else if (method == "independent-nostrand"){
    message("This method is only applicable when the mutation signatures are provided as row names to fit$theta")

    topic_break <- theta_breakdown(theta)
    table <- t(apply(reads_data, 1, function(x){

      symbol_vec <- rownames(theta)[which(x == 1)]

      if(length(symbol_vec) == 0){
        sym_prob <- rep(1/dim(theta)[2], dim(theta)[2])
      }else{
        sym_prob <- array(1, length(topic_break))
        for(v in 1:length(symbol_vec)){
          sym_prob <- sym_prob * get_prob_symbol_nostrand(topic_break, as.character(symbol_vec[v]))
        }
        sym_prob <- sym_prob/sum(sym_prob)
        sym_prob <- round(sym_prob, round_off)
        sym_prob <- sym_prob/sum(sym_prob)
      }

      if(add_sig){
        sym_prob_2 <- c(sym_prob, paste0(rownames[which(x == 1)], collapse = " ; "))
      }else{
        sym_prob_2 <- sym_prob
      }
      return(sym_prob_2)
    }))
    table <- as.data.frame(table)
  }else{
    stop("method must be either of independent, independent-nostrand, lik or map")
  }

  if(!is.null(rownames(reads_data))){
    rownames(table) <- rownames(reads_data)
  }
  colnames(table) <- c(paste("cluster-", 1:dim(theta)[2]), "signature")
  return(table)
}


get_prob_symbol_2 <- function(broken_topics, symbol){

  left_flank <- substring(as.character(symbol), 1,1)
  mismatch <- substring(as.character(symbol), 2, 5)
  right_flank <- substring(as.character(symbol), 6, 6)
  strand <- substring(as.character(symbol), 8, 8)
  if(strand == "+"){strand <- "plus"}else{strand = "minus"}
  strand_break <- substring(as.character(symbol), 10, 10)
  pos <- as.numeric(unlist(strsplit(as.character(symbol), "[_]")[[1]])[4])

  if(pos > 20){
    pos <- 20
  }

  net_prob <- array(0, length(broken_topics))

  for(m in 1:length(broken_topics)){
    ll <- broken_topics[[m]]
    net_prob[m] <- as.numeric(ll$`mismatch-flank`[paste0(left_flank),1] *
                               ll$`mismatch-flank`[paste0(right_flank),3] *
                               ll$`mismatch-flank`[paste0(mismatch),2] *
                               ll$`strand-break-composition`[paste0(strand_break),1] *
                               as.numeric(ll$`strand-orientation`[paste0(strand)]) *
                               ll$`mismatch-trend-read`[pos])

  }

  return(net_prob)
}


get_prob_symbol_nostrand_2 <- function(broken_topics, symbol){

  left_flank <- substring(as.character(symbol), 1,1)
  mismatch <- substring(as.character(symbol), 2, 5)
  right_flank <- substring(as.character(symbol), 6, 6)
  pos <- as.numeric(unlist(strsplit(as.character(symbol), "[_]")[[1]])[4])

  if(pos >= 20){
    pos <- 20
  }

  net_prob <- array(0, length(broken_topics))

  for(m in 1:length(broken_topics)){
    ll <- broken_topics[[m]]
    net_prob[m] <- as.numeric(ll$`mismatch-flank`[paste0(left_flank),1] *
                                ll$`mismatch-flank`[paste0(right_flank),3] *
                                ll$`mismatch-flank`[paste0(mismatch),2] *
                                ll$`mismatch-trend-read`[pos])

  }

  return(net_prob)
}
