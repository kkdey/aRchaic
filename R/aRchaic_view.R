#' @title Logo plot representation of an MFF file.
#'
#' @description Aggregate the signatures into bins for an MFF file and then based on the counts
#' of mutational signatures, creates a logo plot representation of the mutational patterns in the
#' sample represented by the MFF file.
#'
#' @param file a MFF file.
#' @param breaks The breaks used for binning the distance from ends of the reads of the mutations, when
#'               processing the MFF files.
#' @param flanking_bases The numbe rof flanking bases to the mutation considered. Defaults to 1.
#' @param logo.control The control parameters for the logo plot representation of the sample.
#'                     Check the input parameters of \code{damageLogo5()} for the control parameters.
#' @param title The name of the sample. If not provided, defaults to the filename.
#' @param output_dir The directory where the logo plot is saved.
#' @param filename The name of the saved logo image file. Defaults to "logo".
#' @return The function creates a logo plot representation of the mutational features for a single
#'         MFF file.
#' @keywords aRchaic_view
#' @import Logolas
#' @import gridBase
#' @import grid
#' @import ggplot2
#' @export



aRchaic_view = function(file,
                        breaks = c(-1, seq(1,20,1)),
                        flanking_bases =1,
                        logo.control = list(),
                        title = NULL,
                        output_dir = NULL,
                        filename = "logo"){

  header <- head(strsplit(rev((as.vector(strsplit(file, "/" )[[1]])))[1], ".csv")[[1]],1)
  if(is.null(title)){
    title <- header
  }
  if(length(strsplit(file, "/")[[1]]) == 1){
    dir <- getwd()
  }else{
    ss <- strsplit(file, "/")[[1]]
    ll <- length(ss)
    dir <- paste0(paste0(ss[1:(ll-1)],collapse="/"), "/")
  }

  if(!file.exists(paste0(dir, header, ".csv"))){
    stop("The file is not found in the given directory")
  }


  logo.control.default <- list(sig_names = NULL,  max_pos = 20,
                               flanking_bases=1,mutlogo.control = list(),
                               breaklogo.control = list(), base_probs_list = NULL,
                               clip = 0, mut_ranges  = c(0, 0), break_ranges = c(0, 0),
                               logoport_x = 0.7, logoport_y= 0.5, logoport_width= 1.2,
                               logoport_height= 1.2, breaklogoport_x = 0.6,
                               breaklogoport_y = 0.45, breaklogoport_width=0.4,
                               breaklogoport_height=1, lineport_x = 0.8,
                               lineport_y=0.5, lineport_width=1,lineport_height=1.4,
                               output_width = 1200, output_height = 700)

  logo.control <- modifyList(logo.control.default, logo.control)


  if(file.exists(paste0(dir, tail(strsplit(dir, "/")[[1]],1), ".rda"))){
    message("Aggregated MFF (.rda) file present: skipping the feature aggregation step")
    mff_dat <- get(load(paste0(dir, tail(strsplit(dir, "/")[[1]],1), ".rda")))
    index <- grep(paste0(header), rownames(mff_dat))
    clubbed_counts <- mff_dat[index, ]
    clubbed_counts_norm <- clubbed_counts/ sum(clubbed_counts)

    temp <- as.matrix(clubbed_counts_norm)
    rownames(temp) <- names(clubbed_counts_norm)

  }else{
    message("Aggregated MFF (.rda) file not present: performing the feature aggregation step ")
    pattern = header
    out <- aggregate_signature_counts(dir = dir,
                                      pattern = paste0(pattern, ".csv"),
                                      breaks = breaks,
                                      flanking_bases = flanking_bases)
    clubbed_counts <- club_signature_counts(out, flanking_bases = 1)
    clubbed_counts_norm <- clubbed_counts/ sum(clubbed_counts)

    temp <- t(clubbed_counts_norm)
  }

  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/")
  }


  do.call(damageLogo_six, append(list(theta_pool = temp,
                                       output_dir = output_dir,
                                       filename = filename),
                                  logo.control))
  }
