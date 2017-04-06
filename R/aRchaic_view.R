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
#' @return The function creates a logo plot representation of the mutational features for a single
#'         MFF file.
#' @keywords aRchaic_view
#' @import Logolas, gridBase, grid, ggplot2
#' @export



aRchaic_view = function(file,
                        breaks = c(-1, seq(1,20,1)),
                        flanking_bases =1,
                        logo.control = list(),
                        title = NULL,
                        output_dir = NULL){

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

  logo.control.default <- list(sig_names = NULL, ic.scale=TRUE,
                               max_pos = 20, flanking_bases=1,
                               yscale_change = TRUE, xaxis=TRUE,
                               yaxis=TRUE, xlab = " ", xaxis_fontsize=20,
                               xlab_fontsize=10, title_aligner = 11,
                               y_fontsize=20, title_fontsize = 35,
                               mut_width=2, start=0.0001,
                               renyi_alpha = 5, inflation_factor = c(3,1,3),
                               pop_names=title,
                               logoport_x = 0.25, logoport_y= 0.50, logoport_width= 0.25, logoport_height= 0.50,
                               lineport_x = 0.9, lineport_y=0.40, lineport_width=0.32, lineport_height=0.28,
                               breaklogoport_x = 0.94, breaklogoport_y = 0.40, breaklogoport_width=0.30, breaklogoport_height=0.45,
                               barport_x = 0.60, barport_y=0.60, barport_width=0.25, barport_height=0.25,
                               output_width = 1200, output_height = 700)

  logo.control <- modifyList(logo.control, logo.control.default)


  if(file.exists(paste0(dir, tail(strsplit(dir, "/")[[1]],1), ".rda"))){
    message("MutationFeatureFormat file present: skipping the signature aggregation step")
    mff_dat <- get(load(paste0(dir, tail(strsplit(dir, "/")[[1]],1), ".rda")))
    index <- grep(paste0(header), rownames(mff_dat))
    clubbed_counts <- mff_dat[index, ]
    clubbed_counts_norm <- clubbed_counts/ sum(clubbed_counts)

    temp <- as.matrix(clubbed_counts_norm)
    rownames(temp) <- names(clubbed_counts_norm)

  }else{
    message("MutationFeatureFormat file not present: performing the signature aggregation step ")
    pattern = header
    out <- aggregate_signature_counts(dir = dir,
                                      pattern = pattern,
                                      breaks = breaks,
                                      flanking_bases = flanking_bases)
    clubbed_counts <- club_signature_counts(out, flanking_bases = 1)
    clubbed_counts_norm <- clubbed_counts/ sum(clubbed_counts)

    temp <- t(clubbed_counts_norm)
  }

  if(is.null(output_dir)){
    output_dir <- paste0(getwd(), "/")
  }


  do.call(damageLogo_five, append(list(theta_pool = temp,
                                       output_dir = output_dir),
                                  logo.control))
  }
