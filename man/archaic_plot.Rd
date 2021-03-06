% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/archaic_plot.R
\name{archaic_plot}
\alias{archaic_plot}
\title{STRUCTURE plot and logo plot representation of the clusterring
model from archaic_fit}
\usage{
archaic_plot(
  model,
  topic_cols = c("red", "blue", "darkgoldenrod1", "cyan", "firebrick", "green",
    "hotpink", "burlywood", "yellow", "darkgray", "deepskyblue", "darkkhaki", "brown4",
    "darkorchid", "magenta", "yellow", "azure1", "azure4"),
  background = "modern",
  structure.control = list(),
  logo.control = list(),
  output_dir = NULL
)
}
\arguments{
\item{model}{Fitted model from \code{archaic_fit}.}

\item{topic_cols}{A vector of color assignment to the clusters/topics used
for cluster representation in the STRUCTURE (Rosenberg2002)
representation in \code{archaic_plot()}.}

\item{background}{if equals "modern", as in the default, compares enrichment
of mismatch features against a modern background - else
uses a background with equal probability of all mismatch
features.}

\item{structure.control}{The control or tuning parameters for the STRUCTURE
plot representation \code{structure.pdf} output
of \code{archaic_plot} (Dey2017)}

\item{logo.control}{The control or tuning parameters for the EDLogo plots
representation \code{logo_clus_*.pdf} output
of \code{archaic_plot} (Dey2018)}

\item{output_dir}{The path/directory where to save the output plots}
}
\value{
Returns a \code{structure.pdf} and as many logo plots of the form
        \code{logo_clus_k.pdf} for each cluster k, in the output path
        provided \code{output_dir}.
}
\description{
Takes the clustering model fit  from \code{archaic_fit}
as an input and plots the clusters using as EDLogo plots (Dey et al 2018)
and the proportional mixing of the clusters for each sample using a
STRUCTURE plot (Pritchard et al 2000, Dey2017) representation.
}
\references{
Rosenberg2002.
    Rosenberg, N.A., Pritchard, J.K., Weber, J.L., Cann, H.M., Kidd, K.K.,
    Zhivotovsky, L.A. and Feldman, M.W., 2002. Genetic structure of
    human populations. science, 298(5602), pp.2381-2385.

    Dey2017.
    Dey, K.K., Hsiao, C.J. and Stephens, M., 2017. Visualizing the structure
    of RNA-seq expression data using grade of membership models.
    PLoS genetics, 13(3), p.e1006599.

    Dey2018.
    Dey, K.K., Xie, D. and Stephens, M., 2017. A new sequence logo plot to
    highlight enrichment and depletion. bioRxiv, p.226597.

    Pritchard2002.
    Pritchard, J.K., Stephens, M. and Donnelly, P., 2000.
    Inference of population structure using multilocus genotype data.
    Genetics, 155(2), pp.945-959.
}
\keyword{EDLogo}
\keyword{structure}
