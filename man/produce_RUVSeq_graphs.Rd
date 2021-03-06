% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RUVg.R
\name{produce_RUVSeq_graphs}
\alias{produce_RUVSeq_graphs}
\title{Produce RUVSeq recommended graph to evaluate normalization}
\usage{
produce_RUVSeq_graphs(txi, output = "RUVg_graphs.pdf")
}
\arguments{
\item{txi}{The \code{txi} object returned by the \code{ruvg_normalization}
function.}

\item{output}{The name of the pdf that will be produced.
Default: RUVg_graphs.pdf}
}
\value{
Silently return the \code{txi} object used as input.
}
\description{
Produce RUVSeq recommended graph to evaluate normalization
}
\examples{
txi <- get_demo_txi()
# For this demo, we use a subset of housekeeping genes. In a real example,
# It is recommended to use the default value or a complete set of
# housekeeping genes.
txi_ruv <- ruvg_normalization(txi, housekeeping_genes = c("RPL9", "RPL24"))
\dontrun{
   produce_RUVSeq_graphs(txi_ruv)
}

}
