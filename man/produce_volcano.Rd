% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcano.R
\name{produce_volcano}
\alias{produce_volcano}
\title{produce a volcano plot from DESeq2 results}
\usage{
produce_volcano(de_res, fc_threshold = 3, graph = TRUE)
}
\arguments{
\item{de_res}{the \code{data.frame} object returned by the
\code{DESeq2::results} function.}

\item{fc_threshold}{The threshold of FC to be considered as significant.
Default: \code{TRUE}.}

\item{graph}{produce the graph. \code{true} or \code{false}. default:
\code{true}.}
}
\value{
produce the volcano plot and silently returns the \code{ggplot}
object and the data.frame used.
}
\description{
produce a volcano plot from DESeq2 results
}
\examples{
txi <- get_demo_txi()
design <- get_demo_design()
dds <- deseq2_analysis(txi, design, ~ group)
de <- DESeq2::results(dds, contrast = c("group", "A", "B"))
volcano <- produce_volcano(de, graph = FALSE)

}
