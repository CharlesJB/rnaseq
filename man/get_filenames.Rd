% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import.R
\name{get_filenames}
\alias{get_filenames}
\title{Parse kallisto quant directory to retrieve filenames}
\usage{
get_filenames(dir_kallisto, file_extension = "h5")
}
\arguments{
\item{dir_kallisto}{Kallisto quantification base directory}

\item{file_extension}{The extension to parse (h5 or tsv). Default: h5}
}
\value{
A named list of the abundance files.
}
\description{
Parse kallisto quant directory to retrieve filenames
}
\examples{
dir_kallisto <- get_demo_kallisto_dir()
filenames <- get_filenames(dir_kallisto, file_extension = "tsv")

}
