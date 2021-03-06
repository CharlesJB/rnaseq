% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anno.R
\name{prepare_anno}
\alias{prepare_anno}
\title{Download and clean ref and prepare anno}
\usage{
prepare_anno(
  prefix,
  org,
  db,
  release,
  removeTxVersion = TRUE,
  ERCC92 = FALSE,
  force_download = FALSE
)
}
\arguments{
\item{prefix}{The prefix to be used for the files that will be produced.}

\item{org}{The organism name. Currently accepted:
* Homo sapiens (Ensembl and Gencode)
* Mus musculus (Ensembl and Gencode)
* Macaca mulata (Ensembl only)
* Rattus norvegicus (Ensembl only)}

\item{db}{The database to use: Ensembl or Gencode}

\item{release}{The version of the database to use. Must be greater than 100
for Ensembl, 35 for Gencode Homo sapiens and 25 for Gencode Mus musculus.}

\item{removeTxVersion}{Remove tx version? Default: TRUE.}

\item{ERCC92}{Add ERCC92 sequence to reference and to anno? Default: TRUE}

\item{force_download}{Re-download raw reference if it is already present?
Default: FALSE}
}
\value{
Invisibly returns a  \code{list} including the reference
transcriptome as a the\code{DNAStringSet} object, the annotation, and the
infos (metadata).
}
\description{
The goal of this function is to download the reference fasta file for a
specific release of Ensembl or Gencode. The reference is then cleaned. We
keep only the transcript id and we remove the transcript version by default.
It is also possible to add ERCC92 sequences.
}
\details{
After calling this function, a <prefix>.raw_ref.fa.gz file will be
downloaded (if not already present) in the current working directory that
corresponds to the raw reference file. There will also be a clean version of
the transcriptome in the <prefix>.fa.gz format that will be different
based on the parameter used to call the function. There will be a
<prefix>.info that will contains metadata about the file download and the
parameters used. Finally, there will be a <prefix>.csv that contains the
annotation formated correctly for the rnaseq packages.

The <prefix>.info file contains the following columns:
   * prefix: The prefix of the file. Must match filename (i.e.: prefix of
             Hs.Gencode38.csv is Hs.Gencode38).
   * org: The organism name (i.e.: Homo sapiens)
   * db: Database where the annotation was downloaded.
   * release: The version of the database.
   * rnaseq_pkg_version: The rnaseq package version, if the rnaseq package
   was used to download the annotation.
   * download_date: The date the annotation was downloaded.
   * download_url: The URL that was used to download the annotation.
   * md5_raw_ref: md5sum of the raw transcriptome file.
   * md5_clean_ref: md5sum of the cleaned transcriptome.
   * md5_anno: md5sum of the annotation file.
}
\examples{
\dontrun{
  prepare_anno("Hs.Ensembl103", org = "Homo sapiens", db = "Ensembl",
               release = 103)
}

}
