#' Get demo kallisto abundance files
#'
#' @return A vector of kallisto abundance filenames
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#'
#' @export
get_demo_abundance_files <- function() {
    filenames <- c(system.file("extdata/quant/a/abundance.tsv", package="rnaseq"),
      system.file("extdata/quant/b/abundance.tsv", package="rnaseq"),
      system.file("extdata/quant/c/abundance.tsv", package="rnaseq"),
      system.file("extdata/quant/d/abundance.tsv", package="rnaseq"))
    names(filenames) <- letters[1:4]
    filenames
}

#' Get demo txi file
#'
#' @return A txi object
#'
#' @examples
#' txi <- get_demo_txi()
#'
#' @export
get_demo_txi <- function() {
    abundances <- get_demo_abundance_files()
    names(abundances) <- basename(dirname(abundances))
    demo_anno <- system.file("extdata/demo_anno.csv", package = "rnaseq")
    import_kallisto(abundances, anno = demo_anno)
}

#' Get demo design
#'
#' @return A data.frame corresponding to the demo txi object. See
#'         ?get_demo_txi.
#'
#' @examples
#' design <- get_demo_design()
#'
#' @export
get_demo_design <- function() {
    data.frame(sample = letters[1:4], group = c("A", "A", "B", "B"))
}

#' Get demo kallisto quant dir
#'
#' @return The path to kallisto results
#'
#' @examples
#' dir_quant <- get_demo_kallisto_dir()
#'
#' @export
get_demo_kallisto_dir <- function() {
    system.file("extdata/quant", package = "rnaseq")
}

#' Get demo annotation file
#'
#' @return The path to the annotation file
#'
#' @examples
#' file_anno <- get_demo_anno_file()
#'
#' @export
get_demo_anno_file <- function() {
    system.file("extdata/demo_anno.csv", package = "rnaseq")
}
