#' Get demo kallisto abundance files
#'
#' @param large If \code{TRUE}, returns 8 filenames. Otherwise returns 4.
#'
#' @return A vector of kallisto abundance filenames
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#'
#' @export
get_demo_abundance_files <- function(large = FALSE) {
    stopifnot(is(large, "logical"))
    path <- system.file("extdata/quant/", package = "rnaseq")
    if (!large) {
        dir_names <- letters[1:4]
    } else {
        dir_names <- letters[1:8]
    }
    filenames <- paste0(path, "/", dir_names, "/abundance.tsv")
    names(filenames) <- dir_names
    filenames
}

#' Get demo txi file
#'
#' @param large If \code{TRUE}, txi matrices will contain 8 samples, otherwise
#' they will contain 4.
#'
#' @return A txi object
#'
#' @examples
#' txi <- get_demo_txi()
#'
#' @export
get_demo_txi <- function(large = FALSE) {
    abundances <- get_demo_abundance_files(large)
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

#' Get demo pca infos file
#'
#' @return The path to the pca infos file
#'
#' @examples
#' pca_infos <- get_demo_pca_infos_file()
#'
#' @export
get_demo_pca_infos_file <- function() {
    system.file("extdata/pca_infos.csv", package = "rnaseq")
}

#' Get demo metadata file
#'
#' @return The path to the metadata file
#'
#' @examples
#' metadata <- get_demo_metadata_file()
#'
#' @export
get_demo_metadata_file <- function() {
    system.file("extdata/metadata.csv", package = "rnaseq")
}
