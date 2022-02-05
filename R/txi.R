filter_txi <- function(txi, samples) {
    filter_matrices <- function(txi, name) {
        stopifnot(all(samples %in% colnames(txi[[name]])))
        txi[[name]] <- txi[[name]][,colnames(txi[[name]]) %in% samples]
        txi
    }
    txi <- filter_matrices(txi, "counts")
    txi <- filter_matrices(txi, "abundance")
    txi <- filter_matrices(txi, "length")
    if (!is.null(txi$fpkm)) {
        txi <- filter_matrices(txi, "fpkm")
    }
    if (!is.null(txi$ruvg_counts)) {
        txi <- filter_matrices(txi, "ruvg_counts")
    }
    if (!is.null(txi$combat_counts)) {
        txi <- filter_matrices(txi, "combat_counts")
    }
    validate_txi(txi)
    txi
}

validate_txi <- function(txi) {
    # Global
    stopifnot(is(txi, "list"))
    stopifnot(all(c("counts", "abundance", "length", "anno") %in% names(txi)))

    # Matrices
    stopifnot(is(txi$counts, "matrix"))
    stopifnot(is(txi$abundance, "matrix"))
    stopifnot(is(txi$length, "matrix"))
    stopifnot(is.numeric(txi$counts))
    stopifnot(is.numeric(txi$abundance))
    stopifnot(is.numeric(txi$length))
    stopifnot(identical(colnames(txi$counts), colnames(txi$abundance)))
    stopifnot(identical(colnames(txi$counts), colnames(txi$length)))
    stopifnot(identical(rownames(txi$counts), rownames(txi$abundance)))
    stopifnot(identical(rownames(txi$counts), rownames(txi$length)))

    if (!is.null(txi$ruvg_counts)) {
        stopifnot(is(txi$ruvg_counts, "matrix"))
        stopifnot(is.numeric(txi$ruvg_counts))
        stopifnot(identical(colnames(txi$ruvg_counts), cownames(txi$counts)))
        stopifnot(identical(rolnames(txi$ruvg_counts), rownames(txi$counts)))
    }

    # Data.frame
    stopifnot(is(txi$anno, "data.frame"))
    expected_col <- c("id", "ensembl_gene", "symbol", "entrez_id", "transcript_type")
    stopifnot(expected_col %in% colnames(txi$anno))
    stopifnot(identical(txi$anno$id, rownames(txi$counts)))
    invisible(TRUE)
}
