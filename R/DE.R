#' DESeq2 analysis
#'
#' @param txi The txi object returned by the import_kallisto function.
#' @param design The experimental design (see ?DESeqDataSetFromTximport).
#' @param formula The design formula in data.frame format (see
#'                ?DESeqDataSetFromTximport).
#' @param filter The minimum number of reads detected for a feature across all
#'               samples. Default: 2
#' @param count_matrix The count matrix to use for the differential analysis.
#' Will use the \code{DESeq2::DESeqDataSetFromMatrix} instead of the
#' \code{DESeq2::DESeqDataSetFromTximport} function, so will work even if txi
#' object is incomplete (i.e.: length matrix is missing). Default: NA
#' @param ... Extra param for the DESeq2::DESeq function
#'
#' @return A DESeqDataSet object.
#'
#' @examples
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de <- DESeq2::results(dds, contrast = c("group", "A", "B"))
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeqDataSetFromTximport
#' @importFrom DESeq2 DESeq
#'
#' @export
deseq2_analysis <- function(txi, design, formula, filter = 2,
                            count_matrix = NULL, ...) {
    validate_txi(txi)
    stopifnot(all(c("sample") %in% colnames(design)))
    stopifnot(identical(colnames(txi$counts), as.character(design$sample)))
    if (!is.null(count_matrix)) {
        stopifnot(is(count_matrix, "character"))
        stopifnot(count_matrix %in% names(txi))
        stopifnot(is(txi[[count_matrix]], "matrix"))
    }
    if (!is.null(txi$dummy)) {
        stopifnot(is(txi$dummy, "character"))
        if (is.null(count_matrix)) {
            stopifnot("count" %in% txi$dummy)
            stopifnot("length" %in% txi$dummy)
        } else {
            stopifnot(count_matrix %in% txi$dummy)
        }
    }

    if (!is.null(count_matrix)) {
        dds <- DESeq2::DESeqDataSetFromMatrix(txi[[count_matrix]], design,
                                              formula)
    } else {
        dds <- DESeq2::DESeqDataSetFromTximport(txi, design, formula)
    }
    dds <- dds[rowSums(counts(dds)) >= filter]
    dds <- DESeq2::DESeq(dds, ...)
    dds
}

#' Prepare formated DE table.
#'
#' The table contains annotation and the DE results.
#'
#' @param dds The DESeqDataSet object returned by deseq2_analysis.
#' @param txi The txi object returned by the import_kallisto function.
#' @param contrast The contrast for the comparison (see ?DESeq2::results).
#' @param ignoreTxVersion Should the transcript version be ignored for anno
#' mapping. Default: \code{FALSE}.
#' @param digits Integer indicating the number of decimal places
#'
#' @return A data.frame with the anno and the merged counts values.
#'
#' @examples
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de <- format_de(dds, txi, c("group", "A", "B"))
#'
#' @importFrom magrittr %>%
#' @importFrom DESeq2 results
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom stringr str_replace
#'
#' @export
format_de <- function(dds, txi, contrast, ignoreTxVersion = FALSE, digits = 4) {
    if (!is.null(txi$dummy)) {
        stopifnot(is(txi$dummy, "character"))
        stopifnot("abundance" %in% txi$dummy)
    }
    res <- DESeq2::results(dds, contrast = contrast) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("id")
    if (!ignoreTxVersion) {
        res <- dplyr::left_join(res, txi$anno, by = "id")
    } else {
        res$id <- stringr::str_replace(res$id, "\\..*$", "")
        res <- dplyr::left_join(res,
                                dplyr::mutate(txi$anno,
                                              id = str_replace(id, "\\..*$", "")),
                                by = "id")
    }
    res <- dplyr::mutate(res,
               mean_TPM_grp1 = get_mean_tpm(dds, txi, contrast[2]),
               mean_TPM_grp2 = get_mean_tpm(dds, txi, contrast[3]),
               ratio = 2^log2FoldChange,
               fold_change = if_else(ratio < 1, -1 * (1/ratio), ratio)) %>%
               splicing_analysis(txi)

    res <- dplyr::select(res, id, ensembl_gene, symbol, entrez_id, transcript_type,
                  mean_TPM_grp1, mean_TPM_grp2, pV = pvalue, qV = padj,
                  percent_grp1, percent_grp2, main_isoform_grp1,
                  main_isoform_grp2, baseMean, lfcSE, log2FoldChange,
                  fold_change, ratio, stat)

    res <- dplyr::mutate(res,
           mean_TPM_grp1 = round_values(as.numeric(mean_TPM_grp1), digits),
           mean_TPM_grp2 = round_values(as.numeric(mean_TPM_grp2), digits),
           pV = round_values(as.numeric(pV), digits),
           qV = round_values(as.numeric(qV), digits),
           percent_grp1 = round_values(as.numeric(percent_grp1), digits),
           percent_grp2 = round_values(as.numeric(percent_grp2), digits),
           baseMean = round_values(as.numeric(baseMean), digits),
           lfcSE = round_values(as.numeric(lfcSE), digits),
           fold_change = round_values(as.numeric(fold_change), digits),
           log2FoldChange = round_values(as.numeric(log2FoldChange), digits),
           stat = round_values(as.numeric(stat), digits))
    as.data.frame(res)
}

get_mean_tpm <- function(dds, txi, group) {
    samples <- dds@colData[,"sample", drop = TRUE]
    samples <- samples[dds@colData[,"group", drop = TRUE] == group]
    mean_tpm <- rowMeans(txi$abundance[,samples])
    mean_tpm[names(dds)]
}

splicing_analysis <- function(res, txi) {
    is.max <- function(x) seq_along(x) == which.max(x)
    if (txi$txOut) {
        sums <- dplyr::group_by(res, ensembl_gene) %>%
            dplyr::summarize(sum_g1 = sum(mean_TPM_grp1),
            sum_g2 = sum(mean_TPM_grp2))
        dplyr::left_join(res, sums, by = "ensembl_gene") %>%
            dplyr::group_by(ensembl_gene) %>%
            dplyr::mutate(percent_grp1 = mean_TPM_grp1 / sum_g1,
                   percent_grp2 = mean_TPM_grp2 / sum_g2) %>%
            dplyr::mutate(percent_grp1 = if_else(is.nan(percent_grp1),
                                          0, percent_grp1),
                   percent_grp2 = if_else(is.nan(percent_grp2),
                                          0, percent_grp2)) %>%
            dplyr::mutate(main_isoform_grp1 = is.max(mean_TPM_grp1),
                   main_isoform_grp2 = is.max(mean_TPM_grp2)) %>%
            dplyr::select(-sum_g1, -sum_g2)
    } else {
        dplyr::mutate(res,
               percent_grp1 = NA,
               percent_grp2 = NA,
               main_isoform_grp1 = NA,
               main_isoform_grp2 = NA)
    }
}

#' Split differential expression results
#'
#' This function will split the DE results table into 4 elements:
#'     1) de: The original DE table
#'     2) signif: All the genes considered statistically differentially
#'                expressed
#'     3) up: The significant genes that are up-regulated
#'     4) down: The significant genes that are down-regulated
#'
#' To be considered significant, a gene must have a padj (qV) value lower or
#' equal to the \code{p_threshold} param, an absolute fold change greater or
#' equal to the \code{fc_threshold} param. Also, if \code{tpm_threshold} is not
#' \code{NULL}, the average TPM across all samples in the comparison must be
#' greater or equal to the \code{tpm_threshold} value.
#'
#' It is important to use the results produced with the \code{format_de}
#' function.
#'
#' @param de_res the \code{data.frame} object returned by the
#' \code{format_de} function.
#' @param fc_threshold The threshold of FC to be considered as significant.
#' Default: 0.05
#' @param p_threshold The threshold of p stat to be considered as significant.
#' Default: 1.5
#' @param tpm_threshold The threshold of the mean TPM of the current samples to
#' be considered as significant
#'
#' @return A list of DE tables.
#'
#' @importFrom dplyr filter
#'
#' @examples
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de <- format_de(dds, txi, c("group", "A", "B"))
#' split_de <- split_de_results(de)
#'
#' @export
split_de_results <- function(de_res, p_threshold = 0.05, fc_threshold = 1.5,
                             tpm_threshold = NULL) {
    stopifnot(is(de_res, "data.frame"))
    expected_cols <- c("qV", "fold_change", "mean_TPM_grp1", "mean_TPM_grp2")
    stopifnot(all(expected_cols %in% colnames(de_res)))
    stopifnot(is(p_threshold, "numeric"))
    stopifnot(p_threshold >= 0 & p_threshold <= 1)
    stopifnot(is(fc_threshold, "numeric"))
    stopifnot(fc_threshold >= 0)
    if (!is.null(tpm_threshold)) {
        stopifnot(is(tpm_threshold, "numeric"))
        stopifnot(tpm_threshold >= 0)
    }

    idx_p <- de_res$qV <= p_threshold
    idx_fc <- abs(de_res$fold_change) >= fc_threshold
    de_res$p_threshold <- FALSE
    de_res$p_threshold[idx_p] <- TRUE
    de_res$fc_threshold <- FALSE
    de_res$fc_threshold[idx_fc] <- TRUE

    if (!is.null(tpm_threshold)) {
        mean_tpm_1 <- de_res$mean_TPM_grp1
        mean_tpm_2 <- de_res$mean_TPM_grp2
        idx_tpm <- mean(mean_tpm_1, mean_tpm_2) >= tpm_threshold
        de_res$tpm_threshold <- FALSE
        de_res$tpm_threshold[idx_tpm] <- TRUE
        signif <- dplyr::filter(de_res, p_threshold, fc_threshold,
                                tpm_threshold)
    } else {
        signif <- dplyr::filter(de_res, p_threshold, fc_threshold)
    }

    up <- dplyr::filter(signif, fold_change > 0)
    down <- dplyr::filter(signif, fold_change < 0)

    list(de_res = de_res, signif = signif, up = up, down = down)
}

round_values <- function(values, digits = 4) {
    stopifnot(is.numeric(values))
    i <- is.na(values)
    new_values <- round(values, digits) %>% format(scientific = FALSE)
    new_values[i] <- NA
    as.numeric(new_values)
}
