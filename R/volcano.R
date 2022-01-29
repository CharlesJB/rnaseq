#' produce a volcano plot from DESeq2 results
#'
#' @param de_res the \code{data.frame} object returned by the
#' \code{DESeq2::results} function.
#' @param fc_threshold The threshold of FC to be considered as significant.
#' Default: \code{TRUE}.
#' @param graph produce the graph. \code{true} or \code{false}. default:
#' \code{true}.
#'
#' @return produce the volcano plot and silently returns the \code{ggplot}
#' object and the data.frame used.
#'
#' @examples
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de <- DESeq2::results(dds, contrast = c("group", "A", "B"))
#' volcano <- produce_volcano(de, graph = FALSE)
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_colour_identity
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 annotate
#' @importFrom utils head
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom dplyr if_else
#' @importFrom stringr str_detect
#'
#' @export
produce_volcano <- function(de_res, fc_threshold = 3, graph = TRUE) {
    stopifnot(is.numeric(fc_threshold))
    stopifnot(fc_threshold > 0)

    red <- "#E73426"
    blue <- "#0020F5"
    grey <- "#7C7C7C"

    # Rename qV to padj
    i <- stringr::str_detect(colnames(de_res), "qV")
    stopifnot(sum(i) %in% c(0,1))
    if (sum(i) == 1) {
        colnames(de_res)[i] <- "padj"
    }
    de_res <- dplyr::mutate(de_res, padj = as.numeric(padj),
                            log2FoldChange = as.numeric(log2FoldChange))

    # Remove NA padj
    de_res <- dplyr::filter(de_res, !is.na(padj))

    # Remove inf -log10(padj), i.e.: padj == 0
    min_padj <- dplyr::filter(de_res, !is.na(padj), padj != 0) %>%
        dplyr::arrange(padj) %>% head(1) %>% dplyr::pull(padj)
    i <- de_res$padj == 0
    de_res$padj[i] <- min_padj

    # Add color
    de_res <- dplyr::mutate(de_res, color = grey) %>%
        dplyr::mutate(color = dplyr::if_else(log2FoldChange <= -log2(fc_threshold) & padj <= 0.05,
                                             blue, color)) %>%
        dplyr::mutate(color = dplyr::if_else(log2FoldChange >= log2(fc_threshold) & padj <= 0.05,
                                      red, color))

    count_blue <- dplyr::filter(de_res, color == blue) %>% nrow
    count_red <- dplyr::filter(de_res, color == red) %>% nrow
    lbl <- c(count_blue, count_red) %>% as.character
    count_y <- round(max(-log10(de_res$padj), na.rm = TRUE))
    count_y <- count_y * 0.925
    min_x <- round(min(de_res$log2FoldChange, na.rm = TRUE))
    min_x <- min_x * 0.925
    max_x <- round(max(de_res$log2FoldChange, na.rm = TRUE))
    max_x <- max_x * 0.925
    p <- ggplot2::ggplot(de_res, ggplot2::aes(x = log2FoldChange,
                                              y = -log10(padj),
                                              color = color)) +
        ggplot2::geom_point(size = 3, alpha = 0.8) +
        ggplot2::scale_colour_identity() +
        ggplot2::theme_minimal() +
        ggplot2::geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)), linetype = "dashed") +
        ggplot2::annotate("text",
                 x = c(min_x, max_x),
                 y = count_y,
                 label = lbl,
                 size = 8,
                 fontface = 2,
                 color = c(blue, red))
    if (isTRUE(graph)) {
        print(p)
    }
    invisible(list(p = p, df = de_res))
}
