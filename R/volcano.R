#' produce a volcano plot from DESeq2 results
#'
#' @param de_res the \code{data.frame} object returned by the
#' \code{DESeq2::results} function.
#' @param fc_threshold The threshold of FC to be considered as significant.
#' Default: \code{TRUE}.
#' @param graph produce the graph. \code{true} or \code{false}. default:
#' \code{true}.
#'
#' @return produce the volcano plot and silently returns the \code{ggplot} object.
#'
#' @examples
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de <- DESeq2::results(dds, contrast = c("group", "A", "B"))
#' volcano <- produce_volcano(de, graph = FALSE)
#'
#' @import ggplot2
#'
#' @export
produce_volcano <- function(de_res, fc_threshold = 3, graph = TRUE) {
    red <- "#E73426"
    blue <- "#0020F5"
    grey <- "#7C7C7C"

    # Rename qV to padj
    i <- str_detect(colnames(de_res), "qV")
    stopifnot(sum(i) %in% c(0,1))
    if (sum(i) == 1) {
        colnames(de_res)[i] <- "padj"
    }
    de_res <- mutate(de_res, padj = as.numeric(padj), log2FoldChange = as.numeric(log2FoldChange))

    # Remove inf -log10(padj), i.e.: padj == 0
    min_padj <- filter(de_res, !is.na(padj), padj != 0) %>%
        arrange(padj) %>% head(1) %>% pull(padj)
    i <- de_res$padj == 0
    de_res$padj[i] <- min_padj

    # Add color
    de_res <- mutate(de_res, color = grey) %>%
        mutate(color = if_else(log2FoldChange <= log2(1/fc_threshold) & padj <= 0.05, blue, color)) %>%
        mutate(color = if_else(log2FoldChange >= log2(fc_threshold) & padj <= 0.05, red, color))

    count_blue <- filter(de_res, color == blue) %>% nrow
    count_red <- filter(de_res, color == red) %>% nrow
    lbl <- c(count_blue, count_red) %>% as.character
    count_y <- round(max(-log10(de_res$padj), na.rm = TRUE))
    count_y <- count_y * 0.925
    min_x <- round(min(de_res$log2FoldChange, na.rm = TRUE))
    min_x <- min_x * 0.925
    max_x <- round(max(de_res$log2FoldChange, na.rm = TRUE))
    max_x <- max_x * 0.925
    p <- ggplot(de_res, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
        geom_point(size = 3, alpha = 0.8) +
        scale_colour_identity() +
        theme_minimal() +
        geom_vline(xintercept = c(log2(1/fc_threshold), log2(fc_threshold)), linetype = "dashed") +
        annotate("text",
                 x = c(min_x, max_x),
                 y = count_y,
                 label = lbl,
                 size = 8,
                 fontface = 2,
                 color = c(blue, red))
    if (isTRUE(graph)) {
        print(p)
    }
    invisible(p)
}