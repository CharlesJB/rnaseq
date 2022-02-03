#' Produce report in Rmd format
#'
#' The function will build a rmarkdown file based on the the
#' \code{report_infos} table. The table must contain at least the following
#' columns: \code{add} and \code{value}. The \code{add} must be one of the
#' following: text, file or plot. The value will vary depending on the type of
#' add.
#'
#' For text, it must be a character string corresponding to the text that
#' will be directly added to the file specified with \code{report_filename}
#' file.
#'
#' For file, it must be a valid path to a filename
#'
#' For plot, it must be a valid path to a \code{rds} file that contains the
#' plot to show. The plot must be into a format that can be shown using the
#' \code{print} function
#'
#' The content of the rmarkdown will be created by parsing each line of the
#' \code{report_infos} table and adding the requested elements in the same
#' order as they are found in the table.
#'
#' There will be no validation that the created rmarkdown file is in the
#' correct format!
#'
#' @param report_infos A \code{data.frame} or the path to a csv file that
#' describes the report to produce.
#' @param report_filename The name of the
#' rmarkdown file to create. Default: report.Rmd
#'
#' @return Invisibly returns the \code{report_infos} table.
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_csv
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom tools file_path_sans_ext
#' @importFrom stringr str_detect
#'
#' @export
produce_report <- function(report_infos, report_filename = "report.Rmd") {

    # 1. Data validation
    stopifnot(is(report_infos, "data.frame") | is(report_infos, "character"))
    if (is.character(report_infos)) {
        stopifnot(file.exists(report_infos))
        report_infos <- readr::read_csv(report_infos, show_col_types = FALSE)
    }
    stopifnot(all(c("add", "value") %in% colnames(report_infos)))
    stopifnot(nrow(report_infos) > 0)
    stopifnot(all(report_infos$add %in% c("text", "file", "plot")))
    all_ids <- dplyr::filter(report_infos, add == "plot") %>%
            dplyr::pull(value) %>%
            basename %>%
            tools::file_path_sans_ext()
    stopifnot(!any(duplicated(all_ids)))
    for (i in 1:nrow(report_infos)) {
        current_add <- report_infos$add[i]
        current_value <- report_infos$value[i]
        if (current_add == "text") {
            stopifnot(is(current_value, "character"))
        } else if (current_add == "file") {
            stopifnot(file.exists(current_value))
        } else if (current_add == "plot") {
            stopifnot(file.exists(current_value))
            stopifnot(stringr::str_detect(current_value, "\\.rds$"))
        }
    }

    # 2. Parse report_infos
    lines <- character()
    for (i in 1:nrow(report_infos)) {
        current_add <- report_infos$add[i]
        current_value <- report_infos$value[i]

        if (current_add == "text") {
            lines <- c(lines, current_value, "\n")
        } else if (current_add == "file") {
            lines <- c(lines, readLines(current_value), "\n")
        } else if (current_add == "plot") {
            current_id <- basename(current_value) %>%
                tools::file_path_sans_ext()
            current_line <- paste0("```{r ", current_id, "}\n")
            current_line <- paste0(current_line, "print(readRDS('", current_value, "'))", "\n")
            current_line <- paste0(current_line, "```")
            lines <- c(lines, current_line, "\n")
        }
    }
    file_conn<-file(report_filename)
    writeLines(lines, file_conn)
    close(file_conn)
    invisible(lines)
}
