#' Parse metadata LO
#'
#'`r lifecycle::badge('experimental')`
#' Parse metadata and generate a draft for pca_info, volcano_info, report_info
#' files.
#'
#'
#' @importFrom dplyr left_join mutate case_when select pull
#' @importFrom purrr imap_dfr
#'
#' @export


parse_metadata_for_LO_report <- function(metadata){

    stopifnot(is(metadata, "data.frame"))
    stopifnot(all(c("ID", "Compound", "Cell", "Dose", "Time", "Vehicule") %in%
                      colnames(metadata)))

    ## start report info
    ########################
    counter_report <- 1
    report_info_df = list()
    report_info_df[[counter_report]] <- list(id = counter_report,
                                             add = "text",
                                             value = "# PCA")


    ## PCA
    ########################

    # report section: General PCA
    counter_report <- counter_report + 1
    report_info_df[[counter_report]] <- list(id = counter_report,
                                             add = "text",
                                             value = "## General")

    # general PCA
    # colored by cell
    counter_obj = 0
    pca_info_df <- list()

    counter_obj <- counter_obj + 1
    id_pca <- paste0("pca_",counter_obj)
    pca_info_df[[id_pca]] <- list(
        id_plot=id_pca,
        id_metadata = "ID",
        group = NA,
        group_val = NA,
        use_normalisation = "none",
        min_counts = 5,
        size = 3,
        shape = NA,
        color = "Cell",
        title = "Cells (tpm)",
        legend.position = "right",
        legend.box = "vertical",
        show_names = TRUE)


    # split by cells
    cell_types <- unique(metadata$Cell)
    for(ct in cell_types){

        # report subsection: PCA/Cell
        counter_report <- counter_report + 1
        report_info_df[[counter_report]] <- list(id = counter_report,
                                                 add = "text",
                                                 value = paste("##", ct))

        for(met in c("Compound", "Dose", "Time", "Vehicule")){

            # increment counter
            counter_obj <- counter_obj + 1
            id_pca <- paste0("pca_",counter_obj)

            pca_info_df[[id_pca]] <- list(id_plot=id_pca, id_metadata = "ID",
                                           group = "Cell",  # filtered by cell types
                                           group_val = ct,
                                           use_normalisation = "none",
                                           min_counts = 5, size = 3,
                                           shape = NA,
                                           color = met,  # met here
                                           title = paste(ct, met, sep = " - "),
                                           legend.position = "right",
                                           legend.box = "vertical",
                                           show_names = TRUE)

            # report add figure
            counter_report <- counter_report + 1
            report_info_df[[counter_report]] <- list(id = counter_report,
                                                     add = "plot",
                                                     value = id_pca)
        }
    }
    pca_info_df <- purrr::imap_dfr(pca_info_df, ~.x)

    # report section: Volcano
    counter_report <- counter_report + 1
    report_info_df[[counter_report]] <- list(id = counter_report,
                                             add = "text",
                                             value = "# Differential Expression")

    ## DE and volcano
    ########################
    de_info_df <- list()
    volcano_info_df <- list()
    counter_obj <- 0

    metadata_design <- metadata %>%
        dplyr::select(Cell, Compound, Time, Dose, Vehicule) %>%
        unique
    design_df <- list("sample" = metadata$ID)

    for(line in 1:nrow(metadata_design)){
        # DE: Compound vs Control (Vehicule), split by Cell, Time, Dose

        current_contrast_1 <- metadata_design[line, "Compound", drop = TRUE]
        current_contrast_2 <- metadata_design[line, "Vehicule", drop = TRUE]

        if(current_contrast_1 == current_contrast_2){next}

        # increment counter
        counter_obj <- counter_obj + 1

        # DE
        id_de <- paste0("de_",counter_obj)

        de_info_df[[id_de]] <- list(id_de = id_de,
                                    group = id_de, # always Compound vs Control
                                    contrast_1 = current_contrast_1,
                                    contrast_2 = current_contrast_2,
                                    formula = paste("~", id_de),
                                    filter = 2,
                                    count_matrix = "extra_count_matrix")

        samples_contrast_1 <- metadata_design[line,] %>%
            dplyr::left_join(metadata, by = c("Cell", "Compound", "Time", "Dose", "Vehicule")) %>%
            dplyr::pull(ID)
        samples_contrast_2 <- metadata_design[line,] %>%
            dplyr::mutate(Compound = metadata_design[line, "Vehicule", drop = TRUE])  %>%
            dplyr::mutate(Dose = "Control") %>%
            dplyr::left_join(metadata, by = c("Cell", "Compound", "Time", "Dose", "Vehicule")) %>%
            dplyr::pull(ID)

        design_df[[id_de]] <- dplyr::case_when(
            design_df$sample %in% samples_contrast_1 ~ current_contrast_1,
            design_df$sample %in% samples_contrast_2 ~ current_contrast_2,
            TRUE ~ "-")

        # volcano
        id_volcano <- paste0("de_",counter_obj)

        volcano_info_df[[id_volcano]] <- list(id_plot = id_volcano,
                                         id_de = id_de,
                                         y_axis = "padj",
                                         p_threshold = 0.05,
                                         fc_threshold = 1.5,
                                         show_signif_counts = TRUE,
                                         show_signif_lines = "vertical",
                                         show_signif_color = TRUE,
                                         col_up = "#E73426",
                                         col_down = "#0020F5",
                                         size = 3)

        # report add figure: volcano
        counter_report <- counter_report + 1
        report_info_df[[counter_report]] <- list(id = counter_report,
                                                 add = "plot",
                                                 value = id_volcano)
    }
    de_info_df <- purrr::imap_dfr(de_info_df, ~.x)
    design_df <- as.data.frame(design_df)
    volcano_info_df <- purrr::imap_dfr(volcano_info_df, ~.x)
    report_info_df <- purrr::imap_dfr(report_info_df, ~.x)

    return(list(pca_info = pca_info_df,
                de_info = de_info_df,
                design_info = design_df,
                volcano_info = volcano_info_df,
                report_info = report_info_df))
}



