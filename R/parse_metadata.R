#' Parse metadata LO
#'
#' Parse metadata and generate a draft for pca_info, volcano_info, report_info
#' files.
#'
#' @param metadata a dataframe
#'
#' @return a list of data.frame
#'
#' @importFrom dplyr left_join mutate case_when select pull
#' @importFrom purrr imap_dfr
#'
#' TODO: export and documentation

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

    # report add figure
    counter_report <- counter_report + 1
    report_info_df[[counter_report]] <- list(id = counter_report,
                                             add = "plot",
                                             value = id_pca)


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
        id_volcano <- paste0("volcano_",counter_obj)

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


#' wrapper
#' @import checkmate checkPathForOutput
#' @import dplyr mutate
#' TODO: export documentation
wrapper_report_LO <- function(metadata, txi, outdir){

    # check metadata
    stopifnot(is(metadata, "data.frame"))
    stopifnot(colnames(metadata) %in% c("ID", "Compound", "Time", "Vehicule", "Dose",  "Cell"))

    # check txi
    validate_txi(txi)

    # check outdir and create out folders
    checkmate::checkPathForOutput(outdir, overwrite = TRUE)
    r_objects <- paste0(outdir,"/r_objects/")
    path_png <- paste0(outdir,"/png/")
    dir.create(r_objects, showWarnings = TRUE, recursive = TRUE)
    dir.create(path_png, showWarnings = TRUE, recursive = TRUE)

    results <- list()

    # 1) parse metadata
    parse_res <- parse_metadata_for_LO_report(metadata.file)

    # 2) from metadata, do batch pca
    results[["pca"]] <- batch_pca(pca_infos = parse_res$pca_info,
                                  txi = txi, metadata = metadata.file,
                                  r_objects = r_objects, outdir = path_png)


    # 3) from metadata, do batch de
    results[["de"]] <- batch_de(de_infos = parse_res$de_info,
                                txi = txi,
                                design = parse_res$design_info,
                                r_objects = r_objects,
                                outdir = path_png)

    # 4) from metadata and batch_de results, do batch volcano
    results[["volcano"]] <- batch_volcano(volcano_infos = parse_res$volcano_info,
                                          de_results = r_objects, # unique ids
                                          r_objects = r_objects, # unique ids
                                          outdir = path_png)

    # 5) produce report
    # TODO: need to change parse_res$report_info to include path of object
    parse_res$report_info <- parse_res$report_info %>%
        dplyr::mutate(value = ifelse(add == "plot",
                                     paste0(r_objects, value, ".rds"),
                                     value))
    results[["parse_metadata"]] <- parse_res

    results[["report"]] <- produce_report(report_infos = parse_res$report_info,
                   report_filename = paste0(outdir, "/report.rmd"))

    return(invisible(results))
}


