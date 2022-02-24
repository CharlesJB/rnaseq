#' Parse metadata LO
#'
#' Parse metadata and generate a draft for pca_info, volcano_info, report_info
#' files.
#'
#' @param metadata (dataframe)
#' @param pca_subset (character) column of metadata, filter for automated PCA (ex.: Cell type)
#' @param pca_batch_metadata (character) extra columns for pca coloring (biological or batch effects)
#' @param extra_count_matrix
#'
#' @return a list of data.frame
#'
#' @importFrom dplyr left_join mutate case_when select pull
#' @importFrom purrr imap_dfr
#'
parse_metadata_for_LO_report <- function(metadata,
                                         pca_subset = "Cell",
                                         pca_batch_metadata = c("Cell", "Compound"),
                                         extra_count_matrix = NULL){
    # TODO: export and documentation

    stopifnot(is(metadata, "data.frame"))
    stopifnot(all(c("ID", "Compound", "Cell", "Dose", "Time", "Vehicule") %in%
                      colnames(metadata)))
    stopifnot(pca_subset %in% colnames(metadata.file))
    stopifnot(pca_batch_metadata %in% colnames(metadata.file))

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
    # colored by metadata
    counter_obj <-  0
    pca_info_df <- list()
    for(i in pca_batch_metadata){

        # report section: General PCA
        counter_report <- counter_report + 1
        report_info_df[[counter_report]] <- list(id = counter_report,
                                                 add = "text",
                                                 value = paste0("### ", pca_batch_metadata))

        counter_obj <- counter_obj + 1
        id_pca <- paste(c("pca",counter_obj,i), collapse = '_')
        pca_info_df[[id_pca]] <- list(id_plot=id_pca, id_metadata = "ID",
                                      group = NA, group_val = NA,
                                      use_normalisation = "none",
                                      min_counts = 5,
                                      size = 3,
                                      shape = NA,
                                      color = i,
                                      title = paste0("General ",i),
                                      legend.position = "right",
                                      legend.box = "vertical",
                                      show_names = TRUE)
        # report add figure
        counter_report <- counter_report + 1
        report_info_df[[counter_report]] <- list(id = counter_report,
                                                 add = "plot",
                                                 value = id_pca)
    }

    sub <- unique(metadata[,pca_subset, drop = TRUE])
    if(length(sub)>1){ # if ==1 -> same as general
        for(i in sub){
            # report section:
            counter_report <- counter_report + 1
            report_info_df[[counter_report]] <- list(id = counter_report,
                                                     add = "text",
                                                     value = paste0("## ", sub))
            for(j in pca_batch_metadata){
                # report section:
                counter_report <- counter_report + 1
                report_info_df[[counter_report]] <- list(id = counter_report,
                                                         add = "text",
                                                         value = paste0("### ", pca_batch_metadata))

                counter_obj <- counter_obj + 1
                id_pca <- paste(c("pca",counter_obj, pca_subset, i, j), collapse = "_")
                pca_info_df[[id_pca]] <- list(id_plot=id_pca, id_metadata = "ID",
                                              group = pca_subset, group_val = i,
                                              use_normalisation = "none",
                                              min_counts = 5,
                                              size = 3,
                                              shape = NA,
                                              color = j,
                                              title = paste(pca_subset, i, j),
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
        id_de <- paste(c("de",counter_obj, current_contrast_1, current_contrast_2), collapse = "_")

        de_info_df[[id_de]] <- list(id_de = id_de,
                                    group = id_de, # always Compound vs Control
                                    contrast_1 = current_contrast_1,
                                    contrast_2 = current_contrast_2,
                                    formula = paste("~", id_de),
                                    filter = 2)
        if(!is.null(extra_count_matrix)){
            de_info_df[[id_de]][["count_matrix"]] <- extra_count_matrix
        }


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
        id_volcano <- paste(c("volcano",counter_obj, current_contrast_1, current_contrast_2), collapse = "_")

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


#' wrapper report LO
#'
#' From metadata and txi, this wrapper parse the metadata and perform all the PCA and volcano (DE).
#' It produces a Rmd file which can be rendered.
#'
#'
#' @param metadata (data.frame) sample sheet with at least the columns: ID, Compound, Time, Vehicule, Dose, Cell
#' @param txi (list) txi object
#' @param pca_subset (character) column of metadata, filter for automated PCA (ex.: Cell type)
#' @param pca_batch_metadata (character) extra columns for pca coloring (biological or batch effects)
#' @param do_pca (logical) do PCA and included the results in report (default = TRUE)
#' @param do_DE (logical) do DE and included the results in report (default = TRUE)
#' @param render_repport (logical) render Rmd report
#' @param custom_parsed_metadata (list of data.frame) when automatic parsing can be problematic
#'
#' @return a list containing pca, volcano and report
#'
#' @importFrom checkmate checkPathForOutput assert_logical assert_character
#' @importFrom dplyr mutate filter
#' @importFrom stringr str_detect
#' @importFrom rmarkdown render
#'
#' @export
wrapper_report_LO <- function(metadata, txi, outdir, pca_subset, pca_batch_metadata,
                              do_pca = TRUE, do_DE = TRUE, render_repport = TRUE,
                              extra_count_matrix = NULL,
                              custom_parsed_metadata = NULL,
                              rmd_out_filepath = "./report.Rmd"){

    # check metadata
    stopifnot(is(metadata, "data.frame"))
    stopifnot(c("ID", "Compound", "Time", "Vehicule", "Dose",  "Cell") %in% colnames(metadata))

    stopifnot(is(metadata, "data.frame"))
    stopifnot(all(c("ID", "Compound", "Cell", "Dose", "Time", "Vehicule") %in%
                      colnames(metadata)))

    stopifnot(pca_subset %in% colnames(metadata))
    stopifnot(pca_batch_metadata %in% colnames(metadata))

    checkmate::assert_character(extra_count_matrix, null.ok = TRUE)

    checkmate::assert_logical(do_pca)
    checkmate::assert_logical(do_DE)
    checkmate::assert_logical(render_repport)

    # check txi
    validate_txi(txi)

    # check outdir and create out folders
    checkmate::checkPathForOutput(outdir, overwrite = TRUE)

    r_objects <- paste0(outdir,"/r_objects/")
    path_png <- paste0(outdir,"/pdf/")
    path_pca <- paste0(path_png,"/pca/")
    path_csv <- paste0(outdir,"/de_csv/")
    path_volcano <- paste0(path_png,"/volcano/")
    rmd_out_filepath <- paste0(outdir, "/", rmd_out_filepath)
    checkmate::checkPathForOutput(rmd_out_filepath, overwrite = TRUE)


    dir.create(r_objects, showWarnings = TRUE, recursive = TRUE)
    dir.create(path_png, showWarnings = TRUE, recursive = TRUE)
    dir.create(path_pca, showWarnings = TRUE, recursive = TRUE)
    dir.create(path_volcano, showWarnings = TRUE, recursive = TRUE)
    dir.create(path_csv, showWarnings = TRUE, recursive = TRUE)

    results <- list()

    # 1) parse metadata
    if(!is.null(custom_parsed_metadata)){
        stopifnot(all(sapply(custom_parsed_metadata, function(x) is(x, "data.frame"))))
        stopifnot(all(names(custom_parsed_metadata) %in% c("pca_info", "de_info", "design_info", "volcano_info", "report_info")))
        parse_res <- custom_parsed_metadata
    } else {
        parse_res <- parse_metadata_for_LO_report(metadata.file, pca_subset = pca_subset, pca_batch_metadata = pca_batch_metadata,
                                                  extra_count_matrix = extra_count_matrix)
    }
    report_info <- parse_res$report_info

    # check infos sheets: only pca and de
    # volcano needs DE files
    # report need pca, volcano files
    tmp_pca_info <- complete_pca_infos(parse_res$pca_info)
    stopifnot(length(validate_pca_infos(tmp_pca_info, metadata, txi)) == 0)
    tmp_de_infos <- complete_de_infos(parse_res$de_info)
    stopifnot(length(validate_de_infos(tmp_de_infos, parse_res$design_info, txi)) == 0)

    if(do_pca){
        # 2) from metadata, do batch pca
        results[["pca"]] <- batch_pca(pca_infos = parse_res$pca_info,
                                      txi = txi, metadata = metadata.file,
                                      r_objects = r_objects, outdir = path_pca)
        results[["parse_metadata"]] <- parse_res
    } else {
        # remove pca from report
        report_info <- report_info %>%
            dplyr::filter(!stringr::str_detect(value, "(?i)PCA"))
    }

    if(do_DE){
        # 3) from metadata, do batch de
        results[["de"]] <- batch_de(de_infos = parse_res$de_info,
                                    txi = txi,
                                    design = parse_res$design_info,
                                    r_objects = r_objects,  # DDS
                                    outdir = path_csv)  # DE res as .csv file

        # 4) from metadata and batch_de results, do batch volcano
        results[["volcano"]] <- batch_volcano(volcano_infos = parse_res$volcano_info,
                                              de_results = path_csv, # unique ids
                                              r_objects = r_objects, # unique ids
                                              outdir = path_volcano)

    } else {
        # remove volcano from report info
        report_info <- report_info %>% dplyr::filter(!stringr::str_detect(value, "volcano")) %>%
            dplyr::filter(!str_detect(value, "# Differential Expression"))
    }

    # 5) produce report anyway
    # TODO: need to change parse_res$report_info to include path of object
    # path should be relative to rmd file
    results[["parse_metadata"]][["report_info"]] <- report_info <- report_info %>%
        dplyr::mutate(value = ifelse(add == "plot",
                                     paste0("./r_objects/", value, ".rds"),
                                     value))

    results[["report"]] <- produce_report(report_infos = report_info,
                                          report_filename = rmd_out_filepath)

    if(render_repport){
        ## rmarkdown::render(...)
        rmarkdown::render(rmd_out_filepath)  # to htmls
    }

    return(invisible(results))
}


