get_pca_df <- function(X, metadata, ncp = 2){
    # X = object from prepare_data
    pca.res <- FactoMineR::PCA(X, ncp = ncp, graph = FALSE)
    pca.coord <- pca.res$ind$coord[,c(1:ncp)]  %>% # only 2 comp
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        left_join(metadata, by = c("sample" = "ID"))

    xlab <- paste0("Dim1 (", pca.res$eig[1,2] %>% round(2), "%)")
    ylab <- paste0("Dim2 (", pca.res$eig[2,2] %>% round(2), "%)")

    # pca.obj
    return(list(coord = pca.coord,
                xlab = xlab,
                ylab = ylab))
}

plot_pca <- function(pca.obj, color = NULL, shape = NULL, show_name = TRUE, title = NULL){

    # check pca.obj: list with pca.coord, xlab, ylab
    # check color shape: NULL or "string", must be present in colnames pca.obj$coord
    # check show_name: T/F
    # check_title: string or NULL

    # to avoid 6 shapes limitation
    shape_vec <- c(15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14)

    # check if data has Dim.1, Dim.2
    gg <- ggplot(pca.obj$coord, aes(x = Dim.1, y = Dim.2)) +
        theme_bw() +
        labs(x = pca.obj$xlab, y = pca.obj$ylab)

    # add color and shape
    if(!is.null(color) & !is.null(shape)){
        gg <- gg +  geom_point(aes_string(col = color, shape = shape), size = 3) +
            scale_shape_manual(values = shape_vec[1:length(unique(pca.obj$coord$Dose))])

    } else if(!is.null(color) & is.null(shape)){
        gg <- gg +  geom_point(aes_string(col = color), size = 3)

    } else if(is.null(color) & !is.null(shape)){
        gg <- gg +  geom_point(aes_string(shape = shape), size = 3) +
            scale_shape_manual(values = shape_vec[1:length(unique(pca.obj$coord$Dose))])

    } else{
        gg <- gg +  geom_point(size = 3)
    }

    # show_name
    if(!is.null(show_name)){
        gg <- gg + ggrepel::geom_text_repel(ggplot2::aes(label = sample),
                                            color = "black", force = 10)
    }

    # title
    if(!is.null(title)){
        gg <- gg + ggtitle(title)
    }

    return(gg)
}
