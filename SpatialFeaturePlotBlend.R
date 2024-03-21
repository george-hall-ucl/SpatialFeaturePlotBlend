SpatialFeaturePlotBlend <- function(object, features, combine = TRUE,
                                    feature_1_alt_name = NULL,
                                    feature_2_alt_name = NULL, assay = NULL,
                                    ...)  {
    # Convert decimal number to hexadecimal. Pad with 0s if only a single
    # character following conversion.
    as_hex <- function(num) {
        hex_str <- as.character(as.hexmode(num))
        if (nchar(hex_str) == 1) {
            hex_str <- paste0("0", hex_str)
        }

        return(hex_str)
    }

    metadata_to_hexadecimal <- function(in_dat) {
        apply(in_dat, 2,
              function(x) {
                  # Make minimum 0
                  x - min(x)
              }) %>%
        apply(2,
              function(x) {
                  # Constrain to range [0, 255]
                  round(255 * (x / max(x)))
              }) %>%
        apply(1,
              function(x) {
                  # Convert to hexadecimal codes
                  toupper(paste0("#", as_hex(x[1]), as_hex(x[2]), "00"))
              })
    }

    if (length(features) != 2) {
        stop(paste(c("Incorrect number of features. ",
                     "Requires two features, received ",
                     length(features))))
    }

    if (!is.null(assay)) {
        DefaultAssay(object) <- assay
    }

    blend_plot_theme <- theme(legend.position = "none",
                              plot.title = element_text(hjust = 0.5))

    plot_list_outer <- list()

    object$image_id <- unlist(lapply(Images(object),
                                        function(x) {
                                            rep(x, nrow(object[[x]]@coordinates))
                                        }))


    for (i in Images(object)) {
        plot_list <- lapply(features,
                            function(feature) {
                                max_color <- ifelse(feature == features[1],
                                                    "#FF0000", "#00FF00")
                                SpatialFeaturePlot(object, feature,
                                                   images = i, ...) +
                                    scale_fill_gradient(low = "#000000",
                                                        high = max_color) +
                                    ggtitle(feature) +
                                    blend_plot_theme
                            })

        cells_obj_sub <- subset(object, image_id == i)
        dat <- FetchData(cells_obj_sub, features)
        colors <- as.matrix(dat) %>% metadata_to_hexadecimal()

        new_md_column <- paste0(features[1], "_vs_", features[2])
        cells_obj_sub[[new_md_column]] <- colors
        names(colors) <- as.character(colors)

        plot_list[[3]] <- SpatialDimPlot(cells_obj_sub, new_md_column,
                                         cols = colors, images = i, ...) +
                            ggtitle(paste0(features[1], "_", features[2])) +
                            blend_plot_theme

        side_length <- 100
        legend_grid <- expand.grid(seq(from = min(dat[, features[1]]),
                                       to = max(dat[, features[1]]),
                                       length.out = side_length),
                                   seq(from = min(dat[, features[2]]),
                                       to = max(dat[, features[2]]),
                                       length.out = side_length))
        colnames(legend_grid) <- features
        legend_colors <- metadata_to_hexadecimal(legend_grid)
        legend_grid$color <- legend_colors
        names(legend_colors) <- legend_colors

        legend <- ggplot(legend_grid,
                         aes(x = .data[[features[1]]], y = .data[[features[2]]],
                             color = color)) +
                    geom_point(shape = 15, size = 1.9) +
                    scale_color_manual(values = legend_colors) +
                    coord_cartesian(expand = FALSE) +
                    theme(legend.position = "none", aspect.ratio = 1,
                          panel.background = element_blank(),
                          axis.text.x = element_text(angle = 45, hjust = 1)) +
                    xlab(ifelse(is.null(feature_1_alt_name),
                                features[1], feature_1_alt_name)) +
                    ylab(ifelse(is.null(feature_2_alt_name),
                                features[2], feature_2_alt_name))

        plot_list[[4]] <- wrap_plots(ggplot() + theme_void(), legend,
                                     ggplot() + theme_void(), ncol = 1,
                                     heights = c(0.2, 0.6, 0.2))

        plot_list_outer[[i]] <- plot_list
    }

    if (combine == FALSE) {
        return(plot_list_outer)
    } else {
        plot_list_outer <- lapply(plot_list_outer,
                                  function(p) {
                                      wrap_plots(p, nrow = 1,
                                                 widths = c(0.28, 0.28,
                                                            0.28, 0.16))
                                  })
        p <- wrap_plots(plot_list_outer, ncol = 1)

        return(p)
    }
}
