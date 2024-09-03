SpatialFeaturePlotBlend <- function(object, features, combine = TRUE,
                                    feature_1_alt_name = NULL,
                                    feature_2_alt_name = NULL, assay = NULL,
                                    bottom_left = "#000000",
                                    bottom_right = "#FF0000",
                                    top_left = "#00FF00",
                                    top_right = "#FFFF00", ...)  {

    # Generate a grid of RGB color values given the requested corner colours.
    gen_color_grid <- function(side_length, bottom_left, bottom_right,
                               top_left, top_right) {

        grad_gen <- function(start, end, n = side_length) {
            colfunc <- colorRampPalette(c(start, end))
            return(colfunc(n))
        }

        # x_y = "x to y"; "bl" = "bottom left", etc
        bl_tl <- grad_gen(bottom_left, bottom_right)
        br_tr <- grad_gen(top_left, top_right)

        l <- lapply(seq_len(length(bl_tl)),
               function(i) {
                   start <- bl_tl[i]
                   end <- br_tr[i]
                   new_grad <- grad_gen(start, end)
               })

        return(t(matrix(unlist(l), ncol = side_length, nrow = side_length)))
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

    for (i in Images(object)) {
        plot_list <- lapply(features,
                            function(feature) {
                                max_color <- ifelse(feature == features[1],
                                                    bottom_right, top_left)
                                SpatialFeaturePlot(object, feature,
                                                   images = i, ...) +
                                    scale_fill_gradient(low = bottom_left,
                                                        high = max_color) +
                                    ggtitle(feature) +
                                    blend_plot_theme
                            })

        cell_barcodes <- Seurat:::CellsByImage(object, images = i,
                                               unlist = TRUE)
        cells_obj_sub <- object[, cell_barcodes]
        images_sub_list <- list(object[[i]])
        names(images_sub_list) <- i
        cells_obj_sub@images <- images_sub_list
        dat <- FetchData(cells_obj_sub, features)
        side_length <- 100
        col_grid <- gen_color_grid(side_length, bottom_left, bottom_right,
                                   top_left, top_right)
        dat_norm <- apply(dat, 2,
                          function(x) {
                              round((side_length - 1) * x / max(x)) + 1
                          })
        colors <- sapply(seq_len(nrow(dat_norm)),
                       function(x) {
                           col_grid[dat_norm[x, 1], dat_norm[x, 2]]
                       })

        new_md_column <- paste0(features[1], "_vs_", features[2])
        cells_obj_sub[[new_md_column]] <- colors
        names(colors) <- as.character(colors)

        plot_list[[3]] <- SpatialDimPlot(cells_obj_sub, new_md_column,
                                         cols = colors, images = i, ...) +
                            ggtitle(paste0(features[1], "_", features[2])) +
                            blend_plot_theme

        legend_grid <- expand.grid(seq(from = min(dat[, features[1]]),
                                       to = max(dat[, features[1]]),
                                       length.out = side_length),
                                   seq(from = min(dat[, features[2]]),
                                       to = max(dat[, features[2]]),
                                       length.out = side_length))
        colnames(legend_grid) <- features
        legend_colors <- c(col_grid)
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
