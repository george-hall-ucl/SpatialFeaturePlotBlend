SpatialFeaturePlotBlend <- function(cells_obj, column_1, column_2,
                                    combine = TRUE, column_1_alt_name = NULL,
                                    column_2_alt_name = NULL, assay = NULL,
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

    if (!is.null(assay)) {
        DefaultAssay(cells_obj) <- assay
    }

    blend_plot_theme <- theme(legend.position = "none",
                              plot.title = element_text(hjust = 0.5))

    plot_list <- lapply(c(column_1, column_2),
                        function(column) {
                            max_color <- ifelse(column == column_1,
                                                "#FF0000", "#00FF00")
                            SpatialFeaturePlot(cells_obj, column, ...) +
                                scale_fill_gradient(low = "#000000",
                                                    high = max_color) +
                                ggtitle(column) +
                                blend_plot_theme
                        })

    dat <- FetchData(cells_obj, c(column_1, column_2))
    colors <- as.matrix(dat) %>% metadata_to_hexadecimal()

    new_md_column <- paste0(column_1, "_vs_", column_2)
    cells_obj[[new_md_column]] <- colors
    names(colors) <- as.character(colors)

    plot_list[[3]] <- SpatialDimPlot(cells_obj, new_md_column,
                                     cols = colors, ...) +
                        ggtitle(paste0(column_1, "_", column_2)) +
                        blend_plot_theme

    side_length <- 100
    legend_grid <- expand.grid(seq(from = min(dat[, column_1]),
                                   to = max(dat[, column_1]),
                                   length.out = side_length),
                               seq(from = min(dat[, column_2]),
                                   to = max(dat[, column_2]),
                                   length.out = side_length))
    colnames(legend_grid) <- c(column_1, column_2)
    legend_colors <- metadata_to_hexadecimal(legend_grid)
    legend_grid$color <- legend_colors
    names(legend_colors) <- legend_colors

    legend <- ggplot(legend_grid,
                     aes(x = .data[[column_1]], y = .data[[column_2]],
                         color = color)) +
                geom_point(shape = 15, size = 1.9) +
                scale_color_manual(values = legend_colors) +
                coord_cartesian(expand = FALSE) +
                theme(legend.position = "none", aspect.ratio = 1,
                      panel.background = element_blank()) +
                xlab(ifelse(is.null(column_1_alt_name), column_1, column_1_alt_name)) +
                ylab(ifelse(is.null(column_2_alt_name), column_2, column_2_alt_name))

    plot_list[[4]] <- wrap_plots(ggplot() + theme_void(), legend,
                                 ggplot() + theme_void(), ncol = 1,
                                 heights = c(0.2, 0.6, 0.2))

    if (combine == FALSE) {
        return(plot_list)
    } else {
        p <- wrap_plots(plot_list, nrow = 1,
                        widths = c(0.28, 0.28, 0.28, 0.16))
        return(p)
    }
}
