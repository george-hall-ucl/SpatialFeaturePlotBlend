SpatialFeaturePlotBlend <- function(object, features, combine = TRUE,
                                    feature_1_alt_name = NULL,
                                    feature_2_alt_name = NULL, assay = NULL,
                                    bottom_left = "#000000",
                                    bottom_right = "#FF0000",
                                    top_left = "#00FF00",
                                    top_right = "#FFFF00",
                                    use_seurat_backend = FALSE,
                                    fp_extra_arguments = list(),
                                    sfp_extra_arguments = list())  {

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

    custom_color_SpatialDimPlot <- function(cells_obj, image_name,
                                            new_md_column_name,
                                            colors_per_spot, ...) {
        cells_obj[[new_md_column_name]] <- colors_per_spot
        names(colors_per_spot) <- as.character(colors_per_spot)

        p <- SpatialDimPlot(cells_obj, new_md_column_name,
                            cols = colors_per_spot, images = image_name, ...) +
                            ggtitle(new_md_column_name) +
                            blend_plot_theme
        return(p)
    }

    extract_colors_from_ggplot <- function(p) {
        built <- ggplot_build(p)$data[[1]]
        if (!is.na(built[1, "fill"])) {
            col_to_use <- "fill"
        } else {
            col_to_use <- "colour"
        }
        return(built[, col_to_use])
    }

    if (length(features) != 2) {
        stop(paste(c("Incorrect number of features. ",
                     "Requires two features, received ",
                     length(features))))
    }

    if (!is.null(assay)) {
        DefaultAssay(object) <- assay
    }

    if (length(fp_extra_arguments) > 0) {
        use_seurat_backend <- TRUE
    }

    blend_plot_theme <- theme(legend.position = "none",
                              plot.title = element_text(hjust = 0.5))

    plot_list_outer <- list()

    for (i in Images(object)) {
        cell_barcodes <- Seurat:::CellsByImage(object, images = i,
                                               unlist = TRUE)
        cells_obj_sub <- object[, cell_barcodes]
        images_sub_list <- list(object[[i]])
        names(images_sub_list) <- i
        cells_obj_sub@images <- images_sub_list
        if (!use_seurat_backend) {
            plot_list <- lapply(features,
                                function(feature) {
                                    max_color <- ifelse(feature == features[1],
                                                        bottom_right, top_left)
                                    SpatialFeaturePlot(object, feature,
                                                       images = i,
                                                       sfp_extra_arguments) +
                                        scale_fill_gradient(low = bottom_left,
                                                            high = max_color) +
                                        ggtitle(feature) +
                                        blend_plot_theme
                                })
            colors_list <- lapply(plot_list, extract_colors_from_ggplot)

            # Now construct the blended plot
            dat <- FetchData(cells_obj_sub, features)
            side_length <- 100
            col_grid <- gen_color_grid(side_length, bottom_left, bottom_right,
                                       top_left, top_right)
            dat_norm <- apply(dat, 2,
                              function(x) {
                                  round((side_length - 1) * x / max(x)) + 1
                              })
            colors_list[[3]] <- sapply(seq_len(nrow(dat_norm)),
                                       function(x) {
                                           col_grid[dat_norm[x, 1],
                                                    dat_norm[x, 2]]
                                       })
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
                             aes(x = .data[[features[1]]],
                                 y = .data[[features[2]]],
                                 color = color)) +
                        geom_point(shape = 15, size = 1.9) +
                        scale_color_manual(values = legend_colors) +
                        coord_cartesian(expand = FALSE) +
                        theme(legend.position = "none", aspect.ratio = 1,
                              panel.background = element_blank(),
                              axis.text.x = element_text(angle = 45,
                                                         hjust = 1)) +
                        xlab(ifelse(is.null(feature_1_alt_name),
                                    features[1], feature_1_alt_name)) +
                        ylab(ifelse(is.null(feature_2_alt_name),
                                    features[2], feature_2_alt_name))
        } else {
            if (top_right != "#FFFF00") {
                warning(paste("Cannot alter color in top right corner when",
                              "use_seurat_backend is TRUE"))
            }
            vis_reduc <- cells_obj_sub@images[[i]]@coordinates[, c(3, 2)]
            colnames(vis_reduc) <- c("vis_1", "vis_2")
            vis_reduc$vis_2 <- -1 * vis_reduc$vis_2
            vis_reduc_mat <- as.matrix(vis_reduc)
            vis_reduc_obj <- CreateDimReducObject(embeddings = vis_reduc_mat,
                                                  key = "vis_")
            cells_obj_sub@reductions$vis <- vis_reduc_obj
            seurat_fp <- do.call(FeaturePlot, c(list(object = cells_obj_sub,
                                                     features = features,
                                                     reduction = "vis",
                                                     blend = TRUE,
                                                     cols = c(bottom_left,
                                                              bottom_right,
                                                              top_left),
                                                     combine = FALSE),
                                                fp_extra_arguments))
            colors_list <- lapply(seurat_fp[1:3], extract_colors_from_ggplot)
            legend <- seurat_fp[[4]]
        }

        names(colors_list) <- c(features, paste0(features[1], "_", features[2]))
        plot_list <- lapply(names(colors_list),
                            function(x) {
                                do.call(custom_color_SpatialDimPlot,
                                        c(list(cells_obj = cells_obj_sub, i,
                                               x, colors_list[[x]]),
                                          sfp_extra_arguments))
                            })

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
