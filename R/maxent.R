#' Fits ecological niche models using Maxent.
#'
#' @inheritParams do_bioclim
#' @return A data frame with the evaluation statistics (TSS, AUC, etc.)
#' @export
do_maxent <- function(sp,
                      coordinates,
                      partitions,
                      buffer = FALSE,
                      seed = 512,
                      predictors,
                      models.dir = "./models",
                      project.model = FALSE,
                      projections = NULL,
                      mask = NULL,
                      write_png = FALSE,
                      n.back,
                      bootstrap = T,
                      boot_proportion = 0.8,
                      n_boot = 10,
                      crossvalidation = F,
                      n_cv = 10) {

    cat(paste("maxent", "\n"))

    #sdmdatasetup
    partition.folder <- paste0(models.dir, "/", sp, "/present", "/partitions")
    if (file.exists(paste0(partition.folder, "/sdmdata.txt"))) {
        sdmdata <- read.table(paste0(partition.folder, "/sdmdata.txt"))
        } else {
            sdmdata <- setup_sdmdata(
                sp = sp,
                coordinates = coordinates,
                partitions = partitions,
                buffer = buffer,
                seed = seed,
                predictors = predictors,
                models.dir = models.dir,
                plot_sdmdata = plot_sdmdata,
                n.back = n.back,
                bootstrap = bootstrap,
                boot_proportion = boot_proportion,
                n_boot = n_boot,
                crossvalidation = crossvalidation,
                n_cv = n_cv
            )
        }

  ##### Hace los modelos
    runs <- which(names(sdmdata) == "pa") - 1
      #para cada columna de la matriz de diseño
    for (i in seq_along(1:runs)) {
        group.all <- sdmdata[, i]
        group <- group.all[sdmdata$pa == 1]
        bg.grp <- group.all[sdmdata$pa == 0]
        backgr <- sdmdata[sdmdata$pa == 0, c("lon", "lat")]
      #para cada grupo
        for (g in unique(group)) {
            cat(paste(sp,"run number", i, "partition number", g, "\n"))
            pres_train <- coordinates[group != g, ]
            if (nrow(coordinates) == 1)
                pres_train <- coordinates[group == g,]
            pres_test <- coordinates[group == g,]
            backg_test <- backgr[bg.grp == g,]

            mx <- dismo::maxent(predictors, pres_train)
            emx <- dismo::evaluate(pres_test, backg_test, mx, predictors)
            thresholdmx <- emx@t[which.max(emx@TPR + emx@TNR)]
            thmx <- dismo::threshold(emx)
            mx_TSS <- max(emx@TPR + emx@TNR) - 1
            mx_cont <- dismo::predict(predictors, mx, progress = "text")
            mx_bin <- mx_cont > thresholdmx
            mx_cut <- mx_cont * mx_bin
            thmx$AUC <- emx@auc
            thmx$TSS <- mx_TSS
            thmx$algoritmo <- "maxent"
            thmx$run <- i
            thmx$partition <- g
            row.names(thmx) <- paste(sp, i, g, "maxent")

            write.table(thmx, file = paste0(partition.folder, "/evaluate",
                                            sp, "_", i, "_", g, "_maxent.txt"))

            if (class(mask) == "SpatialPolygonsDataFrame") {
                mx_cont <- crop_model(mx_cont, mask)
                mx_bin <- crop_model(mx_bin, mask)
                mx_cut <- crop_model(mx_cut, mask)
            }
            raster::writeRaster(
                x = mx_cont,
                filename = paste0(partition.folder, "/maxent_cont_",
                                  sp, "_", i, "_", g, ".tif"), overwrite = T)
            raster::writeRaster(
                x = mx_bin,
                filename = paste0(partition.folder, "/maxent_bin_",
                                  sp, "_", i, "_", g, ".tif"), overwrite = T)
            raster::writeRaster(
                x = mx_cut,
                filename = paste0(partition.folder, "/maxent_cut_",
                                  sp, "_", i, "_", g, ".tif"), overwrite = T)

            if (write_png == T) {
                png(paste0(partition.folder,
                           "/maxent_cont", sp, "_", i, "_", g, ".png"))
                raster::plot(mx_cont,
                             main = paste("maxent raw", "\n",
                                          "AUC =", round(emx@auc, 2), "-",
                                          "TSS =", round(mx_TSS, 2)))
                dev.off()
                png(paste0(partition.folder,
                           "/maxent_bin", sp, "_", i, "_", g, ".png"))
                raster::plot(mx_bin,
                             main = paste("maxent P/A", "\n",
                                          "AUC =", round(emx@auc, 2), "-",
                                          "TSS =", round(mx_TSS, 2)))
                dev.off()
                png(paste0(partition.folder,
                           "/maxent_cut", sp, "_", i, "_", g, ".png"))
                raster::plot(mx_cut,
                             main = paste("maxent cut", "\n",
                                          "AUC =", round(emx@auc, 2), "-",
                                          "TSS =", round(mx_TSS, 2)))
                dev.off()
                }

            if (project.model == T) {
                for (proj in projections) {
                    projection.folder <- paste0(models.dir, "/", sp, "/", proj)
                    if (file.exists(projection.folder) == FALSE)
                        dir.create(paste0(projection.folder), recursive = T)
                    data <- list.files(paste0("./env/", proj), pattern = proj)
                    data2 <- stack(data)
                    mx_proj <- predict(data2, mx, progress = "text")
                    mx_proj_bin <- mx_proj > thresholdmx
                    mx_proj_cut <- mx_proj_bin * mx_proj
                    # Normaliza o modelo cut
                    mx_proj_cut <- mx_proj_cut / maxValue(mx_proj_cut)
                    if (class(mask) == "SpatialPolygonsDataFrame") {
                        mx_proj <- crop_model(mx_proj, mask)
                        mx_proj_bin <- crop_model(mx_proj_bin, mask)
                        mx_proj_cut <- crop_model(mx_proj_cut, mask)
                        }
                    writeRaster(
                        x = mx_proj,
                        filename = paste0(projection.folder,
                                          "/maxent_cont_", sp, "_", i, "_", g, ".tif"),
                        overwrite = T)
                    writeRaster(
                        x = mx_proj_bin,
                        filename = paste0(projection.folder,
                                          "/maxent_bin_", sp, "_", i, "_", g, ".tif"),
                        overwrite = T)
                    writeRaster(
                        x = mx_proj_cut,
                        filename = paste0(projection.folder, "/maxent_cut_", sp, "_", i, "_", g, ".tif"),
                        overwrite = T)
                    rm(data2)
                }
            }
        }
    }
    return(thmx)
}

