#' Fits ecological niche models using bioclim.
#'
#' @param sp A character string with the species name
#' @param coordinates A two-column data frame with the occurrence points
#' @param partitions The number of partitions for a cross validation
#' @param buffer Defines if a buffer will be used to sample pseudo-absences
#'        (F, "mean", "median", "max")
#' @param seed For reproducibility purposes
#' @param predictors A RasterStack of predictor variables
#' @param models.dir Folder path to save the output files
#' @param project.model Logical, whether to perform a projection
#' @param projections The RasterStack of projeciton variables
#' @param mask A SpatialPolygonsDataFrame to be used to mask the final models
#' @param write_png Logical, whether png files will be written
#' @param n.back Number of pseudoabsence points
#' @return A data frame with the evaluation statistics (TSS, AUC, etc.)
#' @import grDevices
#' @importFrom utils write.table
#' @export
do_bioclim <- function(sp,
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
  cat(paste("bioclim", "\n"))

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
      group.all <- sdmdata[,i]
      group <- group.all[sdmdata$pa == 1]
      bg.grp <- group.all[sdmdata$pa == 0]
      backgr <- sdmdata[sdmdata$pa == 0, c("lon", "lat")]
      #para cada grupo
      for (g in unique(group)) {
          cat(paste(sp,"run number", i, "partition number", g, "\n"))
          pres_train <- coordinates[group != g, ]
          if (nrow(coordinates) == 1)
              pres_train <- coordinates[group == g, ]
          pres_test <- coordinates[group == g, ]
          backg_test <- backgr[bg.grp == g, ]
          
          bc <- dismo::bioclim(predictors, pres_train)
          ebc <- dismo::evaluate(pres_test, backg_test, bc, predictors)
          thresholdbc <- ebc@t[which.max(ebc@TPR + ebc@TNR)]
          thbc <- dismo::threshold(ebc)
          bc_TSS <- max(ebc@TPR + ebc@TNR) - 1
          bc_cont <- dismo::predict(predictors, bc, progress = "text")
          bc_bin <- bc_cont > thresholdbc
          bc_cut <- bc_cont * bc_bin
          thbc$AUC <- ebc@auc
          thbc$TSS <- bc_TSS
          thbc$algoritmo <- "bioclim"
          thbc$run <- i
          thbc$partition <- g
          row.names(thbc) <- paste(sp, i, g, "bioclim")
          
          write.table(thbc, file = paste0(partition.folder, "/evaluate",
      sp, "_", i, "_", g, "_bioclim.txt"))
      
          if (class(mask) == "SpatialPolygonsDataFrame") {
              bc_cont <- crop_model(bc_cont, mask)
              bc_bin <- crop_model(bc_bin, mask)
              bc_cut <- crop_model(bc_cut, mask)
              }
    raster::writeRaster(x = bc_cont,
                        filename = paste0(partition.folder, "/bioclim_cont_",
                                          sp, "_", i,"_", g, ".tif"), overwrite = T)
    raster::writeRaster(x = bc_bin,
                        filename = paste0(partition.folder, "/bioclim_bin_",
                                          sp, "_", i, "_", g, ".tif"), overwrite = T)
    raster::writeRaster(x = bc_cut,
                        filename = paste0(partition.folder, "/bioclim_cut_",
                                          sp, "_", i, "_", g, ".tif"), overwrite = T)

  if (write_png == T) {
      png(paste0(partition.folder, "/bioclim_cont", sp, "_", i, "_", g, ".png"))
      raster::plot(bc_cont,
                   main = paste("bioclim raw", "\n",
                                "AUC =", round(ebc@auc, 2), "-",
                                "TSS =", round(bc_TSS, 2)))
      dev.off()
      png(paste0(partition.folder, "/bioclim_bin", sp, "_", i, "_", g, ".png"))
      raster::plot(bc_bin, main = paste("bioclim P/A", "\n",
                                        "AUC =", round(ebc@auc, 2), "-",
                                        "TSS =", round(bc_TSS, 2)))
      dev.off()
      png(paste0(partition.folder, "/bioclim_cut", sp, "_", i, "_", g, ".png"))
      raster::plot(bc_cut, main = paste("bioclim cut", "\n",
                                        "AUC =", round(ebc@auc, 2), "-",
                                        "TSS =", round(bc_TSS, 2)))
      dev.off()
      }

    if (project.model == T) {
      for (proj in projections) {
          projection.folder <- paste0(models.dir, "/", sp, "/", proj)
          if (file.exists(projection.folder) == FALSE)
              dir.create(paste0(projection.folder), recursive = T)
          data <- list.files(paste0("./env/", proj), pattern = proj)
          data2 <- stack(data)
          bc_proj <- predict(data2, bc, progress = "text")
          bc_proj_bin <- bc_proj > thresholdbc
          bc_proj_cut <- bc_proj_bin * bc_proj
        # Normaliza o modelo cut
         bc_proj_cut <- bc_proj_cut / maxValue(bc_proj_cut)
        if (class(mask) == "SpatialPolygonsDataFrame") {
          bc_proj     <- crop_model(bc_proj, mask)
          bc_proj_bin <- crop_model(bc_proj_bin, mask)
          bc_proj_cut <- crop_model(bc_proj_cut, mask)
          }
         writeRaster(x = bc_proj,
                     filename = paste0(projection.folder, "/bioclim_cont_", sp,
                                       "_", i, "_", g, ".tif"), overwrite = T)
         writeRaster(x = bc_proj_bin,
                     filename = paste0(projection.folder, "/bioclim_bin_", sp,
                                       "_", i, "_", g,".tif"), overwrite = T)
         writeRaster(x = bc_proj_cut,
                     filename = paste0(projection.folder, "/bioclim_cut_", sp,
                                       "_", i, "_", g, ".tif"), overwrite = T)
         rm(data2)
         }
        }
      }
    }

  return(thbc)
}
