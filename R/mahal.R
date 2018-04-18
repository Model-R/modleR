#' Fits ecological niche models using Mahalanobis distance.
#'
#' @inheritParams do_bioclim
#' @return A data frame with the evaluation statistics (TSS, AUC, etc.)
#' @export
do_mahal <- function(sp,
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
                     n.back) {

  cat(paste("Mahalanobis distance", "\n"))


  if (file.exists(paste0(models.dir)) == FALSE)
       dir.create(paste0(models.dir))
    if (file.exists(paste0(models.dir, "/", sp)) == FALSE)
     dir.create(paste0(models.dir, "/", sp))
    partition.folder <- paste0(models.dir, "/", sp, "/present", "/partitions")
    if (file.exists(partition.folder) == FALSE)
        dir.create(partition.folder, recursive = T)

  # tabela de valores
  presvals <- raster::extract(predictors, coordinates)

  if (buffer %in% c("mean", "max", "median")) {
      backgr <- create_buffer(coord = coordinates,
                              n.back = n.back,
                              buffer.type = buffer,
                              seed = seed,
                              predictors = predictors)
  } else {
    set.seed(seed + 2)
    backgr <- dismo::randomPoints(mask = predictors,
                                  n = n.back,
                                  p = coordinates,
                                  excludep = T)
  }

  colnames(backgr) <- c("lon", "lat")

  # Extraindo dados ambientais dos bckgr
  backvals <- raster::extract(predictors, backgr)
  pa <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals)))

  # Data partition
  if (nrow(coordinates) < 11)
    partitions <- nrow(coordinates)
  set.seed(seed)  #reproducibility
  group <- dismo::kfold(coordinates, partitions)
  set.seed(seed + 1)
  bg.grp <- dismo::kfold(backgr, partitions)
  group.all <- c(group, bg.grp)

  pres <- cbind(coordinates, presvals)
  back <- cbind(backgr, backvals)
  rbind_1 <- rbind(pres, back)
  sdmdata <- data.frame(cbind(group.all, pa, rbind_1))
  rm(rbind_1)
  rm(pres)
  rm(back)
  gc()
  write.table(sdmdata, file = paste0(partition.folder, "/sdmdata.txt"))


  ##### Hace los modelos
  for (i in unique(group)) {
    cat(paste(sp, "partition number", i, "\n"))
    pres_train <- coordinates[group != i, ]
    if (nrow(coordinates) == 1)
      pres_train <- coordinates[group == i, ]
    pres_test <- coordinates[group == i, ]

    backg_train <- backgr[bg.grp != i, ]  #not used?
    backg_test <- backgr[bg.grp == i, ]  #new

    sdmdata_train <- subset(sdmdata, group != i)  #new
    sdmdata_test <- subset(sdmdata, group == i)  #new

    envtrain <- subset(sdmdata_train, select = c(-group, -lon, -lat))  #new
    envtest <- subset(sdmdata_test, select = c(-group, -lon, -lat))
    envtest_pre <- subset(sdmdata_test, pa == 1, select = c(-group, -lon, -lat,
      -pa))  #new
    envtest_back <- subset(sdmdata_test, pa == 0, select = c(-group, -lon, -lat,
      -pa))  #new

    ma <- dismo::mahal(predictors, pres_train)
    if (exists("ma")) {
        ema <- dismo::evaluate(pres_test, backg_test, ma, predictors)
        thresholdma <- ema@t[which.max(ema@TPR + ema@TNR)]
        thma <- dismo::threshold(ema)
        ma_TSS <- max(ema@TPR + ema@TNR) - 1
        ma_cont <- dismo::predict(ma, predictors, progress = "text")
        ma_cont[ma_cont < dismo::threshold(ema, "no_omission")] <-
            dismo::threshold(ema, "no_omission")
        ma_bin <- ma_cont > thresholdma
        ma_cut <- ma_cont
        ma_cut[ma_cut < thresholdma] <- thresholdma
        if (raster::minValue(ma_cut) < 0) {
            ma_cut <-
                (ma_cut - raster::minValue(ma_cut)) /
                raster::maxValue(ma_cut - raster::minValue(ma_cut))
        }

    thma$AUC <- ema@auc
    thma$TSS <- ma_TSS
    thma$algoritmo <- "mahal"
    thma$partition <- i
    row.names(thma) <- paste(sp, i, "mahal")

    write.table(thma, file = paste0(partition.folder, "/evaluate",
      sp, "_", i, "_mahal.txt"))

    if (class(mask) == "SpatialPolygonsDataFrame") {
      ma_cont <- crop_model(ma_cont, mask)
      ma_bin <- crop_model(ma_bin, mask)
      ma_cut <- crop_model(ma_cut, mask)
    }
    raster::writeRaster(x = ma_cont,
                        filename = paste0(partition.folder, "/mahal_cont_",
                                          sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = ma_bin,
                        filename = paste0(partition.folder, "/mahal_bin_",
                                          sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = ma_cut,
                        filename = paste0(partition.folder, "/mahal_cut_",
                                          sp, "_", i, ".tif"), overwrite = T)

       if (write_png == T) {
           png(filename = paste0(partition.folder,
                                 "/mahal_cont_", sp, "_", i, ".png"))
           plot(ma_cont,
                main = paste("mahal raw", "\n",
                             "AUC =", round(ema@auc, 2), "-",
                             "TSS =", round(ma_TSS, 2)))
         dev.off()
         png(filename = paste0(partition.folder,
                                 "/mahal_bin_", sp, "_", i, ".png"))
           plot(ma_bin, main = paste("mahal P/A", "\n",
                                     "AUC =", round(ema@auc, 2), "-",
                                     "TSS =", round(ma_TSS, 2)))
         dev.off()
         png(filename = paste0(partition.folder,
                                 "/mahal_cut_", sp, "_", i, ".png"))
           plot(ma_cut, main = paste("mahal cut", "\n",
                                     "AUC =", round(ema@auc, 2), "-",
                                     "TSS =", round(ma_TSS, 2)))
           dev.off()
           }

    if (project.model == T) {
      for (proj in projections) {
      projection.folder <- paste0(models.dir, "/", sp, "/", proj)
            if (file.exists(projection.folder) == FALSE)
                dir.create(paste0(projection.folder), recursive = T)

        data <- list.files(paste0("./env/", proj), pattern = proj)
        data2 <- stack(data)
        ma_proj <- predict(data2, ma, progress = "text")
        ma_proj_bin <- ma_proj > thresholdma
        ma_proj_cut <- ma_proj_bin * ma_proj
        # Normaliza o modelo cut
        ma_proj_cut <- ma_proj_cut / maxValue(ma_proj_cut)
        if (class(mask) == "SpatialPolygonsDataFrame") {
          ma_proj <- crop_model(ma_proj, mask)
          ma_proj_bin <- crop_model(ma_proj_bin, mask)
          ma_proj_cut <- crop_model(ma_proj_cut, mask)
        }
          writeRaster(x = ma_proj,
                      filename = paste0(projection.folder, "/mahal_cont_", sp,
                                        "_", i, ".tif"), overwrite = T)
          writeRaster(x = ma_proj_bin,
                      filename = paste0(projection.folder, "/mahal_bin_", sp,
                                        "_", i, ".tif"), overwrite = T)
          writeRaster(x = ma_proj_cut,
                      filename = paste0(projection.folder, "/mahal_cut_", sp,
                                        "_", i, ".tif"), overwrite = T)
          rm(data2)
        }
      }
    } else cat("mahalanobis distance did not run")
  }
  return(thma)
}
