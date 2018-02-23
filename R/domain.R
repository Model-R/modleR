#' Fits ecological niche models using Domain.
#'
#' @inheritParams do_bioclim
#' @return A data frame with the evaluation statistics (TSS, AUC, etc.)
#' @export
do_domain <- function(sp,
                      coordinates,
                      partitions,
                      buffer,
                      seed,
                      predictors,
                      models.dir,
                      project.model,
                      projections,
                      mask,
                      write_png,
                      n.back) {
  cat(paste("Domain", "\n"))

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

    backg_test <- backgr[bg.grp == i, ]  #new

    do <- dismo::domain(predictors, pres_train)
    edo <- dismo::evaluate(pres_test, backg_test, do, predictors)
    thresholddo <- edo@t[which.max(edo@TPR + edo@TNR)]
    thdo <- dismo::threshold(edo)
    do_TSS <- max(edo@TPR + edo@TNR) - 1
    do_cont <- dismo::predict(predictors, do, progress = "text")
    do_bin <- do_cont > thresholddo
    do_cut <- do_cont * do_bin
    thdo$AUC <- edo@auc
    thdo$TSS <- do_TSS
    thdo$algoritmo <- "Domain"
    thdo$partition <- i
    row.names(thdo) <- paste(sp, i, "Domain")

    write.table(thdo, file = paste0(partition.folder, "/evaluate",
      sp, "_", i, "_domain.txt"))

    if (class(mask) == "SpatialPolygonsDataFrame") {
      do_cont <- crop_model(do_cont, mask)
      do_bin  <- crop_model(do_bin, mask)
      do_cut  <- crop_model(do_cut, mask)
    }
    raster::writeRaster(x = do_cont,
                        filename = paste0(partition.folder, "/Domain_cont_",
      sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = do_bin,
                        filename = paste0(partition.folder, "/Domain_bin_",
      sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = do_cut,
                        filename = paste0(partition.folder, "/Domain_cut_",
      sp, "_", i, ".tif"), overwrite = T)

    if (write_png == T) {
        png(filename = paste0(partition.folder,
                              "/Domain", sp, "_", i, "%03d.png"))
        plot(do_cont,
             main = paste("Domain raw", "\n", "AUC =", round(edo@auc, 2), "-",
                          "TSS =", round(do_TSS, 2)))
        plot(do_bin,
             main = paste("Domain P/A", "\n", "AUC =", round(edo@auc, 2), "-",
                          "TSS =", round(do_TSS, 2)))
        plot(do_cut,
             main = paste("Domain cut", "\n", "AUC =", round(edo@auc, 2), "-",
                          "TSS =", round(do_TSS, 2)))
        dev.off()
        }

    if (project.model == T) {
      for (proj in projections) {
      projection.folder <- paste0(models.dir, "/", sp, "/", proj)
            if (file.exists(projection.folder) == FALSE)
                dir.create(paste0(projection.folder), recursive = T)

        data <- list.files(paste0("./env/", proj), pattern = proj)
        data2 <- stack(data)
        do_proj <- predict(data2, do, progress = "text")
        do_proj_bin <- do_proj > thresholddo
        do_proj_cut <- do_proj_bin * do_proj
        # Normaliza o modelo cut
        do_proj_cut <- do_proj_cut / maxValue(do_proj_cut)
        if (class(mask) == "SpatialPolygonsDataFrame") {
          source("./fct/cropModel.R")
          do_proj     <- crop_model(do_proj, mask)
          do_proj_bin <- crop_model(do_proj_bin, mask)
          do_proj_cut <- crop_model(do_proj_cut, mask)
        }
        writeRaster(x = do_proj,
                    filename = paste0(projection.folder, "/Domain_cont_", sp,
                                      "_", i, ".tif"), overwrite = T)
        writeRaster(x = do_proj_bin,
                    filename = paste0(projection.folder, "/Domain_bin_", sp,
                                      "_", i, ".tif"), overwrite = T)
        writeRaster(x = do_proj_cut,
                    filename = paste0(projection.folder, "/Domain_cut_", sp,
                                      "_", i, ".tif"), overwrite = T)
        rm(data2)
      }
    }
  }
  return(thdo)
}
