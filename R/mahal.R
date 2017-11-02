#' Faz modelagem de distribuição de espécies com distância de Mahalanobis
#'
#' @inheritParams do_bioclim
#' @return Um data.frame com metadados da modelagem (TSS, AUC, algoritmo etc.)
#' @export
do_mahal <- function(sp,
		     coordinates,
		     partitions,
		     buffer = FALSE,
		     seed = 512,
		     predictors,
		     models.dir,
		     project.model,
		     projections,
		     mask,
		     n.back = 500) {
  cat(paste("Mahalanobis distance", "\n"))


  if (file.exists(paste0(models.dir)) == FALSE)
    dir.create(paste0(models.dir))
  if (file.exists(paste0(models.dir, "/", sp)) == FALSE) 
    dir.create(paste0(models.dir, "/", sp))
  if (project.model == T) {
    for (proj in projections) {
      if (file.exists(paste0(models.dir, "/", sp, "/", proj)) == FALSE) 
        dir.create(paste0(models.dir, "/", sp, "/", proj))
    }
  }

  # tabela de valores
  presvals <- raster::extract(predictors, coordinates)

  if (buffer %in% c("mean", "max")) {
    backgr <- createBuffer(coord = coordinates, n.back = n.back, buffer.type = buffer,
                           sp = sp, seed = seed, predictors = predictors)
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
  write.table(sdmdata, file = paste0(models.dir, "/", sp, "/sdmdata.txt"))

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
    ma_cont[ma_cont < dismo::threshold(ema, "no_omission")] <- dismo::threshold(ema, "no_omission")
    ma_bin <- ma_cont > thresholdma
    ma_cut <- ma_cont
    ma_cut[ma_cut < thresholdma] <- thresholdma
    if (raster::minValue(ma_cut) < 0) {
      ma_cut <- (ma_cut - raster::minValue(ma_cut))/raster::maxValue(ma_cut - raster::minValue(ma_cut))
    }
    
    thma$AUC <- ema@auc
    thma$TSS <- ma_TSS
    thma$algoritmo <- "Mahal"
    thma$partition <- i
    row.names(thma) <- paste(sp, i, "Mahal")
    
    write.table(thma, file = paste0(models.dir, "/", sp, "/evaluate", 
      sp, "_", i, "_mahal.txt"))

    if (class(mask) == "SpatialPolygonsDataFrame") {
      ma_cont <- cropModel(ma_cont, mask)
      ma_bin <- cropModel(ma_bin, mask)
      ma_cut <- cropModel(ma_cut, mask)
    }
    raster::writeRaster(x = ma_cont, filename = paste0(models.dir, "/", sp, "/Mahal_cont_", 
      sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = ma_bin, filename = paste0(models.dir, "/", sp, "/Mahal_bin_", 
      sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = ma_cut, filename = paste0(models.dir, "/", sp, "/Mahal_cut_", 
      sp, "_", i, ".tif"), overwrite = T)
  
  
    if (project.model == T) {
      for (proj in projections) {
        data <- list.files(paste0("./env/", proj), pattern = proj)
        data2 <- stack(data)
        ma_proj <- predict(data2, ma, progress = "text")
        ma_proj_bin <- ma_proj > thresholdma
        ma_proj_cut <- ma_proj_bin * ma_proj
        # Normaliza o modelo cut 
        ma_proj_cut <- ma_proj_cut/maxValue(ma_proj_cut)        
        if (class(mask) == "SpatialPolygonsDataFrame") {
          source("./fct/cropModel.R")
          ma_proj <- cropModel(ma_proj, mask)
          ma_proj_bin <- cropModel(ma_proj_bin, mask)
          ma_proj_cut <- cropModel(ma_proj_cut, mask)
        }
          writeRaster(x = ma_proj, filename = paste0(models.dir, "/", sp, "/", 
            proj, "/mahal_cont_", sp, "_", i, ".tif"), overwrite = T)
          writeRaster(x = ma_proj_bin, filename = paste0(models.dir, "/", sp, 
            "/", proj, "/mahal_bin_", sp, "_", i, ".tif"), overwrite = T)
          writeRaster(x = ma_proj_cut, filename = paste0(models.dir, "/", sp, 
            "/", proj, "/mahal_cut_", sp, "_", i, ".tif"), overwrite = T)
          rm(data2)
        }
      }
    } else cat("Mahalanobis distance did not run")
  }
  return(thma)
}
