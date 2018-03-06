#' Fits ecological niche models using GLM.
#'
#' @inheritParams do_bioclim
#' @return A data frame with the evaluation statistics (TSS, AUC, etc.)
#' @importFrom stats formula
#' @importFrom stats step
#' @export
do_GLM <- function(sp,
                   coordinates,
                   partitions,
                   buffer = FALSE,
                   seed = 512,
                   predictors,
                   models.dir = "./models",
                   project.model = FALSE,
                   projections = NULL,
                   mask,
                   write_png = FALSE,
                   n.back) {
    cat(paste("GLM", "\n"))

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

    null.model <- glm(sdmdata_train$pa ~ 1, data = envtrain,
                      family = "binomial")
    full.model <- glm(sdmdata_train$pa ~ ., data = envtrain,
                      family = "binomial")
    glm <- step(object = null.model, scope = formula(full.model),
                direction = "both", trace = F)
    eglm <- dismo::evaluate(envtest_pre, envtest_back, model = glm,
                            type = "response")  #####
    # eglm <- evaluate(pres_test, backg_test, glm, predictors, type="response")
    glm_TSS <- max(eglm@TPR + eglm@TNR) - 1
    thresholdglm <- eglm@t[which.max(eglm@TPR + eglm@TNR)]
    thglm <- dismo::threshold(eglm)
    thglm$AUC <- eglm@auc
    thglm$TSS <- glm_TSS
    thglm$algoritmo <- "glm"
    thglm$partition <- i
    row.names(thglm) <- paste(sp, i, "glm")

    glm_cont <- dismo::predict(predictors, glm, progress = "text",
                               type = "response")
    glm_bin <- glm_cont > thresholdglm
    glm_cut <- glm_bin * glm_cont

    write.table(thglm, file = paste0(partition.folder, "/evaluate",
      sp, "_", i, "_glm.txt"))

    if (class(mask) == "SpatialPolygonsDataFrame") {
        glm_cont <- crop_model(glm_cont, mask)
        glm_bin <- crop_model(glm_bin, mask)
        glm_cut <- crop_model(glm_cut, mask)
    }
    raster::writeRaster(x = glm_cont,
                        filename = paste0(partition.folder, "/glm_cont_",
      sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = glm_bin,
                        filename = paste0(partition.folder, "/glm_bin_", sp,
      "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = glm_cut,
                        filename = paste0(partition.folder, "/glm_cut_", sp,
      "_", i, ".tif"), overwrite = T)

    if (write_png == T) {
        png(filename = paste0(partition.folder, "/glm", sp, "_", i, "%03d.png"))
        plot(glm_cont, main = paste("GLM raw", "\n",
                                    "AUC =", round(eglm@auc, 2), "-",
                                    "TSS =", round(glm_TSS, 2)))
        plot(glm_bin, main = paste("GLM P/A", "\n",
                                   "AUC =", round(eglm@auc, 2), "-",
                                   "TSS =", round(glm_TSS, 2)))
        plot(glm_cut, main = paste("GLM cut", "\n",
                                   "AUC =", round(eglm@auc, 2), "-",
                                   "TSS =", round(glm_TSS, 2)))
        dev.off()
        }

    if (project.model == T) {
      for (proj in projections) {
      projection.folder <- paste0(models.dir, "/", sp, "/", proj)
            if (file.exists(projection.folder) == FALSE)
                dir.create(paste0(projection.folder), recursive = T)

        data <- list.files(paste0("./env/", proj), pattern = proj)
        data2 <- stack(data)
        glm_proj <- predict(data2, glm, progress = "text")
        glm_proj_bin <- glm_proj > thresholdglm
        glm_proj_cut <- glm_proj_bin * glm_proj
        # Normaliza o modelo cut
        glm_proj_cut <- glm_proj_cut / maxValue(glm_proj_cut)
        if (class(mask) == "SpatialPolygonsDataFrame") {
          glm_proj <- crop_model(glm_proj, mask)
          glm_proj_bin <- crop_model(glm_proj_bin, mask)
          glm_proj_cut <- crop_model(glm_proj_cut, mask)
        }
        writeRaster(x = glm_proj,
                    filename = paste0(projection.folder, "/glm_cont_", sp,
                                      "_", i, ".tif"), overwrite = T)
        writeRaster(x = glm_proj_bin, filename = paste0(projection.folder,
                                                        "/glm_bin_", sp, "_", i,
                                                        ".tif"), overwrite = T)
        writeRaster(x = glm_proj_cut, filename = paste0(projection.folder,
                                                        "/glm_cut_", sp, "_", i,
                                                        ".tif"), overwrite = T)
        rm(data2)
      }
    }
  }
  return(thglm)
}
