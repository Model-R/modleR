#' Prepara os dados para modelagem
#'
#' @param species.name Um nome de espécie
#' @param coordinates Uma tabela com pontos de ocorrência
#' @param buffer Define se será usado buffer e de que tipo ("mean" ou "max")
#' @param seed Para reprodutibilidade
#' @param predictors Objeto do tipo RasterStack com variáveis preditoras
#' @param models.dir Path do diretório onde serão escritos os arquivos de saída
#' @param plot_sdmdata Se o png vai ser criado - defaults to F
#' @param n.back Número de pontos de background
#' @return Um data.frame com metadados da modelagem (TSS, AUC, algoritmo etc.)
#' @export
#'
#'
# tabela de valores
setup_sdmdata <- function(species.name = species.name,
                          coordinates = coordinates,
                          buffer = FALSE,
                          seed = 512,
                          predictors = predictors,
                          models.dir = models.dir,
                          plot_sdmdata = T,
                          n.back = 1000,
                          bootstrap = F,
                          boot_proportion = 0.8,
                          n_boot = 10,
                          crossvalidation = F,
                          partitions = partitions,
                          n_cv = 10) {
    if (file.exists(paste0(models.dir)) == FALSE)
        dir.create(paste0(models.dir))
    if (file.exists(paste0(models.dir, "/", species.name)) == FALSE)
        dir.create(paste0(models.dir, "/", species.name))
    partition.folder <- paste0(models.dir, "/", species.name, "/present", "/partitions")
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

    pres <- cbind(coordinates, presvals)
    back <- cbind(backgr, backvals)
    coord_env_all <- rbind(pres, back)
    sdmdata <- cbind(pa,coord_env_all)
    # Data partition-----
    #Crossvalidation, repetated crossvalidation and jacknife
    if (crossvalidation == TRUE) {
        if (nrow(coordinates) < 11) {
            #forces jacknife
            partitions <- nrow(coordinates)
            n_cv <- 1
        }
        if (n_cv == 1) {
            #Crossvalidation
            set.seed(seed)  #reproducibility
            group <- dismo::kfold(coordinates, partitions)
            set.seed(seed + 1)
            bg.grp <- dismo::kfold(backgr, partitions)
            group.all <- c(group, bg.grp)
        }
        if (n_cv > 1) {
            # Repeated CV
            cv.pres <- replicate(n = n_cv,
                                 dismo::kfold(coordinates, partitions))
            dimnames(cv.pres) <- list(NULL, paste0("cv", 1:n_cv))
            cv.back <- replicate(n = n_cv,
                                 dismo::kfold(backgr, partitions))
            dimnames(cv.back) <- list(NULL, paste0("cv", 1:n_cv))
            cv.matrix <- rbind(cv.pres, cv.back)
        }
    }
    # Bootstrap
    if (bootstrap == TRUE) {
    boot.pres <- replicate(n = n_boot,
                           sample(
                               x = seq_along(1:nrow(coordinates)),
                               size = nrow(coordinates) * boot_proportion,
                               replace = FALSE
                           ))
    boot.back <- replicate(n = n_boot,
                           sample(
                               x = seq_along(1:nrow(backgr)),
                               size = nrow(backgr) * boot_proportion,
                               replace = FALSE
                           ))
    boot_p <- matrix(data = "test",
                     nrow = nrow(coordinates),
                     ncol = n_boot,
                     dimnames = list(NULL,paste0("boot",1:n_boot)))
    boot_a <- matrix(data = "test",
                     nrow = nrow(backgr),
                     ncol = n_boot,
                     dimnames = list(NULL,paste0("boot",1:n_boot)))
    for (i in seq_along(1:n_boot)) {
        boot_p[, i][boot.pres[, i]] <- "train"
        }
    for (i in seq_along(1:n_boot)) {
        boot_a[, i][boot.back[, i]] <- "train"
    }
    boot.matrix <- rbind(boot_p, boot_a)
}

    if (exists("group.all"))   sdmdata <- data.frame(group.all, sdmdata)
    if (exists("cv.matrix"))   sdmdata <- data.frame(cv.matrix, sdmdata)
    if (exists("boot.matrix")) sdmdata <- data.frame(boot.matrix, sdmdata)
    write.table(sdmdata, file = paste0(partition.folder, "/sdmdata.txt"))



    if (plot_sdmdata) {
        #Creates a .png plot of the initial dataset
        cat(paste("Plotting the dataset...",'\n'))
        png(filename = paste0(partition.folder, "/sdmdata_", species.name,".png"))
        par(mfrow = c(1, 1), mar = c(5, 4, 3, 0))
        raster::plot(predictors[[1]], legend = F, col = "grey90", colNA = NA)
        points(back, pch = ".", col = "black")
        points(pres, pch = 3, col = "grey50")
        legend("topleft", pch = c("+","."),
               col = c("grey50", "black"), legend = c("Occs","Back"))
        dev.off()
    }
    rm(coord_env_all)
    rm(pres)
    rm(back)
    gc()
    return(sdmdata)
}
