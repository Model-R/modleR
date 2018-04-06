#' Prepara os dados para modelagem
#'
#' @param sp Um nome de espécie
#' @param coordinates Uma tabela com pontos de ocorrência
#' @param buffer Define se será usado buffer e de que tipo ("mean" ou "max")
#' @param seed Para reprodutibilidade
#' @param predictors Objeto do tipo RasterStack com variáveis preditoras
#' @param models.dir Path do diretório onde serão escritos os arquivos de saída
#' @param write_png Se o png vai ser criado - defaults to F
#' @param n.back Número de pontos de background
#' @return Um data.frame com metadados da modelagem (TSS, AUC, algoritmo etc.)
#' @export
#'
#'
# tabela de valores
setup_sdmdata <- function(sp = sp,
                          coordinates = coordinates,
                          partitions = partitions,
                          buffer = FALSE,
                          seed = 512,
                          predictors = predictors,
                          models.dir = models.dir,
                          plot_sdmdata = T,
                          n.back = 500) {
    if (buffer %in% c("mean", "max", "median")) {
        backgr <- create_buffer(coord = coordinates,
         n.back = n.back,
                               buffer.type = buffer, seed = seed,
                               predictors = predictors)
        } else {
            set.seed(seed + 2)
            backgr <- dismo::randomPoints(mask = predictors,
                                          n = n.back,
                                          p = coordinates,
                                          excludep = T)
        }

    colnames(backgr) <- c("lon", "lat")

    # Extraindo dados ambientais
    presvals <- raster::extract(predictors, coordinates)
    backvals <- raster::extract(predictors, backgr)

    pa <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals)))

    # Data partition
#Jacknife
if (crossvalidation == TRUE) {
    if (nrow(coordinates) < 11)
        partitions <- nrow(coordinates)
    #Crossvalidation
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
    return(sdmdata)
}
