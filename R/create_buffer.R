create_buffer <- function(coord,
                          n.back,
                          buffer.type,
                          seed,
                          predictors) {

# Transformando em spatial points
sp::coordinates(coord) <- ~lon + lat
raster::crs(coord) <- raster::crs(predictors)
    if (buffer.type == "mean")
        dist.buf <- mean(sp::spDists(x = coord, segments = FALSE))
    if (buffer.type == "max")
        dist.buf <- max(sp::spDists(x = coord, segments = FALSE))
    if (buffer.type == "median")
        dist.buf <- stats::median(sp::spDists(x = coord, segments = FALSE))

    #cria o buffer - é um shape
    buffer.shape <- raster::buffer(coord,
                                   width = dist.buf * 100,
                                   dissolve = TRUE)

    #Rasterizando o buffer p/ geração dos ptos aleatorios
    r_buffer <- raster::rasterize(buffer.shape, predictors,
                                  field = buffer.shape@plotOrder)
    #o mask é o mínimo necessário para que o buffer fique sem NAs nos preditores
    r_buffer <- raster::mask(r_buffer, predictors[[1]])


    # Gerando pontos aleatorios no buffer
    set.seed(seed + 2)
    backgr <- dismo::randomPoints(mask = r_buffer,
                                  n = n.back,
                                  p = coord,
                                  excludep = T)
    rm(buffer.shape)
    rm(r_buffer)
    gc()
    return(backgr)
}
