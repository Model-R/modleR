createBuffer <- function(coord, n.back, buffer.type, sp, seed, predictors) {

    if (buffer.type == "mean")
        dist.buf <- mean(sp::spDists(x = coord, longlat = FALSE, segments = TRUE))
    if (buffer.type == "max")
        dist.buf <- max(sp::spDists(x = coord, longlat = FALSE, segments = TRUE))
    
    buffer <- raster::buffer(sp::SpatialPoints(coord), width = dist.buf, dissolve = TRUE)
    
    #Rasterizando o buffer p/ geração dos ptos aleatorios
    r_buffer <- raster::rasterize(buffer, predictors, field = buffer_dist@plotOrder)
    r_buffer <- raster::mask(r_buffer, predictors[[1]])
    
    
    # Gerando pontos aleatorios no buffer
    set.seed(seed + 2)
    backgr <- dismo::randomPoints(mask = r_buffer,
                                  n = n.back,
                                  p = coord, 
                                  excludep = T)
    rm(buffer)
    rm(r_buffer)
    gc()
    return(backgr)
}
