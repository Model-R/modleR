#this function calculates the mean distance to the centroid of a distribution
#x is a stack
# occs are the occurrence points
# centroid.val are the values at the centroid of the distribution

euclidean <- function(x, occs, algo, filename = '', ...) {
    x.st <- raster::scale(x)
    pres.vals <- raster::extract(x.st, occs)
    centroid.val <- apply(pres.vals, 2, mean, na.rm = TRUE)

    out <- raster(x)
    big <- !canProcessInMemory(out, 3)
    filename <- trim(filename)
    if (big & filename ==  '') {
        filename <- rasterTmpFile()
    }
    if (filename != '') {
        out <- writeStart(out, filename, ...)
        todisk <- TRUE
    } else {
        vv <- matrix(ncol = nrow(out), nrow = ncol(out))
        todisk <- FALSE
    }

    bs <- blockSize(x)
    pb <- pbCreate(bs$n, ...)

    if (todisk) {
        for (i in 1:bs$n) {
            v <- getValues(x.st, row = bs$row[i], nrows = bs$nrows[i])
            if (algo == "centroid") {
            dist.vals <- apply(v, 1, FUN = function(x) {
                d <- dist(rbind(centroid.val, x))
                vals <- -d
            return(vals)
            })
            }
            if (algo == "mindist") {
            dist.vals <- apply(v, 1, FUN = function(x) {
                mindata <- rbind(x, pres.vals)
                d <- min(as.matrix(dist(mindata))[1,][-1], na.rm = T)
                vals <- -d
            return(vals)
            })
            }
            out <- writeValues(out, dist.vals, bs$row[i])
            pbStep(pb, i)
        }
        out <- writeStop(out)
    } else {
        for (i in 1:bs$n) {
            v <- getValues(x.st, row = bs$row[i], nrows = bs$nrows[i])
            if (algo == "centroid") {
                dist.vals <- apply(v, 1, FUN = function(x) {
                    d <- dist(rbind(centroid.val, x))
                    vals <- -d
                    return(vals)
                })
            }
            if (algo == "mindist") {
                dist.vals <- apply(v, 1, FUN = function(x) {
                    if (complete.cases(x) == T) { #aqui dà um warning porque só está vendo se o primeiro raster do stack
                        mindata <- rbind(x, pres.vals)
                        d <- min(as.matrix(dist(mindata))[1,][-1], na.rm = T)
                        } else {
                        d <- NA
                    }
                    vals <- -d
                    return(vals)
                })
            }
            cols <- bs$row[i]:(bs$row[i] + bs$nrows[i] - 1)
            vv[, cols] <- matrix(dist.vals, nrow = out@ncols)
            pbStep(pb, i)
        }
        out <- setValues(out, as.vector(vv))
    }
    pbClose(pb)
    return(out)
}
