#' This function calculates the mean distance to the centroid of a distribution
#' @param predictors is an environmental variables stack
#' @param occurrences are the occurrence points
#' @param algo is "minimum" or "centroid"
#' @param ... other parameters in raster::writeRaster()
#' @param filename Optional. The raster that will be created in disk
#' @importFrom stats median
#' @examples
#' centr <- euclidean(example_vars, occurrences = coordenadas[[1]][,c(2,3)])
#' raster::plot(centr)
#' @export

euclidean <- function(predictors,
                      occurrences,
                      algo = "centroid",
                      #probs,
                      filename = '',
                      ...) {
    x.st <- raster::scale(predictors)
    pres.vals <- raster::extract(x.st, occurrences)
    if (!is.null(dim(pres.vals))) {
    centroid.val <- apply(pres.vals, 2, mean, na.rm = TRUE)
} else if (is.vector(pres.vals))
    centroid.val <- median(pres.vals, na.rm = T)
    #ou mensagem de erro?
    #stop("only one raster provided")


    out <- raster(predictors)
    big <- !canProcessInMemory(out, 2)
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

    bs <- blockSize(predictors)
    pb <- pbCreate(bs$n, ...)

    if (todisk) {
        for (i in 1:bs$n) {
            v <- getValues(x.st, row = bs$row[i], nrows = bs$nrows[i])
            if (algo == "centroid") {
            dist.vals <- apply(v, 1, FUN = function(x) {
                d <- dist(rbind(centroid.val, x))
                vals <- (-d) + 1
            #cortar
                #q <- quantile(vals, probs = probs, na.rm = T, names = F)
                #vals[vals < q] <- q
                return(vals)
            })
            }
            if (algo == "mindist") {
            dist.vals <- apply(v, 1, FUN = function(x) {
                mindata <- rbind(x, pres.vals)
                d <- min(as.matrix(dist(mindata))[1,][-1], na.rm = T)
                vals <- (-d) + 1
                # cortar
                #q <- quantile(vals, probs = probs, na.rm = T, names = F)
                #vals[vals < q] <- q
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
                    vals <- (-d) + 1
                    #cortar
                    #q <- quantile(vals, probs = probs, na.rm = T, names = F)
                    #vals[vals < q] <- q
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
                    vals <- (-d) + 1
                    #cortar
                    #q <- quantile(vals, probs = probs, na.rm = T, names = F)
                    #vals[vals < q] <- q
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
