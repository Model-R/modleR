#' Calculates environmental euclidean distances to the occurrences of a species
#'
#' This function calculates the euclidean distance from each pixel in a RasterStack of predictor variables to the environmental centroid of a distribution or the minimum distance from each pixel to any occurrence point
#' @inheritParams setup_sdmdata
#' @param algo Character. Either "centroid" or "mindist"
#' @param ... other parameters in raster::writeRaster()
#' @param filename Optional. The raster that will be created in disk
#' @importFrom stats median
#' @examples
#' centr <- euclidean(example_vars, occurrences = coordenadas[[1]][, c(2,3)])
#' raster::plot(centr)
#' @export

euclidean <- function(predictors,
                      occurrences,
                      algo = "centroid",
                      #probs,
                      filename = '',
                      ...) {
  if (!algo %in% c("centroid", "mindist")) {
    stop('Algorithm must be either "centroid" or "mindist"')
  }
    x.st <- raster::scale(predictors)
    pres.vals <- raster::extract(x.st, occurrences)
    if (!is.null(dim(pres.vals))) {
    centroid.val <- apply(pres.vals, 2, median, na.rm = TRUE)
    } else if (is.vector(pres.vals))
    centroid.val <- median(pres.vals, na.rm = TRUE)
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
                #q <- quantile(vals, probs = probs, na.rm = TRUE, names = FALSE)
                #vals[vals < q] <- q
                return(vals)
            })
            }
            if (algo == "mindist") {
            dist.vals <- apply(v, 1, FUN = function(x) {
                mindata <- rbind(x, pres.vals)
                d <- min(as.matrix(dist(mindata))[1,][-1], na.rm = TRUE)
                vals <- (-d) + 1
                # cortar
                #q <- quantile(vals, probs = probs, na.rm = TRUE, names = FALSE)
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
                    #q <- quantile(vals, probs = probs, na.rm = TRUE, names = FALSE)
                    #vals[vals < q] <- q
                    return(vals)
                })
            }
            if (algo == "mindist") {
                dist.vals <- apply(v, 1, FUN = function(x) {
                    if (complete.cases(x) == TRUE) { #aqui dá um warning porque só está vendo se o primeiro raster do stack
                        mindata <- rbind(x, pres.vals)
                        d <- min(as.matrix(dist(mindata))[1,][-1], na.rm = TRUE)
                        } else {
                        d <- NA
                    }
                    vals <- (-d) + 1
                    #cortar
                    #q <- quantile(vals, probs = probs, na.rm = TRUE, names = FALSE)
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
    raster::pbClose(pb)
    return(out)
}
