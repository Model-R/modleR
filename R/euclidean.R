euclidean <- function(predictors,
                      occurrences,
                      env_dist = "centroid",
                      #probs,
                      filename = "",
                      ...) {
  if (!env_dist %in% c("centroid", "mindist")) {
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
            if (env_dist == "centroid") {
            dist.vals <- apply(v, 1, FUN = function(x) {
                d <- dist(rbind(centroid.val, x))
                vals <- (-d) + 1
            #cortar
                #q <- quantile(vals, probs = probs, na.rm = TRUE, names = FALSE)
                #vals[vals < q] <- q
                return(vals)
            })
            }
            if (env_dist == "mindist") {
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
            if (env_dist == "centroid") {
                dist.vals <- apply(v, 1, FUN = function(x) {
                    d <- dist(rbind(centroid.val, x))
                    vals <- (-d) + 1
                    #cortar
                    #q <- quantile(vals, probs = probs, na.rm = TRUE, names = FALSE)
                    #vals[vals < q] <- q
                    return(vals)
                })
            }
            if (env_dist == "mindist") {
                dist.vals <- apply(v, 1, FUN = function(x) {
                  #aqui dá um warning porque só está vendo se o primeiro raster
                    if (complete.cases(x) == TRUE) {
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
