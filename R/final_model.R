#' Joins ENM from several partitions, creating a model per algorithm.
#'
#' This function reads the output from dismo.mod and creates a model per species
#' @param species.name A character string with the species name
#' @param algorithms Which algorithms will be processed. If no name is given it
#' will process all algorithms present in the evaluation files
#' @param weight.par Which performance statistic should be used to weight the
#'  partitions defaults to NULL but either \code{c("AUC", "TSS")} can be used.
#' @param select.partitions TRUE ou FALSE
#' @param TSS.value Threshold to select models from TSS values
#' @param threshold Which selecting threshold will be used to cut the mean
#'                  models in final_model_3 approach (see vignettes), it
#'                  defaults to "spec_sens" but any dismo threshold
#'                  can be used: "kappa", "no_omission", "prevalence",
#'                  "equal_sens_spec", "sensitivity".
#' @param consensus.level Which proportion of models will be kept when creating
#'                   \code{final_model_8} (binary)
#' @param models.dir Character. Folder path where the input files are located
#' @param final.dir Character. Name of the folder to save the output files.
#'                  A subfolder will be created.
#' @return A set of ecological niche models and figures (optional) written in the
#'          \code{final.dir} subfolder
#' @import raster
#' @importFrom utils read.table
#' @export
final_model <- function(species.name,
                        algorithms = NULL,
                        weight.par = NULL,
                        select.partitions = TRUE,
                        TSS.value = 0.7,
                        threshold = c("spec_sens"),
                        consensus.level = 0.5,
                        models.dir = "./models",
                        final.dir = "final_models") {

    if (file.exists(paste0(models.dir, "/", species.name, "/present/",
                           final.dir)) == FALSE)
        dir.create(paste0(models.dir, "/", species.name, "/present/", final.dir),
                   recursive = TRUE)
    print(date())

    cat(paste(species.name, "\n"))
    cat(paste("Reading the evaluation files", "\n"))
    evall <- list.files(
        path = paste0(models.dir, "/", species.name, "/present/partitions"),
        pattern = "evaluate", full.names = T)
    lista <- list()
    for (i in 1:length(evall)) {
        lista[[i]] <- read.table(file = evall[i],
                                 header = T,
                                 row.names = 1)
    }
    stats <- data.table::rbindlist(lista)
    stats <- as.data.frame(stats)
# ö escribir los evaluate de todos?
    # Extracts only for the selected algorithm
    if (is.null(algorithms)) {
        algorithms <- unique(stats$algoritmo)
    }
    algorithms <- as.factor(algorithms)
    n.runs <-  length(unique(stats$run)) #How many runs were there
    for (run in seq_along(1:n.runs)) {
        stats.run <- stats[stats$run == run, ]

    for (algo in algorithms) {
        cat(paste("Extracting data for", algo, "\n"))
        stats.algo <- stats.run[stats.run$algoritmo == algo, ]
        #part <- nrow(stats.algo)  #How many partitions were there
        n.part <-  length(unique(stats.algo$partition)) #How many partitions were there
        cat(paste("Reading models from .tif files", "\n"))
        modelos.cont <-
            list.files(
                path = paste0(models.dir, "/", species.name, "/present/partitions"),
                full.names = T,
                pattern = paste0(algo, "_cont_",sp,"_",run,"_")
            )

        modelos.bin <-
            list.files(
                path = paste0(models.dir, "/", species.name, "/present/partitions"),
                full.names = T,
                pattern = paste0(algo, "_bin_",sp,"_",run,"_")
            )

        mod.cont <- raster::stack(modelos.cont)  #(0)
        mod.bin <- raster::stack(modelos.bin)  #(0)
        #names(mod.cont) <- paste0(algo, "_cont_", species.name, "_Run_", run, "_Partition_", 1:n.part)
        #names(mod.bin) <- names(mod.cont)

        if (select.partitions == T) {
            cat(paste("selecting partitions for", species.name, algo, "run",run, "\n"))
            sel.index <- which(stats.algo[, "TSS"] >= TSS.value)
        } else {
            #it will use everything
            sel.index <- 1:n.part
            }
            cont.sel.1  <- mod.cont[[sel.index]]  #(1)
            bin.sel.5   <- mod.bin[[sel.index]]  #(5)
            th.mean <- mean(stats.algo[, names(stats.algo) == threshold][sel.index])

            if (length(sel.index) == 0) {
                cat(paste("No partition was selected for", species.name, algo,"run",run, "\n"))
                }
            # if length(sel.index) == 1 many of the final models are the same
            # 1 and 2 = continuous
            # 3, 5, 7, 8 = binary
            # 4, 9 = cut
            #I build the stack by repeating those, for homogeneity
            if (length(sel.index) == 1) {
                cat(paste(length(sel.index), "partition was selected for",
                    species.name, algo, "run",run,"\n"))

                # bin.sel #[3] bin.sel #[7]
                final <- raster::stack(bin.sel.5, bin.sel.5, bin.sel.5)
                names(final) <- c("final_model_3", "final_model_7", "final_model_8")
            }

            # When the selected models are more than one, refer to the map in the vignette
            if (length(sel.index) > 1) {
                cat(paste(length(sel.index), "partitions were selected for",
                          species.name, algo,"run",run, "\n"))

                final.cont.mean.2 <- raster::mean(cont.sel.1)  #(2)
                final.bin.mean.3 <- (final.cont.mean.2 > th.mean)  #(3)
                final.sel.bin.7 <- raster::mean(bin.sel.5)  #(7)
                final.sel.bin.8 <- final.sel.bin.7 > consensus.level  #(8)

                final <- raster::stack(final.bin.mean.3, final.sel.bin.7, final.sel.bin.8)
                names(final) <- c("final_model_3", "final_model_7", "final_model_8")
                cat(paste("selected final models for", species.name, algo, "run", run, "DONE", "\n"))
            }

            if (!is.null(weight.par)) {
                final.w <- stack()
                for (wpar in unique(weight.par)) {
                    pond.stats <- stats.algo[, wpar]
                    if (wpar == "TSS")
                        pond.stats <- (pond.stats + 1) / 2
                    #pond <- mod[[1:part]] * pond.stats
                    cat(paste("Calculating the weighted mean for", species.name, wpar, algo, "run", run))
                    final.w.cont <- raster::weighted.mean(mod.cont, w = pond.stats)
                    names(final.w.cont) <- paste0("final_model_weighted_",wpar)
                    final.w <- addLayer(final.w, final.w.cont)
                    # png(filename = paste0(models.dir, "/", species.name, "/present/",
            #                       final.dir, "/final_", par, "_Wmean_",
            #                       species.name,"_", algo, ".png"))
            # par(mar = c(4, 4, 3, 3), mfrow = c(1, 1))
            # plot(final.w, main = paste(par, "Wmean", species.name, algo))
            # dev.off()
            # ##

            #names(final.w)[length(names(final))] <- paste("Weighted",par)
                    cat("weighted final models for", species.name, algo,"run" ,run,"DONE", "\n")

                    }
    }

            if (exists("final")) {
                if (exists("final.w")) {
            final <- addLayer(final, final.w)
                }
                # plot(final) Escribe final
                raster::rasterOptions(setfileext = F)
                raster::writeRaster(
                    x = final,
                    filename = paste0(models.dir, "/", species.name, "/present/",
                                      final.dir, "/", names(final), species.name,
                                      algo, "_run",run, ".tif"),
                    bylayer = T,
                    overwrite = T,
                    format = "GTiff"
                )

                #        for (i in 1:dim(final)[[3]]) {
                #          png(filename = paste0(models.dir, "/", species.name, "/", final.dir, "/",
                #          names(final)[i], species.name, algo, ".png"))
                #          plot(final[[i]], main = names(final)[i])
                #          dev.off()
                #        }
            }


        }

    }
    print(paste("DONE", algo, "\n"))
    print(date())
}
