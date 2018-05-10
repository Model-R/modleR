#' Joins ENM from several partitions, creating a model per algorithm.
#'
#' This function reads the output from dismo.mod and creates a model per species
#' @param species_name A character string with the species name
#' @param algorithms Which algorithms will be processed. If no name is given it
#' will process all algorithms present in the evaluation files
#' @param weight_par Which performance statistic should be used to weight the
#'  partitions. Defaults to NULL but either \code{c("AUC", "TSS")} can be used.
#' @param select_partitions TRUE ou FALSE
#' @param threshold Which selecting threshold will be used to cut the mean
#'                  models in final_model_3 approach (see vignettes), it
#'                  defaults to "spec_sens" but any dismo threshold
#'                  can be used: "kappa", "no_omission", "prevalence",
#'                  "equal_sens_spec", "sensitivity".
#' @param select_par Which performance statistic should be used to select the
#'  partitions- Defaults to NULL but either \code{c("AUC", "TSS")} can be used.
#' @param select_par_val Threshold to select models from TSS values
#' @param which_models Which final_model() will be used? Currently it can be:
#' @param consensus_level Which proportion of models will be kept when creating
#'                   \code{final_model_8} (binary)
#' @param models_dir Character. Folder path where the input files are located
#' @param final_dir Character. Name of the folder to save the output files.
#'                  A subfolder will be created.
#' \describe{
#'   \item{\code{weighted_AUC} or \code{weighted_TSS}}{the models weighted
#'   by TSS or AUC}
#'   \item{\code{final_model_3}}{the binary model created by selecting or not the
#'    partitions, taking their mean and cutting by the mean threshold that
#'    maximizes TSS (or other dismo thresholds)}
#'   \item{\code{final_model_7}}{the mean of the selected binary models}
#'   \item{\code{final_model_8}}{the binary consensus from \code{final_model_7}}
#' }
#' @param write_png Writes png files of the final models
#' @return A set of ecological niche models and figures (optional) written in the
#'          \code{final_dir} subfolder
#' @import raster
#' @importFrom utils read.table write.csv
#' @export
final_model <- function(species_name,
                        algorithms = NULL,
                        weight_par = NULL,
                        select_partitions = TRUE,
                        threshold = c("spec_sens"),
                        select_par = "TSS",
                        select_par_val = 0.7,
                        consensus_level = 0.5,
                        models_dir = "./models",
                        final_dir = "final_models",
                        which_models = c("final_model_3"),
                        write_png = T,
                        ...) {

    if (file.exists(paste0(models_dir, "/", species_name, "/present/",
                           final_dir)) == FALSE)
        dir.create(paste0(models_dir, "/", species_name, "/present/", final_dir),
                   recursive = TRUE)
    print(date())

    cat(paste(species_name, "\n"))
    cat(paste("Reading the evaluation files", "\n"))
    evall <- list.files(
        path = paste0(models_dir, "/", species_name, "/present/partitions"),
        pattern = "evaluate", full.names = T)
    lista <- list()
    for (i in 1:length(evall)) {
        lista[[i]] <- read.table(file = evall[i],
                                 header = T,
                                 row.names = 1)
    }
    stats <- data.table::rbindlist(lista)
    stats <- as.data.frame(stats)
    write.csv(stats, file = paste0(models_dir,"/", species_name, "/present/",
                                   final_dir,"/",species_name,
                                   "final_statistics.csv"))

    # Extracts only for the selected algorithm
    if (is.null(algorithms)) {
        algorithms <- unique(stats$algoritmo)
    }
    algorithms <- as.factor(algorithms)

    for (algo in algorithms) {
        cat(paste("Extracting data for", species_name, algo, "\n"))
        stats.algo <- stats[stats$algoritmo == algo, ]
        #stats.algo <- stats.run[stats.run$algoritmo == algo, ]
        n.part <- nrow(stats.algo)  #How many partitions were there
        #n.part <-  length(unique(stats.algo$partition)) #How many partitions were there
        cat(paste("Reading models from .tif files", "\n"))
        modelos.cont <-
            list.files(
                path = paste0(models_dir, "/", species_name, "/present/partitions"),
                full.names = T,
                #pattern = paste0(algo, "_cont_",species_name,"_",run,"_")
                pattern = paste0(algo, "_cont_",".*tif$")
            )

        modelos.bin <-
            list.files(
                path = paste0(models_dir, "/", species_name, "/present/partitions"),
                full.names = T,
                #pattern = paste0(algo, "_bin_",species_name,"_",run,"_")
                pattern = paste0(algo, "_bin_",".*tif$")
            )
        modelos.cut <-
            list.files(
                path = paste0(models_dir, "/", species_name, "/present/partitions"),
                full.names = T,
                #pattern = paste0(algo, "_bin_",species_name,"_",run,"_")
                pattern = paste0(algo, "_cut_",".*tif$")
            )

        mod.cont <- raster::stack(modelos.cont)  #(0)
        mod.bin <- raster::stack(modelos.bin)  #(0)
        mod.cut <- raster::stack(modelos.cut)  #(0)
        #names(mod.cont) <- paste0(algo, "_cont_", species_name, "_Run_", run, "_Partition_", 1:n.part)
        #names(mod.bin) <- names(mod.cont)

        if (select_partitions == T) {
            cat(paste("selecting partitions for", species_name, algo, "\n"))
            sel.index <- which(stats.algo[, select_par] >= select_par_val)
        } else {
            #it will use everything
            sel.index <- 1:n.part
            }
            cont.sel.1  <- mod.cont[[sel.index]]  #(1)
            bin.sel.5   <- mod.bin[[sel.index]]  #(5)
            cut.sel     <- mod.cut[[sel.index]]  #(5)
            th.mean <- mean(stats.algo[, names(stats.algo) == threshold][sel.index])

            if (length(sel.index) == 0) {
                cat(paste("No partition was selected for", species_name, algo, "\n"))
                }
            # if length(sel.index) == 1 many of the final models are the same
            # 1 and 2 = continuous
            # 3, 5, 7, 8 = binary because th.mean is TSSth
            # 4, 9 = cut
            #I build the stack by repeating those, for homogeneity
            if (length(sel.index) == 1) {
                #cat(paste(length(sel.index), "partition was selected for",
                 #   species_name, algo, "run",run,"\n"))
                cat(paste(length(sel.index), "partition was selected for",
                    species_name, algo, "\n"))

                # bin.sel #[3] bin.sel #[7]
                final <- raster::stack(cont.sel.1,#2
                                       bin.sel.5, #3
                                       cut.sel, #4
                                       bin.sel.5, #7
                                       bin.sel.5, #8
                                       cut.sel #9
                                       )
                names(final) <- c("final_model_2", "final_model_3","final_model_4", "final_model_7", "final_model_8","final_model_9")
                warning("when only one partition is selected some final models are identical")
            }

            # When the selected models are more than one, refer to the map in the vignette
            if (length(sel.index) > 1) {
                cat(paste(length(sel.index), "partitions were selected for",
                          species_name, algo, "\n"))

                final.cont.mean.2 <- raster::mean(cont.sel.1)  #(2)
                final.bin.mean.3 <- (final.cont.mean.2 > th.mean)  #(3)
                final_model_4 <- final.bin.mean.3 * final.cont.mean.2 #(4)
                final.sel.bin.7 <- raster::mean(bin.sel.5)  #(7)
                final.sel.bin.8 <- final.sel.bin.7 > consensus_level  #(8)
                final_model_9 <- final.sel.bin.8 * final.sel.bin.7
                final <- raster::stack(final.cont.mean.2, final.bin.mean.3,
                                       final_model_4,
                                       final.sel.bin.7, final.sel.bin.8,
                                       final_model_9)
                names(final) <- c("final_model_2","final_model_3",
                                  "final_model_4", "final_model_7",
                                  "final_model_8",
                                  "final_model_9")
                #cat(paste("selected final models for", species_name, algo, "run", run, "DONE", "\n"))
                cat(paste("selected final models for", species_name, algo, "DONE", "\n"))
            }

            if (!is.null(weight_par)) {
                final.w <- stack()
                for (wpar in unique(weight_par)) {
                    pond.stats <- stats.algo[, wpar]
                    if (wpar == "TSS")
                        pond.stats <- (pond.stats + 1) / 2
                    #pond <- mod[[1:part]] * pond.stats
                    cat(paste("Calculating the weighted mean for", species_name, wpar, algo, "\n"))
                    final.w.cont <- raster::weighted.mean(mod.cont, w = pond.stats)
                    names(final.w.cont) <- paste0("final_model_weighted_",wpar)
                    final.w <- addLayer(final.w, final.w.cont)

            #names(final.w)[length(names(final))] <- paste("Weighted",par)
                    cat(paste("weighted final models for", species_name, algo, "DONE", "\n"))

                    }
            }

            if (!exists("final")) {
                final <- stack()
                } else {
                    if (exists("final.w")) {
                        final <- addLayer(final, final.w)
                    }
                # Escribe final
                    if (length(final) != 0) {

                        #pero solo los que sean pedidos en which_model
                        which_final <- final[[which_models]]
                        for (i in 1:dim(which_final)[[3]]) {
                        raster::writeRaster(x = which_final[[i]],
                                            filename = paste0(models_dir, "/",
                                                species_name, "/present/",
                                                final_dir, "/", species_name,
                                                "_", algo, "_",
                                                names(which_final)[i], ".tif"),
                                            overwrite = T,
                                            format = "GTiff")
                        }
                        if (write_png == T) {
                            for (i in 1:dim(which_final)[[3]]) {
                                png(filename = paste0(models_dir, "/", species_name,
                                                      "/present/", final_dir, "/",
                                                      species_name,"_", algo, "_",
                                                      names(which_final)[i], ".png"))
                                plot(which_final[[i]], main = names(which_final)[i])
                                dev.off()
                            }
                        }
                    } else {
                        warning(paste("no models were selected for", species_name, algo))
                    }
            }

    }
    print(paste("DONE", algo, "\n"))
    print(date())

}
