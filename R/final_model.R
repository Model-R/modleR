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
#' @param scale_model Logical. Whether input models should be scaled between 0
#' and 1
#' @param select_par Which performance statistic should be used to select the
#'  partitions- Defaults to NULL but either \code{c("AUC", "TSS")} can be used.
#' @param select_par_val Threshold to select models from TSS values
#' @param consensus_level Which proportion of models will be kept when creating
#'                   \code{bin_consensus} (binary)
#' @param models_dir Character. Folder path where the input files are located
#' @param final_dir Character. Name of the folder to save the output files.
#'                  A subfolder will be created.
#' @param proj_dir Character. The name of the subfolder with the projection.
#' Defaults to "present" but can be set according to the other projections (i.e.
#' to execute the function in projected models)
#' @param which_models Which final_model() will be used? Currently it can be:
#' \describe{
#'   \item{\code{weighted_AUC} or \code{weighted_TSS}}{the models weighted
#'   by TSS or AUC}
#'   \item{\code{raw_mean}}{the mean of the selected raw models}
#'   \item{\code{bin_mean_th}}{the binary model created by cutting
#'   \code{raw_mean} by the mean of the thresholds that
#'    maximize the selected evaluation metric (e.g. TSS (\code{spec_sens}) or
#'    other dismo thresholds)}
#'    \item{\code{cut_mean_th}}{the cut model created by recovering
#'    \code{raw_mean} values above the mean threshold that
#'    maximizes the selected evaluation metric (e.g. TSS (\code{spec_sens}) or
#'    other dismo thresholds)}
#'   \item{\code{bin_mean}}{the mean of the selected binary models}
#'   \item{\code{bin_consensus}}{the binary consensus from \code{bin_mean}.
#'   \code{consensus_level} must be defined, 0.5 means a majority consensus}
#'   \item{\code{cut_mean}}{the mean of the selected cut models}
#' }
#' @param write_png Writes png files of the final models
#' @return A set of ecological niche models and figures (optional) written in the
#'          \code{final_dir} subfolder
#' @import raster
#' @importFrom utils read.table write.csv read.csv
#' @export
final_model <- function(species_name,
                        algorithms = NULL,
                        weight_par = NULL,
                        select_partitions = TRUE,
                        threshold = c("spec_sens"),
						scale_models = TRUE,
                        select_par = "TSS",
                        select_par_val = 0.7,
                        consensus_level = 0.5,
                        models_dir = "./models",
                        final_dir = "final_models",
						proj_dir = "present",
                        which_models = c("raw_mean"),
                        write_png = T) {

    if (file.exists(paste0(models_dir, "/", species_name, "/", proj_dir, "/",
                           final_dir)) == FALSE)
        dir.create(paste0(models_dir, "/", species_name, "/", proj_dir, "/",
                          final_dir),
                   recursive = TRUE)
    print(date())

    cat(paste(species_name, "\n"))
    cat(paste("Reading the evaluation files for",species_name,"in", proj_dir, "\n"))
    evall <- list.files(
        path = paste0(models_dir, "/", species_name, "/present/partitions"),
        pattern = "^evaluate.+.csv$", full.names = T)
    lista_eval <- lapply(evall, read.csv, header = T)
    stats <- data.table::rbindlist(lista_eval)
    stats <- as.data.frame(stats)
    names(stats)[1] <- "species"
    write.csv(stats, file = paste0(models_dir,"/", species_name, "/present/",
                                   final_dir,"/",species_name,
                                   "_final_statistics.csv"))
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
                path = paste0(models_dir, "/", species_name, "/", proj_dir, "/partitions"),
                full.names = T,
                #pattern = paste0(algo, "_cont_",species_name,"_",run,"_")
                pattern = paste0(algo, "_cont_",".*tif$")
            )
		mod.cont <- raster::stack(modelos.cont)  #(0)

        mod.bin <- mod.cont > stats.algo[,threshold] #(0)
        mod.cut <- mod.cont * mod.bin #(0)


		if (scale_models == T) {
		mod.cont <- rescale.layer(mod.cont)
		mod.cut <- rescale.layer(mod.cut)
		}

		#names(mod.cont) <- paste0(algo, "_cont_", species_name, "_Run_", run, "_Partition_", 1:n.part)
        #names(mod.bin) <- names(mod.cont)
#select partitions----
        if (select_partitions == T) {
            cat(paste("selecting partitions for", species_name, algo, "\n"))
            sel.index <- which(stats.algo[, select_par] >= select_par_val)
        } else {
            #it will use everything
            sel.index <- 1:n.part
            }
            cont.sel.1  <- mod.cont[[sel.index]]  #(1)
            bin.sel.2   <- mod.bin[[sel.index]]  #(2)
            cut.sel.3     <- mod.cut[[sel.index]]  #(3)
            th.mean <- mean(stats.algo[, threshold][sel.index])

            if (length(sel.index) == 0) {
                cat(paste("NO partition selected", species_name, algo,proj_dir, "\n"))
                }
            # if length(sel.index) == 1 the mean models are equal to the originals
            # 1 raw and 4 rawmean = continuous
            # 2 bin and 5 binmean = binary
            # 3 cut and 6 cutmean = cut
            # 7 and 2 and 5 because th.mean is TSSth
            # 8 se parece a 7 pero está cortado por 0.5. no tiene sentido porque es cortar un modelo binario por 0.5

            #I build the stack by repeating those, for homogeneity
            if (length(sel.index) == 1) {
                #cat(paste(length(sel.index), "partition was selected for",
                 #   species_name, algo, "run",run,"\n"))
                message(paste(length(sel.index), "partition was selected for",
                    species_name, algo, proj_dir,"\n"))

                final <- raster::stack(cont.sel.1,#4
                                       bin.sel.2, #5
                                       cut.sel.3, #6
                                       bin.sel.2, #7
                                       bin.sel.2 > consensus_level, #8
                                       cut.sel.3 #9
                                       )
                names(final) <- c("raw_mean",
                                  "bin_mean",
                                  "cut_mean",
                                  "bin_mean_th",
                                  "bin_consensus",
                                  "cut_mean_th")
                warning("when only one partition is selected some final models are identical")
            }

            # When the selected models are more than one, refer to the map in the vignette
            if (length(sel.index) > 1) {
                message(paste(length(sel.index), "partitions were selected for",
                          species_name, algo, "\n"))

                raw_mean_4 <- raster::mean(cont.sel.1)  #(4)
                bin_mean_5 <- raster::mean(bin.sel.2)  #(5)
                cut_mean_6 <- raster::mean(cut.sel.3)  #(6)

                mean_TSS_7 <- (raw_mean_4 > th.mean)  #(7)
                bin_consensus_8 <- (bin_mean_5 > consensus_level)  #(8)
                cut_tss_9 <- mean_TSS_7 * raw_mean_4 #(9)

                final <- raster::stack(raw_mean_4,
                                       bin_mean_5,
                                       cut_mean_6,
                                       mean_TSS_7,
                                       bin_consensus_8,
                                       cut_tss_9)
                names(final) <- c("raw_mean",
                                  "bin_mean",
                                  "cut_mean",
                                  "bin_mean_th",
                                  "bin_consensus",
                                  "cut_mean_th")
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
                                                species_name, "/", proj_dir, "/",
                                                final_dir, "/", species_name,
                                                "_", algo, "_",
                                                names(which_final)[i], ".tif"),
                                            overwrite = T,
                                            format = "GTiff")
                        }
                        if (write_png == T) {
                            for (i in 1:dim(which_final)[[3]]) {
                                png(filename = paste0(models_dir, "/", species_name,
                                                      "/", proj_dir, "/", final_dir, "/",
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
