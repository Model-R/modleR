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
#' @param scale_models Logical. Whether input models should be scaled between 0
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
#' @param uncertainty Whether an uncertainty map, measured as range (max-min)
#' should be calculated
#' @param write_png Writes png files of the final models
#' @param ... Other parameters from writeRaster
#' @return A set of ecological niche models and figures (optional) written in
#' the \code{final_dir} subfolder
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
                        uncertainty = F,
                        write_png = T,
                        ...) {
    # Escribe final
    final_path <- paste(models_dir, species_name, proj_dir,
                        final_dir, sep = "/")
    if (file.exists(final_path) == FALSE) {
        dir.create(final_path)
    }
    print(date())
    cat(paste(species_name, "\n"))

    cat(paste("Reading evaluation files for", species_name, "in", proj_dir, "\n"))
    evall <- list.files(
        path = paste0(models_dir, "/", species_name, "/present/partitions"),
        pattern = "^evaluate.+.csv$", full.names = T)
    lista_eval <- lapply(evall, read.csv, header = T)
    stats <- data.table::rbindlist(lista_eval)
    stats <- as.data.frame(stats)
    names(stats)[1] <- "species"

    write.csv(stats, file = paste0(models_dir, "/", species_name, "/present/",
                                   final_dir, "/", species_name,
                                   "_final_statistics.csv"))

    # Extracts only for the selected algorithm
    # if the user doesnt specify, it will take all of them
    if (is.null(algorithms)) {
        algorithms <- unique(stats$algoritmo)
    }
    algorithms <- as.factor(algorithms)

    for (algo in algorithms) {
        final_algo <- raster::stack()
        cat(paste("Extracting data for", species_name, algo, "\n"))
        stats.algo <- stats[stats$algoritmo == algo, ]
        #stats.algo <- stats.run[stats.run$algoritmo == algo, ]
        n.part <- nrow(stats.algo)  #How many partitions were there
        #n.part <-  length(unique(stats.algo$partition)) #How many partitions were there
        cat(paste("Reading models from .tif files", "\n"))
        modelos.cont <-
            list.files(
                path = paste0(models_dir, "/", species_name, "/", proj_dir,
                              "/partitions"),
                full.names = T,
                #pattern = paste0(algo, "_cont_", species_name, "_", run, "_")
                pattern = paste0(algo, "_cont_", ".*tif$")
            )
        mod.cont <- raster::stack(modelos.cont)  #(0)

        #select partitions----
        sel.index <- 1:n.part
        if (select_partitions == T) {
            cat(paste("selecting partitions for", species_name, algo, "\n"))
            sel.index <- which(stats.algo[, select_par] >= select_par_val)
        }
        if (!is.null(weight_par)) {
            pond.stats <- stats.algo[, weight_par][sel.index]
            if ("TSS" %in% weight_par)
                pond.stats <- (pond.stats + 1) / 2
        } else {
            pond.stats <- rep(1, length(sel.index))#either selected or not
        }

        if (length(sel.index) == 0) {
            cat(paste("No partition selected", species_name, algo, proj_dir, "\n"))
        } else if (length(sel.index) != 0) {
            message(paste(length(sel.index), "/", n.part,
                          "partitions will be used for", species_name, algo, "\n"))
            if (length(sel.index) == 1) {
                warning(paste("when only one partition is selected some final models
                          are identical", "\n"))
                cont.sel.1  <- mod.cont[[c(sel.index, sel.index)]]
                pond.stats <- c(pond.stats, pond.stats)#(1)
                }
            if (length(sel.index) > 1) {
                cont.sel.1  <- mod.cont[[sel.index]]  #(1)
                }
            #first column of the map. takes raw means and makes them binary or cut by a single mean threshold
            raw_mean <- raster::weighted.mean(cont.sel.1, w = pond.stats)
            if ("raw_mean" %in% which_models) {
                names(raw_mean) <- "raw_mean"#(4)
                final_algo <- raster::addLayer(final_algo, raw_mean)####layerz#
            }
            if (any(c("raw_mean_th", "raw_mean_cut") %in% which_models)) {
                if (is.numeric(threshold)) {#este threshold se repite na outra coluna, verificar que seja equivalente ¬¬ [ö]
                    th.mean <- threshold
                    } else {
                        th.mean <- mean(stats.algo[, threshold][sel.index])
                        }
                raw_mean_th <- (raw_mean > th.mean)  #(7)
                if ("raw_mean_th" %in% which_models) {
                names(raw_mean_th) <- "raw_mean_th"
                final_algo <- raster::addLayer(final_algo, raw_mean_th)
                }
                if ("raw_mean_cut" %in% which_models) {
                    raw_mean_cut <- raw_mean * raw_mean_th #(9)
                    names(raw_mean_cut) <- "raw_mean_cut"
                    final_algo <- raster::addLayer(final_algo, raw_mean_cut)####layerz#
                    }
                }
             #second column of the figure. creates binary selected
             if (any(c("bin_mean", "cut_mean", "bin_consensus") %in% which_models)) {
                if (is.numeric(threshold)) {#este aqui se repete, linha 145, é equivalente cortar aqui e lá?
                    cont.sel.1_scaled <- rescale_layer(cont.sel.1)
                    mod.sel.bin <- cont.sel.1_scaled > threshold #(0)
                } else {
                    mod.sel.bin <- cont.sel.1 > (stats.algo[, threshold][sel.index]) #(0)
                }
                if (any(c("bin_mean", "bin_consensus") %in% which_models)) {
                    bin_mean <- raster::weighted.mean(mod.sel.bin, w = pond.stats)  #(5)
                    names(bin_mean) <- "bin_mean"
                    final_algo <- raster::addLayer(final_algo, bin_mean)####layerz#
                    if ("bin_consensus" %in% which_models) {
                        if (is.null(consensus_level)) {
                            stop( "consensus_level must be specified")
                        }
                        bin_consensus <- (bin_mean > consensus_level)  #(8)
                        names(bin_consensus) <- "bin_consensus"
                        final_algo <- raster::addLayer(final_algo, bin_consensus)####layerz#
                    }
                }
                #third column of the figure depends on mod.sel.bin
                 if ("cut_mean" %in% which_models) {
                     mod.cut.sel <- mod.sel.bin * cont.sel.1
                     cut_mean <- raster::weighted.mean(mod.cut.sel, w = pond.stats)  #(6)
                     names(cut_mean) <- "cut_mean"
                     final_algo <- raster::addLayer(final_algo, cut_mean)####layerz#
                 }
             }

            if (scale_models == T) {
             final_algo <- rescale_layer(final_algo)
            }



            #incerteza #ö está criando esta camada duplicada com cada algoritmo
            if (uncertainty == T) {
                raw_inctz <- raster::calc(cont.sel.1,
                                          fun = function(x) {max(x) - min(x)})
                names(raw_inctz) <- "raw_uncertainty"
                final_algo <- raster::addLayer(final_algo, raw_inctz)####layerz#
                }
            #creation ok
                #cat(paste("selected final models for", species_name, algo, "run", run, "DONE", "\n"))
                cat(paste("selected final models for", species_name, algo, "DONE", "\n"))
        }
#################

        if (raster::nlayers(final_algo) != 0) {
            if (uncertainty == T) {
                which_f <- c(which_models, "raw_uncertainty")
                } else {
                    which_f <- which_models
                }
            which_final <- final_algo[[which_f]]

           message(paste("writing models", algo, names(which_final), "\n"))
           if (raster::nlayers(which_final) > 1 ) {
           raster::writeRaster(which_final,
                                filename = paste0(final_path,
                                                  "/", species_name, "_", algo),
                                suffix = "names",
                                bylayer = T,
                                format = "GTiff", ...)
               }
           if (raster::nlayers(which_final) == 1 ) {
           raster::writeRaster(which_final,
                                filename = paste0(final_path,
                                                  "/", species_name, "_", algo,
                                                  "_", names(which_final)),
                                format = "GTiff", ...)
               }

            if (write_png == T) {
                for (i in 1:raster::nlayers(which_final)) {
                    png(filename = paste0(final_path, "/",
                                          species_name, "_", algo, "_",
                                          names(which_final)[i], ".png"))
                    raster::plot(which_final[[i]], main = names(which_final)[i])
                    dev.off()
                }
            }

        }

    } #else {
      #  warning(paste("no models were selected for", species_name, algo, "\n"))
    #}
      #  }
    print(paste("DONE", algo, "\n"))
    print(date())
}
