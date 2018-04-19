# ModelR: a workflow for ecological niche models based on dismo"

ModelR is a workflow based on package dismo, designed to automatize some of the common steps when performing ecological niche models. Given the occurrence records and a set of environmental predictors, it prepares the data to crossvalidation or jacknife procedures depending on the number of occurrence points, then it performs ecological niche models using the several algorithms implemented in the `dismo` package.

## Installing 

```
devtools::install_github("Model-R/modelr_pkg", ref = "sdmdata", build_vignettes = TRUE)
```

(`build_vignettes` will include this vignette on the installation, `ref = sdmdata` refers to this branch)

## The workflow

The workflow consists of mainly three functions that should be used sequentially.

1. `do_enm()` and each independent `do_[algorithm]` function make the ENM for each partition of the crossvalidation
2. `final_model()` selects and joins the partition models into a final model per species per algorithm, 
3. `ensemble_model()` joins the final models per algorithm into an ensemble model.

# Folder structure created by this package

`modelr` writes the outputs in the hard disk, according to the following folder structure:   

    `models.dir/projection/partitions`  
    `models.dir/projection/final_models`  
    `models.dir/projection/ensemble_models`  

+ We define a partition as the individual modeling round that takes part of the data to train the algorithms and the other parts to test them. Currently we follow k-fold crossvalidations but _bootstrap is being implemented_.
+ We define the final models as joining together the partitions and obtaining __one model per species per algorithm__.
+ "Ensemble" models join together the results obtained by different algorithms.
+ When projecting models into the present, the projection folder is called `present`. _The projection unto other areas and/or climate scenarios is being implemented_.
+ You can set `models.dir` wherever you want in the hard disk, but if you do not modify the default value, it will create the output under the working directory (its default value is `./models.dir`, where the period points to the working directory)
+ However, the nested subfolder structure will remain the same, and only the _names_ of the final and ensemble folders can be modified. If you change `final_models` default value ("final_model") you will need to include the new value when calling ensemble() (`final_dir = "[new name]"`), to indicate the function where to look for models. This flexibility allows for experimenting with final model and ensemble construction. 


# Fitting a model per partition

All `do_[algorithm]` functions and `do_enm()` create a *model per partition, per algorithm*. The available algorithms are `bioclim`, `domain`, `maxent`, mahalanobis distances (`mahal`), as implemented in `dismo`, and support vector machines (SVM), as implemented by packages `kernlab` and `e1071`. GLM and Random Forest (from package `randomForest`) are also implemented.


## Setting up the data

`modelr` comes with example data, a data frame called `coordenadas`, with occurrence data for four species, and predictor variables called `variaveis_preditoras`


```{r lib, echo = T, eval = T}
library(modelr)
library(rJava) 
library(raster)
head(coordenadas)
species <- unique(coordenadas$sp)
species
```



```{r dataset, fig.width= 5, fig.height=5, fig.cap= "Figure 1. The example dataset: predictor variables and occurrence for four species."}
raster::plot(variaveis_preditoras[[1]])
points(sp::SpatialPoints(coordenadas[,c(2,3)]), 
       bg = as.numeric(unclass(coordenadas$sp)), pch = 21)
```


We'll filter the `coordenadas` file to select only the data for the first species: 

```{r occs, message = F}
library(dplyr)
species[1]
occs <- filter(coordenadas, sp == species[1]) %>% select(lon, lat)
head(occs)
```

## Fitting one algorithm: `do_[algorithm]`


```{r do_bioclim}
args(do_bioclim)
```
The arguments in each individual `do_[algorithm]` are:

+ `partitions`: It implements a k-fold cross-validation (argument `part`, defaults to 3) but overwrites part when n < 10, setting part to the number of occurrence records (a jacknife partition).  
+ `buffer`: can build a distance buffer around the occurrence points, by taking either the maximal, median or mean distance between points. Pseudoabsence points will be sampled (using `dismo::randomPoints()`) within this buffer.
+ `seed`: for reproducilibity purposes 
+ `mask`: will crop and mask the partition models into a ShapeFile


```{r, eval = F}
do_bioclim(sp = species[1],
           coordinates = occs,
           partitions = 3,
           buffer = "mean",
           predictors = variaveis_preditoras,
           mask = NULL,
           models.dir = "~/temp",
           write_png = T,
           n.back = 500)
```



You can explore the list of files created at this phase, for example:

```{r partfiles, echo=FALSE, eval=FALSE}
partitions.folder <- 
    list.files("~/temp", recursive = T, pattern = "partitions", include.dirs = T, full.names = T)
head(list.files(partitions.folder, recursive = T))
```
 
At the end of a modeling round, the partition folder containts: 

+ A file called `sdmdata.txt` with the data used for each partition. 
+ A .tif file for each partition, continuous, binary and cut by the threshold that maximizes its TSS.
+ Figures in .png to explore the results readily, without reloading them into R or opening them in a SIG program. The creation of these .png can be controlled with the `write_png` parameter. 
+ A .txt table with the evaluation data for each partition: `evaluate_[Species name ]_[partition number]_[algorithm].txt`. These files will be read by the `final_model()` function, to generate the final model per species.


### Fitting several algorithms per species: `do_enm()`

The same modeling procedure can be performed by using `do_enm()`, that receives the same parameters as `do_bioclim()` but allows the user to call the algorithms using TRUE or FALSE statements (just as BIOMOD2 functions do). 

```{r do_enm1, eval = F}
args(do_enm)
do_enm(sp = species[1],
       coordinates = occs,
       partitions = 3,
       buffer = "mean",
       predictors = variaveis_preditoras,
       mask = NULL,
       models.dir = "~/temp",
       write_png = T,
       n.back = 500,
       bioclim = T)
```

`do_enm()` calls several independent functions for each algorithm, called `do_[algorithm]`. This allows to parallelize by species and/or species algorithms in multi-cluster environments. The following lines call for bioclim, GLM, maxent, random forests and smv.

```{r do_enm2, eval = F, cache = F}
args(do_enm)
do_enm(sp = species[1],
       coordinates = occs,
       partitions = 3,
       buffer = "mean",
       predictors = variaveis_preditoras,
       mask = NULL,
       models.dir = "~/temp",
       write_png = T,
       n.back = 500,
       bioclim = T, 
       glm = T,
       maxent = T,
       rf = T,
       svm = T)
```


## Joining partitions: `final_model()`

There are many ways to create a final model per algorithm per species. `final_model()` follows the following logic:

<!--![](final_model_english.png)-->


+ It can weigh the partitions by a performance metric `weigh.partitions = TRUE` and `weight.par = "spec_sens"`, and give larger weights to partitions with better performance. This results in a continuous, uncut surface. 
+ It can select the best partitions if the parameter `select.partitions = TRUE`, selecting only those who obtained a TSS value above `TSS.value` (TSS varies between -1 and 1, defaults to 0.7). If `select.partitions` is set to FALSE, it will create the final model using all partitions. 
+ The final models can be done using a subset of the algorithms avaliable on the hard disk, using the parameter `algorithms`. If left unspecified, all algorithms listed in the `evaluate` files will be used.
+ The selected models/algorithms form a `raster::rasterStack()` object. Their mean can be calculated (step 2) and the binary model can be obtained by cutting by the mean threshold (meanTSSth) that maximizes the individual partition's TSS (step 3)
+ The selected binary models (step 5) can also be joined by a mean (step 7) and a binary (step 8) or cut (step 9) model can be obtained through levels of consensus (defaults to 0.5: majority consensus approach).


```{r final_model}
args(final_model)
```


```{r, cache = T, eval = F}
final_model(species.name = species[1],
            select.partitions = T,
            algorithms = c("bioclim", "maxent"),
            TSS.value = 0.5,
            models.dir = "~/temp")
```

`final_model()` creates a .tif file for each final.model (one per algorithm) under the specified folder (default: `final_models`)
 
We can explore these models from the files:

```{r final_folder, eval = F}
final.folder <- list.files("~/temp", 
                           recursive = T,
                           pattern = "final_models",
                           include.dirs = T,
                           full.names = T)
final_mods <- list.files(final.folder, full.names = T)
head(final_mods)
```

```{r plot_final, fig.width = 7, fig.height = 6, eval = F}
library(raster)
final_models <- stack(final_mods)
plot(final_models)
```

## ensemble()

The third step of the workflow is joining the models for each algorithm into a final ensemble model. `ensemble()` calculates the mean of the final models and saves them under the folder specified by `ensemble_dir`. It can also cut this mean using a consensus rule (what proportion of final models predict a presence in each pixel, 0.5 is a majority rule).

`ensemble()` uses a `which.model` parameter to specify which final model should be assembled together (the default is
`which.models = c("final_model_3", "final_model_7", "final_model_8")`) referring to step 3, 7 and 8 of the `final_model()` approach.

```{r ensemble, cache = T, eval = F}
ensemble(species[1],
         occs = occs,
         which.models = "final_model_7",
         models.dir = "~/temp")
```

At any point we can explore the outputs in the folders: 

```{r check_ensemble, fig.width = 5, fig.height = 5,eval = F}
ensemble_files <-  list.files("~/temp/", recursive = T,
                              pattern = "final_model_7_ensemble.+tif",
                               full.names = T, )

ensemble_files
ens.mod <- stack(ensemble_files)
plot(ens.mod)
plot(ens.mod[[2]])
maps::map(,,add = T)
```


# Workflows with multiple species

Our `coordenadas` dataset has data for four species. 
An option to do the several models is to use a `for` loop

```{r, eval = F}
especies <- unique(coordenadas$sp)
for (especie in especies) {
    occs <- coordenadas[coordenadas$sp == especie, c("lon", "lat")]
    do_enm(sp = especie,
           coordinates = occs,
           partitions = 3,
           buffer = "mean",
           predictors = variaveis_preditoras,
           models.dir = "~/tempforloop",
           n.back = 500,
           write_png = T,
           bioclim = T,
           maxent = T,
           rf = T,
           svm = T)
}

```

Another option is to use the `purrr` package:

```{r purrr example, eval = F}
library(purrr)
coordenadas %>% split(.$sp) %>%
    purrr::map(~ do_enm(sp = unique(.$sp), 
                        coordinates = .[, c("lon", "lat")],
                        partitions = 3,
                        buffer = "mean",
                        predictors = variaveis_preditoras,
                        models.dir = "~/temp_purrr",
                        n.back = 500,
                        write_png = T,
                        bioclim = T,
                        maxent = T,
                        rf = T,
                        svm = T))
```

```{r purrr_final, eval = F}
coordenadas %>% 
    split(.$sp) %>%
    purrr::map(~ final_model(species.name = unique(.$sp),
                             select.partitions = TRUE,
                             TSS.value = 0.5,
                             consensus.level = 0.5,
                             models.dir = "~/temp_purrr"))
```

```{r purrr_ensemble, eval = F}
coordenadas %>% 
    split(.$sp) %>%
    purrr::map(~ ensemble(species.name = unique(.$sp),
                          occs = .[, c("lon", "lat")],
                          which.models = "final_model_7",
                          write_png = T,
                          models.dir = "~/temp_purrr"))

```

# NEWS
## 2018-04-20
+ sdmdata ficou separado da geração do modelo. agora quando cada algoritmo vai começar simplesmente procura o sdmdata.txt, se já tem ele não o gera de novo. 
+ sdmdata pode criar uma tabela com o desneho experimental para bootstrap (n repetições de um sampling proporcional), k-fold crosvalidação (separa o data set em k grupos) e crosvalidação repetida (faz cv varias vezes)
+ como agora tem runs rodadas, os arquivos são escritos acorde com isto

+ Implementei uma única função de modelagem, `do_any`com um parâmetro `algo` que seleciona entre as opções, mas chama `sdmdata` uma vez só. 

        setupsdmdata
        #seleção
        if algo == (...) {
          bc <- dismo::bioclim(predictors, pres_train)
          bc <- dismo::maxent(predictors, pres_train)
          bc <- dismo::mahal(predictors, pres_train)
          bc <- dismo::domain(predictors, pres_train)
          bc <- dismo::domain(predictors, pres_train)
        }
        
        avaliação, threshold etc...
        returns(th_table)

Isto requer de um tratamento especial para glm, rf e svm pois eles não vêm de um chamado direto a `dismo`, mas a avaliação com `evaluate()` é igual
