## This script loops through sets of hyper-parameters that specify a
## boosted regression tree algorithm. Specifically, this script passes
## a set of parameters and data-set to the train.pathogen.learners
## function, which, in turn, performs the boosting algorithm and
## generates output.

## Require packages used for gis, parallel,
## boosted tree modeling, plotting, and
## calculating the AUC statistic
require(raster)
require(dismo)
require(gbm)
require(parallel)
require(sf)
require(fields)
require(verification)
require(ggplot2)
require(plyr)
require(rgdal)
require(viridis)

## Directory name that will contain output for all hyperparameter sets
prefix <- 'pathogen_v2'

## set version of random number generator to use to ensure
## reproducibility across R versions
RNGversion('3.6.0')

## Storage locations
storage.fold <- '../Storage'
raster.data.fold <- paste(storage.fold,'/Raster_Data', sep = '')

## Load in predictor stack. Used for predicting across West Africa
all.stack <- stack(paste(raster.data.fold, "/predictor_stack.grd", sep = ''))

## Load in shape files for plotting
foc.shp.ogr <- readOGR(dsn = paste(storage.fold, '/Shapefiles/West_Africa', sep = ''), layer = 'foc', verbose = FALSE)
masto.rangemap <-readOGR(dsn = paste(storage.fold, '/Shapefiles/Masto_Range', sep = ''), layer = 'data_0', verbose = FALSE)

## Additional functions used by Train_Pathogen_Learners
source("Tools/Functions.r")
source('Prep_Pathogen_Data.r', local = TRUE)
source("Calc_Sig_Pathogen_Preds.r")

## Load function that trains learners with a given set of hyperparameters
source("Train_Pathogen_Learners.r")

## Set hyperparameters.
## nboots - number of times the data is sampled and the model fit
## tree.complexity - depth of trees that are fit at each model iteration
## mllr - minus log10 of the learning rate
## lmt - log10 of the max.trees parameter
grid <- expand.grid(nboots = 25,
                    tree.complexity = c(1),
                    mllr = 3,
                    lmt = 7
                    )

cat(paste0('\n\n\n\n' , '---- ', prefix, ' ----', '\n\n\n\n'))

starttime <- Sys.time()

for(gi in 1:nrow(grid)){

    ## Prep the Mastomys natalensis LASV occurrence dataset
    rodlsv.survey.dat <- Prep.Pathogen.Data()

    ## Predictors that are deemed significant by the Wilcox test
    var.names <- Calc.Sig.Pathogen.Preds(rodlsv.survey.dat)

    ## Train learners with the given set of hyperparameters
    mod.out <- train.pathogen.learners(rodlsv.survey.dat, grid[gi,])

    ## Extract AUC statistics
    fold <- generate.lassa.name(grid[gi,])
    tree.dat <- read.table(paste('Figures_Fits/', prefix, '/', fold,'/tree.dat',sep=''),
                           header = TRUE)
    grid[gi,'auc'] <- mean(unlist(tree.dat[,'model.auc']), na.rm = TRUE)
    grid[gi,'sd.auc'] <- sd(unlist(tree.dat[,'model.auc']), na.rm = TRUE)

    ## Write AUC statistics to file
    write.table(grid, file = paste('Figures_Fits/', prefix,'/grid_lsv_data', sep = ''),
                row.names = FALSE)

    print(paste('--- gi:', gi, ' ---', sep = ''), quote = FALSE)
}

print(paste0('Run time: ', Sys.time() - starttime))
