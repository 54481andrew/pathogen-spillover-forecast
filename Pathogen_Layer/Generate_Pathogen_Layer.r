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
options("rgdal_show_exportToProj4_warnings"="none") ## Turn off gdal warnings
require(rgdal)
require(viridis)

## Directory name that will contain output for all hyperparameter sets
prefix <- 'pathogen_v3'

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
hypers.dat <- expand.grid(nboots = 25,
                          set.ambiguous.to = c(1,NA),
                          min.test = 5,
                          tree.complexity = c(1),
                          mllr = 3,
                          lmt = 7
                          )

cat(paste0('\n\n\n\n' , '---- ', prefix, ' ----', '\n\n\n\n'))

starttime <- Sys.time()

## Create [prefix] directory
dirpath <- paste('Figures_Fits/',prefix,sep='')
if(!dir.exists(dirpath)){dir.create(dirpath, showWarnings = FALSE)}

for(ii in 1:nrow(hypers.dat)){

    ## Set up [prefix] subdirectories
    fold <- generate.lassa.name(hypers.dat[ii,])
    dirpath <- paste('Figures_Fits/', prefix, '/', fold,sep='')
    models.folder = paste0('Figures_Fits/', prefix, '/', fold, '/Models')
    if(!dir.exists(dirpath)){
        dir.create(dirpath, showWarnings = FALSE)
    }else{
        unlink(models.folder, recursive = TRUE)
    }
    dir.create(models.folder,showWarnings = FALSE)
    
    ## Prep the Mastomys natalensis LASV occurrence dataset
    rodlsv.survey.dat <- Prep.Pathogen.Data(hypers.dat[ii,])[[1]]

    ## Maps the prepared dataset
    source("Map_Cases.r")
    
    ## Predictors that are deemed significant by the Wilcox test
    var.names <- Calc.Sig.Pathogen.Preds(rodlsv.survey.dat)

    ## Train learners with the given set of hyperparameters
    mod.out <- train.pathogen.learners(rodlsv.survey.dat, hypers.dat[ii,])

    ## Extract AUC statistics
    tree.dat <- read.table(paste('Figures_Fits/', prefix, '/', fold,'/tree.dat',sep=''),
                           header = TRUE)
    hypers.dat[ii,'auc'] <- mean(unlist(tree.dat[,'model.auc']), na.rm = TRUE)
    hypers.dat[ii,'sd.auc'] <- sd(unlist(tree.dat[,'model.auc']), na.rm = TRUE)

    ## Write AUC statistics to file
    write.table(hypers.dat, file = paste('Figures_Fits/', prefix,'/hypers_lsv_data', sep = ''),
                row.names = FALSE)

    print(paste('--- ii:', ii, ' ---', sep = ''), quote = FALSE)
}

print(paste0('Run time: ', Sys.time() - starttime))
