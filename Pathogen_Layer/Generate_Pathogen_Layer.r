## This script loops through sets of hyper-parameters that specify a
## boosted regression tree algorithm and fits the corresponding
## model. Specifically, this script passes a set of parameters and
## data-set to the train.pathogen.learners function, which, in turn,
## performs the boosting algorithm and generates output.

require(raster)
require(dismo)
require(gbm)
require(parallel)
require(sf)
require(fields)
require(verification)
require(ggplot2)
require(ggthemes)
require(plyr)
require(rgdal)
require(viridis)
require(gridExtra)

## Directory name that will contain output for all hyperparameter sets
prefix <- 'pathogen_v7'

## Set version of random number generator to use to ensure
## reproducibility across R versions
RNGversion('3.6.0')

## Storage locations
storage.fold <- '../Storage'
raster.data.fold <- paste(storage.fold,'/Raster_Data', sep = '')

## Load in predictor stack. Used for predicting across West Africa
all.stack <- stack(paste(raster.data.fold, "/predictor_stack.grd", sep = ''))

## Load in shape files for plotting West Africa and the
## Mastomys natalensis homerange
foc.shp.ogr <- readOGR(dsn = paste(storage.fold, '/Shapefiles/West_Africa', sep = ''), layer = 'foc', verbose = FALSE)
masto.rangemap <-readOGR(dsn = paste(storage.fold, '/Shapefiles/Masto_Range', sep = ''), layer = 'data_0', verbose = FALSE)

## Additional functions used by Train_Pathogen_Learners
source("Tools/Functions.r")
source('Prep_Pathogen_Data.r') 
source("Calc_Sig_Pathogen_Preds.r")

## Load function that trains learners with a given set of hyperparameters
source("Train_Pathogen_Learners.r")

## Set hyperparameters.
## nboots - number of times that the data is sampled and a model is fit
## set.ambiguous.to - NA or 1; what to do with pixels with no LASV virus detected but + rodent serology
## min.test - minimum number of tested rodents required for (LASV-) pixel status. 
## tree.complexity - depth of trees that are fit at each model iteration
## mllr - minus log10 of the learning rate
## lmt - log10 of the max.trees parameter
hypers.dat <- expand.grid(nboots = 25,
                          set.ambiguous.to = c(NA),
                          min.test = c(5),
                          tree.complexity = c(1),
                          mllr = 3,
                          lmt = 7
                          )

writeLines(paste0('\n\n\n\n' , '---- Prefix: ', prefix, ' ----', '\n\n\n\n'))

starttime <- Sys.time()

## Create [prefix] directory
dirpath <- paste('Figures_Fits/',prefix,sep='')
if(!dir.exists(dirpath)){dir.create(dirpath, showWarnings = FALSE, recursive = TRUE)}

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

    ## Print model name
    writeLines(paste('\n\n\n --------- Model name:', fold, '------------- \n\n'))
    
    ## Prep the Mastomys natalensis LASV occurrence dataset and human serosurvey data
    rodlsv.survey.dat <- Prep.Pathogen.Data(hypers.dat[ii,])[[1]]

    ## Reviewer requested to find the mean LASV seroprevalence in LASV+ rodent sites
    ## LASV.sites.sero = subset(rodlsv.survey.dat, ArenaStat==1 & NumTestAb > 0,
    ##                          select = c('NumTestAb','NumPosAb', 'PropAb'))
    ## LASV.sites.sero$PropAb <- a$NumPosAb/a$NumTestAb
    ## mean(LASV.sites.sero$PropAb)
    ## ## Mean Ab is 0.24; 24%
    
    ## Map the prepared dataset
    source("Map_Cases.r")
    
    ## Predictors that are deemed significant by the Wilcox test are
    ## stored in var.names; this is a global variable that will be
    ## used by the function Train.Reservoir.Learners
    var.names <- Calc.Sig.Pathogen.Preds(rodlsv.survey.dat)

    ## Train learners with the given set of hyperparameters
    mod.out <- train.pathogen.learners(rodlsv.survey.dat, hypers.dat[ii,])

    ## Extract AUC statistics
    tree.dat <- read.table(paste('Figures_Fits/', prefix, '/', fold,'/tree_metrics.dat',sep=''),
                           header = TRUE)
    hypers.dat[ii,'auc'] <- mean(unlist(tree.dat[,'model.oob.auc']), na.rm = TRUE)
    hypers.dat[ii,'sd.auc'] <- sd(unlist(tree.dat[,'model.oob.auc']), na.rm = TRUE)

    ## Write model metrics and hyperparameter set to file
    write.table(hypers.dat, file = paste('Figures_Fits/', prefix,'/hypers_lsv_data', sep = ''),
                row.names = FALSE)

    writeLines(paste0('\n \n --- Finished hyper-parameter set:', ii, ' ---'))

}

print(paste0('Run time: ', Sys.time() - starttime))
