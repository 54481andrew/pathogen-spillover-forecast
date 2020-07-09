## This script loops through sets of hyper-parameters that specify a
## boosted regression tree algorithm. Specifically, this script passes a
## set of parameters and data-set to the train.reservoir.learners
## function, which, in turn, performs the boosting algorithm and
## generates output. Most of the work in the boosting algorithm for this
## layer is contained in the Train_Reservoir_Learners.r script, which
## defines the train.reservoir.learners.r function (see below).

## Require packages used for gis, parallel,
## boosted tree modeling, plotting, and
## calculating the AUC statistic
require(raster)
require(rgdal)
require(sf)
require(parallel)
require(dismo)
require(gbm)
require(fields)
require(ggplot2)
require(gridExtra)
require(plyr)
require(verification)
require(viridis)

## Directory name that will contain output for all hyperparameter sets
prefix <- 'reservoir_v2'

## Set version of random number generator to use to ensure
## reproducibility across R versions
RNGversion('3.6.0')

## Storage locations
storage.fold <- '../Storage'
raster.data.fold <- paste(storage.fold,'/Raster_Data', sep = '')

## Load in predictor stack. Used for predicting across West Africa
all.stack <- stack(paste(raster.data.fold, "/predictor_stack.grd", sep = ''))

## Load in shape files for plotting West Africa, and the
## Mastomys natalensis homerange
foc.shp.ogr <- readOGR(dsn = paste(storage.fold, '/Shapefiles/West_Africa', sep = ''), layer = 'foc', verbose = FALSE)
africa.shp.ogr <- readOGR(dsn = paste(storage.fold, '/Shapefiles/Africa', sep = ''), layer = 'Africa')
masto.rangemap <-readOGR(dsn = paste(storage.fold, '/Shapefiles/Masto_Range', sep = ''), layer = 'data_0', verbose = FALSE)

## The raw GBIF data provides background points
raw.dat <- read.csv('Data/GBIF_Background_Murid_Data.csv', sep = "\t")

## Additional functions that prep the data and predictors
source("Tools/Functions.r")
source('Prep_Reservoir_Data.r')
source("Calc_Sig_Reservoir_Preds.r")

## Load function that trains learners with a given set of hyperparameters
source("Train_Reservoir_Learners.r")

## Set hyperparameters.
## nboots - number of times the data is sampled and the model fit
## num.bg.points - number of background points. Character entry sets this to the number of presences
## tree.complexity - depth of trees that are fit at each model iteration
## mllr - minus log10 of the learning rate
## lmt - log10 of the max.trees parameter
grid <- expand.grid(Species = c('Mastomys natalensis'),
                    nboots = 25,
                    num.bg.points = c('Same'),
                    tree.complexity = c(1),
                    mllr = 2,
                    lmt = 7
                    )

cat(paste0('\n\n\n\n' , '---- ', prefix, ' ----', '\n\n\n\n'))

starttime <- Sys.time()

for(gi in 1:nrow(grid)){

    Species = paste(grid$Species[gi])

    ## Prep in the presence/absence data
    classi.dat.rod <- Prep.Reservoir.Data(Species)

    ## Predictors that are deemed significant by the Wilcox test
    var.names <- Calc.Sig.Reservoir.Preds(Species, classi.dat.rod)

    ## Train learners with the given set of hyperparameters
    out <- Train.Reservoir.Learners(dataset = classi.dat.rod,
                                    grid[gi,])

    ## Extract AUC statistics with (auc.pwd) and without (auc) pairwise distance sampling
    fold <- generate.res.name(grid[gi,])
    tree.dat <- read.table(paste('Figures_Fits/', prefix, '/', fold,'/tree.dat',sep=''),
                           header = TRUE)
    grid[gi,'auc'] <- mean(unlist(tree.dat[,'model.oob.auc']), na.rm = TRUE)
    grid[gi,'sd.auc'] <- sd(unlist(tree.dat[,'model.oob.auc']), na.rm = TRUE)

    grid[gi,'auc.pwd'] <- mean(unlist(tree.dat[,'model.oob.auc.pwd']), na.rm = TRUE)
    grid[gi,'sd.auc.pwd'] <- sd(unlist(tree.dat[,'model.oob.auc.pwd']), na.rm = TRUE)

    ## Write AUC statistics to file
    write.table(grid, file = paste('Figures_Fits/', prefix, '/grid_res_data', sep = ''),
                row.names = FALSE)

    print(paste('--- gi:', gi, ' ---', sep = ''), quote = FALSE)
}

print(paste0('Run time: ', Sys.time() - starttime))



a = read.table(file = paste('Figures_Fits/', prefix, '/', fold, '/assess.dat', sep = ''))
spec.names = unique(paste(classi.dat.rod$Species))
meana = a[,1:length(spec.names)]
names(meana) = spec.names

sda <- a[,-c(1:length(spec.names))]
names(sda) = spec.names

dat <- data.frame(mean = colMeans(meana, na.rm = TRUE),
                  sd = colMeans(sda, na.rm = TRUE))
dat[order(dat$mean),]

spec.names[order(rowMeans(meana, na.rm = TRUE))]
