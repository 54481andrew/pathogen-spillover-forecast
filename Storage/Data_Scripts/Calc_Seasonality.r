## This script calculates the measures of seasonality from the original
## MODIS data-set. Specifically, this code will calculate the means, Colwell's
## indices, coefficient of variation, duration, and minimum/maximum measures
## of precipitation, temperature, and NDVI.

## These boolean variables determine which seasonal attributes are calculated.
calc.mean <- TRUE ## mean rasters
calc.colwell <- TRUE ## colwell indices
calc.cv <- TRUE ## coefficient of variation
calc.duration <- TRUE ## duration
calc.minmax <- TRUE ## average min/max

## Packages for handling rasters and shapefiles
require(raster)
require(rgdal)
require(sf)

## set version of random number generator to use to ensure
## reproducibility across R versions
RNGversion('3.6.0')
set.seed(1)

## Locations of raster and shapefile
storage.fold <- '../../Storage'
shapefile.storage.fold <- paste(storage.fold,'/Shapefiles', sep = '')
raster.data.storage.fold = paste(storage.fold, '/Raster_Data', sep = '')

## Load in West Africa shapefile for plotting
foc.shp <- st_read(dsn = '../Shapefiles/West_Africa', layer = 'foc')

starttime <- Sys.time()

## These are the directories that contained the original MODIS data from 2001 - 2019
prec.fold <- paste(raster.data.storage.fold,'/Original_Data/MODIS_Precip/', sep = '')
temp.fold <- paste(raster.data.storage.fold,'/Original_Data/MODIS_Temp/', sep = '')
ndvi.fold <- paste(raster.data.storage.fold,'/Original_Data/MODIS_NDVI/', sep = '')

## Load in raster locations, then load in rasters as stacks
prec.names <- list.files(prec.fold)
temp.names <- list.files(temp.fold)
ndvi.names <- list.files(ndvi.fold)
prec.stack <- stack(paste(prec.fold, prec.names, sep = ''))
temp.stack <- stack(paste(temp.fold, temp.names, sep = ''))
ndvi.stack <- stack(paste(ndvi.fold, ndvi.names, sep = ''))

## Calculate mean rasters

if(calc.mean){
    print('Calculating means', quote = FALSE)
    ## Precipitation mean and standard deviation
    stat = calc(prec.stack, fun = function(x){mean(x, na.rm = TRUE)})
    writeRaster(x = stat, filename = paste(storage.fold, '/Raster_Data','/Rast_Pmu.tif',sep=''),
                overwrite = TRUE)

    ## Temperature mean and standard deviation
    stat = calc(temp.stack, fun = function(x){mean(x, na.rm = TRUE)})
    writeRaster(x = stat, filename = paste(storage.fold, '/Raster_Data','/Rast_Tmu.tif',sep=''),
                overwrite = TRUE)

    ## NDVI mean
    stat = calc(ndvi.stack, fun = function(x){mean(x, na.rm = TRUE)})
    writeRaster(x = stat, filename = paste(storage.fold, '/Raster_Data','/Rast_Nmu.tif',sep=''),
                overwrite = TRUE)
}

## The mean precipitation raster serves as a template for all other raster data-sets
template.rast <- raster('../Raster_Data/Rast_Pmu.tif')

## Calculate minimum/maximum rasters
if(calc.minmax){
    calc.minmax <- function(stack.rasts){
        stack.names <- names(stack.rasts)
        min.stack <- stack()
        max.stack <- stack()
        for(year in 2001:2019){
            foc.stack <- stack.rasts[[grep(paste(year), stack.names) ]]
            min.year <- min(foc.stack, na.rm = TRUE)
            max.year <- max(foc.stack, na.rm = TRUE)
            min.stack <- stack(min.year, min.stack)
            max.stack <- stack(max.year, max.stack)
            print(year)
        }
        min.mu = mean(min.stack, na.rm = TRUE)
        max.mu = mean(max.stack, na.rm = TRUE)
        return(c(min.mu,max.mu))
    }

    ## Precipitation
    print('Calculating precipitation minmax', quote = FALSE)
    minmax.prec <- calc.minmax(prec.stack)
    minmax.reproj <- stack(resample(minmax.prec[[1]], template.rast, method = 'bilinear'),resample(minmax.prec[[2]], template.rast, method = 'bilinear'))
    minmax.clamp  <- clamp(minmax.reproj, lower = 0, upper = Inf, useValue = TRUE)
    names(minmax.clamp) <- c('Pmin','Pmax')
    writeRaster(minmax.clamp, filename = paste(storage.fold, '/Raster_Data','/Rast_Pminmax.grd',sep=''),
                overwrite = TRUE)
    jpeg('Output/minmax_precip.jpg')
    plot(minmax.clamp)
    plot(foc.shp[,2], add = TRUE)
    dev.off()

    ## Temperature
    print('Calculating temperature minmax', quote = FALSE)
    minmax.temp <- calc.minmax(temp.stack)
    minmax.reproj <- stack(resample(minmax.temp[[1]], template.rast, method = 'bilinear'),
                           resample(minmax.temp[[2]], template.rast, method = 'bilinear'))
    minmax.clamp  <- clamp(minmax.reproj, lower = 0, upper = Inf, useValue = TRUE)
    names(minmax.clamp) <- c('Tmin','Tmax')
    writeRaster(minmax.clamp, filename = paste(storage.fold, '/Raster_Data','/Rast_Tminmax.grd',sep=''),
                overwrite = TRUE)
    jpeg('Output/minmax_temp.jpg')
    plot(minmax.clamp)
    plot(foc.shp[,2], add = TRUE)
    dev.off()

    ## NDVI
    print('Calculating NDVI minmax', quote = FALSE)
    minmax.ndvi <- calc.minmax(ndvi.stack)
    minmax.reproj <- stack(resample(minmax.ndvi[[1]], template.rast, method = 'bilinear'),
                           resample(minmax.ndvi[[2]], template.rast, method = 'bilinear'))
    minmax.clamp  <- clamp(minmax.reproj, lower = 0, upper = Inf, useValue = TRUE)
    names(minmax.clamp) <- c('Nmin','Nmax')
    writeRaster(minmax.clamp, filename = paste(storage.fold, '/Raster_Data','/Rast_Nminmax.grd',sep=''),
                overwrite = TRUE)
    jpeg('Output/minmax_ndvi.jpg')
    plot(minmax.clamp)
    plot(foc.shp[,2], add = TRUE)
    dev.off()

}## End if calc.minmax


## Calculate measures of duration of small NDVI, small Precipitation, and low Temperature
if(calc.duration){

    ## Calculate coefficient of variation
    calc.duration <- function(stack.rasts){
        stack.names <- names(stack.rasts)
        cv.stack <- stack()
        for(year in 2001:2019){
            foc.stack <- stack.rasts[[grep(paste(year), stack.names) ]]
            duration.year <- sum(foc.stack < thresh) ## 1 cm/mo is dry (avg is 12.5 cm/mo)
            cv.stack <- stack(duration.year, cv.stack)
            print(year)
        }
        cv.mu = mean(cv.stack, na.rm = TRUE)
        return(cv.mu)
    }

    ## Precipitation
    print('Calculating precipitation duration', quote = FALSE)
    thresh = 1
    dur.prec <- calc.duration(prec.stack)
    dur.reproj <- resample(dur.prec, template.rast, method = 'bilinear')
    dur.clamp  <- clamp(dur.reproj, lower = 0, upper = 12, useValue = TRUE)
    writeRaster(dur.clamp, filename = paste(storage.fold, '/Raster_Data','/Rast_Pdur.tif',sep=''),
                overwrite = TRUE)
    jpeg('Output/duration_precip.jpg')
    plot(dur.clamp)
    plot(foc.shp[,2], add = TRUE)
    dev.off()

    ## NDVI
    print('Calculating NDVI duration', quote = FALSE)
    thresh <- 0.5
    dur.ndvi <- calc.duration(ndvi.stack)
    dur.reproj <- resample(dur.ndvi, template.rast, method = 'bilinear')
    dur.clamp  <- clamp(dur.reproj, lower = 0, upper = 12, useValue = TRUE)
    writeRaster(dur.clamp, filename = paste(storage.fold, '/Raster_Data','/Rast_Ndur.tif',sep=''),
                overwrite = TRUE)
    jpeg('Output/duration_ndvi.jpg')
    plot(dur.clamp)
    plot(foc.shp[,2], add = TRUE)
    dev.off()

}## End if calc.duration


## Calculate coefficient of variation for Temperature, Precipitation, and NDVI
if(calc.cv){

    ## Calculate coefficient of variation
    calc.cv <- function(stack.rasts){
        stack.names <- names(stack.rasts)
        cv.stack <- stack()
        for(year in 2001:2019){
            foc.stack <- stack.rasts[[grep(paste(year), stack.names) ]]
            cv.year <- cv(foc.stack, na.rm = TRUE)/100
            cv.stack <- stack(cv.year, cv.stack)
            print(year)
        }
        cv.mu = mean(cv.stack, na.rm = TRUE)
        return(cv.mu)
    }

    ## Precipitation
    print('Calculating precipitation coefficient of variation', quote = FALSE)
    cv.prec <- calc.cv(prec.stack)
    cv.prec[is.na(cv.prec)] <- 0
    cv.reproj <- resample(cv.prec, template.rast, method = 'bilinear')
    cv.clamp  <- clamp(cv.reproj, lower = 0, upper = Inf, useValue = TRUE)
    writeRaster(cv.clamp, filename = paste(storage.fold, '/Raster_Data','/Rast_Pcv.tif',sep=''),
                overwrite = TRUE)
    jpeg('Output/coef_var_precip.jpg')
    plot(cv.clamp)
    plot(foc.shp[,2], add = TRUE)
    dev.off()

    ## Temperature
    print('Calculating temperature coefficient of variation', quote = FALSE)
    cv.temp <- calc.cv(temp.stack)
    cv.temp[is.na(cv.temp)] <- 0
    cv.reproj <- resample(cv.temp, template.rast, method = 'bilinear')
    cv.clamp  <- clamp(cv.reproj, lower = 0, upper = Inf, useValue = TRUE)
    writeRaster(cv.clamp, filename = paste(storage.fold, '/Raster_Data','/Rast_Tcv.tif',sep=''),
                overwrite = TRUE)
    jpeg('Output/coef_var_temp.jpg')
    plot(cv.clamp)
    plot(foc.shp[,2], add = TRUE)
    dev.off()

    ## NDVI
    print('Calculating NDVI coefficient of variation', quote = FALSE)
    cv.ndvi <- calc.cv(ndvi.stack)
    cv.ndvi[is.na(cv.ndvi)] <- 0
    cv.reproj <- resample(cv.ndvi, template.rast, method = 'bilinear')
    cv.clamp  <- clamp(cv.reproj, lower = 0, upper = Inf, useValue = TRUE)
    writeRaster(cv.clamp, filename = paste(storage.fold, '/Raster_Data','/Rast_Ncv.tif',sep=''),
                overwrite = TRUE)
    jpeg('Output/coef_var_ndvi.jpg')
    plot(cv.clamp)
    plot(foc.shp[,2], add = TRUE)
    dev.off()

}## End if calc.cv


## Calculate colwell's indices
if(calc.colwell){

    print('Calculating Colwells indices', quote = FALSE)

    colwell.fxn <- function(vec){
        for(mo in 1:12){
            wi.mo <- seq(mo,length(vec), by = 12)
            prec.hist = hist(vec[wi.mo], breaks, plot = FALSE)$counts
            mat[,mo] <- prec.hist
        }
        Xj = colSums(mat)
        Yi = rowSums(mat)
        Z = sum(mat)

        Hx = -sum(Xj/Z*log(Xj/Z), na.rm = TRUE)
        Hy = -sum(Yi/Z*log(Yi/Z), na.rm = TRUE)
        Hxy = -sum(mat/Z*log(mat/Z), na.rm = TRUE)

        P = 1 - (Hxy - Hx)/log(s) ## predictability
        C = 1 - Hy/log(s) ## Constancy
        M = (Hx + Hy - Hxy)/log(s) ## Contingency

        return(c(C,M,P))
    }

    ## Precipitation
    print("Calculating precipitation Colwell's index", quote = FALSE)
    prec.max <- max(values(prec.stack), na.rm = TRUE)
    prec.max <- ceiling(prec.max/10)*10 + 10
    breaks <- seq(-0.5,prec.max, by = 1 )
    mat <- matrix(NA, nrow = length(breaks) - 1, ncol = 12)
    s = nrow(mat)
    t = ncol(mat)
    stat = calc(prec.stack, fun = colwell.fxn, progress = 'text')
    names(stat) = c('Pcolc','Pcolm','Pcolp')
    stat <- resample(stat, template.rast, method = 'bilinear')
    stat <- clamp(stat, lower = 0, upper = 1, useValue = TRUE)
    writeRaster(stat, filename = paste(storage.fold, '/Raster_Data/Rast_Pcol.tif', sep = ''),
                format = "raster",
                overwrite = TRUE)

    ## Temperature
    print("Calculating temperature Colwell's index", quote = FALSE)
    temp.max <- max(values(temp.stack), na.rm = TRUE)
    temp.max <- ceiling(temp.max/10)*10 + 10
    breaks <- seq(-0.5,temp.max, by = 1 )
    mat <- matrix(NA, nrow = length(breaks) - 1, ncol = 12)
    s = nrow(mat)
    t = ncol(mat)
    stat = calc(temp.stack, fun = colwell.fxn, progress = 'text')
    names(stat) = c('Tcolc','Tcolm','Tcolp')
    stat <- resample(stat, template.rast, method = 'bilinear')
    stat <- clamp(stat, lower = 0, upper = 1, useValue = TRUE)
    writeRaster(stat, filename = paste(storage.fold, '/Raster_Data/Rast_Tcol.tif', sep = ''),
                format = "raster",
                overwrite = TRUE)

    ## NDVI
    print("Calculating NDVI Colwell's index", quote = FALSE)
    ndvi.max <- 1.02
    breaks <- seq(0,ndvi.max, by = 0.02 )
    mat <- matrix(NA, nrow = length(breaks) - 1, ncol = 12)
    s = nrow(mat)
    t = ncol(mat)

    stat = calc(ndvi.stack, fun = colwell.fxn, progress = 'text')
    names(stat) = c('Ncolc','Ncolm','Ncolp')
    stat <- resample(stat, template.rast, method = 'bilinear')
    stat <- clamp(stat, lower = 0, upper = 1, useValue = TRUE)
    writeRaster(stat, filename = paste(storage.fold, '/Raster_Data/Rast_Ncol.tif', sep = ''),
                format = "raster",
                overwrite = TRUE)

}## End if calc.colwell
