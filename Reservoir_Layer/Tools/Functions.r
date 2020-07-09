
## Return interpretable predictor names for plots
        pretty.labels <- function(x){
            x.out <- c()
            for(xi in x){
                xi.name <- switch(xi,
                                  TreeCover = 'Tree Cover',
                                  ShrubCover = 'Shrub Cover',
                                  Grasslands = 'Grassland',
                                  Cropland = 'Cropland',
                                  AquaticVeg = 'Aquatic Veg.',
                                  SparseVeg = 'Sparse Veg.',
                                  Bare = 'Bare',
                                  BuiltUp = 'Built Up',
                                  SnowIce = 'Snow and Ice',
                                  OpenWater = 'Open Water',
                                  Pmu = 'Mean Precipitation',
                                  Tmu = 'Mean Temperature',
                                  Nmu = 'NDVI',
                                  Pmin = 'Minimum Precip',
                                  Pmax = 'Maximum Precip',
                                  Nmin = 'Minimum NDVI',
                                  Nmax = 'Maximum NDVI',
                                  Pcv = 'Precip. Coef. of Variation',
                                  Ncv = 'NDVI Coef. of Variation',
                                  Pc = 'Precip. Constancy',
                                  Pm = 'Precip. Contingency',
                                  Nc = 'NDVI Constancy',
                                  Nm = 'NDVI Contingency',
                                  Pdur = 'Dry Duration',
                                  Ndur = 'Brown Duration',
                                  Elev = 'Elevation',
                                  Pop = 'Population',
                                  TotCrop = 'Croplands',
                                  Evergreen_Needleleaf_Forest = 'Ev. Needle Forest',
                                  Evergreen_Broadleaf_Forest = 'Ev. Broad Forest',
                                  Deciduous_Needleleaf_Forest = 'De. Needle Forest',
                                  Deciduous_Broadleaf_Forest = 'De. Broad Forest',
                                  Mixed_Forest= 'Mixed Forest',
                                  Closed_Shrubland = 'Cl. Shrubland',
                                  Open_Shrubland = 'Op. Shrubland',
                                  Woody_Savanna = 'Woody Savanna',
                                  Savannas = 'Savanna',
                                  Grasslands = 'Grassland',
                                  Permanent_Wetlands = 'Wetland',
                                  Croplands = 'Cropland',
                                  Urban_BuiltUp = 'Urban',
                                  Cropland_Natural_Mosaic = 'Cropland Mosaic',
                                  Permanent_Snow_Ice = 'Snow/Ice',
                                  Barren = 'Barren',
                                  Water_Bodies = 'Water',
                                  Unclassified = 'NA')
                x.out <- c(x.out, xi.name)
            }
            return(x.out)
        }

## Impute missing values
impute <- function(dataset, var.names, impute.fun = 'median'){
    fun <- get(impute.fun)
    for(var in var.names){
        dataset[is.na(dataset[,var]),var] <- fun(na.omit(dataset[,var]))
    }
    return(dataset)
}

## Function that removes doubles (points that fall in the same 5x5 grid square
purge.repeats <- function(dat, template){
    dat.with.repeats.removed <- c()
    rem.dat <- c()
    ## Step 1: get cell numbers of each point
    points <- cbind(dat$Longitude, dat$Latitude)
    cells <- raster::extract(template, points, cellnumbers=TRUE)[,'cells']

    ## Loop through rows of dat, retain points that fall in unique
    ## cells, omit repeats based on some priority. Here, priority is
    ## given to Presence==1.
    omit.points <- c()
    keep.points <- c()
    for(jj in 1:nrow(dat)){
        if(jj %in% omit.points){
            ## Do nothing if jj is already omitted
        }else{
            repeat.set <- which(cells[jj] == cells) ## Set of repeats for point jj
            keep.point <- jj
            ## Keep first point with presence==1, if any presence points exist in repeat.set
            if(sum(dat$Presence[repeat.set]) > 0){
                keep.point <- repeat.set[which(dat$Presence[repeat.set]==1)[1]]
            }
            keep.points <- c(keep.points, keep.point)
            omit.points <- c(omit.points, repeat.set)
        } ## End if that checks if jj is in omit.points
    } ## Loop through jj
    dat.with.repeats.removed <- dat[keep.points,]
    rem.dat <- data.frame(orig.len = nrow(dat),
                          purg.len = length(keep.points))
    rem.dat$n.rem <- with(rem.dat, orig.len - purg.len)
    return(list(dat.with.repeats.removed, rem.dat))
}## End Function

## Generate name for model run directories
generate.res.name <- function(gridgi){
    fold.name <- with(gridgi, paste(paste0(substr(paste(Species),1,1),
                                           substr(strsplit(paste(Species),' ')[[1]][2],1,1)),
                                    'pa_nboots',nboots,
                                    'nbg', num.bg.points,
                                    'tc', tree.complexity,
                                    'mllr', mllr,
                                    'lmt', lmt,
                                    sep = '_'))
    return(fold.name)
}

