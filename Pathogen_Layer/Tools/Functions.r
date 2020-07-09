## Add columns describing total individuals tested, total
## number of individuals that were positive, and
## general positive / absence of LASV.
add.total.columns <- function(dat){
    num.pos.ag <- with(dat, ifelse(is.na(NumPosAg), 0, NumPosAg))
    num.test.ag <- with(dat, ifelse(is.na(NumTestAg), 0, NumTestAg))
    num.pos.ab <- with(dat, ifelse(is.na(NumPosAb), 0, NumPosAb))
    num.test.ab <- with(dat, ifelse(is.na(NumTestAb), 0, NumTestAb))
    tot.test <- num.test.ag + num.test.ab
    dat$TotTest <- tot.test
    dat$TotPos <- num.pos.ag + num.pos.ab

    ## Add column indicating general arenavirus status. Note that for
    ## an absence, the survey must have tested at least 5 individuals.
    dat$ArenaStat <- NA
    pos.result <- with(dat, NumPosAg > 0)
    neg.result <- with(dat, (TotTest >= 5 & TotPos==0))

    dat[which(pos.result), 'ArenaStat'] <- 1
    dat[which(neg.result), 'ArenaStat'] <- 0

    return(dat)
}

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

## Function that removes doubles (points that fall in the same 0.05x0.05 degree grid square
purge.repeats <- function(dat, template){
    dat.with.repeats.removed <- c()
    rem.dat <- c()
    ## Get cell numbers of each point
    points <- cbind(dat$Longitude, dat$Latitude)
    cells <- raster::extract(template, points, cellnumbers=TRUE)[,'cells']

    ## Loop through rows of dat dataframe, retain points that fall in unique
    ## cells, omit repeats based on some ranking
    omit.points <- c()
    keep.points <- c()
    for(jj in 1:nrow(dat)){
        if(jj %in% omit.points){
            ## Do nothing if jj is already omitted
        }else{
            ## Set of repeats for point jj
            repeat.set <- which(cells[jj] == cells)
            keep.points <- c(keep.points, jj)
            omit.points <- c(omit.points, repeat.set)
            keep.repeat <- jj

            ## --- Aggregate information on repeats
            dat.with.repeats.removed <- rbind(dat.with.repeats.removed,
                                              dat[keep.repeat,])

            ## Total tested/positive
            ldat <- nrow(dat.with.repeats.removed)
            dat.with.repeats.removed$TotTest[ldat] <- sum(dat[repeat.set, 'TotTest'], na.rm = TRUE)
            dat.with.repeats.removed$TotPos[ldat] <- sum(dat[repeat.set, 'TotPos'], na.rm = TRUE)
            ## Serology
            dat.with.repeats.removed$NumTestAb[ldat] <- sum(dat[repeat.set, 'NumTestAb'], na.rm = TRUE)
            dat.with.repeats.removed$NumPosAb[ldat] <- sum(dat[repeat.set, 'NumPosAb'], na.rm = TRUE)
            ## PCR
            dat.with.repeats.removed$NumTestAg[ldat] <- sum(dat[repeat.set, 'NumTestAg'], na.rm = TRUE)
            dat.with.repeats.removed$NumPosAg[ldat] <- sum(dat[repeat.set, 'NumPosAg'], na.rm = TRUE)
            ## Other aggregated stats
            dat.with.repeats.removed$ArenaStat[ldat] <- max(dat[repeat.set, 'ArenaStat'])
        } ## End if checking for repeats
    } ## Loop through jj
    ## Statistics on the aggregation
    rem.dat <- data.frame(orig.len = nrow(dat),
                          purg.len = length(keep.points))
    rem.dat$n.rem <- with(rem.dat, orig.len - purg.len)
    return(list(dat.with.repeats.removed, rem.dat))
}## End Function

## Generate a folder name for tree with specified properties
generate.lassa.name <- function(gridgi){
    fold.name <- with(gridgi, paste('pa_nboots', nboots,
                                    tree.complexity,
                                    'mllr', mllr,
                                    'lmt', lmt,
                                    sep='_')) ## Folder suffix for data storage and results
    return(fold.name)
}


