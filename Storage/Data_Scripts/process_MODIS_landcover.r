## This script processes the MODIS 2001 - 2018 landcover classification data.
##

## Load in packages for GIS, parallel processing, and plotting
require(raster)
require(sf)
require(parallel)
require(fields)
require(ggplot2)
require(expm)

## set version of random number generator to use to ensure
## reproducibility across R versions
RNGversion('3.6.0')
set.seed(1)

## Calculate landcover transition probabilities
calc.lc.bool <- TRUE

## Key that describes all land cover types in the original data-set
hab.key <- data.frame(value = c(0:16, 255),
                      class_name = c(
                          'Water_Bodies',
                          'Evergreen_Needleleaf_Forest',
                          'Evergreen_Broadleaf_Forest',
                          'Deciduous_Needleleaf_Forest',
                          'Deciduous_Broadleaf_Forest',
                          'Mixed_Forest',
                          'Closed_Shrubland',
                          'Open_Shrubland',
                          'Woody_Savanna',
                          'Savannas',
                          'Grasslands',
                          'Permanent_Wetlands',
                          'Croplands',
                          'Urban_BuiltUp',
                          'Cropland_Natural_Mosaic',
                          'Permanent_Snow_Ice',
                          'Barren',
                          'Unclassified'))

## Breaks used to bin landcover types. Should be chosen so that each
## bin will contain one and only one land cover value
breaks = c(-0.5 + 0:16, 254, 256)

## Make list of each year's landcover classifications, load into a stack
lc.files <- list.files("../Raster_Data/Original_Data/MODIS_LandCover/", full.names = TRUE)
lc.stack <- stack(lc.files)

## Rainfall mean raster (generated in Calc_Seasonal_Weather.r) is used as template
template.rast <- raster('../Raster_Data/Rast_Pmu.tif')

## Get list of years for which we have MODIS land cover data
lc.file.list <- strsplit(lc.files, '_')
lc.years <- sapply(lc.file.list, `[`, 6)
nyears <- length(lc.years)

## This function takes year as input, and, for that year's land cover data, loops through each
## land cover type. The function saves and returns a raster stack that describes the proportion of
## surrounding pixels that are of a particular type, for the given year.

mcl.fun = function(year){
    year.fold <- paste('../Raster_Data/Processed_MODIS_LandCover/',year,sep='')
    dir.create(path = year.fold,
               showWarnings = FALSE)
    wi.stack <- grepl(year, lc.file.list)
    lcc <- lc.stack[[which(wi.stack)]]
    for(ii in sort(hab.key$value)){

        starttime <-Sys.time()
        print(paste('Started aggregating land type', ii), quote = FALSE)

        lcc.agg <- aggregate(lcc, fact = 3,
                             fun = function(x,...){x = na.omit(x); sum(x==ii)/length(x)})
        print(paste('Finished aggregating land type', ii), quote = FALSE)
        print(Sys.time() - starttime)

        if(ii==min(hab.key$value)){
            plc.stack <- stack(lcc.agg)
        }else{
            plc.stack <- stack(plc.stack, lcc.agg)
        }
    }
    fn <- paste(year.fold, '/LC_', year, '.grd',sep='')
    pc.stack.reproj <- resample(plc.stack, template.rast, method = 'bilinear')
    pc.stack.reproj <- clamp(pc.stack.reproj, lower = 0, upper = 1, useValue = TRUE)
    names(pc.stack.reproj) <- hab.key$class_name

    writeRaster(pc.stack.reproj, filename = fn, format="raster",
                overwrite = TRUE)
    return(pc.stack.reproj)
}


if(calc.lc.bool){
    dir.create('../Raster_Data/Processed_MODIS_LandCover/', showWarnings = FALSE)
    ## The line below runs the mcl.fun function on multiple cores
    out <- stack(mclapply(X=lc.years, FUN = mcl.fun, mc.cores = detectCores() - 2))

    ## Next, form a transition matrix that describes the switching time
    ## for each land cover type, and what land cover transitions occur.
    ## Note that this uses the original data in lc.stack.

    ## count.mat tracks (counts) the number of different types of transitions,
    ## storing them in matrix form.
    count.mat = matrix(0, nrow = nrow(hab.key), ncol = nrow(hab.key))
    nl <- nlayers(lc.stack)
    for(ti in 1:(nl-1)){
        mat = as.matrix(lc.stack[[ti]])
        next.mat = as.matrix(lc.stack[[ti + 1]])
        for(ii in hab.key$value){
            ## Which pixels are of type ii
            wi.ii <- which(mat==ii, arr.ind = TRUE)
            ## Find the distribution of land cover types for wi.ii locations in the next year
            a = hist(as.vector(next.mat[wi.ii]), breaks = breaks, plot=FALSE)
            wi.hab <- which(ii == hab.key$value)
            count.mat[wi.hab,] <- count.mat[wi.hab,] + a$counts
            print(paste('ti: ', ti, 'ii: ', ii))
        }
        print(paste('Finished ', ti))
    }

    ## Remove habitats that do not occur in West Africa
    omit <- which(rowSums(count.mat) == 0 & colSums(count.mat) == 0)
    count.mat.omit <- count.mat[ -omit, -omit ]

    ## Rescale these transition counts to find the matrix describing probabilities of
    ## yearly transitions
    trans.mat <- t(apply(count.mat.omit, MARGIN = 1,
                         FUN = function(x){x/sum(x, na.rm=TRUE)}))
    write.table(trans.mat, file = 'LandCover_Figures/transition_matrix_one_year.txt')

    ## Raise transition matrix to 10th power to find 10-year transition probabilities
    trans.mat.10 <- t( t( trans.mat )%^%10 )
    write.table(trans.mat.10, file = 'LandCover_Figures/transition_matrix_ten_year.txt')
}else{
    trans.mat <- read.table(file = 'LandCover_Figures/transition_matrix_one_year.txt',
                            header = TRUE)
    trans.mat.10 <- read.table(file = 'LandCover_Figures/transition_matrix_ten_year.txt',
                               header = TRUE)

}

## Land cover classification names that show up on the figures
pretty.labels <- function(x){
    x.out <- c()
    for(xi in x){
    xi.name <- switch(xi,
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

## Make heat map of transition matrix
n.habs <- nrow(trans.mat)
png(file = 'LandCover_Figures/transition_prob_oneyear.png', width = 10, height = 8, units = 'in',
    res = 400)
par(las = 2, mar = c(10,10,2,5))
image.plot(x = 1:n.habs, y = 1:n.habs, z = t(trans.mat), xaxt = 'n', yaxt = 'n',
           xlab = '', ylab = '', main = 'Yearly Transition Probabilities',
           col = rev(heat.colors(12)))
axis(side = 1, labels = pretty.labels(paste(hab.key$class_name))[-omit], at = 1:n.habs)
axis(side = 2, labels = pretty.labels(paste(hab.key$class_name))[-omit], at = 1:n.habs)
dev.off()

## Make heat map of 10-year transition matrix
png(file = 'LandCover_Figures/transition_prob_tenyear.png', width = 10, height = 8, units = 'in',
    res = 400)
par(las = 2, mar = c(10,10,2,5))
image.plot(x = 1:n.habs, y = 1:n.habs, z = t(trans.mat.10), xaxt = 'n', yaxt = 'n',
           xlab = '', ylab = '', main = 'Decadel Transition Probabilities',
           col = rev(heat.colors(12)))
axis(side = 1, labels = pretty.labels(paste(hab.key$class_name))[-omit], at = 1:n.habs)
axis(side = 2, labels = pretty.labels(paste(hab.key$class_name))[-omit], at = 1:n.habs)
dev.off()

## Make a plot of expected transition time, E[transition time] = 1/p, where
## p is the vector of diagonal elements (non-transitions)
p.same = diag(trans.mat)
p.trans = 1 - p.same
trans.time = 1/p.trans
barplot(trans.time, names.arg = pretty.labels(paste(hab.key$class_name))[-omit], ylim = c(0,100),
        horiz = FALSE, ylab = 'Years to Transition')
dat = data.frame(LandType = pretty.labels(paste(hab.key$class_name))[-omit], Time = trans.time)
ord = order(dat$Time)
dat$LandType = factor(dat$LandType, levels = dat$LandType[ord], ordered = TRUE)
ggplot(data = dat[ord,], aes(x = LandType, y = Time)) + geom_bar(stat = 'identity') +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab('Land Type') + ylab('Years') + ggtitle('Duration of Land Types in West Africa') + coord_cartesian(ylim = c(0, 100))
ggsave(filename = 'LandCover_Figures/LC_Duration.png', width = 5, height = 5, units = 'in')

## Print list of land cover types that take, on average, at least 20 years to transition
print(paste(dat[which(dat$Time > 20),'LandType']))

## Save the list to file, omitting any types that don't occur in West Africa
dat$class_name <- hab.key$class_name[-omit]
dat.20 <- dat[which(dat$Time > 20), c('class_name', 'Time')]
print(dat.20)
write.table(x = dat.20,
            file = 'Output/landcover_variable_names_to_use',
            row.names = FALSE)

