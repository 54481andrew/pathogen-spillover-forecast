## This script processes Mastomys natalensis presence data. The output is a
## dataframe that contains Latitude, Longitude, and associated environmental
## predictors of pixels that contain Mastomys natalensis captures and
## background captures of other Murids.

Prep.Reservoir.Data <- function(Species){

    ## Load in [Species] occurrence dataset. Keep only those columns that are needed,
    ## and discard any entries without Lat/Lon values. The data-set is called classi.dat, short
    ## for classification data-set.
    Species.name = gsub(' ', '_', Species)
    data = read.table(file = paste0('Data/Dataset_', Species.name, '.csv'), sep = ";", header = TRUE)
    keep <- c('Species', 'Latitude', 'Longitude', 'Country', 'Confidence')
    dat.1 <- data[,keep]
    classi.dat <- dat.1[which(!is.na(dat.1$Longitude) & !is.na(dat.1$Latitude)),]

    ## Only keep entries that are high confidence (species id is based on
    ## PCR, rather than based on morphology).
    keep <- which(classi.dat$Confidence=='high')
    classi.dat <- classi.dat[keep,]
    
    ## Add column indicating that these are Presence data
    classi.dat$Presence = 1

    ## -- Incorporate GenBank dataset into classi.dat
    
    ## Load GenBank data
    genbank.dat <- read.csv('Data/Rodents_Genbank_June16_2020.csv', stringsAsFactors = FALSE,
                         sep = '\t')

    ## Prune the dataset to those entries that are of Mastomys natalensis, and
    ## for which latitude and longitude are available.
    wi.genbank <- (genbank.dat$gbHost=='Mastomys natalensis') &
        !is.na(genbank.dat$Lat) & !is.na(genbank.dat$Long) &
        genbank.dat$ID_method %in% c('CytB','DNA')
    genbank.dat <- genbank.dat[wi.genbank, ]
    genbank.dat$Species <- 'Mastomys natalensis'

    ## Reformat to be compatible with survey data
    genbank.dat.form <- with(genbank.dat,
                             data.frame(Longitude = Long, Latitude = Lat,
                                        Species = Species, Presence = 1))
    classi.dat <- rbind.fill(classi.dat, genbank.dat.form)

    ## Load in GBIF background points of Muridae occurrences. Only keep
    ## entries that 1) contain non-NA lat/lon, and 2) are not [Species].
    background.dat <- data.frame(Longitude = raw.dat$decimalLongitude,
                                 Latitude = raw.dat$decimalLatitude,
                                 Species = raw.dat$species)

    keep <- !is.na(background.dat$Latitude) & !is.na(background.dat$Longitude) &
        !is.na(background.dat$Species) & paste(background.dat$Species) != Species
    background.dat <- background.dat[keep,]
    background.dat$Presence = 0

    ## Join the Presence and Background data-sets together
    classi.dat <- rbind.fill(classi.dat, background.dat)

    ## Restrict dataset to those points that occur in the focal countries of West Africa (This
    ## step does not omit any of the Mastomys natalensis presence data)
    ## First, load in shapefile of Africa and define focal countries.
    foc.countries <- c('Mali', "Guinea", "Ivory Coast", "Sierra Leone",
                       "Nigeria", 'Liberia', "Cote d`Ivoire",
                       'Ghana', 'Togo', 'Benin',
                       'Burkina Faso', 'Mauritania', 'Niger', 'Senegal',
                       'Gambia', 'Guinea-Bissau')

    ## Next, make a temporary spatial points object (dat) out of rodent coordinates.
    ## Extract feature info from locations contained in dat from the africa shape file.
    dat <- data.frame(Longitude = classi.dat$Longitude,
                      Latitude = classi.dat$Latitude)
    coordinates(dat) <- ~ Longitude + Latitude
    proj4string(dat) <- proj4string(africa.shp.ogr)
    over.WA <- over(dat, africa.shp.ogr)

    ## Remove any entries with an NA value of COUNTRY feature, as well as
    ## entries that are not in the focal countries.
    keep <- !is.na(over.WA$COUNTRY) & over.WA$COUNTRY %in% foc.countries
    classi.dat <- classi.dat[keep,]

    ## Next, only keep data that falls within the IUCN Mastomys natalensis range map.
    ## This step does not exclude any M natalensis presence points.
    dat <- classi.dat ## Make a copy that will become a spatialpointsdataframe
    coordinates(dat) <- ~ Longitude + Latitude
    proj4string(dat) <- proj4string(masto.rangemap)
    overMastobg <- over(dat, masto.rangemap)
    wi.range <- which(!is.na(overMastobg$PRESENCE))
    classi.dat <- classi.dat[wi.range,]

    ## Finally, only keep points that fall within unique pixels. This is
    ## achieved with the purge.repeats function. The pixelation of
    ## West Africa is determined by the 0.05x0.05 degree MODIS precipitation
    ## data-set, which is a raster inside of the all.stack raster-stack.
    template.rast <- all.stack[['Pmu']]
    purged.out <- purge.repeats(classi.dat, template.rast)
    classi.dat <- purged.out[[1]]
    purge.info <- purged.out[[2]] ## diagnostics

    ## Change coordinates to be at cell center
    orig.points <- cbind(classi.dat$Longitude, classi.dat$Latitude)
    cells <- raster::extract(template.rast, orig.points, cellnumbers=TRUE)[,'cells']
    cell.points <- as.data.frame(xyFromCell(object = template.rast, cell = cells))
    classi.dat$Longitude <- cell.points$x
    classi.dat$Latitude <- cell.points$y

    ## Next, incorporate the environmental predictors into the
    ## dataset. All predictors are extracted from all.stack.
    points <- cbind(classi.dat$Longitude, classi.dat$Latitude)
    classi.dat <- cbind(classi.dat, raster::extract( all.stack, points, sp = TRUE))

    ## Impute any missing values with median
    var.names <- names(all.stack)
    classi.dat <- impute(classi.dat, var.names)

    ## Write the data-set to the Data directory
    write.table(classi.dat, file = paste0('Data/Prepped_PresAbs_', Species.name, '_Data'),
                col.names = TRUE, row.names = FALSE)

    ## Print number of pixels classified as presences for each species
    print('-- Count Statistics --', quote = FALSE)
    tab <- table(classi.dat$Species)
    print(addmargins(tab, FUN = list(Total = sum), quiet = TRUE))

    ## Quick map figure showing pixel [Species] capture locations. Load shapefile for plotting.
    png(file = paste0('Figures_Fits/', Species.name,'_captures.png'),
        width = 6, height = 4, units = 'in', res = 400)
    par(mai = 1*c(0.2,0.2,0.2,0.2))
    plot(foc.shp.ogr, col = 'cornsilk')
    points(Latitude~Longitude, classi.dat[classi.dat$Presence==1,], bg = 'green', col = 'black',
           pch = 21, cex = 0.5)
    legend(x = 'topleft', legend = 'Confirmed captures', pch = 21, pt.bg = 'green', col = 'black')
    dev.off()

    return(classi.dat)
} ## End function
