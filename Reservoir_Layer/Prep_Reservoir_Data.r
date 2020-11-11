## This script processes Mastomys natalensis presence data. The output is a
## dataframe that contains Latitude, Longitude, and associated environmental
## predictors of pixels that contain Mastomys natalensis captures and
## background captures of other Murids.

Prep.Reservoir.Data <- function(Species){
    ## Load in [Species] occurrence dataset. Discard any entries without Lat/Lon values.
    ## The data-set is called classi.dat, short for classification data-set.
    Species.name = gsub(' ', '_', Species)
    presence.data = read.table(file = paste0('Data/', Species.name, '_presences.csv'),
                               sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    presence.data$Latitude = as.numeric(presence.data$Latitude)
    presence.data$Longitude = as.numeric(presence.data$Longitude)

    ## Make "Year" entry a range with minYear and maxYear. Remove entries without year
    presence.data$Year = as.character(presence.data$Year)
    keep <- presence.data$Year!=""
    presence.data <- presence.data[keep,]
    for(yi in 1:nrow(presence.data)){
        yr <- presence.data$Year[yi]
        if(grepl(',', yr)){
            yr.split <- strsplit(yr,',')[[1]]
        }else{
            yr.split <- strsplit(yr,'-')[[1]]
        }        
        yr.split <- as.numeric(yr.split)
        presence.data$min.Year[yi] <- min(yr.split)
        presence.data$max.Year[yi] <- max(yr.split)
        presence.data$Year[yi] <- with(presence.data[yi,], (min.Year + max.Year)/2)
    }
    presence.data$Year <- as.numeric(presence.data$Year)

    ## Keep only those columns that are needed
    keep <- c('Latitude', 'Longitude', 'Country', 'Species', 'Year', 'min.Year', 'max.Year')
    dat.1 <- presence.data[,keep]
    classi.dat <- dat.1[which(!is.na(dat.1$Longitude) & !is.na(dat.1$Latitude)),]

    ## Add column indicating that these are Presence data
    classi.dat$Presence = 1

    ## Add column indicating that these data are not directly from genbank
    classi.dat$genbank = FALSE
    
    ## -- Incorporate GenBank dataset into classi.dat
    
    ## Load GenBank data
    genbank.dat <- read.csv('Data/Rodents_Genbank_Oct_2020.csv', stringsAsFactors = FALSE,
                         sep = '\t')

    ## Prune the dataset to those entries that are of Mastomys natalensis, and
    ## for which latitude and longitude are available.
    wi.genbank <- (genbank.dat$gbHost=='Mastomys natalensis') &
        !is.na(genbank.dat$Lat) & !is.na(genbank.dat$Long) &
        genbank.dat$ID_method %in% c('CytB','DNA') &
        !is.na(genbank.dat$gbCollectYear)
    genbank.dat <- genbank.dat[wi.genbank, ]
    genbank.dat$Species <- 'Mastomys natalensis'
    genbank.dat$min.Year <- genbank.dat$gbCollectYear
    genbank.dat$max.Year <- genbank.dat$gbCollectYear
    
    ## Reformat to be compatible with survey data
    genbank.dat.form <- with(genbank.dat,
                             data.frame(Longitude = Long, Latitude = Lat,
                                        Country = Country,
                                        Genus = 'Mastomys',
                                        Species = Species, Presence = 1,
                                        Year = gbCollectYear, min.Year = min.Year,
                                        max.Year = max.Year,
                                        genbank = TRUE,
                                        uID = UniqueID))
    classi.dat <- rbind.fill(classi.dat, genbank.dat.form)
    minimum.year = min(classi.dat$min.Year)

    ###
    ## Incorporate GBIF background points of Muridae occurrences.
    background.dat <- data.frame(countryCode = raw.dat$countryCode,
                                 Longitude = raw.dat$decimalLongitude,
                                 Latitude = raw.dat$decimalLatitude,
                                 Species = raw.dat$species,
                                 Genus = raw.dat$genus,
                                 Year = raw.dat$year,
                                 min.Year = raw.dat$year,
                                 max.Year = raw.dat$year,
                                 basis = raw.dat$basisOfRecord
                                 )
    background.dat$countryCode = paste(background.dat$countryCode)
    
    ## Replace country codes with country names. First, define the focal countries
    ## that make up the study region. 
    foc.countries <- c('Mali', "Guinea", "Ivory Coast", "Sierra Leone",
                       "Nigeria", 'Liberia', "Cote d`Ivoire",
                       'Ghana', 'Togo', 'Benin',
                       'Burkina Faso', 'Mauritania', 'Niger', 'Senegal',
                       'Gambia', 'Guinea-Bissau')

    ## Define a table that associates country codes with countries
    country.codes = data.frame(code = c('ML', 'GN', 'CI', 'SL',
                                        'NG','LR','CI',
                                        'GH','TG','BJ',
                                        'BF','MR','NE','SN',
                                        'GM','GW'),
                               country = foc.countries)
    ## Remove repeated Ivory Coast (kept in foc.countries for later)    
    country.codes <- country.codes[-which(country.codes$country=="Cote d`Ivoire"),]
    country.codes$code = paste(country.codes$code)
    country.codes$country = paste(country.codes$country)
    
    ##Only keep entries that
    ## 1) contain non-NA lat/lon
    ## 2) have a species entry
    ## 3) have a species entry other than [Species]
    ## 4) have year of capture
    ## 5) have year of capture that is in-line with presence data
    ## 6) document a rodent that was preserved in some way
    ## 7) have a country code in the study region
    ## 8) remove other entries from the genus Mastomys
    keep <- !is.na(background.dat$Latitude) & !is.na(background.dat$Longitude) &
        !is.na(background.dat$Species) & paste(background.dat$Species)!="" &
        paste(background.dat$Species) != Species & !is.na(background.dat$Year) &
        background.dat$Year >= minimum.year &
        background.dat$basis %in% c('PRESERVED_SPECIMEN', 'MATERIAL_SAMPLE') &
        background.dat$countryCode %in% country.codes$code &
        !is.na(background.dat$Genus) &
        !(background.dat$Genus %in% c('Mastomys'))

    
    background.dat <- background.dat[keep,]
    background.dat$Presence = 0
    background.dat$genbank = FALSE
    
    ## Join the Presence and Background data-sets together
    classi.dat <- rbind.fill(classi.dat, background.dat)
    
    ## Double check that all points occur in the focal countries of West Africa (This
    ## step does not omit any of the Mastomys natalensis presence data)
    
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

    ## Now store information on the amount of each type of data
    mask.pres <- classi.dat$Presence==1
    mask.gen <- classi.dat$genbank
    background.dat <- classi.dat[!mask.pres,]
    gen.dat <- classi.dat[mask.pres & mask.gen,]
    lit.dat <- classi.dat[mask.pres & !mask.gen,]
    info <- paste0('-- Count Statistics -- \n \n',
                   'Before Purge: \n\n', 
                   'Number Background: ', nrow(background.dat), '\n', 
                   paste0('Collected Between: ', min(background.dat$min.Year),
                          ' - ', max(background.dat$max.Year)), '\n', 
                   paste0('----\n'),
                   'Literature Presences: ', nrow(lit.dat), '\n',
                   paste0('Collected Between: ', min(lit.dat$min.Year),
                          ' - ', max(lit.dat$max.Year)), '\n',
                   paste0('----\n'),
                   'Genbank Presences: ', nrow(gen.dat), '\n',
                   paste0('Collected Between: ', min(gen.dat$min.Year),
                          ' - ', max(gen.dat$max.Year)), '\n')
    
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


    ## Finally, transform country codes into country names, remove unnecessary columns
    mask.bg <- classi.dat$Presence==0
    Background.countries = sapply(X = which(mask.bg),
                                  FUN = function(X){country.codes$country[which(classi.dat$countryCode[X]==country.codes$code)]})
    classi.dat[mask.bg, 'Country'] = Background.countries
    classi.dat <- classi.dat[,-which(names(classi.dat) %in% c('countryCode','basis'))]

    
    ## Write the data-set to the Data directory
    write.csv(classi.dat, file = paste0('Figures_Fits/', prefix, '/',fold,'/',
                                        'Prepped_PresAbs_', Species.name, '_Data.csv'),
              row.names = FALSE)
    
    ## Print number of pixels classified as presences for each species
    print('-- Count Statistics --', quote = FALSE)
    nback <- sum(classi.dat$Presence==0)
    npres <- sum(classi.dat$Presence==1)
    yrange <- paste0(min(classi.dat$min.Year), ' - ', max(classi.dat$max.Year))
    print(paste0('Number Background: ', nback), quote = FALSE)
    print(paste0('Number Presence: ', npres), quote = FALSE)
    print(paste0('Collected Between: ', yrange), quote = FALSE)

    ## Lines below can be used during debug to see which genbank points were used. Helps avoid duplicates. 
    ## uids = classi.dat[!is.na(classi.dat$uID),'uID']    
    ## genbank.dat[genbank.dat$UniqueID%in%uids, c('Reference', 'LocVillage', 'Lat', 'Long', 'min.Year')]
    
    ## Write information on data preparation to file
    info <- paste0(info, '\n\n', 'After Purge: \n\n', 
                   'Number Background: ', nback, '\n', 
                   'Number Presences: ', npres, '\n',
                   'All Collected Between: ', yrange)
    write(info, paste0('Figures_Fits/', prefix, '/',fold,'/','data_prep_info.txt'))

    
    ## Map figure showing pixel [Species] capture locations. Load shapefile for plotting.
    png(file = paste0('Figures_Fits/',prefix, '/',fold,'/', Species.name,'_captures.png'),
        width = 6, height = 4, units = 'in', res = 400)
    par(mai = 1*c(0.2,0.2,0.2,0.2))
    plot(foc.shp.ogr, col = 'cornsilk')
    points(Latitude~Longitude, classi.dat[classi.dat$Presence==1,], bg = 'green', col = 'black',
           pch = 21, cex = 0.5)
    points(Latitude~Longitude, classi.dat[classi.dat$Presence==0,], col = 'black',
           pch = 4, cex = 0.5)
    legend(x = 'topleft', legend = c('Confirmed captures', 'Background'), pch = c(21,4),
           pt.bg = c('green',NA), col = 'black')
    dev.off()
    
    return(classi.dat)
} ## End function



