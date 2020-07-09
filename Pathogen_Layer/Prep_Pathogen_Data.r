## This script processes Lassa presence/absence data from published surveys. The
## outputs are two dataframes, one of human serosurvey data, and another
## of rodent Lassa presence/absence data. The rodent dataframe contains
## environmental predictors that are used to train the pathogen layer of
## the model.

Prep.Pathogen.Data <- function(){

    ## Load in Lassa survey data-set
    full.dat <- read.csv(file = 'Data/Lassa_Occurrence_Data_v10_WB.csv', stringsAsFactors = FALSE,
                         sep = '\t')
    
    ## combine Genus and Species columns into single column Species (v10)
    full.dat$Species = with(full.dat, paste(Genus, Species, sep = ' '))
        
    ## Calculate proportion infected and proportion seropositive for all entries
    full.dat$PropAg <- full.dat$NumPosAg/full.dat$NumTestAg
    full.dat$PropAb <- full.dat$NumPosAb/full.dat$NumTestAb
    
    ## Keep only those columns that are needed, and discard any entries
    ## without Lat/Lon values. This dataset is called classi.dat, short
    ## for classification dataset.
    full.dat <- full.dat[,c('Longitude', 'Latitude', 'Country', 'Village',
                            'Species', 'NumPosAg', 'NumTestAg', 'PropAg',
                            'NumPosAb', 'NumTestAb', 'PropAb', 'Source',
                            'DiagnosticMethod', 'Target', 'Human_Random_Survey')]
    keep <- !(is.na(full.dat$Latitude) | is.na(full.dat$Longitude))
    classi.dat <- full.dat[keep,]
    
    ## Add columns describing total individuals tested (TotTest), total
    ## number of individuals that were arenavirus positive (TotPos),
    ## and binary ArenaStat column indicating whether Lassa virus
    ## was found (1) or not (0)
    classi.dat = add.total.columns(classi.dat)
    
    ## If the ArenaStat column is neither positive or negative, there is
    ## no information. Remove rows with no information.
    classi.dat <- classi.dat[which(!is.na(classi.dat$ArenaStat |
                                          classi.dat$Species=='Homo sapiens')),]
    
    ## -- Incorporate GenBank dataset into classi.dat
    
    ## Load GenBank data
    genbank.dat <- read.csv('Data/Rodents_Genbank_June16_2020.csv', stringsAsFactors = FALSE,
                         sep = '\t')
    
    ## Prune the dataset to those entries that are of Mastomys natalensis, and
    ## for which latitude and longitude are available.
    wi.genbank <- (genbank.dat$gbHost=='Mastomys natalensis') &
        !is.na(genbank.dat$Lat) & !is.na(genbank.dat$Long)
    genbank.dat <- genbank.dat[wi.genbank, ]
    genbank.dat$Species <- 'Mastomys natalensis'
    
    ## Reformat to be compatible with historical survey data
    genbank.dat.form <- with(genbank.dat,
                             data.frame(Longitude = Long, Latitude = Lat,
                                        Country = Country, Village = LocVillage,
                                        Species = Species, NumPosAg = 1,
                                        NumTestAg = 1, PropAg = NA,
                                        NumPosAb = 0, NumTestAb = 0,
                                        PropAb = NA, Source = Reference,
                                        DiagnosticMethod = 'Various', Target = 'Ag',
                                        Human_Random_Survey= FALSE,
                                        ArenaStat = 1, TotPos = 1, TotTest = 1
                                        ))
    
    ## Add into classi.dat
    classi.dat <- rbind(classi.dat, genbank.dat.form)

    ## Restrict dataset to those points that occur in the focal countries of West Africa
    foc.countries <- c('Mali', "Guinea", "Ivory Coast", "Sierra Leone",
                       "Nigeria", 'Liberia', "Cote d`Ivoire",
                       'Ghana', 'Togo', 'Benin',
                       'Burkina Faso', 'Mauritania', 'Niger', 'Senegal',
                       'Gambia', 'Guinea-Bissau')
    
    ## Remove entries that lie outside the study region
    classi.dat <- classi.dat[classi.dat$Country %in% foc.countries,]

    
    ## Split the dataset into two:
    ## 1. Just rodent data (used for presence/absence)
    ## 2. Just human data (used for seroprevalence data)
    classi.dat.jrod <- classi.dat[classi.dat$Species=='Mastomys natalensis', ]
    human.test.dat <-  classi.dat[classi.dat$Human_Random_Survey==TRUE &
                                  classi.dat$NumTestAb > 0, ]
    human.test.dat$ArenaStat <- 1*(human.test.dat$NumPosAb > 0)
    
    
    ## Aggregate data that fall within unique pixels. This is
    ## achieved with the purge.repeats function. The pixelation of
    ## West Africa is determined by the 0.05x0.05 degree MODIS precipitation
    ## data-set, which is a raster inside of the all.stack raster-stack.
    template.rast <- all.stack[['Pmu']]
    purged.out <- purge.repeats(classi.dat.jrod, template.rast)
    jrod.purged <- purged.out[[1]]
    
    ## Change coordinates to be at cell center
    orig.points <- cbind(jrod.purged$Longitude, jrod.purged$Latitude)
    cells <- raster::extract(template.rast, orig.points, cellnumbers=TRUE)[,'cells']
    cell.points <- data.frame(xyFromCell(object = template.rast, cell = cells))
    jrod.purged$Longitude <- cell.points$x
    jrod.purged$Latitude <- cell.points$y
    
    ## Incorporate the environmental predictors into the
    ## dataset. All predictors are extracted from all.stack.
    points <- cbind(jrod.purged$Longitude, jrod.purged$Latitude)
    jrod.purged <- cbind(jrod.purged, raster::extract( all.stack, points, sp = TRUE))
    
    ## Impute any missing values with median
    var.names <- names(all.stack)
    jrod.purged <- impute(jrod.purged, var.names)
    
    ## Save a copy of rodent presence/absence dataset, as well as the human seroprevalence data
    write.csv(jrod.purged, file = 'Data/Prepped_Pathogen_PresAbs_Data.csv', row.names = FALSE)
    write.csv(human.test.dat, file = 'Data/Prepped_Human_Seroprevalence_Data.csv', row.names = FALSE)
    
    cat('--Lassa Rodent Count Statistics--')
    tab <- table(jrod.purged$ArenaStat)
    print(addmargins(tab, FUN = list(Total = sum), quiet = TRUE))
    
    ## Number of seroprevalence studies of humans
    print(paste('Number of Human Sources: ', length(unique(human.test.dat$Source))), quote = FALSE)
    
    ## Number of countries
    print(paste('Number of Countries: ', length(unique(human.test.dat$Country))), quote = FALSE)
    
    ## Time span
    temp <- with(human.test.dat, strsplit(paste(Source), '_'))
    years = sapply(temp, FUN = function(x){x[length(x)]})
    years = as.numeric(years)
    years = years[!is.na(years)] ## Filter out entries that are not years
    print(paste0('Human seroprevalence collected between ', paste(range(years), collapse = '-')), quote = FALSE)
    
    png(file = paste0('Figures_Fits/Human_Test_Dat.png'), width = 6, height = 4, units = 'in', res = 400)
    par(mai = 1*c(0.2,0.2,0.2,0.2))
    plot(foc.shp.ogr, col = 'cornsilk')
    points(Latitude~Longitude, human.test.dat, bg = 'green', col = 'black', pch = 21, cex = 0.5)
    legend(x = 'topleft', legend = 'Human LASV sero-survey', pch = 21, pt.bg = 'green', col = 'black')
    dev.off()
    
    ## --Lassa Rodent Count Statistics--
    ##     0     1 Total
    ##    36    38    74
    
    
    ##"Human seroprevalence collected between 1973-2019"

    return(jrod.purged)
} ## End function    
