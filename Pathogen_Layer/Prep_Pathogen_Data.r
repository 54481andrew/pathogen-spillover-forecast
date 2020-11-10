## This script processes Lassa presence/absence data from published surveys. The
## outputs are two dataframes, one of human serosurvey data, and another
## of rodent Lassa presence/absence data. The rodent dataframe contains
## environmental predictors that are used to train the pathogen layer of
## the model.
## Certain hyper-parameters defined in Generate_Pathogen_Layer.r determine
## how the data will be processed. Specifically, [min.test] determines the
## number of rodents that need to be tested to call a site LASV-; if [min.test]
## rodents are tested and none have LASV exposure, the pixel is classified as
## LASV-. [set.ambiguous.to] determines how pixels without virus detection, but
## with seropositive rodents, are handled. The natural choices are NA (omitted)
## or 1 (counted as LASV+). 

Prep.Pathogen.Data <- function(hypers.i){
    
    ## Load in Lassa survey data-set
    full.dat <- read.csv(file = 'Data/Cleaned_Lassa_Literature.csv', stringsAsFactors = FALSE,
                         sep = '\t')
    
    ## Combine Genus and Species columns into single column Species
    full.dat$Species = with(full.dat, paste(Genus, Species, sep = ' '))
    
    ## Calculate proportion infected and proportion seropositive for all entries
    full.dat$PropVirus <- full.dat$NumPosVirus/full.dat$NumTestVirus
    full.dat$PropAb <- full.dat$NumPosAb/full.dat$NumTestAb

    ## The human entries in the data-set all come from sero-surveys that
    ## sample a random human population
    full.dat$Human_Random_Survey = full.dat$Species== 'Homo sapiens'

    ## Keep only those columns that are needed, and discard any entries
    ## without Lat/Lon values. This dataset is called classi.dat, short
    ## for classification dataset.
    full.dat <- full.dat[,c('Longitude', 'Latitude', 'Country', 'Village',
                            'Species', 'NumPosVirus', 'NumTestVirus', 'PropVirus',
                            'Virus_Diagnostic_Method',
                            'NumPosAb', 'NumTestAb', 'PropAb',
                            'Ab_Diagnostic_Method', 'Antibody_Target',
                            'Source',
                            'Human_Random_Survey',
                            'Year', 'Bibtex')]
                            
    keep <- !(is.na(full.dat$Latitude) | is.na(full.dat$Longitude))
    classi.dat <- full.dat[keep,]
    
    ## Add columns describing total individuals tested (TotTest), and total
    ## number of individuals that were arenavirus positive (TotPos).
    classi.dat = add.total.columns(classi.dat, hypers.i$min.test)

    ## Specify that these data are not from genbank
    classi.dat$genbank <- FALSE
    
    ## -- Incorporate GenBank dataset into classi.dat
    
    ## Load GenBank data
    genbank.dat <- read.csv('Data/Rodents_Genbank_Oct_2020.csv', stringsAsFactors = FALSE,
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
                                        Species = Species, NumPosVirus = 1,
                                        NumTestVirus = 1, PropVirus = NA,
                                        Virus_Diagnostic_Method = 'PCR',
                                        NumPosAb = 0, NumTestAb = 0,
                                        PropAb = NA,
                                        Ab_Diagnostic_Method = NA,
                                        Antibody_Target = NA,
                                        Source = Reference,
                                        Human_Random_Survey= FALSE,
                                        Year = gbCollectYear,
                                        Bibtex = NA,
                                        TotPos = 1, TotTest = 1,
                                        genbank = TRUE
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
    ## Add ArenaStat column based on serology. This is used to define LASV status
    ## in the rodents. For the human data, this is only used for graphing. 
    human.test.dat$ArenaStat <- 1*(human.test.dat$NumPosAb > 0)

    
    ## Aggregate data that fall within unique pixels. This is
    ## achieved with the purge.repeats function. The pixelation of
    ## West Africa is determined by the 0.05x0.05 degree MODIS precipitation
    ## data-set, which is a raster inside of the all.stack raster-stack.
    template.rast <- all.stack[['Pmu']]
    purged.out <- purge.repeats(classi.dat.jrod, template.rast)
    jrod.purged <- purged.out[[1]]

    ## -----------------
    
    ## Keep track of the number of pixels that are 
    ## ambiguous given our definition of LASV +. Column names refer to the
    ## literature data except for those with suffix y refers to genbank.
    ## These are pixels in which not enough rodents were tested, or were serology-only.
    classi.dat.jrod$Cell <- raster::extract(template.rast,
                                            classi.dat.jrod[,c('Longitude','Latitude')],
                                            cellnumbers=TRUE)[,'cells']
    rod.lit <- classi.dat.jrod[classi.dat.jrod$genbank==FALSE,]
    rod.gen <- classi.dat.jrod[classi.dat.jrod$genbank==TRUE,]
    rod.gen <- aggregate(NumPosVirus~Cell,rod.gen, FUN = sum)
    rod.lit.ambiguous <- rod.lit[(rod.lit$NumPosAb > 0 &
                                            rod.lit$NumPosVirus==0),]

    ## Merge genbank information into those entries that are ambiguous
    merged.ambiguous <- merge(rod.lit.ambiguous, rod.gen[,c('NumPosVirus','Cell')], by = 'Cell',
                              all.x = TRUE, all.y = FALSE)
    merged.ambiguous <- merged.ambiguous[,c('Cell','Source', 'Village', 'Country',
                                            'NumPosVirus.x', 'NumTestVirus',
                                            'NumPosAb', 'NumTestAb',
                                            'NumPosVirus.y')]
    merged.ambiguous.data <- merged.ambiguous[order(merged.ambiguous$NumPosVirus.y),]

    ## Pixels for which no genbank data was found have NumPosVirus.y= NA; set these to 0
    merged.ambiguous.data[is.na(merged.ambiguous.data$NumPosVirus.y),'NumPosVirus.y'] <- 0
    ## Points with confirmed LASV virus from genbank are not ambiguous, so remove them
    ambiguous.data <- merged.ambiguous.data[merged.ambiguous.data$NumPosVirus.y==0,]

    ## -----------------
    
    ## Define sites as positive or negative or NA for arenavirus (indicated in ArenaStat column).
    
    ## Three kinds of sites are negative (meet minimum # tested, all tested are negative),
    ## positive (1 or more rodent with Lassa virus detected), and ambiguous (positive serology
    ## with no virus explicitly detected). Ambiguous sites are set according to
    ## set.ambiguous.sites option in hypers.dat
    jrod.purged$ArenaStat <- NA
    mask.neg <- jrod.purged$TotTest >= hypers.i$min.test & jrod.purged$NumPosVirus ==0 &
        jrod.purged$NumPosAb ==0
    mask.pos <- jrod.purged$NumPosVirus > 0
    mask.ambi <- jrod.purged$NumPosAb > 0 & jrod.purged$NumPosVirus==0
        
    jrod.purged[mask.pos, 'ArenaStat'] <- 1
    jrod.purged[mask.neg, 'ArenaStat'] <- 0
    jrod.purged[mask.ambi, 'ArenaStat'] <- hypers.i$set.ambiguous.to

    ## Save ambiguous pixels for debugging purposes
    ambiguous.pixels <- jrod.purged[mask.ambi,c()]
    
    ## Now that ambiguous entries have been set to what the user wants, exclude any entries
    ## from jrod.purged that still have NA ArenaStat
    mask.na <- is.na(jrod.purged$ArenaStat)
    jrod.purged <- jrod.purged[!mask.na,]
    
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
    write.csv(jrod.purged, file = paste0('Figures_Fits/', prefix, '/',fold,'/',
                                         'Prepped_Pathogen_PresAbs_Data.csv'),
              row.names = FALSE)
    write.csv(ambiguous.data, file = paste0('Figures_Fits/', prefix, '/',fold,'/',
                                            'Prepped_Pathogen_Ambiguous_Data.csv'),
              row.names = FALSE)
    write.csv(ambiguous.pixels, file = paste0('Figures_Fits/', prefix, '/',fold,'/',
                                              'Prepped_Pathogen_Ambiguous_Pixels.csv'),
              row.names = FALSE)
    write.csv(human.test.dat, file = paste0('Figures_Fits/', prefix, '/',fold,'/',
                                            'Prepped_Human_Seroprevalence_Data.csv'),
              row.names = FALSE)
    
    cat('--Lassa Rodent Count Statistics--')
    tab <- table(jrod.purged$ArenaStat)
    print(addmargins(tab, FUN = list(Total = sum), quiet = TRUE))
    
    ## Number of seroprevalence studies of humans
    print(paste('Number of Human Sources: ', length(unique(human.test.dat$Source))), quote = FALSE)
    
    ## Number of countries
    print(paste('Number of Countries: ', length(unique(human.test.dat$Country))), quote = FALSE)
    
    ## Time span
    temp <- with(human.test.dat, strsplit(paste(Year), '-'))
    years = sapply(temp, FUN = function(x){range(as.numeric(x))})
    hum.yrange <- paste(range(years), collapse = '-')
    print(paste0('Human seroprevalence collected between ', hum.yrange),
          quote = FALSE)

    ## Store some of the printed output to file
    nneg <- sum(jrod.purged$ArenaStat==0)
    npos <- sum(jrod.purged$ArenaStat==1)
    info <- paste0('Rodent: \n', 
                   'Number Absences: ', nneg, '\n', 
                   'Number Presences: ', npos, '\n\n',
                   'Human: \n',
                   'Number sites: ', nrow(human.test.dat), '\n',
                   'Number studies: ', length(unique(human.test.dat$Source)), '\n',
                   'Year range: ', hum.yrange)
    write(info, paste0('Figures_Fits/', prefix, '/',fold,'/','data_prep_info.txt'))

    png(file = paste0('Figures_Fits/', prefix, '/',fold,'/',
                      'Human_Test_Dat.png'), width = 6, height = 4, units = 'in', res = 400)
    par(mai = 1*c(0.2,0.2,0.2,0.2))
    plot(foc.shp.ogr, col = 'cornsilk')
    points(Latitude~Longitude, human.test.dat, bg = 'green', col = 'black', pch = 21, cex = 0.5)
    legend(x = 'topleft', legend = 'Human LASV sero-survey', pch = 21, pt.bg = 'green', col = 'black')
    dev.off()
    
    return(list(jrod.purged, ambiguous.data, ambiguous.pixels))
} ## End function    
