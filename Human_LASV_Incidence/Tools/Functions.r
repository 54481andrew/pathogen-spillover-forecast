## Generate a folder name for pathogen fit
generate.lassa.name <- function(gridgi){
    fold.name <- with(gridgi, paste('pa_nboots', nboots,
                                    'tc', tree.complexity,
                                    'mllr', mllr,
                                    'lmt', lmt,
                                    sep='_')) ## Folder suffix for data storage and results
    return(fold.name)
}

## Generate a folder name for reservoir fit
generate.masto.name <- function(gridgi){
    fold.name <- with(gridgi, paste('pa_nboots',nboots,
                                    'tc', tree.complexity,
                                    'mllr', mllr,
                                    'lmt', lmt,
                                    sep = '_'))
    fold.name <- paste('Mn_', fold.name, sep = '') ## Data
    return(fold.name)
}


## Given raster input that describes the number of new human
## LASV cases in each pixel, this function returns the summed total
## cases in each West African country.
get.country.sums <- function(human.lasv.cases){

    ## Subset africa.ogr to West African countries
    ss.africa.ogr <- africa.ogr[africa.ogr$COUNTRY %in% foc.countries,]
    
    ## Extract and sum human LASV cases within each region
    totcase.dat <- extract(human.lasv.cases, ss.africa.ogr)
    totcase <- lapply(totcase.dat, FUN=sum, na.rm = TRUE)
    totcase <- unlist(totcase)

    ## Organize sums of each region into a dataframe, with country column. Then
    ## sum all the regions within each country.
    case.df <- data.frame(Country = ss.africa.ogr$COUNTRY, Cases = totcase)
    case.df.agg <- aggregate(. ~  Country, data = case.df, sum)

    ## Extract and sum total population within each region
    totpop.dat <- extract(pop.rast, ss.africa.ogr)
    totpop <- lapply(totpop.dat, FUN=sum, na.rm = TRUE)
    totpop <- unlist(totpop)

    ## Organize sums of each region into a dataframe
    pop.df <- data.frame(Country = ss.africa.ogr$COUNTRY, Pop = totpop)
    pop.df.agg <- aggregate(. ~  Country, data = pop.df, sum)

    ## Combine dataframes into one, organized by country row. Pare this list down to focal countries
    wi.foc.countries <- case.df.agg$Country%in% foc.countries
    case.df.foc <- case.df.agg[wi.foc.countries,]
    case.df.foc$Pop = pop.df.agg[wi.foc.countries,'Pop']
    case.df.foc$Country = paste(case.df.foc$Country) ## Recast country name as string

    ## Calculate total cases
    total.cases <- cellStats(human.lasv.cases, sum)

    ## Add a row with totals
    case.df.foc <- rbind(case.df.foc,
                         data.frame(Country = 'Total', Cases = total.cases,
                                    Pop = sum(case.df.foc$Pop)))
    case.df.foc$Pop <- round(case.df.foc$Pop)
    case.df.foc$Cases <- case.df.foc$Cases

    ## Keep only certain columns, and order rows by number of cases
    case.df.foc <- case.df.foc[,c('Country', 'Cases', 'Pop')]
    case.df.foc <- case.df.foc[order(case.df.foc$Cases, decreasing = TRUE),]

    
    ## Return case data
    return(case.df.foc)
}## End Function

