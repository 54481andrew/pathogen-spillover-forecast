
## Calculates deviance-based pseudo r-squared
## from Cameron, Windmeijer, 1996
pseudor2.fun <- function(obj){
    a = summary(obj)
    1 - a$deviance/a$null.deviance
}


## Generate a folder name for tree with specified properties
generate.lassa.name <- function(gridgi){
    fold.name <- with(gridgi, paste('pa_nboots', nboots,
                                    'tc', tree.complexity,
                                    'mllr', mllr,
                                    'lmt', lmt,
                                    sep='_')) ## Folder suffix for data storage and results
    return(fold.name)
}
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

    ## Countries over which sums will be calculated
    foc.countries <- c('Mali', "Guinea", "Ivory Coast", "Sierra Leone",
                       "Nigeria", 'Liberia', "Cote d`Ivoire",
                       'Ghana', 'Togo', 'Benin',
                       'Burkina Faso', 'Mauritania', 'Niger',
                       'Senegal', 'Guinea-Bissau')

    ## Extract and sum human LASV cases within each region
    totcase.dat <- extract(human.lasv.cases, africa.ogr)
    totcase <- lapply(totcase.dat, FUN=sum, na.rm = TRUE)
    totcase <- unlist(totcase)

    ## Organize sums of each region into a dataframe, with country column. Then
    ## sum all the regions within each country.
    case.df <- data.frame(Country = africa.ogr$COUNTRY, Cases = totcase)
    case.df.agg <- aggregate(. ~  Country, data = case.df, sum)

    ## Extract and sum total population within each region
    totpop.dat <- extract(pop.rast, africa.ogr)
    totpop <- lapply(totpop.dat, FUN=sum, na.rm = TRUE)
    totpop <- unlist(totpop)

    ## Organize sums of each region into a dataframe
    pop.df <- data.frame(Country = africa.ogr$COUNTRY, Pop = totpop)
    pop.df.agg <- aggregate(. ~  Country, data = pop.df, sum)

    ## Combine dataframes into one, organized by country row. Pare this list down to focal countries
    wi.foc.countries <- case.df.agg$Country%in% foc.countries
    case.df.foc <- case.df.agg[wi.foc.countries,]
    case.df.foc$Pop = pop.df.agg[wi.foc.countries,'Pop']
    case.df.foc$Country = paste(case.df.foc$Country) ## Recast country name as string

    ## Add new columns of country population size and rate of infection
    case.df.foc$Rate = with(case.df.foc, Cases/Pop*1000)

    ## Calculate total cases
    total.cases <- cellStats(human.lasv.cases, sum)

    ## Add a row with totals
    case.df.foc <- rbind(case.df.foc,
                         data.frame(Country = 'Total', Cases = total.cases,
                                    Pop = sum(case.df.foc$Pop),
                                    Rate = NA))
    case.df.foc$Pop <- round(case.df.foc$Pop)
    case.df.foc$Cases <- case.df.foc$Cases/1000
    case.df.foc$Rate <- round(case.df.foc$Rate, 1)

    ## Keep only certain columns, and order rows by number of cases
    case.df.foc <- case.df.foc[,c('Country', 'Cases', 'Rate', 'Pop')]
    case.df.foc <- case.df.foc[order(case.df.foc$Cases, decreasing = TRUE),]

    ## Return case data
    return(case.df.foc)
}## End Function

## ---
## Log-likelihood
LL = function(P,D, tol = 1e-8){
  mask.notna = which(!is.na(P) & !is.na(D))
  P = P[mask.notna]
  P[P < tol] = tol
  D = D[mask.notna]
  ll = 0
    for(i in 1:length(P)){
    ll = ll + -P[i] + D[i]*log(P[i]) - sum(log(1:max(1,D[i])))
      }
  return(ll)
  }

## ---
## Deviance
Dev = function(P,D, tol = 1e-8){
  mask.notna = which(!is.na(P) & !is.na(D))
  P = P[mask.notna]
  P[P < tol] = tol
  D = D[mask.notna]
  dev = 0
    for(i in 1:length(P)){
      devi = ifelse(D[i]==0, 2*P[i],  2*(D[i]*log(D[i]/P[i]) - (D[i]-P[i])))
      dev = dev + devi
  }
  return(dev)
  }
## ---
