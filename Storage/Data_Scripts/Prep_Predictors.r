## This script joins all predictor rasters into a raster stack, and writes the
## stack to a file.

## Packages for handling rasters and shapefiles
require(raster)
require(rgdal)
require(sf)

## set version of random number generator to use to ensure
## reproducibility across R versions
RNGversion('3.6.0')
set.seed(1)

## Storage locations
shapefile.storage.fold <- '../Shapefiles'
raster.data.storage.fold ='../Raster_Data'

## Load in environmental rasters that describe seasonality

## - Read in means of Precipitation, temperature, and NDVI
stat <- raster(paste(raster.data.storage.fold,'/Rast_Pmu.tif',sep=''))
Pmu.rast <- stat
stat <- raster(paste(raster.data.storage.fold,'/Rast_Tmu.tif',sep=''))
Tmu.rast <- stat
stat <- raster(paste(raster.data.storage.fold,'/Rast_Nmu.tif',sep=''))
Nmu.rast <- stat

## - Colwell's index for precipitation (Pcol) and NDVI (Ncol).
## Here, Colwell's indices data is a stack of 3 layers:
## Layers 1,2, and 3 refer to constancy, contingency, and predictability,
## respectively. We only use constancy and contingency. Temperature
## seasonality predictors are omitted due to NA values.
stat <- stack(paste(raster.data.storage.fold,'/Rast_Pcol.grd',sep=''))
Pcol.stack <- stat
stat <- stack(paste(raster.data.storage.fold,'/Rast_Ncol.grd',sep=''))
Ncol.stack <- stat

## - Coefficient of variation (Pcv and Ncv).
stat <- raster(paste(raster.data.storage.fold,'/Rast_Pcv.tif',sep=''))
Pcv.rast <- stat
stat <- raster(paste(raster.data.storage.fold,'/Rast_Ncv.tif',sep=''))
Ncv.rast <- stat

## - Duration of dry (Pdur) and brown (Ndur).
## Defined as the duration of precipitation that is less than a daily
## average of 1 mm/day, and of NDVI that is less than 0.5. These
## thresholds are somewhat arbitrary, but capture the spirit of the
## measurement.
stat <- raster(paste(raster.data.storage.fold,'/Rast_Pdur.tif',sep=''))
Pdur.rast <- stat
stat <- raster(paste(raster.data.storage.fold,'/Rast_Ndur.tif',sep=''))
Ndur.rast <- stat

## - Minimum and maximum precipitation and NDVI for each region.
## Precipitation min and max are in a stack (Pmm) that has layers
## Pmin and Pmax. NDVI is similar.
stat <- stack(paste(raster.data.storage.fold,'/Rast_Pminmax.grd',sep=''))
Pmm.stack <- stat
stat <- stack(paste(raster.data.storage.fold,'/Rast_Nminmax.grd',sep=''))
Nmm.stack <- stat

## -----

## - Elevation data
elev.rast <- raster(paste(raster.data.storage.fold,"Rast_Elevation_gt30w020n40.tif", sep = '/'))
elev.rast <- projectRaster(elev.rast, Pmu.rast)

## - WorldPop 2020 population data
## For the predictor dataset, we first recast this as density (per km2), then
## project to the appropriate coordinate system.
pop.rast <- raster(paste(raster.data.storage.fold, 'Rast_Population_AFR_PPP_2020_adj_v2.tif', sep = '/'))
pop.rast <- crop(pop.rast, Pmu.rast)
## Convert to density, then reproject
pop.rast <- pop.rast / raster::area(pop.rast)
pop.rast <- projectRaster(pop.rast, to = Pmu.rast)

## - MODIS Landcover data-set
## Form stack for landcover layers. Only use those that exist for an average
## time of at least 20 years (see process_MODIS_landcover.r script for details
file <- paste(raster.data.storage.fold, '/Processed_MODIS_LandCover/2001/LC_2001.grd', sep='')
hab.stack <- stack(file)
use.hab.names <- read.table('Output/landcover_variable_names_to_use',
                            header = TRUE)
hab.stack <- hab.stack[[paste(use.hab.names$class_name)]]


## Combine all prediction rasters into a stack
all.stack <- stack(list(elev.rast, Tmu.rast, Pmu.rast, Nmu.rast,
                        Pcv.rast, Ncv.rast,
                        Pcol.stack[[1:2]],
                        Ncol.stack[[1:2]],
                        Pmm.stack,
                        Nmm.stack,
                        Pdur.rast, Ndur.rast,
                        pop.rast,  hab.stack
                        ))

## Assign names to each layer in the stack
names(all.stack) <- c('Elev', 'Tmu', 'Pmu', 'Nmu',
                      'Pcv', 'Ncv',
                      'Pc', 'Pm',
                      'Nc', 'Nm',
                      names(Pmm.stack),
                      names(Nmm.stack),
                      'Pdur', 'Ndur',
                      'Pop', names(hab.stack)
                      )

## Save all.stack to Storage directory
writeRaster(all.stack,paste(raster.data.storage.fold,
                            "/predictor_stack.grd", sep = ''),
            format="raster",
            overwrite = TRUE)

## Write shapefile containing only countries within the study region. This will be
## created from a larger shapefile that contains all countries in Africa.
foc.countries <- c('Mali', "Guinea", "Ivory Coast", "Sierra Leone",
                   "Nigeria", 'Liberia', "Cote d`Ivoire",
                   'Ghana', 'Togo', 'Benin',
                   'Burkina Faso', 'Mauritania', 'Niger', 'Senegal',
                   'Gambia', 'Guinea-Bissau')

## Africa shapefile source: http://www.maplibrary.org/library/stacks/Africa/index.htm
all.africa.shp = st_read(dsn = paste(shapefile.storage.fold,'/Africa/',sep=''),
                          layer = 'Africa')
## Resave with crs info added
st_crs(all.africa.shp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
st_write(all.africa.shp, dsn = paste(shapefile.storage.fold,'/Africa/',sep=''),
         layer = 'Africa', driver = 'ESRI Shapefile', delete_dsn = TRUE)

adf <- as.data.frame(all.africa.shp)
country.ord = adf$COUNTRY
foc.shp <- all.africa.shp[country.ord %in% foc.countries,4]

## Write West Africa shapefile
wa.shapefile.fold <- paste(shapefile.storage.fold, '/West_Africa', sep = '')
dir.create(wa.shapefile.fold, showWarnings = FALSE)
st_write(foc.shp, dsn = wa.shapefile.fold, layer = 'foc',
         driver = 'ESRI Shapefile',
         delete_dsn = TRUE)

