## This script reads in the tif outputs from Generate_Reservoir_Layer.r
## and Generate_Pathogen_Layer.r, and, using data from human LASV
## serosurveys, outputs a prediction of human LASV seroprevalence and
## incidence across West Africa. This prediction is made in three steps:

## 1. Let $D_M$ and $D_L$ denote the reservoir and pathogen risk layers
## outputed Generate_*_Layer.r scripts above. Calculate the combined risk
## raster, $D_X$ as
## $$
## D_X = D_M \times D_L.
## $$
## 2. Regress the human seroprevalence data onto the combined risk layer
##    $D_X$. Use this regression to predict the prevalence of LASV across
##    West Africa.
## 3. Use equations derived in the Appendix to convert the predicted
##    seroprevalence into predictions of LASV incidence.

## Load in packages for GIS, regression, and plotting
require(raster) ## GIS
require(sf) ## GIS
require(rgdal) ## rgeos
require(MASS) ## glm.nb
require(DescTools) ## PseudoR2
require(pscl) ## odTest
require(fields) ## image.plot
require(viridis) ## color-blind friendly color scheme

## set version of random number generator to use to ensure
## reproducibility across R versions
RNGversion('3.6.0')
set.seed(1)

## grid.name specifies the location of the directory to which the regression
## results are written. Prefix names specify where the target .tif files that
## describe reservoir and pathogen risk are located.
grid.name <- 'human_v2' ## name for model.dat dataframe
prefix.masto <- 'reservoir_v2' ## prefix used in Generate_Reservoir_Layer
prefix.lassa <-  'pathogen_v2' ## prefix used in Generate_Lassa_Layer

## Load additional functions
source("Tools/Functions.r")

## --- Load in human data
human.test.dat <- read.csv('../Pathogen_Layer/Data/Prepped_Human_Seroprevalence_Data.csv')

## Parameters for calculating incidence of human LASV from human LASV seroprevalence
## Time units are years.
d = 1/50 ## Life span of 50 years
del = 12 ## LASV recovery is 1 month
lam.mcc87 = 0.064 ## from McCormick, Webb, Krebs, 1987

## Load in Africa shapefiles needed for plotting
storage.fold <- '../Storage'
foc.shp.ogr <- readOGR(dsn = paste(storage.fold, '/Shapefiles/West_Africa', sep = ''),
                       layer = 'foc', verbose = FALSE)
africa.ogr <- readOGR(dsn = paste(storage.fold, '/Shapefiles/Africa', sep = ''),
                       layer = 'Africa', verbose = FALSE)
masto.rangemap <-readOGR(dsn = paste(storage.fold, '/Shapefiles/Masto_Range',sep=''),
                         layer = 'data_0', verbose = FALSE)


## ----- Null prediction of human LASV incidence

## Calculate null prediction of new human cases per year. This treats West Africa
## as one population, with seroprevalence equal to average seroprevalence in the
## human.test.dat.
## Calculate the weighted mean seroprevalence across all sero-surveys
wmean <- sum(human.test.dat$NumPosAb)/sum(human.test.dat$NumTestAb)

## Next, use equations derived in the mathematica file to convert seroprevalence
## into human incidence. First step is to load in human population data,
## which is in compressed format.
unzip('../Storage/Raster_Data/Rast_Population_AFR_PPP_2020_adj_v2.tif.zip',
      exdir = '../Storage/Raster_Data')
pop.rast <- raster('../Storage/Raster_Data/Rast_Population_AFR_PPP_2020_adj_v2.tif')
pop.rast <- crop(pop.rast, foc.shp.ogr)
pop.rast <- mask(pop.rast, foc.shp.ogr)
pop.rast <- mask(pop.rast, masto.rangemap, updatevalue = 0)
unlink('../Storage/Raster_Data/Rast_Population_AFR_PPP_2020_adj_v2.tif')
unlink('../Storage/Raster_Data/__MACOSX', recursive = TRUE)

## Get the total population of the study region in West Africa
tot.pop = cellStats(pop.rast, sum)
## 374,279,977

## Estimate of number of human LASV cases without seroreversion
lam = 0
null.low.estimate = wmean*tot.pop*(d + del)*(d + lam)*del^-1
## Estimate, but with seroreversion
lam = lam.mcc87
null.high.estimate = wmean*tot.pop*(d + del)*(d + lam)*del^-1
## If every location has seroprevalence at 0.184 (wmean), total cases
## is 1,380,375 - 5,797,574, depending on lam in the interval [0, 0.064]
##
## -----

## Load Mastomys natalensis prediction
masto.fold <- 'Mn_pa_nboots_25_nbg_Same_tc_1_mllr_2_lmt_7'
filename.masto <- paste0('../Reservoir_Layer/Figures_Fits/', prefix.masto, '/',
                         masto.fold, '/', "Reservoir_Layer_",masto.fold, '.tif')
Masto.Risk <- raster(filename.masto)
Masto.Risk <- mask(Masto.Risk, masto.rangemap, updatevalue = 0)

## Load in Lassa prediction
lassa.fold <- 'pa_nboots_25_1_mllr_3_lmt_7'
filename.lassa <- paste0('../Pathogen_Layer/Figures_Fits/', prefix.lassa, '/',
                         lassa.fold, '/', 'Lassa_Layer_', lassa.fold,".tif")
Lassa.Risk <- raster(filename.lassa)
Lassa.Risk <- mask(Lassa.Risk, masto.rangemap, updatevalue = 0)


## Plot product, D_X = D_L \times D_M, that describes the combined risk of LASV
Dx = Lassa.Risk*Masto.Risk
heat.cols <- viridis(120, begin = 0.1, end = 1, option = 'D')
xlims = c(-18,16)
ylims = c(16, 16.5)
png(file = paste('Figures_Fits/Product_Risk_Layer.png', sep = ''),
    width = 6, height = 4, units = 'in', res = 400)
par(mai = 1*c(0.2,0.2,0.2,0.6))
image.plot(Dx, col = heat.cols, zlim = c(0, 1),
           bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',
           xlim = xlims, ylim = ylims,
           asp = 1, legend.lab = 'Combined score', legend.line = 2.5,
           main = '')
mtext(text = as.expression('Combined Risk ('~'D'['X']~' = D'['M'] %.% 'D'['L']~')'), side = 3, line = -1)
plot(foc.shp.ogr, add = TRUE, bty = 'n', asp = 1)
plot(rgeos::gIntersection(foc.shp.ogr, masto.rangemap), add = TRUE, bty = 'n', asp = 1, lwd = 3)
dev.off()

## ---- Regress human seroprevalence on D_X layer

## Extract predictions from each layer at each human sero-survey location
human.test.dat$Dx <- extract(Dx, human.test.dat[,c('Longitude', 'Latitude')])

## Fit a standard binomial as a baseline
bin.mod <- glm(PropAb~ 1 + Dx, human.test.dat,
              family = 'binomial',
              weights = NumTestAb)
summary(bin.mod)
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -2.01988    0.04950  -40.81   <2e-16 ***
## Dx           0.99637    0.08419   11.83   <2e-16 ***


## Fit quasi-binomial model
qb.mod <- glm(PropAb~ 1 + Dx, human.test.dat,
              family = 'quasibinomial',
              weights = NumTestAb)
summary(qb.mod)
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  -2.0199     0.2014 -10.029  < 2e-16 ***
## Dx            0.9964     0.3426   2.908  0.00449 ** 


## Calculate correlation
cor(qb.mod$fitted.values, qb.mod$y)
## 0.23
cov.wt(cbind(qb.mod$fitted.values, qb.mod$y), wt = human.test.dat$NumTestAb,
                        cor = TRUE)$cor
## 0.29


## --- Create location in which images will be saved

fold.name <- paste('Figures_Fits/', grid.name, sep = '')
dir.create(fold.name, showWarnings = FALSE)


## --- Use regression model to predict human LASV seroprevalence across West Africa
pred.stack <- stack(Dx)
names(pred.stack) = 'Dx'
pred.rast <- raster::predict(object = pred.stack, model = qb.mod, type = 'response')

## --- Plot predicted seroprevalence vs actual seroprevalence, and also plot
## predicted seroprevalence across West Africa
human.test.dat$fitted <- extract(x = pred.rast,
                          y = human.test.dat[,c('Longitude', 'Latitude')])

sc = 0.01 ## Scaling of points
png(filename = paste(fold.name,'/Pred_vs_Fit.png',sep=''),
    width = 4, height = 4, units = 'in', res = 400)
par(mar = c(3.25,3.25,0.5,0.5))
plot(PropAb~Dx, data = human.test.dat, cex = sc*NumTestAb,
     xlab = '', ylab = '',
     main = '', ylim = c(0,0.6), xlim = c(0,0.85))
mtext(side = 1, text = as.expression(bquote('D'['x']~'Layer')), line = 2.4)
mtext(side = 2, text = 'Seroprevalence', line = 2.4)
legend(x = 'topleft', legend = c('10', '100', '400'), pch = c(1,1,1), pt.cex = sc*c(10,100,400))
xseq <- seq(0,1,by = 0.01)
newdata <- data.frame(Dx = xseq)
yseq <- predict(object = qb.mod, newdata = newdata, type = 'response')
lines(xseq, yseq)
dev.off()

## Seroprevalence across West Africa
heat.cols <- viridis(120, begin = 0.1, end = 1, option = 'D')
png(file = paste(fold.name,'/Lassa_Seroprevalence.png', sep = ''),
    height = 4, width = 6, units = 'in', res = 400)
par(mai = 1*c(0.2,0.2,0.2,1))
zmax <- maxValue(pred.rast)
zmin <- minValue(pred.rast)
image.plot(pred.rast, col = heat.cols, zlim = c(zmin, zmax),
           bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',
               xlim = xlims, ylim = ylims,
           asp = 1, legend.lab = 'Predicted Human Seroprevalence', legend.line = 4,
           bigplot = c(.0,.75,0.1,0.9),
           smallplot = c(.765,.8,0.1,0.9), legend.mar = 0)
plot(foc.shp.ogr, add = TRUE, bty = 'n', asp = 1)
plot(rgeos::gIntersection(foc.shp.ogr, masto.rangemap), add = TRUE, bty = 'n', asp = 1, lwd = 3)
dev.off()

#########################

## Make plot of residuals
human.test.dat$Resid <- with(human.test.dat, PropAb - fitted)
human.test.dat$ResidSign <- with(human.test.dat, ifelse(abs(Resid) < 0.1, 0, sign(Resid)))

## Plot residuals
xlims = c(-18,18)
ylims = c(16, 16.5)
cex.int <- 1
cex.sc <- 2

## First plot has discrete point coloration
jpeg(file = paste(fold.name,'/resid_v1.jpg', sep=''), width = 6.5,
     height = 4, units = 'in', res = 400)
par(mai = 1*c(0.2,0.2,0.2,0.2))
image.plot(pred.rast, col = heat.cols, zlim = c(zmin, zmax),
           bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',
           main = '', xlim = xlims, ylim = ylims,
           asp = 1,
           legend.lab = 'Predicted Seroprevalence',
           legend.line = 2.75)
cexvals <- cex.int + cex.sc*abs(human.test.dat$Resid)
points(Latitude~Longitude, data = human.test.dat[human.test.dat$NumTestAb > 50,],
       asp = 1, pch = 21, lwd = 0.8,
       bg = c('deepskyblue1', 'white', 'red')[ResidSign + 2],
       col = 'black', cex = 0.75)##cexvals)
plot(foc.shp.ogr, main = '', add = TRUE, bty = 'n', asp = 1)
plot(rgeos::gIntersection(foc.shp.ogr, masto.rangemap), add = TRUE, bty = 'n', asp = 1, lwd = 3)
dev.off()


## 2nd Plot has continuous color scheme
n.cols <- 100
pt.cols <- colorRampPalette(c('darkblue', 'white', 'red'), bias=1)(n.cols)
seq.vals <- seq(-0.4, 0.4, length.out = 100)
col.fun <- function(x){
    cols <- c()
    for(xi in 1:length(x)){
        wi <- which.min(abs(x[xi] - seq.vals))
        cols <- c(cols, pt.cols[wi])
    }
    return(cols)
}
cex.fun <- function(x){
        cex.vals <- c()
        for(xi in 1:length(x)){
            cex.vals <- c(cex.vals, min(x[xi]/100, 1))
        }
        return(cex.vals)
}

jpeg(file = paste(fold.name,'/resid_v2.jpg', sep=''), width = 6,
     height = 5, units = 'in', res = 400)
par(mai = 1*c(0.0,0.2,0.0,0.6))
image.plot(pred.rast, col = heat.cols, zlim = c(0, zmax),
           bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',
           main = '', xlim = xlims, ylim = ylims,
           asp = 1, legend.lab = 'Predicted Seroprevalence',
           legend.line = 2.75)
image.plot(col = pt.cols,legend.only = TRUE, horizontal = TRUE,
           zlim = c(-0.4, 0.4), smallplot = c(0.1, 0.8, 0.93, 0.95),
           legend.lab = 'Actual - Predicted Seroprevalence',
           legend.line = -1.7, bg = 'white')
cexvals <- cex.fun(human.test.dat$NumTestAb)
points(Latitude~Longitude, data = human.test.dat, asp = 1, pch = 21, lwd = 0.8,
       bg = col.fun(Resid),
       col = 'black', cex = cexvals)
plot(foc.shp.ogr, main = '', add = TRUE, bty = 'n', asp = 1)
plot(rgeos::gIntersection(foc.shp.ogr, masto.rangemap), add = TRUE, bty = 'n', asp = 1, lwd = 3)
legend(x = -15, y = 3.5, legend = c('10', '100', '1000'),
       pch = 19, pt.cex = cex.fun(c(10,100,1000)),
       pt.bg = 'gray', horiz = TRUE, bty = 'n',
       title = 'Number Tested')
dev.off()

#########################

## --- Use predicted human seroprevalence raster (pred.rast) to predict the
## incidence of Lassa across West Africa

## reproject to coordinate reference system of human population data
pred.rast <- resample(pred.rast, pop.rast)
human.lasv.cases <- pred.rast*pop.rast*(d + del)*(d)*del^-1

## Plot human.lasv.cases
xlims = c(-18,16)
ylims = c(16, 16.5)
heat.cols <- viridis(120, begin = 0.1, end = 1, option = 'D')
png(file = paste(fold.name,'/New_Cases.png', sep = ''),
    height = 4, width = 6, units = 'in', res = 400)
par(mai = 1*c(0.2,0.2,0.2,1))
agg.cases <- aggregate(human.lasv.cases, fact = 10, fun = sum)
case.density = agg.cases/raster::area(agg.cases)
case.density = mask(case.density, masto.rangemap, updatevalue = 0)

## Set arbitrary lower bound on case.density (so its not -Inf)
ymin <- 1e-3
case.density[case.density < ymin] <- ymin
yaxlab = c(0,sapply(paste(-2:2), function(x){as.expression(bquote(10^ .(x)))}))
image.plot(log10(case.density), col = heat.cols, zlim = c(-3, 2.5),
           bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',
           xlim = xlims, ylim = ylims,
           asp = 1, legend.lab = as.expression(bquote('Predicted Cases / ( Year x km'^2*' )')),
           legend.line = 4,
           bigplot = c(.0,.75,0.1,0.9),
           smallplot = c(.765,.8,0.1,0.9), legend.mar = 0,
           axis.args = list(at = c(-3,-2,-1,0,1,2), labels=yaxlab))
plot(foc.shp.ogr, add = TRUE, bty = 'n', asp = 1)
plot(rgeos::gIntersection(foc.shp.ogr, masto.rangemap), add = TRUE, bty = 'n', asp = 1, lwd = 3)
dev.off()


## Calculate the number of cases over each country. This step uses a function
## defined in the Tools/Functions.r  script. This function takes ~10 minutes to
## sum all cases within West African countries.
print("Calculating incidence of LASV in humans", quote = FALSE)
starttime <- Sys.time()
cases.dat <- get.country.sums(human.lasv.cases)
write.table(cases.dat, file = paste(fold.name, '/foc_case', sep = ''))
print(Sys.time() - starttime)

## Simple dataframe with LASV case count totals
tot.model <- cases.dat[cases.dat$Country=='Total','Cases']
tot.low.round <- round(tot.model,1)
tot.high.round <- round(4.2*tot.model,1)
output <- data.frame(type = 'model', low.est = tot.low.round, high.est = tot.high.round,
                     stringsAsFactors = FALSE)
output <- rbind(output, c('null', round(null.low.estimate/1000,1), round(null.high.estimate/1000,1)))

## For Table 1 of the manuscript
cases.dat.round <- cases.dat
cases.dat.round$Cases <- round(cases.dat$Cases,1)


