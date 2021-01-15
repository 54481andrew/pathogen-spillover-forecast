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
## The pathogen layer was run in two ways that differed on how ambiguous
## pixels were handled. These were set to NA or 1; [set.ambiguous.to]
## allows the user to choose which pathogen layer fit is used for the human
## predictions. 
set.ambiguous.to <- 1
min.test <- 10
grid.name <- paste0('human_v4_ambi',set.ambiguous.to, '_mint_', min.test) ## name for model.dat dataframe
prefix.masto <- 'reservoir_v3' ## prefix used in Generate_Reservoir_Layer
masto.fold <- 'Mn_pa_nboots_25_nbg_Same_tc_1_mllr_2_lmt_7'

prefix.lassa <-  'pathogen_v4' ## prefix used in Generate_Lassa_Layer
lassa.fold <- paste0("pa_nboots_25_1_mllr_3_lmt_7_ambi_", set.ambiguous.to, "_mintest_",min.test)


## Load additional functions
source("Tools/Functions.r")

## Parameters for calculating incidence of human LASV from human LASV seroprevalence.
## All time units are years.
LS.rast <- raster("../Storage/Raster_Data/Lifespan.tif") ## Lifespans in West Africa (WorldBank)
gam = 12 ## LASV recovery is 1 month
lam.mcc87 = 0.064 ## Seroreversion rate from McCormick, Webb, Krebs, 1987
mu = 0.02 ## LASV mortality, Pr(Infected dies) = 0.02, from McCormick, Webb, Krebs, 1987

## Load in Africa shapefiles needed for plotting
storage.fold <- '../Storage'
foc.shp.ogr <- readOGR(dsn = paste(storage.fold, '/Shapefiles/West_Africa', sep = ''),
                       layer = 'foc', verbose = FALSE)
africa.ogr <- readOGR(dsn = paste(storage.fold, '/Shapefiles/Africa', sep = ''),
                       layer = 'Africa', verbose = FALSE)
masto.rangemap <-readOGR(dsn = paste(storage.fold, '/Shapefiles/Masto_Range',sep=''),
                         layer = 'data_0', verbose = FALSE)


## Load in and uncompress human population data.
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

## Load Mastomys natalensis prediction
filename.masto <- paste0('../Reservoir_Layer/Figures_Fits/', prefix.masto, '/',
                         masto.fold, '/', "Reservoir_Layer_",masto.fold, '.tif')
Masto.Risk <- raster(filename.masto)
Masto.Risk <- mask(Masto.Risk, masto.rangemap, updatevalue = 0)

## Load in Lassa prediction
filename.lassa <- paste0('../Pathogen_Layer/Figures_Fits/', prefix.lassa, '/',
                         lassa.fold, '/', 'Lassa_Layer_', lassa.fold,".tif")
Lassa.Risk <- raster(filename.lassa)
Lassa.Risk <- mask(Lassa.Risk, masto.rangemap, updatevalue = 0)

## --- Load in human data from appropriate pathogen prefix
human.test.dat <- read.csv(paste0('../Pathogen_Layer/Figures_Fits/', prefix.lassa, '/',
                         lassa.fold, '/','Prepped_Human_Seroprevalence_Data.csv'))

## Weighted mean seroprevalence across all sero-surveys; good to know
wmean <- sum(human.test.dat$NumPosAb)/sum(human.test.dat$NumTestAb)

## --- Create location in which images will be saved
fold.name <- paste('Figures_Fits/', grid.name, sep = '')
dir.create(fold.name, showWarnings = FALSE, recursive = TRUE)


## Plot product, D_X = D_L \times D_M, that describes the combined risk of LASV
Dx = Lassa.Risk*Masto.Risk
heat.cols <- viridis(120, begin = 0.1, end = 1, option = 'D')
xlims = c(-18,16)
ylims = c(16, 16.5)
png(file = paste(fold.name, '/Product_Risk_Layer.png', sep = ''),
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

## Extract risk prediction Dx from each human sero-survey location
human.test.dat$Dx <- extract(Dx, human.test.dat[,c('Longitude', 'Latitude')])

## Fit a standard binomial as a baseline
bin.mod <- glm(PropAb~ 1 + Dx, human.test.dat,
              family = 'binomial',
              weights = NumTestAb)
bin.summ <- summary(bin.mod)

## Fit quasi-binomial model. Note that the parameter
## estimates of the quasi-binomial are identical to
## those of the standard binomial. 
qb.mod <- glm(PropAb~ 1 + Dx, human.test.dat,
              family = 'quasibinomial',
              weights = NumTestAb)
qbin.summ <- summary(qb.mod)


## --- Use regression model to predict human LASV seroprevalence across West Africa
pred.stack <- stack(Dx)
names(pred.stack) = 'Dx'
pred.rast <- raster::predict(object = pred.stack, model = qb.mod, type = 'response')


## Add some useful columns in human.test.dat data frame
human.test.dat$fitted <- extract(x = pred.rast,
                          y = human.test.dat[,c('Longitude', 'Latitude')])
human.test.dat$res <- residuals(bin.mod, type = 'deviance')
human.test.dat$y.linear <- with(human.test.dat,
                                Dx*qb.mod$coefficients['Dx'] + qb.mod$coefficients['(Intercept)'])


## ---Make a plot of the deviance residuals
png(filename = paste(fold.name,'/Deviance_Residuals.png',sep=''),
    width = 4, height = 4, units = 'in', res = 400)
par(mar = c(3.25,3.25,0.5,0.5))
plot(res~fitted, data = human.test.dat,
     xlab = '', ylab = '',
     main = '')
mtext(side = 1, text = 'Fitted Value', line = 2.4)
mtext(side = 2, text = 'Deviance Residual', line = 2.4)
dev.off()

## ---Plot prediction vs Dx with confidence intervals.
## First compute CI using a bootstrap method.

f <- function(seed){
    ## Set seed for reproducibility
    set.seed(seed) 
    ## Choose a random set of surveys
    boot.set <- sample(nrow(human.test.dat), replace = TRUE)
    boot.dat <- human.test.dat[boot.set,]
    ## The Quasibinomial and binomial GLM give the same mean predictions,
    ## so we use the binomial GLM. 
    bin.mod <- glm(PropAb~ 1 + Dx, boot.dat,
                   family = 'binomial',
                   weights = NumTestAb)
    ## Return predictions
    boot.pred <- predict(bin.mod, newdata = data.frame(Dx = seq(0,1,by = 0.01)), type = 'response')
    return(boot.pred)
}
boot.out <- sapply(X = 1:1000, FUN = f)
## Find 100x[xconf]% confidence intervals
xconf <- 0.95
quantile.fun <- function(boot.row){
    return()
    }
conf.out <- sapply(X = 1:nrow(boot.out),
                   FUN = function(x){as.numeric(quantile(boot.out[x,],
                                                         probs = c((1-xconf)/2, 1-(1-xconf)/2)
                                                         ))})

## Define some variables/functions that are used in graphing
xseq <- seq(0,1,by = 0.01)
y.linear.xseq <- xseq*qb.mod$coefficients['Dx'] + qb.mod$coefficients['(Intercept)']
ilogit <- function(x){exp(x)/(1 + exp(x))} ## inverse logit will come in handy

conf.out <- t(conf.out)
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
lines(xseq, ilogit(y.linear.xseq), col = 'black', lty = 1)
matlines(xseq, conf.out, col = 'red', lty = 2)
dev.off()


## -----Plot seroprevalence prediction across West Africa
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
png(file = paste(fold.name,'/resid_v1.png', sep=''), width = 6.5,
     height = 4, units = 'in', res = 400)
par(mai = 1*c(0.2,0.2,0.2,0.2))
image.plot(pred.rast, col = heat.cols, zlim = c(zmin, zmax),
           bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',
           main = '', xlim = xlims, ylim = ylims,
           asp = 1,
           legend.lab = 'Predicted Seroprevalence',
           legend.line = 2.75)
cexvals <- cex.int + cex.sc*abs(human.test.dat$Resid)
plot(rgeos::gIntersection(foc.shp.ogr, masto.rangemap), add = TRUE, bty = 'n', asp = 1, lwd = 3)
points(Latitude~Longitude, data = human.test.dat,
       asp = 1, pch = 21, lwd = 0.8,
       bg = c('deepskyblue1', 'white', 'red')[ResidSign + 2],
       col = 'black', cex = 0.75)##cexvals)
plot(foc.shp.ogr, main = '', add = TRUE, bty = 'n', asp = 1)
dev.off()

#########################

## --- Use predicted human seroprevalence raster (pred.rast) to predict the
## incidence of Lassa across West Africa

## reproject rasters to coordinate reference system of human population data
pred.rast <- resample(pred.rast, pop.rast)
LS.rast.proj <- resample(LS.rast, pop.rast, method = 'ngb')

## Convert lifespan to country-specific death rate; calculate LASV case rate
## without reinfection. In the mathematica notebook, we show that seroreversion
## rate lam simply multiplies the case rate by a factor (d+lam)/d. We'll
## incorporate this multiplying effect of seroreversion below. 
drast <- 1/LS.rast.proj
human.lasv.cases <- pred.rast*pop.rast*(drast + gam)*(drast) / (gam * (1-mu))

## Aggregate cases 
agg.cases <- aggregate(human.lasv.cases, fact = 10, fun = sum)
case.density = agg.cases/raster::area(agg.cases)
case.density = mask(case.density, masto.rangemap, updatevalue = 0)

## Plot human.lasv.cases
xlims = c(-18,16)
ylims = c(16, 16.5)
heat.cols <- viridis(5, begin = 0.0, end = 1, option = 'D')
png(file = paste(fold.name,'/New_Cases.png', sep = ''),
    height = 4, width = 6, units = 'in', res = 400)
par(mai = 1*c(0.2,0.2,0.2,1))
## Set arbitrary lower bound on case.density (so its not -Inf)
ymin <- 1e-3
case.density[case.density < ymin] <- ymin
yaxlab = c(0,sapply(paste(-2:2), function(x){as.expression(bquote(10^ .(x)))}))
image.plot(log10(case.density), col = heat.cols, breaks = c(-3,-2,-1,0,1,2),
           zlim = c(-3, 2.5),
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
## Countries over which sums will be calculated
foc.countries <- c('Mali', "Guinea", "Ivory Coast", "Sierra Leone",
                   "Nigeria", 'Liberia', "Cote d`Ivoire",
                   'Ghana', 'Togo', 'Benin',
                   'Burkina Faso', 'Mauritania', 'Niger',
                   'Senegal')
org.cases.dat <- get.country.sums(human.lasv.cases)

## Calculate the effect of reinfection and virulence for each country. This is done by
## calculating the multiplying effect of seroreversion on case estimates for each
## country, given that country's average lifespan.

lifespan.dat <- read.csv("../Storage/CSV_Data/foc_Lifespan.csv", stringsAsFactors = FALSE) 
## Add a dummy country "Total" for compatibility with org.cases.dat
lifespan.dat <- rbind(lifespan.dat,
                      data.frame(Country.Name = 'Total', X2018 = mean(lifespan.dat$X2018)))
lifespan.dat$d <- 1/lifespan.dat$X2018
lifespan.dat$lam.mult <- with(lifespan.dat, (d + lam.mcc87)/d) ## Multiplying effect of reinfection

merged.cases.dat = merge(org.cases.dat, lifespan.dat[,c('Country.Name','lam.mult')], by.x = 'Country',
      by.y = 'Country.Name', sort = FALSE)


## Calculate case rates by country with virulence and reinfection incorporated
merged.cases.dat$Rate = with(merged.cases.dat, Cases/Pop)
merged.cases.dat$Cases.high <- with(merged.cases.dat, Cases*lam.mult) 
merged.cases.dat$Rate.high = with(merged.cases.dat, Cases.high/Pop)

## Save original and rounded values to file
## Cases stored in units of 1000's of individual. Rate is per 1000 individuals. 
cases.dat <- merged.cases.dat[,c('Country','Cases', 'Rate', 'Cases.high', 'Rate.high')]
cases.dat$Cases <- round(merged.cases.dat$Cases/1000, 1) 
cases.dat$Rate <- round(1000*merged.cases.dat$Rate, 1)
cases.dat$Cases.high <- round(merged.cases.dat$Cases.high/1000, 1)
cases.dat$Rate.high <- round(1000*merged.cases.dat$Rate.high, 1)

## Write case data and GLM output to file 
write.table(cases.dat, file = paste(fold.name, '/foc_case', sep = ''))
write.table(org.cases.dat, file = paste(fold.name, '/org_foc_case', sep = ''))
capture.output(qbin.summ, file = paste(fold.name, '/glm_summary', sep = ''))


print(Sys.time() - starttime)

## Simple dataframe with LASV case count totals
tot.model <- cases.dat[cases.dat$Country=='Total','Cases']
tot.round <- round(tot.model,1) ## Base estimates
tot.high.round <- round(cases.dat[cases.dat$Country=='Total', 'Cases.high'], 1) ## with Reinfection
output <- data.frame(type = 'model', low.est = tot.round, high.est = tot.high.round,
                     stringsAsFactors = FALSE)
