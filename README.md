
# Overview

This directory contains the scripts and data-sets used to build the
model that is outlined in the paper: 

Basinski, AJ *et al*, 2020. Bridging the Gap: Using reservoir ecology and
human serosurveys to estimate Lassa virus spillover in West
Africa. Submitted to PLOS Computational Biology. 

See the manuscript for full details. The scripts in this directory
describe a machine-learning framework that predicts the number of
spillover infections of an animal-borne pathogen into humans. Almost all
scripts below are written in the statistical language
R. Human_Link_Analysis.cdf contains Mathematica code that derives
expressions described in the manuscript. The model framework is applied
to the Lassa virus (LASV), a rodent-borne pathogen capable of causing
febrile illness (i.e., Lassa Fever) in humans.

The scripts contained in this directory take in data-sets that
describe locations of the LASV reservoir, locations of pathogen
presence/absence from surveys of reservoir populations, as well as
human seroprevalence data, and output a prediction of the number of
spillover infections into human populations each year. 

The model generates predictions of the spillover rate into humans in
two stages.  In the first stage, the model learns to associate the
spatial distribution of the animal reservoir, and the distribution of
the pathogen, with environmental attributes and human population
density.  These associations are learned by boosted classification
trees within the Reservoir Layer and Pathogen Layer directories,
respectively. These fitted models are used to generate a combined risk
of pathogen spillover into humans.  The second stage of the model is
contained in the Human_LASV_Incidence directory, and uses a regression
model to associate seroprevalence in human populations with the
combined risk derived in the first stage. Hypothesis testing that is
associated with this regression automatically provides a rigorous
evaluation of the spillover risk's ability to predict human LASV
seroprevalence.  Lastly, an SIRS model is used to convert the
predicted seroprevalence into estimates of human LASV infections in
West Africa.

The results presented in the manuscript were created with R version
4.0.0 (2020-04-24) "Arbor Day". Different versions of R
will give slightly different results because of differences in the random
number generator that is used. Running the model code requires
installing the following packages:

| Functionality                | Package Names               | 
|------------------------------|----------------------------:|
| GIS capabilities             |  raster, rgdal, sf, rgeos   | 
| Parallel processes           |  parallel                   |
| Boosted models               |  dismo, gbm                 |
| Plotting                     |  ggplot2, gridExtra, fields, ggthemes, viridis | 
| Neg. Binomial                | MASS, pscl, DescTools       |
| Misc                         | verification, dplyr, plyr   |
| <img width=100/>             |<img width=200/>             |


# Predictors

Predictors are included as a pre-processed raster stack
(Storage/Raster_Data/all.stack). While we include the scripts that
were used to process and compile the original raster data into
stack-form, we do not include the original data itself because of
storage considerations. If the user wishes to reproduce the
calculations that produce the data stack, the original data-set must
be re-downloaded from the sources referenced in the manuscript and
placed in the Original_Data directories within the Storage directory.

Predictors that are used in the Reservoir layer, as well as the
Pathogen layer, are chosen from the list of candidate predictors using
a Wilcox significance test (described in the layers below). Candidate
predictors are listed in the tables below.

|Name| Description                                                      | 
|----|------------------------------------------------------------------|
|Elev| Elevation in meters                                              |
|Tmu | Mean daily temperature (C)                                       |
|Pmu | Mean daily precipitation (mm/day)                                |
|Nmu | Mean NDVI (normalized difference vegetation index)               |
|Pcv | Coefficient of variation for precipitation within a year         |
|Ncv | Coefficient of variation for NDVI within a year                  |
|Pc  | Colwell's index of within-year precipitation constancy           |
|Pm  | Colwell's index of within-year precipitation contingency         |
|Nc  | Colwell's index of within-year NDVI constancy                    |
|Nm  | Colwell's index of within-year NDVI contingency                  |
|Pmin| Minimum rainfall within a year                                   |
|Pmax| Maximum rainfall within a year                                   |
|Nmin| Minimum NDVI within a year                                       |
|Nmax| Maximum NDVI within a year                                       |
|Pdur| Duration of rainfall below 1 mm/day                              |
|Ndur| Duration of NDVI below 0.5                                       |
|Pop | Human population density (humans / km^2)                         |

Each environmental data-set was downloaded in its original raster
form, and reprojected into 0.05 x 0.05 degree pixels using the WGS1984
coordinate reference system. Precipitation, temperature, and NDVI
rasters were used to generate measures of seasonality (e.g., Pcv,
Pcol, Ncv, Ncol, etc.). The code used to derive these features is
contained in the **Calc_Seasonality.r** script, contained in the
Storage/Data_Scripts/ directory. Seasonal predictors are calculated
within each year over the time span 2001 - 2019, then averaged across
years.

We also include land cover data derived from the MODIS land cover
classification data-set. Each pixel of the reprocessed data rasters
describes the fraction of surrounding 0.15 x 0.15 degree area that is
of a specific land cover type. Only those land cover predictors that
had an average longevity of 20 or more years were included as
candidate predictors. The **process_MODIS_landcover.r** script,
contained in the Storage/Data_Scripts/ directory, performs these 
calculations, and outputs a list of land cover types that persist 
for an average time of 20 or more years:

| MODIS Land Cover         |
|--------------------------|
|Water_Bodies              |
|Evergreen_Broadleaf_Forest|
|Open_Shrubland            |
|Woody_Savanna             |
|Savannas                  |
|Grasslands                |
|Permanent_Wetlands        |
|Croplands                 |
|Urban_BuiltUp             |
|Cropland_Natural_Mosaic   |
|Barren                    |

The derived seasonality features, as well as the land cover features, 
are then combined into a named raster stack, and saved to the location
Storage/Raster_Data/all.stack. This is performed by the script 
**Prep_Predictors.r**. 


# Reservoir_Layer

This directory contains scripts that take in latitude/longitude
coordinates of *Mastomys natalensis* capture locations, as well as
capture locations of other rodents, and outputs a .tif that predicts
the probability that any given pixel of West Africa contains
*M. natalensis*.

The subdirectories are 

- Data: Contains the data-sets that describe *M. natalensis*
  presences, as well as other rodent presences that are used as
  background data. The Excel file contains 2 spreadsheets. The first
  contains all literature records of genetically (or morphologically)
  confirmed *M. natalensis* rodents. The second spreadsheet is a
  cleaned version of the first; entries that were obtained through
  personal communication, museum specimens only, etc. are filtered
  out.
- Figures_Fits: Contains all figure outputs from the fitting process,
as well as the final tif of the prediction.
- Models: created at run time. Used to temporarily hold model files generated by
Train_Reservoir_Learners.r code, and deleted after the model is fit.
- Tools: Contains a script with additional functions called Functions.r

To re-derive all data-processing and model fitting performed for this
layer, first set Reservoir_Layer as the working directory, then source
the main script by running:

source("Generate_Reservoir_Layer.r").

This script loads and calls three functions, defined in other scripts
detailed below, that prepare the data-set (Prep_Reservoir_Data.r),
determine which predictors significantly vary with the response
variable (Calc_Sig_Reservoir_Preds.r), and uses those predictors to
train a set of 25 boosted classification trees
(Train_Reservoir_Learners.r). Specifically,
Generate_Reservoir_Learners.r loops through sets of hyper-parameters
that specify a boosted classification tree algorithm, and passes said
hyper-parameters to the Train.Reservoir.Learners function for fitting.
All information from the Generate script is stored in a
sub-directory of Figures_Fits, [prefix], that is created when the
code is run. 

The three main functions that Generate_Reservoir_Learners.r calls are
contained in the following scripts: 

###  Prep_Reservoir_Data.r 

This script defines the function Prep.Reservoir.Data. The function
accepts a string that specifies the focal species. For now, this
string can only take on the value "Mastomys natalensis", but other
species will eventually be added. This function processes *M.
natalensis* presence data. The output is a dataframe that contains
Latitude, Longitude, and associated environmental predictors of pixels
that contain Mastomys natalensis captures, as well as background
captures of other Murids.

###  Calc_Sig_Reservoir_Preds.r

This script defines the Calc.Sig.Reservoir.Preds function. The
function accepts a species name (string), as well as the dataframe
output of Prep.Reservoir.Data, and outputs a list of predictors
determined to significantly vary with the *M. natalensis* response
variable. Here, absence refers to background captures of Murid
rodents. Significance is determined by a Wilcox test with a threshold
of p = 0.05.  Significant predictors, in turn, will be used to train
the the model that predicts the reservoir's presence or absence.
 
This function returns a list of the names of predictors that are
deemed significant (also saved to Figures_Fits/Sig_Reservoir_preds),
and also saves a pdf figure showing the distributions of the
significant predictors.

### Train_Reservoir_Learners

This script defines a function, train.reservoir.learners, that takes
as input a data-set (output by Prep.Reservoir.Data) and hyperparameter
list (output of Calc.Sig.Reservoir.Preds), and saves to file a fitted
reservoir risk layer. This script is the workhorse of the
reservoir-model building process. Outputs are saved in the
Figures_Fits/ directory, and include: 1) a .tif file of the averaged
reservoir-risk prediction generated by [nboots] bootstrapped model
fits. [nboots] is a hyperparameter defined in the
Generate_Reservoir_Layer.r script; 2) Figures showing variable
importance, learned relationships, and predicted risk across West
Africa, averaged over all bootstrapped fits; 3) Data file tree.dat
that contains information on the fitting process for each bootstrapped
fit, including "area under the receiver operator curve" (AUC)
calculated on an out-of-bag test set with and without pairwise
distance sampling, and the number of trees-building iterations that
were deemed best by cross-validation; 4) a data file assess.dat that
shows the average classification score assigned to each type of
species in the rodent dataset; 5) a data file grid_res_data that 
summarizes performance across each hyper-parameter set. 

The grid_res_data (stored in prefix directory) gives a performance summary of each hyper-parameter
set tried in Generate_Reservoir_Layer.r. Most columns are
self-explanatory. The auc and auc.pwd columns show the out-of-bag AUC
model performance without, and with, pairwise distance sampling. 

The tree.dat file (stored in subdirectory of prefix) contains contains
information on the out-of-bag model performance for each bootstrapped
fit. The columns are:
- boot.i: ith bootstrapped model, where i is between 1 and [nboots]
- n.tree: number of trees used by the model 
- max.tree: maximum number of trees allowed (set in
  Generate_Reservoir_Layer.r script
- model.oob.auc: out-of-bag AUC on test data without pairwise distance
  sampling
- model.oob.auc.pwd: out-of-bag AUC on test data with pairwise
  distance sampling
- n.oob.auc: number of test points used in calculation of model.oob.auc
- n.oob.auc.pwd: number of test points used in calculation of model.oob.auc.pwd

The assess.dat file (stored in subdirectory of prefix) contains more
information of how the bootstrapped models performed on various
background species from the out of bag, pairwise-distance sampled test
set. The columns are:
- species
- mean: The mean classification score assigned to that species,
  averaged across all [nboots] models
- sd: The standard deviation of the classification score across all models 
- medianCount: The median number of background points that are of each
  species type 




# Pathogen_Layer

This directory contains scripts that take in data that describes the 
presence or absence of LASV in surveys of *Mastomys natalensis*, and 
outputs a .tif file that describes the probability that, given that 
*M. natalensis* is present in a pixel, LASV circulates in the rodent 
population. 

The directory structure is similar to the Reservoir_Layer directory.
Scripts contained in this directory are directly analogous to those 
described in the Reservoir_Layer directory.

- Data: contains the datasets that describe LASV presence and absence
in rodents, as well as data on serosurveys in humans. The Excel file
"Lassa_Occurrence_Data" contains 3 spreadsheets. The "Genbank"
spreadsheet contains LASV records in *M. natalensis* found on GenBank;
the "Raw_Lassa_Literature" spreadsheet contains LASV records in humans
and rodents that were found through a literature search; and the
"Cleaned_Lassa_Literature" spreadsheet contains cleaned version of
LASV records in literature that is used to produce a CSV file that, in
turn, is read into the model pipeline. Cleaning involved removing
human serosurveys that were nonrandom (e.g., tested individuals with
symptoms, individual missionaries), and removing entries that did not
specifically test *M. natalensis* rodents.

- Figures_Fits

To re-derive all data-processing and model fitting for this layer, set
the Pathogen_Layer directory as the working directory, and then run:

source("Generate_Pathogen_Layer.r")

As with the Generate_Reservoir_Layer script, this script loads and
calls three functions that prepare the data-set (Prep_Pathogen_Data.r),
determine which predictors significantly vary with the response
variable (Calc_Sig_Pathogen_Preds.r), and uses those predictors to
train a set of 25 boosted classification trees
(Train_Pathogen_Learners.r).

### Prep_Pathogen_Data.r

This script defines the Prep.Pathogen.Data function. This function
processes Lassa presence/absence data from published surveys of
rodents and humans. The function returns a dataframe containing
*M. natalensis* Lassa presence/absence data. The rodent
dataframe contains environmental predictors that are used to train the
pathogen layer of the model. Another dataset that contains human
serosurvey data is saved to a csv file
(Data/Prepped_Human_Seroprevalence_Data.csv), and used for the
human stage of the framework.

###  Calc_Sig_Pathogen_Preds.r

Analogous to Calc_Sig_Reservoir_Preds.r. Defines a function called
Calc.Sig.Pathogen.Preds that returns a list of those predictors that
are significantly associated with the presence or absence of the
pathogen. These are the predictors the boosted models are trained on.
 
This function returns a list of the names of predictors that
are deemed significant (also saved to Figures_Fits/Sig_LASV_preds), and also
saves a pdf figure showing the distributions of the significant
predictors. 

### Train_Pathogen_Learners.r

This script defines a function, train.pathogen.learners, that takes as
input a data-set and hyperparameter list, and outputs the fitted
pathogen risk layer.  This script is the workhorse of the
pathogen-model building process, and is called by
Generate_Pathogen_Layer.r.  Outputs are saved in the Figures_Fits/
directory, and include: 1) a .tif file of the averaged pathogen-risk
prediction generated by [nboots] bootstrapped model fits. [nboots] is
a hyperparameter defined in the Generate_Pathogen_Layer.r script; 2)
Figures showing variable importance, learned relationships, and
predicted risk across West Africa, averaged over all bootstrapped
fits; 3) Data file tree.dat that contains information on the fitting
process for each bootstrapped fit, including "area under the receiver
operator curve" (AUC) calculated on an out-of-bag test set, and the
number of trees-building iterations that were deemed best by
cross-validation; 4) data file grid_lsv_data (stored in prefix
directory) that summarizes the performance of each hyper-parameter set
run. See Train_Reservoir_Learners.r section for an explanation of 
the data file column names. 


### Map_Cases.r

This script produces a ggplot map that indicates the point locations
of rodent surveys for Lassa that were used to train the pathogen
layer, as well as human sero-surveys that were used to calibrate the
risk layers.


# Human_LASV_Incidence

This directory contains code and output that predicts human 
LASV seroprevalence and incidence from the tif outputs of 
the reservoir and pathogen layers. 

### Human_Link.r

This script reads in the tif outputs from Generate_Reservoir_Layer.r
and Generate_Pathogen_Layer.r, and, using data from human LASV
serosurveys, outputs a prediction of human LASV seroprevalence and
incidence across West Africa. This prediction is made in three steps:

1. Let $D_M$ and $D_L$ denote the reservoir and pathogen risk layers
outputed Generate_*_Layer.r scripts above. Calculate the combined risk
raster, $D_X$ as 
$$
D_X = D_M \times D_L. 
$$
2. Regress the human seroprevalence data onto the combined risk layer
   $D_X$. Use this regression to predict the prevalence of LASV across
   West Africa. 
3. Use equations derived in the manuscript to convert the predicted
   seroprevalence into predictions of LASV incidence. 

Figures are saved in the [grid.name] folder (located within
Figures_Fits), except for the product risk layer. This script also
outputs a data file foc_case, contained in the grid.name folder, 
that summarizes data found in the manuscript table. Column names are
self-explanatory; "Pop" denotes the country's total population size. 

### Human_Link_Analysis.cdf

This is a Mathematica file that derives the expressions found in 
the manuscript. 




