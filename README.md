
# Overview

This directory contains the scripts and datasets used to build the
model that is outlined in the paper: 

Basinski AJ, Fichet-Calvet E, Sjodin AR, Varrelman TJ, Remien CH,
Layman NC, et al. (2021) Bridging the gap: Using reservoir ecology and
human serosurveys to estimate Lassa virus spillover in West
Africa. PLoS Comput Biol 17(3):
e1008811. https://doi.org/10.1371/journal.pcbi.1008811

See the manuscript for full details. The scripts in this directory
describe a machine-learning framework that predicts the number of
spillover infections of an animal-borne pathogen into humans. Almost
all scripts below are written in the statistical language R, however,
Human_Link_Analysis.cdf (in Human_LASV_Incidence) contains Mathematica
code that derives expressions described in the manuscript. The model
framework is applied to the Lassa virus (LASV), a rodent-borne
pathogen capable of causing febrile illness (i.e., Lassa Fever) in
humans.

The scripts contained in this directory take in datasets that
describe locations of the LASV reservoir, locations of pathogen
presence/absence from surveys of reservoir populations, as well as
human seroprevalence data, and output a prediction of the number of
spillover infections into human populations each year. 

The model generates predictions of the spillover rate into humans in
two stages. In the first stage, the model learns to associate the
spatial distribution of the animal reservoir, and the distribution of
the pathogen in the animal reservoir, with environmental attributes
and human population density. These associations are learned by
boosted classification trees within the Reservoir Layer and Pathogen
Layer directories, respectively. These fitted models are used to
generate a combined risk of pathogen spillover into humans. The second
stage of the model is contained in the Human_LASV_Incidence directory,
and uses a regression model to associate seroprevalence in human
populations with the combined risk derived in the first stage. At this
stage, the combined risk, which is fit using data from the animal
reservoir only, is assessed on its ability to correlate with human
seroprevalence, a completely independent test dataset. Correlation
coefficients (weighted by number of humans tested in each serosurvey
or unweighted), as well as the p-value derived from a hypothesis test
in the GLM regression, all provide objective numerical metrics that
describe the predictive ability of the combined risk layer. Lastly, an
SIRS model is used to convert the predictions of the seroprevalence
regression into estimates of human LASV infections per year in West
Africa.

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
| Misc                         | verification, dplyr, plyr   |
| <img width=100/>             |<img width=200/>             |


# Predictors

Predictors are included as a pre-processed raster stack
(Storage/Raster_Data/predictor_stack.grd). While we include the scripts that
were used to process and compile the original raster data into
stack-form, we do not include the original data itself because of
storage considerations. If the user wishes to reproduce the
calculations that produce the predictor stack, the original dataset must
be re-downloaded from the sources referenced in the manuscript and
placed in the Original_Data directories within the Storage directory.

Predictors that are used in the Reservoir layer, as well as the
Pathogen layer, are chosen from the list of candidate predictors using
a Wilcox significance test (described in the layers below). Candidate
predictors are listed in the following tables. More information about
each predictor can be found in the supplemental information of the
manuscript that is cited above. 

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

Each environmental dataset was downloaded in its original raster
form, and reprojected into 0.05 x 0.05 degree pixels using the WGS1984
coordinate reference system. Precipitation, temperature, and NDVI
rasters were used to generate measures of seasonality (e.g., Pcv,
Pcol, Ncv, Ncol, etc.). The code used to derive these features is
contained in the **Calc_Seasonality.r** script, contained in the
Storage/Data_Scripts/ directory. Seasonal predictors are calculated
within each year over the time span 2001 - 2019, then averaged across
years.

We also include land cover features derived from the MODIS land cover
classification dataset. Each pixel of the reprocessed rasters
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
Storage/Raster_Data/predictor_stack.grd. This is performed by the script 
**Prep_Predictors.r**. 


# Reservoir_Layer

This directory contains scripts that take in latitude/longitude
coordinates of *Mastomys natalensis* capture locations, as well as
capture locations of other rodents, and outputs a .tif that predicts
the probability that any given pixel of West Africa contains
*M. natalensis*.

The subdirectories are 

- Data: Contains the datasets that describe *M. natalensis*
  presences, as well as other rodent presences that are used as
  background data. The Excel file contains 2 spreadsheets. The first
  contains all literature records of genetically (or morphologically)
  confirmed *M. natalensis* rodents. The second spreadsheet is a
  cleaned version of the first; entries that were obtained through
  personal communication, museum specimens only, etc. are filtered
  out. A copy of this cleaned version (Mastomys_natalensis_presences.csv)
  is converted into a CSV file that is read into R script files for analysis. 
- Figures_Fits: Contains all figure outputs from the fitting process,
as well as the final tif of the prediction.
- Models: This directory is created at run time and used to
temporarily hold model files generated by Train_Reservoir_Learners.r
code. It can be deleted after the fitting process is complete.
- Tools: Contains a script with additional functions called Functions.r

To re-derive all data-processing and model fitting performed for this
layer, first set Reservoir_Layer as the working directory, then source
the main script by running:

source("Generate_Reservoir_Layer.r").

This script loads and calls three functions, defined in other scripts
detailed below, that prepare the dataset (Prep_Reservoir_Data.r),
determine which predictors significantly vary with the response
variable (Calc_Sig_Reservoir_Preds.r), and uses those predictors to
train a set of 25 boosted classification trees
(Train_Reservoir_Learners.r). Specifically,
Generate_Reservoir_Learners.r loops through sets of hyper-parameters
that specify a boosted classification tree algorithm, and passes said
hyper-parameters to the Train.Reservoir.Learners function for fitting.
All information from the Generate script is stored in a
sub-directory of Figures_Fits, [prefix], that is created when the
code is run. The name [prefix] is a variable that is designated
by the user at the top of the Generate script. 

After the fitting procedure, this script writes a dataframe to
[prefix] / hypers_res_data.csv that describes the performance of
each hyper-parameter set that was simulated. Column names should be
self-explanatory. The auc and auc.pwd columns show the out-of-bag AUC
model performance without, and with, pairwise distance sampling.

In addition, this script creates subdirectories within [prefix] that
contain further information, figures, and metrics, for each
hyper-parameter set that was run. These directories are of the form
[prefix] / [model name], where [model_name] is assigned by the
Generate script as a function of the hyperparameters. See the
script-descriptions below for details. The three main functions that
Generate_Reservoir_Learners.r calls are contained in the following
scripts:

###  Prep_Reservoir_Data.r 

This script defines the function Prep.Reservoir.Data. The function
accepts a string that specifies the focal species. For now, this
string can only take on the value "Mastomys natalensis". This function
processes *M.  natalensis* presence data. The output is a dataframe
that contains Latitude, Longitude, and associated environmental
predictors of pixels that contain *M natalensis* captures, as
well as background captures of other Murids. This script also creates
a PNG map of capture and pseudoabsence locations (saved in
Figures_Fits/[prefix] / [model_name], and named
Mastomys_natalensis_captures.png).


###  Calc_Sig_Reservoir_Preds.r

This script defines the Calc.Sig.Reservoir.Preds function. The
function accepts a species name (string), as well as the dataframe
output of Prep.Reservoir.Data, and outputs a list of predictors that
significantly vary with the *M. natalensis* response
variable. Significance is determined by a Wilcox test with a threshold
of p = 0.05. Significant predictors, in turn, will be used to train
the model that predicts the reservoir's presence or absence.
 
This function returns a list of the names of predictors that are
deemed significant (saved to
Figures_Fits/[prefix]/[model name]/Sig_Reservoir_Predictors_Mastomys_natalensis.csv),
and also saves a png figure showing histograms of the
significant predictors conditioned on response (saved to
Figures_Fits/[prefix]/[model name]/Sig_Reservoir_Predictors_Mastomys_natalensis.png).

### Train_Reservoir_Learners

This script defines a function, train.reservoir.learners, that takes
as input a dataset (output by Prep.Reservoir.Data) and a hyperparameter
list. The
function fits a set of BCT's using the significant variables determined by
Calc.Sig.Reservoir.Preds, and saves the fitted reservoir
risk layer to file. This script is the workhorse of the reservoir-model
building process. Upon running, this function creates several files to
the relevant [prefix]/[model name] directory, including 1) a .tif file
of the averaged reservoir-risk prediction generated by [nboots]
bootstrapped model fits. Here, [nboots] is a hyperparameter defined in the
Generate_Reservoir_Layer.r script; 2) Figures showing variable
importance (Coef_Mn.png), learned feature-response relationships
(Effect_Response_Mn.png), and predicted risk across West Africa
averaged over all bootstrapped fits (Risk_Layer_Mn.png); 3) a data file
tree_metrics.dat that contains information on the fitting process for
each bootstrapped fit, including "area under the receiver operator
curve" (AUC) calculated on an out-of-bag test set with and without
pairwise distance sampling, and the number of boosting iterations that
were deemed best by cross-validation; 4) two other files that describe
performance.

The tree_metrics.dat file (stored in directory [prefix]/[model_name]) contains
detailed information on the out-of-bag model performance for each bootstrapped
fit. The columns are:
- boot.i: ith bootstrapped model, where i is between 1 and [nboots]
- n.tree: number of trees used by the model 
- max.tree: maximum number of trees allowed (set in
  Generate_Reservoir_Layer.r script
- model.oob.auc: out-of-bag (OOB) AUC on test data without pairwise distance
  sampling
- model.oob.auc.pwd: OOB AUC on test data with pairwise
  distance sampling
- model.oob.auc[.pwd]: Accuracy on test data with or without pairwise sampling
- model.oob.adj[.pwd]: Adjusted accuracy on test data with or without pairwise sampling
- model.oob.mcr2[.pwd]: McFadden's R-squared on test data
- model.oob.ll[.pwd]: Log-likelihood of model fit on test data
- null.oob.ll[.pwd]: null model likelihood on test data (used in mcr2 calculation)

Two other files that describe other information on performance:

The species_predictions.dat file (stored in subdirectory of [prefix]/[model_name]) contains more
information on the bootstrapped models' belief that various background points
are Mastomys natalensis, organized around columns that describe the true species
identity of the background point. 

The test_predictions.dat file (stored in subdirectory of [prefix]/[model_name]) contains
the history of the OOB test procedure for each of the bootstrapped models.


# Pathogen_Layer

This directory contains scripts that take in data that describes the 
presence or absence of LASV in surveys of *Mastomys natalensis*, and 
outputs a .tif file that describes the probability that, given that 
*M. natalensis* is present in a pixel, LASV circulates in the rodent 
population. To the extent that is possible, the scripts and the
directory structure are similar to that found in the
Reservoir_Layer directory. 

- Data: contains the datasets that describe LASV presence and absence
in rodents, as well as data on serosurveys in humans. The Excel file
"Lassa_Occurrence_Data" contains 3 spreadsheets. The "Genbank"
spreadsheet contains LASV records in *M. natalensis* found on GenBank;
the "Raw_Lassa_Literature" spreadsheet contains LASV records in humans
and rodents that were found through a literature search; and the
"Cleaned_Lassa_Literature" spreadsheet contains cleaned version of
LASV records from literature. Cleaning involved removing
human serosurveys that were nonrandom (e.g., tested individuals with
symptoms, individual missionaries), and removing entries that did not
specifically test *M. natalensis* rodents. The Genbank and Cleaned
spreadsheets are used to create CSV files that are fed into the
model pipeline. 
- Figures_Fits: Contains all figure outputs from the fitting process,
as well as the final tif of the prediction.
- Models: This directory is created at run time and used to
temporarily hold model files generated by Train_Pathogen_Learners.r
code. It can be deleted after the fitting process is complete.
- Tools: Contains a script with additional functions called Functions.r

To re-derive all data-processing and model fitting for this layer, set
the Pathogen_Layer directory as the working directory, and then run:

source("Generate_Pathogen_Layer.r")

As with the Generate_Reservoir_Layer script, this script loads and
calls three functions that prepare the dataset (Prep_Pathogen_Data.r),
determine which predictors significantly vary with the response
variable (Calc_Sig_Pathogen_Preds.r), and uses those predictors to
train a set of 25 boosted classification trees
(Train_Pathogen_Learners.r).

After the fitting procedure, this script writes a dataframe to
[prefix] / hypers_lsv_data.csv that describes the performance of each
hyper-parameter set that was simulated. In addition, this script
creates subdirectories within [prefix] that contain further
information, figures, and metrics, for each hyper-parameter set that
was run.

### Prep_Pathogen_Data.r

This script defines the Prep.Pathogen.Data function. This function
processes Lassa presence/absence data from published surveys of
rodents and humans. The function returns a dataframe containing
*M. natalensis* Lassa presence/absence data. This rodent dataframe
includes environmental predictors that are used to train the pathogen
layer of the model. In addition, this script processes and writes to
file a dataframe describing human serosurvey data that is used to
assess and calibrate the entire model pipeline (saved to
Figures_Fits/[prefix]/[model
name]/Prepped_Human_Seroprevalence_Data.csv). In addition, this
script writes several CSV files that describe how the data were
processed (Prepped_Pathogen_Rodent_Absence_Data.csv,
Prepped_Pathogen_Ambiguous_Data.csv, and
Prepped_Pathogen_Ambiguous_Pixels.csv). This script
also produces a map figure (Human_Test_Dat.png) showing the locations
of human test data. 

###  Calc_Sig_Pathogen_Preds.r

This script is analogous to Calc_Sig_Reservoir_Preds.r. When run, this
script defines a function called Calc.Sig.Pathogen.Preds that returns
a list of those predictors that are significantly associated with the
presence or absence of the pathogen. These are the predictors that the
boosted models are trained on. This list of predictors is saved to
Figures_Fits/[prefix]/[model name]/Sig_Pathogen_Predictors.csv, along
with a png figure showing histograms of the predictors conditioned on
response outcome.


### Train_Pathogen_Learners.r

This script defines a function, train.pathogen.learners, that takes as
input a dataset and hyperparameter list, and outputs the fitted
pathogen risk layer.  This script is the workhorse of the
pathogen-model building process, and is called by
Generate_Pathogen_Layer.r. As with the Reservoir layer analogue, this
script writes several files to the relevant [prefix]/[model name]
directory, including 1) .tif file of the averaged pathogen-risk
prediction generated by [nboots] bootstrapped model fits (named
Lassa_Layer_[model_name]). Here, [nboots] is a hyperparameter defined
in the Generate_Pathogen_Layer.r script; 2) Figures showing variable
importance (Coef_Lassa.png), learned feature-response relationships
(Effect_Response_Lassa.png), and predicted risk of LASV (given the
presence of *Mastomys natalensis*) for all of West Africa
(Lassa_Risk_Layer.png); 3) Data file tree_metrics.dat that contains
information on the fitting process for each bootstrapped fit,
including "area under the receiver operator curve" (AUC) calculated on
an out-of-bag test set, and the number of boosting iterations that
were deemed best by cross-validation; 4) An additional file
(test_predictions.dat) that describes the OOB test procedure for each
bootstrapped model.


### Map_Cases.r

This script produces a ggplot map that indicates the point locations
of rodent surveys for Lassa that were used to train the pathogen
layer, as well as human sero-surveys that were used to calibrate the
risk layers.


# Human_LASV_Incidence

This directory contains script files that combine the Reservoir and
Pathogen layers above into a prediction of LASV spillover risk to
humans and evaluates the prediction against actual human serosurvey data.
Predictions of human LASV seroprevalence and case rate across West
Africa are also produced. 

### Human_Link.r

This script file 1) tests the combined LASV risk layer's ability to
predict human seroprevalence data, and 2) uses a GLM regression to
translate combined LASV risk predictions into predictions of human
LASV seroprevalence and LASV case rate in West Africa. At the top of
the script, a user specifies a [grid.name] string; figures and text
output are saved in the Figures_Fits/[grid.name] folder.

First, the script calculates the combined risk layer D_X. Let $D_M$
and $D_L$ denote the reservoir and pathogen risk layer outputs of the
Generate_*_Layer.r scripts. Each of these risk layers are the averages
of [nboots] model fits. Calculate the combined risk raster, $D_X$,
as
$$
D_X = D_M \times D_L. 
$$

Next, the script evaluates several metrics that test the combined
risk layer's ability to correlate with human serosurvey data. The
metrics are correlation, weighted correlation, and the p-value that
results when the D_X layer is used in a quasi-binomial regression
with human seroprevalence as a response variable. 

These metrics, in turn, are written to
Figures_Fits/[grid.name]/metrics_output.csv. Columns of this data
frame are correlation of predicted seroprevalence and actual ("corr"),
correlation weighted by number of humans tested ("corr.weighted"),
p-value on Dx predictor in the GLM ("glm.pval"), and fraction of null
deviance explained in GLM ("glm.frac.dev").

Finally, the script use the GLM regression to predict the
seroprevalence and spillover rate of LASV across West Africa.

During processing, several PNG figures are created. These show the
predicted human seroprevalence across West Africa
(Lassa_Seroprevalence.png), the predicted rate of new cases of Lassa
in humans (New_Cases.png), the combined risk layer
(Combined_Risk_Layer.png), the deviance residuals resulting from the
regression of human seroprevalence on the combined risk layer
(Deviance_Residuals.png), and the predicted vs actual seroprevalence
(Pred_vs_Fit.png).

This script also outputs the data files:
- foc_case: summarizes the yearly spillover rate of LASV into humans
by country (Table 4 of the manuscript). 
- glm_summary: text output of GLM summary

### Human_Link_Analysis.cdf

This is a Mathematica file that derives the expressions found in 
the manuscript. 




