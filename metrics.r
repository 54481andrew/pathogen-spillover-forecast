## The goal of this script is to collect assessment metrics for the pathogen,
## reservoir, and human layers. 

## --- Pathogen layer

## Load in the pathogen tree.dat data file that holds various metrics
path.fold <- 'Pathogen_Layer/Figures_Fits/pathogen_v7/pa_nboots_25_1_mllr_3_lmt_7_ambi_NA_mintest_5/'
path.tree <- read.table(paste0(path.fold, 'tree_metrics.dat'), header = TRUE)

## mean accuracy
mean(path.tree$model.oob.acc)
## 0.82

## mean adjusted accuracy
mean(path.tree$model.oob.adj.acc)
## 0.64

## mean auc
mean(path.tree$model.oob.auc)
## 0.85

## ------------------------

## Calculate statistics for reservoir layer

res.fold <- 'Reservoir_Layer/Figures_Fits/reservoir_v7/Mn_pa_nboots_25_nbg_Same_tc_1_mllr_2_lmt_7/'
res.tree <- read.table(paste0(res.fold, 'tree.dat'), header = TRUE)

## mean accuracy
mean(res.tree$model.oob.acc.pwd)
## 0.65

## mean adjusted accuracy
mean(res.tree$model.oob.adj.acc.pwd)
## 0.27

## mean auc calculated with pairwise distance sampling
mean(res.tree$model.oob.auc.pwd)
## 0.68

## ------------------------

## Calculate statistics for the human layer
hum.fold <- 'Human_LASV_Incidence/human_v6_ambiNA_mint_5/metrics_output.csv'
hum.dat <- read.csv(hum.fold)


hum.dat$glm.pval

hum.dat$corr.weighted
