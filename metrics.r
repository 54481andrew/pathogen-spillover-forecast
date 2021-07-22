## The goal of this script is to calculate some additional metric across bootstrap fits. 

## --- Pathogen layer

## Load in the pathogen tree.dat data file that holds various metrics
path.fold <- 'Pathogen_Layer/Figures_Fits/pathogen_v5/pa_nboots_25_1_mllr_3_lmt_7_ambi_NA_mintest_5/'
path.tree <- read.table(paste0(path.fold, 'tree.dat'), header = TRUE)

## mean accuracy
mean(path.tree$model.acc)
## 0.82

## mean adjusted accuracy
mean(path.tree$adj.acc)
## 0.64

## mean auc
mean(path.tree$model.auc)
## 0.85

## ------------------------

## Calculate statistics for reservoir layer

res.fold <- 'Reservoir_Layer/Figures_Fits/reservoir_v5/Mn_pa_nboots_25_nbg_Same_tc_1_mllr_2_lmt_7/'
res.tree <- read.table(paste0(res.fold, 'tree.dat'), header = TRUE)

## mean accuracy
mean(res.tree$model.oob.acc.pwd)
## 0.65

## mean adjusted accuracy
mean(res.tree$model.oob.adj.acc.pwd)
## 0.27

## mean auc
mean(res.tree$model.oob.auc.pwd)
## 0.68
