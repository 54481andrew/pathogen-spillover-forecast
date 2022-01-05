## This script provides a simple means of loading and printing the
## metrics of a reservoir and pathogen layer, as well as the performance
## of the combined risk layer (Dx) on human seroprevalence. 

## --- Pathogen layer

## Load in the pathogen tree.dat data file that holds various metrics
path.fold <- 'Pathogen_Layer/Figures_Fits/pathogen_v6/pa_nboots_25_1_mllr_3_lmt_7_ambi_NA_mintest_5/'
path.tree <- read.table(paste0(path.fold, 'tree_metrics.dat'), header = TRUE)

writeLines('-- Pathogen layer metrics -- \n')

## mean accuracy
writeLines(paste0('Mean accuracy: ', mean(path.tree$model.oob.acc)))
## 0.82

## mean adjusted accuracy
writeLines(paste0('Mean adjusted accuracy: ', mean(path.tree$model.oob.adj.acc)))
## 0.64

## mean auc
writeLines(paste0('Mean AUC: ', mean(path.tree$model.oob.auc)))
## 0.85

## ------------------------

## Calculate statistics for reservoir layer
writeLines('\n\n-- Reservoir layer metrics -- \n')

res.fold <- 'Reservoir_Layer/Figures_Fits/reservoir_v6/Mn_pa_nboots_25_nbg_Same_tc_1_mllr_2_lmt_7/'
res.tree <- read.table(paste0(res.fold, 'tree_metrics.dat'), header = TRUE)

## mean accuracy
writeLines(paste0('Mean accuracy (pairwise distance sampling): ', mean(res.tree$model.oob.acc.pwd)))
## 0.65

## mean adjusted accuracy
writeLines(paste0('Mean adjusted accuracy (pairwise distance sampling): ',
                  mean(res.tree$model.oob.adj.acc.pwd)))
## 0.27

## mean auc calculated with pairwise distance sampling
writeLines(paste0('Mean AUC (pairwise distance sampling): ', mean(res.tree$model.oob.auc.pwd)))
## 0.68

## ------------------------

## Calculate statistics for the human layer

writeLines('\n\n-- Test metrics for human seroprevalence -- \n')

hum.fold <- 'Human_LASV_Incidence/Figures_Fits/human_v6_ambiNA_mint_5/metrics_output.csv'
hum.dat <- read.csv(hum.fold)

## weighted correlation of human seroprevalence with Dx
writeLines(paste0('Weighted correlation with combined risk Dx: ', hum.dat$corr.weighted))
## 0.398

## correlation of human seroprevalence with Dx
writeLines(paste0('Correlation with combined risk Dx: ', hum.dat$corr))
## 0.327

## p-value of regression of human seroprevalence onto combined LASV risk (Dx)
writeLines(paste0('P-value in GLM regression: ', hum.dat$glm.pval))
## 0.000123

