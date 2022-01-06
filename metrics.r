## This script provides a simple means of loading and printing the
## metrics of a reservoir and pathogen layer, as well as the performance
## of the combined risk layer (Dx) on human seroprevalence. 

## See bottom of this script for metric outputs for v6 and v7 model
## sets. 

## --- Pathogen layer

## Load in the pathogen tree.dat data file that holds various metrics
path.fold <- 'Pathogen_Layer/Figures_Fits/pathogen_v7/pa_nboots_25_1_mllr_3_lmt_7_ambi_NA_mintest_5/'
path.tree <- read.table(paste0(path.fold, 'tree_metrics.dat'), header = TRUE)

writeLines('-- Pathogen layer metrics -- \n')

## mean accuracy
writeLines(paste0('Mean accuracy: ', mean(path.tree$model.oob.acc)))

## mean adjusted accuracy
writeLines(paste0('Mean adjusted accuracy: ', mean(path.tree$model.oob.adj.acc)))

## mean auc
writeLines(paste0('Mean AUC: ', mean(path.tree$model.oob.auc)))

## ------------------------

## Calculate statistics for reservoir layer
writeLines('\n\n-- Reservoir layer metrics -- \n')

res.fold <- 'Reservoir_Layer/Figures_Fits/reservoir_v7/Mn_pa_nboots_25_nbg_Same_tc_1_mllr_2_lmt_7/'
res.tree <- read.table(paste0(res.fold, 'tree_metrics.dat'), header = TRUE)

## mean accuracy
writeLines(paste0('Mean accuracy (pairwise distance sampling): ', mean(res.tree$model.oob.acc.pwd)))

## mean adjusted accuracy
writeLines(paste0('Mean adjusted accuracy (pairwise distance sampling): ',
                  mean(res.tree$model.oob.adj.acc.pwd)))

## mean auc calculated with pairwise distance sampling
writeLines(paste0('Mean AUC (pairwise distance sampling): ', mean(res.tree$model.oob.auc.pwd)))

## ------------------------

## Calculate statistics for the human layer

writeLines('\n\n-- Test metrics for human seroprevalence -- \n')

hum.fold <- 'Human_LASV_Incidence/Figures_Fits/human_v7_ambiNA_mint_5/metrics_output.csv'
hum.dat <- read.csv(hum.fold)

## weighted correlation of human seroprevalence with Dx
writeLines(paste0('Weighted correlation with combined risk Dx: ', hum.dat$corr.weighted))

## correlation of human seroprevalence with Dx
writeLines(paste0('Correlation with combined risk Dx: ', hum.dat$corr))

## p-value of regression of human seroprevalence onto combined LASV risk (Dx)
writeLines(paste0('P-value in GLM regression: ', hum.dat$glm.pval))


## ------------------- Output for v6 models (R version 4.0.0 "Arbor Day")

## -- Pathogen layer metrics -- 

## Mean accuracy: 0.820890331890332
## Mean adjusted accuracy: 0.641780663780664
## Mean AUC: 0.850716589779469


## -- Reservoir layer metrics -- 

## Mean accuracy (pairwise distance sampling): 0.648803992947344
## Mean adjusted accuracy (pairwise distance sampling): 0.273635149795386
## Mean AUC (pairwise distance sampling): 0.675479772809565


## -- Test metrics for human seroprevalence -- 

## Weighted correlation with combined risk Dx: 0.397501639879647
## Correlation with combined risk Dx: 0.327204008044843
## P-value in GLM regression: 0.000122538870794965


## ------------------- Output for v7 models (R version 4.1.2 "Bird Hippie")

## -- Pathogen layer metrics -- 

## Mean accuracy: 0.820890331890332
## Mean adjusted accuracy: 0.641780663780664
## Mean AUC: 0.851286002729834


## -- Reservoir layer metrics -- 

## Mean accuracy (pairwise distance sampling): 0.64845968047317
## Mean adjusted accuracy (pairwise distance sampling): 0.273395823630318
## Mean AUC (pairwise distance sampling): 0.675643459068247


## -- Test metrics for human seroprevalence -- 

## Weighted correlation with combined risk Dx: 0.397055147402571
## Correlation with combined risk Dx: 0.326765348490805
## P-value in GLM regression: 0.00012397103506192
