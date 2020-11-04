## This script loops through all candidate predictors, and returns a
## list of those predictors that are significantly associated with
## the presence or absence of the pathogen. Here, significance is
## determined by a Wilcox test with a threshold of p = 0.05.
## Significant predictors, in turn, will be used to train the
## the model that predicts the pathogen's presence or absence.

## This script outputs 1) a list of the names of predictors that
## are deemed significant (Figures_Fits/Sig_LASV_preds), and a
## pdf figure showing the distributions of the significant
## predictors.

Calc.Sig.Pathogen.Preds <- function(dataset){
    
    ## Add column that describes the status of each survey location: lassa present or absent
    dataset$Lsv_Status <- ifelse(dataset$ArenaStat, 'Present', 'Absent')
    dataset$Lsv_Status = factor(dataset$Lsv_Status,
                                        levels = c('Absent', 'Present'), ordered = TRUE)
    
    ## Load in all candidate variables
    var.names <- names(all.stack)
    
    ## Perform Wilcox test on each predictor in the candidate set
    pvalues <- data.frame(var = var.names, p = NA)
    for(i in 1:nrow(pvalues)){
        pvalues[i,'p'] = wilcox.test(get(paste(pvalues[i,'var']))~Lsv_Status, dataset)$p.value
    }
    
    ## Keep predictors that are significant at the p < 0.05 level
    pvalues <- pvalues[order(pvalues$p, decreasing = FALSE),]
    pvalues <- pvalues[pvalues$p < 0.05, ]
    sig.var.names <- paste(na.omit(pvalues$var))
    
    ## Make plot of significant predictors
    graph.list <- list()
    for(i in 1:length(sig.var.names)){
        var = sig.var.names[i]
        graph.list[[which(var==sig.var.names)]] <- local({
            i = i
            p<-ggplot(dataset, aes(x=dataset[,sig.var.names[i]], fill=Lsv_Status)) +
                geom_density(alpha=0.4) + xlab(pretty.labels(sig.var.names[i]))
            print(p)
        })
    }
    pdf(file = paste0('Figures_Fits/', prefix, '/',fold,'/',
                      'Sig_Pathogen_Predictors.pdf'),
        width = 15, height = 10)
    grid.arrange(grobs = graph.list)
    dev.off()
    
    ## Write significant predictor names to file.
    write.csv(sig.var.names, file = paste0('Figures_Fits/', prefix, '/',fold,'/',
                                             'Sig_Pathogen_Preds.csv'),
                row.names = FALSE)

    return(sig.var.names)
    
} ## End Function
    
