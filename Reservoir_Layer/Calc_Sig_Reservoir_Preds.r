## This script loops through all candidate predictors, and outputs a list
## of those predictors that are significantly associated with the
## presence or absence of the reservoir. Here, absence refers to
## background captures of Murid rodents. Significance is determined
## by a Wilcox test with a threshold of p = 0.05.  Significant
## predictors, in turn, will be used to train the the model that predicts
## the reservoir's presence or absence.

## This script outputs 1) a list of the names of predictors that are
## deemed significant (saved to Figures_Fits/Sig_Reservoir_preds), and a
## pdf figure showing the distributions of the significant predictors.

Calc.Sig.Reservoir.Preds <- function(Species, dataset){
    
    ## Variant of species name used for filename
    Species.name = gsub(' ', '_', Species)
    
    ## Add column that describes the type of rodent: Species or other
    Abbrev.name <- paste0(substr(Species,1,1), substr(strsplit(Species,' ')[[1]][2],1,1))
    
    dataset$Rodent <- ifelse(dataset$Species == Species, Abbrev.name, 'Other')
    dataset$Rodent = factor(dataset$Rodent, levels = c('Other', Abbrev.name),
                                   ordered = TRUE)
    
    ## Load in all candidate variables
    var.names <- names(all.stack)
    
    ## Perform Wilcox test on each predictor in the candidate set
    pvalues <- data.frame(var = var.names, p = NA)
    for(i in 1:nrow(pvalues)){
        pvalues[i,'p'] = wilcox.test(get(paste(pvalues[i,'var']))~Rodent, dataset)$p.value
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
            p<-ggplot(dataset, aes(x=dataset[,sig.var.names[i]], fill=Rodent)) +
                geom_density(alpha=0.4) + xlab(pretty.labels(sig.var.names[i]))
            print(p)
        })
    }
    png(file = paste0('Figures_Fits/', prefix, '/',fold,'/',
                      'Sig_Reservoir_Predictors_',Species.name,'.png'),
        width = 25, height = 15, units = 'in', res = 400)
    grid.arrange(grobs = graph.list)
    dev.off()
    
    ## Write significant predictor names to file.
    write.csv(sig.var.names, file = paste0('Figures_Fits/', prefix, '/',fold,'/',
                                             'Sig_Reservoir_Predictors_',Species.name,'.csv'),
                row.names = FALSE)

    return(sig.var.names)

} ## End Function
