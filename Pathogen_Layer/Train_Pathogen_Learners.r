## This script defines a function, train.pathogen.learners, that takes
## as input a data-set and hyperparameter list, and outputs the fitted
## pathogen risk layer.  This script is the workhorse of the
## pathogen-model building process, and is called by
## Generate_Pathogen_Layer.r.  Outputs are saved in the Figures_Fits/
## directory, and include:
## 1. a .tif file of the averaged pathogen-risk prediction generated
## by [nboots] bootstrapped model fits. [nboots] is a hyperparameter
## defined in the Generate_Pathogen_Layer.r script.
## 2. Figures showing variable importance, learned relationships, and
## predicted risk across West Africa, averaged over all bootstrapped
## fits.
## 3. Data file that contains information on the fitting process for
## each bootstrapped fit, including "area under the receiver operator
## curve" (AUC) calculated on an out-of-bag test set, and the number of
## trees-building iterations that were deemed best by cross-validation.

train.pathogen.learners <- function(rodlsv.survey.dat, gridgi = NULL){

    ## Extract variables from gridgi
    nboots = gridgi[,'nboots']
    tree.complexity = gridgi[,'tree.complexity']
    learning.rate <- 10^(-gridgi[,'mllr']) ## Learning rate in the model
    max.trees <- 10^gridgi[,'lmt'] ## Max number of trees in the model: 100000

    ## Generate directory name where all output will be saved
    fold <- generate.lassa.name(gridgi)
    cat(paste0('\n\n\n\n'))
    print(paste('--------- Model fit name:', fold, '-------------'), quote = FALSE)

    ## Create data directories
    dirpath <- paste('Figures_Fits/',prefix,sep='')
    if(!dir.exists(dirpath)){dir.create(dirpath, showWarnings = FALSE)}
    dirpath <- paste('Figures_Fits/', prefix, '/', fold,sep='')
    if(!dir.exists(dirpath)){dir.create(dirpath, showWarnings = FALSE)}
    models.folder = paste0('Figures_Fits/', prefix, '/', fold, '/Models')
    unlink(models.folder, recursive = TRUE)
    dir.create(models.folder,showWarnings = FALSE)

    num.survey <- nrow(rodlsv.survey.dat)

    mcl.fun <- function(boot.i, verbose = FALSE){

        ## Choose a bootstrap sample of training dat
        set.seed(boot.i)

        ## Sample equal number of 1's and 0's for training data
        wi.abs <- which(rodlsv.survey.dat$ArenaStat==0)
        wi.pres <- which(rodlsv.survey.dat$ArenaStat==1)
        num.abs <- length(wi.abs)
        num.pres <- length(wi.pres)
        num.samp.each <- min(num.abs, num.pres)
        wi.train.abs <- sample(wi.abs, size = num.samp.each, replace = TRUE)
        wi.train.pres <- sample(wi.pres, size = num.samp.each, replace = TRUE)
        wi.train.survey <- c(wi.train.abs, wi.train.pres)
        train.survey.data <- rodlsv.survey.dat[wi.train.survey, ]

        ## Sample a test set from out of bag samples
        test.survey.data <- rodlsv.survey.dat[-wi.train.survey,]
        wi.test.abs <- which(test.survey.data$ArenaStat==0)
        wi.test.pres <- which(test.survey.data$ArenaStat==1)
        num.test.abs <- length(wi.test.abs)
        num.test.pres <- length(wi.test.pres)
        num.test.each <- min(num.test.abs, num.test.pres)
        wi.test.abs <- sample(wi.test.abs, size = num.test.each, replace = FALSE)
        wi.test.pres <- sample(wi.test.pres, size = num.test.each, replace = FALSE)
        wi.test.survey <- c(wi.test.abs, wi.test.pres)
        test.survey.data <- test.survey.data[wi.test.survey, ]

        ## Fit the boosted regression tree model
        ## We use a while loop to retry the gbm.step process until a model
        ## is returned. This thread will throw an error after 10 failed
        ## attempts.
        gbm.fit.success <- FALSE
        ntries <- 0
        while(!gbm.fit.success & ntries < 10){
            gbm.mod <- gbm.step(data = train.survey.data,
                                gbm.x = var.names,
                                gbm.y = 'ArenaStat',
                                tree.complexity = tree.complexity,
                                learning.rate = learning.rate, bag.fraction = 0.75,
                                max.trees = max.trees,
                                verbose = verbose, silent = !verbose,
                                plot.main = FALSE)
            if(class(gbm.mod)!='NULL'){
                gbm.fit.success <- TRUE
            }else{
            }
            ntries = ntries + 1
        } ## End While

        ## Test performance on out-of-bag presences and pseudoabsence points

        ## Predict and assess gbm.mod on test.survey.data
        pred.rast <- predict(all.stack[[var.names]], gbm.mod, n.trees = gbm.mod$n.trees,
                             type = 'response')
        model.predictions <- extract(x = pred.rast, y = test.survey.data[,c('Longitude', 'Latitude')])
        rocArea.model = NA
        rocArea.model = roc.area(test.survey.data$ArenaStat,
                                      model.predictions)$A

        ## Save fit statistics to file
        tree.dat <- data.frame(boot.i = boot.i, n.tree = gbm.mod$n.trees, max.tree = max.trees,
                               model.auc = rocArea.model)
        write.table(tree.dat, file = tree.filename,
                    col.names = !file.exists(tree.filename), row.names = FALSE,
                    append = file.exists(tree.filename))

        ## Save model and tif prediction of this bootstrap fit
        mod.filename = paste(models.folder, '/amod_',boot.i,'.rds', sep = '')
        saveRDS(gbm.mod, file = mod.filename)

        tif.filename = paste(models.folder, '/amod_',boot.i,'.tif', sep = '')
        writeRaster(pred.rast, filename = tif.filename, overwrite = TRUE)

        print(paste('-- Finished Bootset: ', boot.i, '; Elapsed Time: ', Sys.time() - starttime), quote = FALSE)
        print(paste('----- gbm.mod.trees', gbm.mod$n.trees, '   max.trees: ', max.trees), quote = FALSE)
        gc() ## helps clear memory
        return(boot.i)
    }## End mcl.fun


    tree.filename = paste0('Figures_Fits/', prefix, '/',
                           fold, '/tree.dat')
    unlink(tree.filename)

    print('Bootstrapping model fits', quote = FALSE)

    starttime <- Sys.time()
    out <- mclapply(1:nboots,
                    mcl.fun, mc.cores = detectCores() - 2)
    print(paste('Bootstrapping finished; Total Time: ', Sys.time() - starttime), quote = FALSE)

    ## Aggregate and analyze the fitted models. Read in all raster predictions, average
    ## them, and save the resulting meta prediction.
    pred.stack.names <- list.files(models.folder , pattern = "*.tif$")
        fullname <- paste(models.folder,pred.stack.names, sep = '/')
    pred.stack <- stack(fullname)
    pred.rast <- overlay(pred.stack, fun = mean)
    writeRaster(pred.rast, file = paste("Figures_Fits/", prefix, '/', fold,"/Lassa_Layer_",
                                        fold,".tif", sep = ''), overwrite = TRUE)

    ## The code below extracts variable importance ranking and the learned relationships
    imp.mat <- matrix(NA, nrow = length(var.names), ncol = nboots)
    imp.dat <- data.frame(coef = NA, imp = NA)
    response.dat <- c() ## Store the effect of each variable in here

    for(boot.i in 1:nboots){
        mod.filename = paste(models.folder,'/amod_',boot.i,'.rds', sep = '')
        if(file.exists(mod.filename)){
                gbm.mod <- readRDS(file = mod.filename)
                summ <- data.frame(summary(gbm.mod))
                var.names.boot <- paste(summ$var)
                imp.mat[as.vector(sapply(var.names.boot,
                                         FUN = function(x){which(x==var.names)})), boot.i] <-
                    summ[var.names.boot,'rel.inf']
                imp.dat <- rbind(imp.dat, data.frame(coef = var.names.boot, imp = summ[var.names.boot,'rel.inf']))
                print(paste('Extracted boot', boot.i, '/', nboots), quote = FALSE)

                ## Save info on the learned relationship for each predictor
                for(pred.i in 1:length(var.names.boot)){
                    gbm.call <- gbm.mod$gbm.call
                    new.dat <- data.frame(gbm::plot.gbm(gbm.mod, pred.i,
                                                        return.grid = TRUE, type = 'response'),
                                          boot = boot.i,
                                          var = gbm.call$predictor.names[pred.i])
                    names(new.dat) <- c('x', 'y', 'boot', 'var')
                    response.dat <- rbind(response.dat,
                                          new.dat)
                }## Loop through predictors
        } ## Check if file exists
    }## End loop through boots

    ## Extract top 10 most important predictors for a boxplot.
    row.names(imp.mat) <- var.names
    med.mat <- apply(imp.mat, 1, median, na.rm = TRUE)
    ord.names <- names(sort(med.mat, TRUE))
    imp.dat.1 <- na.omit(imp.dat[-1,])
    imp.dat.2 <- na.omit(imp.dat.1)
    wi <- which(imp.dat.2$coef %in% ord.names[1:10])
    imp.dat.2 <- imp.dat.2[wi,]
    imp.dat.2$coef <- factor(imp.dat.2$coef, levels = rev(ord.names[1:10]), ordered = TRUE)

    ## Build box plot
    p <- ggplot(imp.dat.2, aes(x=coef, y=imp)) +
        geom_boxplot()
    p + scale_x_discrete(labels=rev(pretty.labels(ord.names[1:10]))) +
        labs(x = 'Predictor', y = 'Importance') + coord_flip()##ylim = c(0,35))
    ggsave(filename = paste('Figures_Fits/', prefix, '/', fold ,'/Coef_Lassa.png', sep = ''), device = 'png')

    ## Plot the learned relationships across all models
    response.dat$var <- factor(response.dat$var, levels = ord.names, ordered = TRUE)
    vars.to.plot <- ord.names[1:6]
    vars.to.plot.pretty <- pretty.labels(vars.to.plot)
    response.dat$var.pretty <- factor(pretty.labels(response.dat$var),
                                      levels = pretty.labels(ord.names),
                                      ordered = TRUE)
    response.dat$Effect = response.dat$y
    response.dat$Value = response.dat$x

    ggplot(data = response.dat[response.dat$var.pretty %in% vars.to.plot.pretty,]) +
        geom_line(aes(x = Value, y = Effect, group = boot), size = 0.05) +
        theme_classic() +
        stat_summary_bin(aes(x = Value, y = Effect, group = var.pretty),
                         geom = 'line', size = 2, fun = mean, fun.args = list(na.rm=TRUE)) +
        facet_wrap(~var.pretty, ncol = 3, nrow = 2, scales = 'free') + xlab('Predictor Value') +
        ylab('Classification Score')
    ggsave(filename = paste('Figures_Fits/', prefix, '/', fold,
                            '/Effect_Response_Lassa.png', sep = ''),
           device = 'png', width = 7, height = 5, units = 'in')

    ggplot(data = response.dat[response.dat$var.pretty %in% vars.to.plot.pretty,]) +
        stat_summary_bin(aes(x = Value, y = Effect, group = var.pretty),
                         fun.data = mean_se, fun.args = list(mult = 2),
                         color = "black", fill = 'blue', linetype = 'dashed',
                         geom = 'ribbon', bins = 30) +
        theme_classic() +
        stat_summary_bin(aes(x = Value, y = Effect, group = var.pretty),
                         geom = 'line', size = 1, fun = mean) +
        facet_wrap(~var.pretty, ncol = 3, nrow = 2, scales = 'free') + xlab('Predictor Value') +
        ylab('Classification Score')
    ggsave(filename = paste('Figures_Fits/', prefix, '/', fold,
                            '/Effect_Response_Lassa_var.png', sep = ''),
           device = 'png', width = 7, height = 5, units = 'in')
    
    
    ## Plot risk map averaged over all boot predictions
    heat.cols <- viridis(120, begin = 0.1, end = 1, option = 'D')
    xlims = c(-18,16)
    ylims = c(16, 16.5)
    points.lasv <- rodlsv.survey.dat[,c('Longitude', 'Latitude')]
    pred.rast = mask(pred.rast, masto.rangemap, updatevalue = 0)
    png(file = paste('Figures_Fits/', prefix, '/', fold, '/Lassa_Risk_Layer.png', sep = ''),
        width = 6, height = 4, units = 'in', res = 400)
    par(mai = 1*c(0.2,0.2,0.2,0.6))
    image.plot(pred.rast, col = heat.cols, zlim = c(0, 1),
               bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',
               xlim = xlims, ylim = ylims,
               asp = 1, legend.lab = 'Occurrence score', legend.line = 2.5,
               main = '')
    mtext(text = expression('Lassa Distribution ('~'D'['L']~')'), side = 3, line = -1)
    points(points.lasv[,1], points.lasv[,2], asp = 1, cex = 1, pch = 21,
           lwd = 1, bg = c('white', 'red')[rodlsv.survey.dat[,'ArenaStat'] + 1], col = 'black')
    plot(foc.shp.ogr, add = TRUE, bty = 'n', asp = 1)
    plot(rgeos::gIntersection(foc.shp.ogr, masto.rangemap), add = TRUE, bty = 'n', asp = 1, lwd = 3)
    legend(x = 'bottomleft', legend = c('LASV -', 'LASV +'),
           pt.bg = c('white', 'red'), pch = 21, pt.lwd = 1, col = 'black', cex = 1,
           bty = 'n')
    dev.off()

    ## Remove fitted models
    unlink('Models', recursive = TRUE)

}## End function
