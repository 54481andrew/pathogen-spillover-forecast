## This script defines a function, train.reservoir.learners, that takes
## as input a data-set and hyperparameter list, and outputs the fitted
## reservoir risk layer.  This script is the workhorse of the
## reservoir-model building process, and is called by
## Generate_Reservoir_Layer.r.  Outputs are saved in the Figures_Fits/
## directory, and include:
## 1. a .tif file of the averaged reservoir-risk prediction generated
## by [nboots] bootstrapped model fits. [nboots] is a hyperparameter
## defined in the Generate_Reservoir_Layer.r script.
## 2. Figures showing variable importance, learned relationships, and
## predicted risk across West Africa, averaged over all bootstrapped
## fits.
## 3. Data file that contains information on the fitting process for
## each bootstrapped fit, including "area under the receiver operator
## curve" (AUC) calculated on an out-of-bag test set, and the number of
## trees-building iterations that were deemed best by cross-validation.

Train.Reservoir.Learners <- function(dataset, gridgi = NULL){

    ## Define gridgi if not provided
    if(length(gridgi)==0){
        gridgi <- data.frame(Species = 'Mastomys natalensis',
                             nboots = 25,
                             num.bg.points = "Same",
                             tree.complexity = 1,
                             mllr = 4,
                             lmt = 7)

    }

    ## Extract variables from gridgi
    Species = paste(gridgi$Species)
    Abbrev.name <- paste0(substr(paste(Species),1,1),
                          substr(strsplit(paste(Species),' ')[[1]][2],1,1))

    nboots = gridgi$nboots
    tree.complexity = gridgi$tree.complexity
    num.bg.points = gridgi$num.bg.points
    learning.rate <- 10^-gridgi[,'mllr']
    max.trees <- 10^gridgi[,'lmt']

    ## Generate directory name where all output will be saved
    fold <- generate.res.name(gridgi)
    cat(paste0('\n\n\n\n'))
    print(paste('--------- Model fit name:', fold, '-------------'), quote = FALSE)

    ## Set up data directories
    dirpath <- paste('Figures_Fits/',prefix,sep='')
    if(!dir.exists(dirpath)){dir.create(dirpath, showWarnings = FALSE)}
    dirpath <- paste('Figures_Fits/', prefix, '/', fold,sep='')
    if(!dir.exists(dirpath)){dir.create(dirpath, showWarnings = FALSE)}
    models.folder = paste0('Figures_Fits/', prefix, '/', fold, '/Models')
    unlink(models.folder, recursive = TRUE)
    dir.create(models.folder,showWarnings = FALSE)

    ## Model fit statistics are stored here
    tree.filename = paste0('Figures_Fits/', prefix, '/',
                           fold, '/tree.dat')
    assess.filename = paste0('Figures_Fits/', prefix, '/',
                           fold, '/assess.dat')

    ## Divide the dataset into two dataframes: the first only contains presences,
    ## the second only background
    presence.data <- dataset[dataset$Presence==1,]
    background.data <- dataset[dataset$Presence==0,]

    mcl.fun <- function(boot.i, verbose = FALSE){
        ## Set random seed so that runs are reproducible
        set.seed(boot.i)

        ## Sample presence data
        wi.train.pres <- sample(nrow(presence.data), replace = TRUE)
        train.presence.data <- presence.data[wi.train.pres,]
        ## If num.bg.points is not a number, set it equal to the number of presence points
        if(!is.numeric(num.bg.points)){num.bg.points = nrow(train.presence.data)}

        ## Choose background points
        wi.train.abs <- sample(nrow(background.data), replace = TRUE, size = num.bg.points)
        train.abs.data <- background.data[wi.train.abs,]

        ## Form full dataset of background and presences
        train.dataset = rbind.fill( train.presence.data, train.abs.data )

        ## Weight presences and absences equally
        num.abs.tot <- sum(train.dataset$Presence==0)
        num.pres <- sum(train.dataset$Presence==1)

        ## Fit the boosted regression tree model on the train.dataset. Occasionally,
        ## the fitting process does not converge. The while loop ensures that the
        ## fitting process will be attempted 10 times, after which an error will
        ## be returned.
        gbm.fit.success <- FALSE
        ntries <- 0
        while(!gbm.fit.success & ntries < 10){
            gbm.mod <- gbm.step(data = train.dataset,
                                gbm.x = var.names,
                                gbm.y = "Presence",
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

        ##---

        ## Test performance on out-of-bag presences and pseudoabsence points.
        ## Assess AUC both with and without pairwise distance sampling.

        ## Obtain out-of-bag presences
        wi.test.pres <- (1:nrow(presence.data))[-wi.train.pres]
        test.presences <- presence.data[wi.test.pres,]

        ## Obtain out-of-bag background points
        wi.test.abs <- (1:nrow(background.data))[-wi.train.abs]
        all.test.background <- background.data[wi.test.abs,]

        ## Form a prediction raster across West Africa
        pred.rast <- predict(all.stack[[var.names]], gbm.mod, n.trees = gbm.mod$n.trees,
                             type = 'response')

        ## Set up a dataframe to store things in
        tree.dat <- data.frame(boot.i = boot.i, n.tree = gbm.mod$n.trees, max.tree = max.trees,
                               model.oob.auc = NA, model.oob.auc.pwd = NA, n.oob.auc = NA,
                               n.oob.auc.pwd = NA)

        ## Method 1: Form the test dataset without pairwise distance sampling
        wi.test.abs.1 <- sample(nrow(all.test.background),
                                length(wi.test.pres),
                                replace = length(wi.test.pres) > nrow(all.test.background))
        test.background.points <- all.test.background[wi.test.abs.1,]
        test.dataset = rbind.fill(presence.data[wi.test.pres,], test.background.points)

        ## Extract and predict on test.dataset
        model.predictions <- extract(x = pred.rast,
                                     y = test.dataset[,c('Longitude', 'Latitude')])
        tree.dat$model.oob.auc <- roc.area(test.dataset$Presence, model.predictions)$A
        tree.dat$n.oob.auc <- nrow(test.dataset)

        ## Method 2: Use pairwise-distance sampling to find appropriate set of test absences
        wi.test.background <- pwdSample(fixed = test.presences[,c('Longitude', 'Latitude')],
                                        sample = all.test.background[,c('Longitude', 'Latitude')],
                                        reference = train.presence.data[ ,c('Longitude',
                                                                            'Latitude')],
                                        lonlat = TRUE)

        ## Omit instances where pairwise background test point wasn't found - this occurs rarely
        wi.test.background <- wi.test.background[!is.na(wi.test.background)]
        test.background.points <- all.test.background[wi.test.background, ]

        ## Form the test dataset by binding oob presences and absences
        test.dataset = rbind.fill(presence.data[wi.test.pres,], test.background.points)

        ## Extract and predict on test.dataset
        model.predictions <- extract(x = pred.rast,
                                     y = test.dataset[,c('Longitude', 'Latitude')])
        tree.dat$model.oob.auc.pwd <- roc.area(test.dataset$Presence, model.predictions)$A
        tree.dat$n.oob.auc.pwd <- nrow(test.dataset)

        ## Store average classification score for each unique species
        spec.names <- unique(paste(dataset$Species))
        spec.dat <- matrix(nrow = 1, ncol = 3*length(spec.names))
        for(spec in spec.names){
            mask.spec <- paste(test.dataset$Species)==spec
            spec.preds <- model.predictions[mask.spec]
            m <- mean(spec.preds, na.rm = TRUE)
            s <- sd(spec.preds, na.rm = TRUE)
            spec.dat[1,which(spec==spec.names)] <- m
            spec.dat[1, length(spec.names) + which(spec==spec.names)] <- s
            spec.dat[1, 2*length(spec.names) + which(spec==spec.names)] <- sum(mask.spec)
        }
        write.table(spec.dat, file = assess.filename,
                    col.names = FALSE, row.names = FALSE,
                    append = file.exists(assess.filename))


        ###

        ## Save fit statistics to file
        write.table(tree.dat, file = tree.filename,
                    col.names = !file.exists(tree.filename), row.names = FALSE,
                    append = file.exists(tree.filename))

        ## Save model and tif prediction of this bootstrap fit
        mod.filename = paste(models.folder,'/amod_',boot.i,'.rds', sep = '')
        saveRDS(gbm.mod, file = mod.filename)

        tif.filename = paste(models.folder,'/amod_',boot.i,'.tif', sep = '')
        writeRaster(pred.rast, filename = tif.filename, overwrite = TRUE)

        print(paste('--Finished Bootset: ', boot.i, '; Elapsed Time: ', Sys.time() - starttime), quote = FALSE)
        print(paste('----- gbm.mod.ntrees', gbm.mod$n.trees, '   max.trees: ', max.trees), quote = FALSE)
        gc() ## Helps clear memory
        return(boot.i)
    } ## End mcl.fun

    unlink(tree.filename)
    unlink(assess.filename)

    print('Bootstrapping model fits', quote = FALSE)

    starttime <- Sys.time()
    out <- mclapply(1:nboots, mcl.fun, mc.cores = detectCores() - 2)
    print(paste('Bootstrapping finished; Total Time: ', Sys.time() - starttime), quote = FALSE)

    ## Aggregate and analyze the fitted models. Read in all raster predictions, average
    ## them, and save the resulting meta prediction.
    print('Averaging bootstrapped fits', quote = FALSE)

    pred.stack.names <- list.files(models.folder, pattern = "*.tif$")
    fullname <- paste(models.folder,pred.stack.names, sep = '/')
    pred.stack <- stack(fullname)
    pred.rast <- mean(pred.stack, na.rm = TRUE)
    writeRaster(pred.rast, file = paste('Figures_Fits/', prefix, '/', fold,'/Reservoir_Layer_', fold,".tif", sep = ''),
                overwrite = TRUE)

    ## The code below investigates the model fits using variable importance ranking and
    ## the learned relationships
    imp.mat <- matrix(NA, nrow = length(var.names), ncol = nboots)
    imp.dat <- data.frame(coef = NA, imp = NA)
    response.dat <- c()
    response.dat.1 <- c()
    for(boot.i in 1:nboots){
        mod.filename = paste(models.folder,'/amod_',boot.i,'.rds', sep = '')
        if(file.exists(mod.filename)){
            try({
                ## Load model, extract and save coefficient importance
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
            })## End try
        } ## Check if the file exists
    } ## Loop through each boot

    ## Extract top 10 most important predictors for a boxplot.
    row.names(imp.mat) <- var.names
    med.mat <- apply(imp.mat, 1, median, na.rm = TRUE)
    ord.names <- names(sort(med.mat, TRUE))
    imp.dat.1 <- imp.dat[-1,]
    imp.dat.2 <- na.omit(imp.dat.1)
    wi <- which(imp.dat.2$coef %in% ord.names[1:10])
    imp.dat.2 <- imp.dat.2[wi,]
    imp.dat.2$coef <- factor(imp.dat.2$coef, levels = rev(ord.names[1:10]), ordered = TRUE)

    ## Build boxplot
    p <- ggplot(imp.dat.2, aes(x=coef, y=imp)) +
        geom_boxplot()
    p + scale_x_discrete(labels=rev(pretty.labels(ord.names[1:10]))) +
        labs(x = 'Predictor', y = 'Importance') + coord_flip()##ylim = c(0,35))
    ggsave(filename = paste('Figures_Fits/', prefix, '/', fold ,'/Coef_', Abbrev.name,'.png', sep = ''), device = 'png')

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
                         geom = 'line', size = 2, fun = mean, fun.args = list(na.rm = TRUE)) +
        facet_wrap(~var.pretty, ncol = 3, nrow = 2, scales = 'free') + xlab('Predictor Value') +
        ylab('Classification Score')
    ggsave(filename = paste('Figures_Fits/', prefix, '/', fold, '/Effect_Response_', Abbrev.name,'.png', sep = ''),
           device = 'png', width = 7, height = 5, units = 'in')


    ggplot(data = response.dat[response.dat$var.pretty %in% vars.to.plot.pretty,]) +
        stat_summary_bin(aes(x = Value, y = Effect, group = var.pretty),
                         fun.data = mean_se, fun.args = list(mult = 2),
                         color = "black", fill = 'blue',
                         geom = 'ribbon', bins = 30) +
    theme_classic() +
        stat_summary_bin(aes(x = Value, y = Effect, group = var.pretty),
                         geom = 'line', size = 1, fun = mean, fun.args = list(na.rm = TRUE)) +
        facet_wrap(~var.pretty, ncol = 3, nrow = 2, scales = 'free') + xlab('Predictor Value') +
        ylab('Classification Score')
    ggsave(filename = paste('Figures_Fits/', prefix, '/', fold, '/Effect_Response_', Abbrev.name,'_var.png', sep = ''),
           device = 'png', width = 7, height = 5, units = 'in')


    ## Plot risk map averaged over all boot predictions
    heat.cols <- viridis(120, begin = 0.1, end = 1, option = 'D')
    xlims = c(-18,16)
    ylims = c(16, 16.5)
    points.sp <- presence.data[,c('Longitude','Latitude')]
    pred.rast = mask(pred.rast, masto.rangemap, updatevalue = 0)
    png(file = paste('Figures_Fits/',prefix,'/',fold,'/', Abbrev.name,'_Risk_Layer.png', sep = ''),
    width = 6, height = 4, units = 'in', res = 400)
    par(mai = 1*c(0.2,0.2,0.2,0.6))
    image.plot(pred.rast, col = heat.cols, zlim = c(0, 1),
               bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',
               xlim = xlims, ylim = ylims, ##main = Species,
               asp = 1, legend.lab = 'Occurrence score', legend.line = 2.5,
               main = '')
    mtext(text = bquote(italic(.(Species))~'Distribution ('~'D'['M']~')'), side = 3, line = -1)
    points(points.sp[,1], points.sp[,2], asp = 1, cex = 1, pch = 21, lwd = 1,
           bg = c('maroon1'), col = 'black')
    plot(foc.shp.ogr, add = TRUE, bty = 'n', asp = 1)
    plot(rgeos::gIntersection(foc.shp.ogr, masto.rangemap), add = TRUE, bty = 'n', asp = 1, lwd = 3)
    legend(x = -18,y = 4.5, legend = 'Confirmed captures',
           pch = 21, pt.bg = 'maroon1', col = 'black',
           cex = 1, bty = 'n')
    plot(rgeos::gIntersection(foc.shp.ogr, masto.rangemap), add = TRUE, bty = 'n', asp = 1, lwd = 3)
    dev.off()

    ## Process and rewrite  data with each species classification score to file
    assess.dat = read.table(file = assess.filename)
    spec.names = unique(paste(dataset$Species))
    assess.mean = assess.dat[,1:length(spec.names)] ## first set of columns store the means
    names(assess.mean) = spec.names
    assess.sd <- assess.dat[,length(spec.names) + 1:length(spec.names)] ## second set of columns store the std dev's
    names(assess.sd) = spec.names
    assess.num <- assess.dat[,2*length(spec.names) + 1:length(spec.names)] ## third set of columns store the number of each species
    dat <- data.frame(species = spec.names,
                      mean = as.vector(colMeans(assess.mean, na.rm = TRUE)),
                      sd = as.vector(colMeans(assess.sd, na.rm = TRUE)),
                      medianCount = as.vector(apply(assess.num, MARGIN = 2, FUN = median)))
    dat <- dat[order(dat$mean, decreasing = TRUE),]
    dat <- na.omit(dat)## omit those for which sd cannot be calculated (fewer than 3 entries)
    write.table(dat, file = assess.filename, row.names = FALSE)

    ## Remove fitted models
    unlink(models.folder, recursive = TRUE)

} ## End function
