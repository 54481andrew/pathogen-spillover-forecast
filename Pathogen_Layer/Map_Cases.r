## This script produces a map containing points that indicate the location of rodent
## LASV surveys, and points that indicate the location of human LASV sero-surveys.

## Load processed rodent dataset
classi.dat.rod <- read.csv(file = paste0('Figures_Fits/', prefix, '/',fold,'/',
                                         'Prepped_Pathogen_PresAbs_Data.csv'))
human.test.dat <- read.csv(file = paste0('Figures_Fits/', prefix, '/',fold,'/',
                                         'Prepped_Human_Seroprevalence_Data.csv'))

## Form graphing dataset
classi.dat <- rbind.fill(human.test.dat, classi.dat.rod)
classi.dat <- classi.dat[,c('Longitude', 'Latitude', 'Species', 'ArenaStat')]
classi.dat[,'Status'] <- as.factor(c('Lassa Absent', 'Lassa Present')[classi.dat$ArenaStat + 1])

which.humans <- classi.dat$Species=='Homo sapiens'

## Read in shape files of focal countries into spatial data frame
storage.fold <- '../Storage'
foc.shp.ogr <- readOGR(dsn = paste(storage.fold, '/Shapefiles/West_Africa', sep = ''),
                       layer = 'foc', verbose = FALSE)

## Read in shapefile of Mastomys natalensis distribution, crop to West Africa
masto.rangemap <-readOGR(dsn = paste(storage.fold, '/Shapefiles/Masto_Range', sep = ''),
                         layer = 'data_0', verbose = FALSE)
masto.rangemap <- crop(masto.rangemap, foc.shp.ogr)

## Organize dataset for plotting
ext = extent(foc.shp.ogr)
xlims <- c(ext@xmin,ext@xmax)
ylims <- c(ext@ymin,ext@ymax)
classi.dat$Species = factor(classi.dat$Sp)
classi.dat$BiStatus <- NA
classi.dat$BiStatus[classi.dat$Status=='Lassa Present' & classi.dat$Sp=='Homo sapiens'] <- 'Human, Arenavirus +'
classi.dat$BiStatus[classi.dat$Status=='Lassa Absent' & classi.dat$Sp=='Homo sapiens'] <- 'Human, Arenavirus -'
classi.dat$BiStatus[classi.dat$Status=='Lassa Present' & classi.dat$Sp=='Mastomys natalensis'] <- 'Rodent, Lassa +'
classi.dat$BiStatus[classi.dat$Status=='Lassa Absent' & classi.dat$Sp=='Mastomys natalensis'] <- 'Rodent, Lassa -'
classi.dat$BiStatus <- factor(classi.dat$BiStatus)

classi.dat <- classi.dat[,c('Longitude', 'Latitude', 'BiStatus')]

classi.dat$BiStatus <- factor(classi.dat$BiStatus,
                              levels = c('Rodent, Lassa +', 'Rodent, Lassa -',
                                         'Human, Arenavirus +', 'Human, Arenavirus -'),
                              ordered = TRUE)

## Plot all data, human and rodent, Lassa + and -
mymap <- fortify(foc.shp.ogr)
g1 <- ggplot() +
      geom_blank(data = mymap, aes(x=long, y=lat)) +
      geom_map(data = mymap, map = mymap,
               aes(group = group, map_id = id),
               fill = 'gray', color = "black", size = 0.3) +
    geom_path(data = masto.rangemap, 
              aes(x = long, y = lat, group = group),
              color = 'blue', size = 1, show.legend = FALSE, linetype = 'dashed') + 
    geom_point(data = classi.dat,  size = ifelse(which.humans, 2, 1.5),
                aes(x = Longitude, y = Latitude,
                    fill = BiStatus, shape = BiStatus),
               stroke = 0.25) +
scale_fill_manual(name = "Population & Status",
                      labels = c('Rodent, Lassa +', 'Rodent, Lassa -', 'Human, Arenavirus +', 'Human, Arenavirus -'),
                      values=c("#F8766D", "deepskyblue1", "#F8766D", "deepskyblue1")) +
  scale_shape_manual(name = "Population & Status",
                      labels = c('Rodent, Lassa +', 'Rodent, Lassa -', 'Human, Arenavirus +', 'Human, Arenavirus -'),
                     values = c(21,21,24,24)) +
    scale_x_continuous(limits = xlims, expand = c(0, 0)) +
    scale_y_continuous(limits = ylims,  expand = c(0, 0)) +
    coord_fixed() +
    theme(panel.background = element_rect(fill = 'brown'),
          legend.position = 'right') +
    xlab('Longitude') + ylab('Latitude')
g1

## These countries are highlighted gray in the inset map
foc.countries <- c('Mali', "Guinea", "Ivory Coast", "Sierra Leone",
                   "Nigeria", 'Liberia',
                   'Ghana', 'Togo', 'Benin',
                   'Burkina Faso', 'Senegal', 'Gambia',
                   'Niger', 'Mauritania', "Guinea-Bissau")
world <- map_data("world")
g2 <- ggplotGrob(
    ggplot() +
    geom_polygon(data = world,
                 aes(x = long, y = lat,
                     fill = c('blue','green')[((region %in% foc.countries) + 1)],
                     group = group),
                 size = 0.3, show.legend = FALSE) +
    scale_x_continuous(limits = c(-25, 50), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-45, 45),  expand = c(0, 0)) +
    coord_quickmap() + coord_fixed() +
theme_map() +
theme(panel.background = element_rect(fill = 'white', color = 'black'),
      panel.border = element_rect(colour = "black", fill=NA, size=5),
      ) +
    scale_fill_manual(values = c("brown", "grey")) + borders()
)

## Put main map and inset together
g3 <- g1 +
      annotation_custom(grob = g2, xmin = 5, xmax = 15,
                        ymin = 17, ymax = 27)
g3
ggsave(filename = paste0('Figures_Fits/', prefix, '/', fold, '/', 'Data_Map.png', sep=''),
       device = 'png', units = 'in', height = 4, width = 6)


