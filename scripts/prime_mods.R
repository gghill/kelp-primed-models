# Setup ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (!file.exists("../outputs")) dir.create("../outputs")
library(terra)
library(geodata)
library(sdmpredictors)
library(fuzzySim)
library(sf)
library(spocc)
library(scales)
library(ggplot2)
library(tidyterra)
library(raster)

my_species <- 'Saccharina latissima'
# my_species <- 'Alaria esculenta'

countries <- world(path = "../outputs/countries")
countries_write <- st_as_sf(countries)
st_write(countries_write, "../outputs/countries/countries.shp")
download_window <- ext(c(-85, 60, 20, 85))

# Data read in ----
load('../inputs/dataBrownAlgaePruned.RData')
pruned <- finalDataBaseBPruned
names(pruned)
head(pruned$acceptedName)
pruned_species <- pruned[pruned$acceptedName==my_species,]

pruned_species <- subset(pruned_species, select = c("decimalLongitude", "decimalLatitude", "country", "occurrenceStatus", "individualCount", "year", "coordinateUncertaintyInMeters", "sourceType", "datasetID", "basisOfRecord"))
pruned_species$decimalLongitude <- as.numeric(as.character(pruned_species$decimalLongitude))
pruned_species$decimalLatitude <- as.numeric(as.character(pruned_species$decimalLatitude))
pruned_species$year <- as.numeric(as.character(pruned_species$year))
pruned_species$coordinateUncertaintyInMeters <- as.numeric(as.character(pruned_species$coordinateUncertaintyInMeters))
pruned_species <- pruned_species[xmin(download_window) <= pruned_species$decimalLongitude & pruned_species$decimalLongitude <= xmax(download_window),]
pruned_species <- pruned_species[ymin(download_window) <= pruned_species$decimalLatitude & pruned_species$decimalLatitude <= ymax(download_window),]

occurrences <- pruned_species
plot(countries, xlim = range(occurrences$decimalLongitude, na.rm = TRUE), ylim = range(occurrences$decimalLatitude, na.rm = TRUE))
points(occurrences[ , c("decimalLongitude", "decimalLatitude")], pch = 20, col = "red")

occurrences_folder <- paste0("../outputs/", my_species, "/occurrences")
if (!file.exists(occurrences_folder)) dir.create(occurrences_folder, recursive = TRUE)
write.csv(occurrences, paste0(occurrences_folder, "/occurrences_raw.csv"), row.names = FALSE)

# from here on, you don't need to download these data again - you can just import them from the .csv:
occurrences <- read.csv(paste0(occurrences_folder, "/occurrences_raw.csv"))
# remove dupes
nrow(occurrences)
occurrences <- unique(occurrences)
nrow(occurrences)
# remove absences
names(occurrences)
sort(unique(occurrences$occurrenceStatus))
nrow(occurrences)
occurrences <- subset(occurrences, occurrenceStatus != "Absent")
nrow(occurrences)
# remove coordinate uncertainty > 10km, but keep NAs
occurrences <- subset(occurrences, is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= 10000)
# add these less uncertain occurrence records with a different colour on top of the previous ones:
points(occurrences[ , c("decimalLongitude", "decimalLatitude")], pch = 20, col = "turquoise")
legend("topright", paste("n =", nrow(occurrences)), bg = "white")
title(my_species)

# remove occurrences with missing or unlikely coordinates (= 0):
occurrences <- subset(occurrences, is.finite(decimalLongitude) & is.finite(decimalLatitude))
nrow(occurrences)
occurrences <- subset(occurrences, decimalLongitude != 0 & decimalLatitude != 0)
nrow(occurrences)
points(occurrences[ , c("decimalLongitude", "decimalLatitude")], pch = 20, col = "blue")  # excluded points are not added in blue colour, so they remain visible in red

# save the cleaned data to disk as a .csv file:
write.csv(occurrences, paste0(occurrences_folder, "/occurrences_Slat_pruned_Euro.csv"), row.names = FALSE)

occurrences <- read.csv(paste0(occurrences_folder, "/occurrences_Slat_pruned_Euro.csv"))

# Layers
layers_ver2 <- list_layers(datasets = "Bio-ORACLE", version = 2)
# convert species occurrences table to a spatial object (like when you import a delimited text file into a GIS, you need to specify the names of the columns that contain the spatial coordinates (point geometry) and the geographic projection / coordinate reference system):
occ_points <- terra::vect(occurrences, geom = c("decimalLongitude", "decimalLatitude"), crs = "epsg:4326")

primed_layers <- c("BO2_tempmax_ss", "BO2_tempmean_ss",	
                   "BO2_chlomin_bdmin", "BO2_salinitymin_ss", "BO2_curvelmean_bdmin")
layers_choice <- subset(layers_ver2,layer_code %in% primed_layers, select = c("name", "layer_code", "version"))
primed_layers <- rast(load_layers(layers_choice$layer_code, datadir = "../outputs/sdmpredictors"))

# Mod region
points_buff <- aggregate(buffer(occ_points, width = 200000))
countries <- disagg(countries)  # disaggregate into separate polygons for mainlands and islands
# may want to restrict to a smaller input model region
include <- c("Albania", "Algeria", "Bosnia and Herzegovina", "Croatia", 
             "Cyprus", "Egypt", "France", "Greece", "Israel", "Italy", "Lebanon", 
             "Libya", "Malta", "Monaco", "Montenegro", "Morocco", "Slovenia", 
             "Spain", "Syria", "Tunisia", "Turkey", "Canada", 
             "Canada", "Germany", "Denmark", 
             "Finland", "Spain", "France", "United Kingdom", "United Kingdom", 
             "United Kingdom", "Isle of Man", 
             "Poland", "Ireland", "Iceland", "Jersey", "Lithuania", "Netherlands", "Norway", "Portugal", "Russia", "Svalbard and Jan Mayen", 
             "Sweden", "Latvia", "Estonia", "United States", "Greenland") # including mediterranean states

pres_buff <- aggregate(buffer(subset(countries, countries$NAME_0 %in% include), width = 250000))
crs(pres_buff) <- crs(points_buff)
mod_region <- pres_buff

# combine point buffer region and country buffer region (if needed)
mod_region <- combineGeoms(mod_region, points_buff)

mod_region <- crop(mod_region, download_window)
primed_layers_cut <- crop(primed_layers, mod_region, mask = T)
plot(download_window, border = "darkred", lwd = 2)
plot(pres_buff, col = "lightgrey", border = NA, add=T)
plot(primed_layers_cut[[1]], add=T)
points(occurrences[ , c("decimalLongitude", "decimalLatitude")], pch = 20)
plot(countries, border = "darkgrey", add = TRUE)
plot(mod_region, border = "blue", lwd = 3, add = TRUE)

primed_layers_mod <- primed_layers_cut
# consider excluding completely Med, Black sea and Baltic
# absences may overwhelming
threshold <- 15 # temperature threshold above which is modified
# model is built using these layers, so temps should be INCREASED by 1C to reflect incr tolerance
add_one <- function(x) {
  ifelse(!is.na(x) & x >= threshold, as.numeric(x + 1), x)
}
subtract_one <- function(x) {
  ifelse(!is.na(x) & x >= threshold, as.numeric(x - 1), x)
}
# primed_layers_mod[] <- unlist(lapply(primed_layers_mod[], subtract_one))

# modify both temp layers
primed_layers_mod$BO2_tempmax_ss[] <- unlist(lapply(primed_layers_mod$BO2_tempmax_ss[], add_one))
primed_layers_mod$BO2_tempmean_ss[] <- unlist(lapply(primed_layers_mod$BO2_tempmean_ss[], add_one))


plot(primed_layers_mod$BO2_tempmax_ss)
plot(primed_layers_cut$BO2_tempmax_ss) # visually inspect for difference

# closely inspect your species data vs. the size of the variables' pixels:
plot(primed_layers_mod[[1]], xlim = c(15, 20), ylim = c(69, 70))
points(occurrences, cex = 0.5)
# looking to have multiple points per pixel, not vice versa (pixels per point grid)

layers_folder <- paste0("../outputs/", my_species, "/layers")
if (!file.exists(layers_folder)) dir.create(layers_folder)

file_deets <- paste0('_',my_species,'TempCurVelChloSal_PrimeEuPruned15C_addOne') # update to reflect layers used

# prime sel is cut, but not forced (unmodified)
gridded_data <- gridRecords(rst = primed_layers_cut, pres.coords = occurrences[ , c("decimalLongitude", "decimalLatitude")])
# should return TRUE
nrow(gridded_data)==sum(!is.na(values(primed_layers_cut[[1]]))) # should be the same number
write.csv(gridded_data, paste0(occurrences_folder, paste0("/occurrences", file_deets, "sel.csv")), row.names = FALSE) # unforced, cut
writeRaster(primed_layers_cut, filename = paste0(layers_folder,"/", file_deets, "layers_cut.tif"), overwrite = TRUE)


# mod is after forcing
gridded_data <- gridRecords(rst = primed_layers_mod, pres.coords = occurrences[ , c("decimalLongitude", "decimalLatitude")])
# should return TRUE
nrow(gridded_data)==sum(!is.na(values(primed_layers_mod[[1]]))) # should be the same number
write.csv(gridded_data, paste0(occurrences_folder, paste0("/occurrences", file_deets, "mod.csv")), row.names = FALSE) # forced based on priming, cut
writeRaster(primed_layers_mod, filename = paste0(layers_folder, "/", file_deets, "layers_forced.tif"), overwrite = TRUE)

# thin occurrences prior to modeling? occurrence_thinner

# Modeling ----
occ_folder <- paste0("../outputs/", my_species, "/occurrences")
dat <- read.csv(paste0(occ_folder, paste0("/occurrences",file_deets,"sel.csv")))
dat_prime <- read.csv(paste0(occ_folder, paste0("/occurrences",file_deets,"mod.csv")))

pres_centroids <- vect(dat[dat$presence == 1, ], geom = c("x", "y"), crs = "epsg:4326")
pres_centroids_prime <- vect(dat_prime[dat_prime$presence == 1, ], geom = c("x", "y"), crs = "epsg:4326")

# format for MaxNet
pres_centroids_spatial <- as(pres_centroids, "Spatial")
pres_centroids_spatial_prime <- as(pres_centroids_prime, 'Spatial')
layers_cut_stack <- raster::stack(crop(primed_layers_cut, download_window, mask=T))
layers_cut_stack_prime <- raster::stack(crop(primed_layers_mod, download_window, mask=T))
# try to speed up by cropping here ^

# Maxent with 'maxnet': ####
library(maxnet)
vars_sel <- names(layers_cut_stack_prime)
vars_sel
?maxnet  # 'maxnet' function in 'maxnet' package, which uses 'glmnet' rather than the Java version of Maxent
# it works on a matrix or data frame, rather than raster maps + presence points

# need both models because dat includes the values for the gridded environmental data
Sys.time()
maxnet_mod <- maxnet(p = dat$presence, data = dat[, vars_sel], maxnet.formula(p = dat$presence, data = dat[, vars_sel]))  # you can add to 'maxnet.formula' e.g. classes="lq", to use only linear ('l') and quadratic ('q') features; read the help file and https://doi.org/10.1111/j.1600-0587.2013.07872.x if you use Maxent for real work!
maxnet_primed_mod <- maxnet(p = dat_prime$presence, data = dat_prime[, vars_sel], maxnet.formula(p = dat_prime$presence, data = dat_prime[, vars_sel]))
Sys.time()
# compute and map maxnet predictions:
# same model, but different enviro layers (primed vs. unprimed)
maxnet_pred <- rast(predict(layers_cut_stack, maxnet_mod, type = "cloglog"))
maxnet_prime_pred <- rast(predict(layers_cut_stack_prime, maxnet_primed_mod, type = "cloglog"))
maxnet_prime_diff <- maxnet_prime_pred - maxnet_pred # positive values should be areas with increased predicted dist

countries_diff <- subset(countries, countries$NAME_0 %in% include)
countries_diff <- crop(countries_diff, download_window)
scenario <- 'present'
plot_deets <- paste0(threshold,'C_',scenario, '_addOne_', my_species)
plot_folder <- paste0('../outputs/',my_species, '/GH_plots/')
if (!file.exists(plot_folder)) dir.create(plot_folder)
png(paste0(plot_folder, plot_deets,'.png'), width = 860, height = 560, units = 'px')
ggplot() +
  geom_spatraster(data = maxnet_prime_diff) +
  scale_fill_gradient2(
    low = "green",
    mid = 'white',
    high = "darkred",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill",
  ) +
  geom_sf(data=pres_centroids,
          size = .1,
          color = alpha('blue', .25),
  ) +
  geom_sf(data=countries_diff) +
  annotate(geom= 'label', label = paste(vars_sel, sep = '\n', collapse = '\n'), x = Inf, y = -Inf, hjust = 1, vjust = 0) +
  ggtitle(paste(my_species, threshold, 'C threshold,', 'Vars:', file_deets)) +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  )
dev.off()
# plot(maxnet_prime_diff, col = clrs, main = "Prime diff")

# Validation and Projection ----
# generate AUC
library(modEvA)
par(mar = c(3, 2, 2, 1), mfrow = c(1, 2))
dat_points <- dat[dat$presence==1, c("x","y")]
dat_primed_points <- dat_prime[dat_prime$presence==1, c("x","y")]
with(dat, AUC(obs = dat[dat$presence==1, c("x","y")], pred = maxnet_pred, main = "Normal"))
with(dat_prime, AUC(obs = dat_prime[dat_prime$presence==1, c("x","y")], pred = maxnet_prime_pred, main = "Primed"))

# attempt to delineate into blocks and test against training dataset, not functional
library(blockCV)
dat_points <- vect(dat, geom = c("x", "y"), crs = crs(layers_cut_stack[[1]]))
dat_points <- (as(dat_points, "Spatial"))
blocks <- cv_spatial(dat_points, r = raster(layers_cut_stack[[1]]),  k = 5, selection = "random", seed = 123)  # you can also use an optional additional argument, species = "presence"; see ?spatialBlock

blocks$folds
blocks$foldID

dat$foldID <- blocks$foldID
head(dat)

dat_points$foldID <- blocks$foldID
par(mfrow = c(1, 1))
plot(dat_points, "foldID", col = hcl.colors(5, "TealRose"), main = my_species)
plot(download_window)
plot(mod_region, add = TRUE)
plot(countries, add = TRUE)
plot(blocks$blocks, add = TRUE)

set_ext <- function(x) {
  final <- rast()
  for (r in x) {
    final <- c(final, crop(rast(r), mod_region), mask=T)
  }
  final
}
future_layers <- list_layers_future()
future_biooracle <- subset(future_layers, dataset_code == "Bio-ORACLE")
scenarios <- c('2050_RCP8-8.5', '2050_RCP2-4.5', '2100_RCP8-8.5', '2100_RCP2-4.5')

# write a function that takes in a scenario and BO2 variables and returns a plot
# ideally vars_sel is the input
# set up for lapply

project_models <- function(scenario, vars) {
    year <- sub("_.*", "", scenario)
    rcp <- sub(".*_", "", scenario)
    if (rcp=="RCP2-4.5") {
      rcp <- "RCP45"
    } else if (rcp=="RCP8-8.5") {
      rcp <- "RCP85"
    }
    vars_scenario <- gsub("BO2", paste0("BO22_", rcp,"_",year), vars)
    vars_scenario <- gsub('bdmin', 'ss', vars_scenario)
    vars_scenario <- subset(vars_scenario, vars_scenario %in% future_biooracle$layer_code)
    
    layers_scenario <- load_layers(layercodes = vars_scenario, datadir = "../outputs/sdmpredictors/future", rasterstack = F)
    layers_scenario <- set_ext(layers_scenario)
    layers_scenario <- crop(layers_scenario, mod_region, mask = T)
    
    names(layers_scenario) <- gsub(pattern = paste0("BO22_", rcp,"_",year), replacement = "BO2", x = names(layers_scenario))
    names(layers_scenario) <- gsub(pattern = "chlomin_ss", replacement = "chlomin_bdmin", x = names(layers_scenario))
    names(layers_scenario) <- gsub(pattern = "curvelmean_ss", replacement = "curvelmean_bdmin", x = names(layers_scenario))
    
    vars_sel_proj <- subset(names(layers_scenario), names(layers_scenario) %in% vars)
    
    maxnet_pred_proj <- rast(predict(raster::stack(layers_scenario), maxnet_mod, type = "cloglog"))
    maxnet_prime_pred_proj <- rast(predict(raster::stack(layers_scenario), maxnet_primed_mod, type = "cloglog"))
    
    maxnet_prime_diff_proj <- maxnet_prime_pred_proj - maxnet_pred_proj # positive values should be areas with increased predicted dist
    
    p <- ggplot() +
      geom_spatraster(data = maxnet_prime_diff_proj) +
      scale_fill_gradient2(
        low = "green",
        mid = 'white',
        high = "darkred",
        midpoint = 0,
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "fill",
      ) +
      geom_sf(data=pres_centroids,
              size = .1,
              color = alpha('blue', .25),
      ) +
      geom_sf(data=countries_diff) +
      ggtitle(paste(my_species, 'projected to', year, rcp, '- diff')) +
      annotate(geom= 'label', label = paste(vars_sel_proj, sep = '\n', collapse = '\n'), x = Inf, y = -Inf, hjust = 1, vjust = 0) +
      theme(axis.title.x=element_blank(), 
            axis.title.y=element_blank(),
            axis.text.x=element_blank(), #remove x axis labels
            axis.ticks.x=element_blank(), #remove x axis ticks
            axis.text.y=element_blank(),  #remove y axis labels
            axis.ticks.y=element_blank()  #remove y axis ticks
      )
}

my_plots <- lapply(scenarios, project_models, vars = vars_sel)
require('dggridR')
require(cowplot)
grid_vis <- plot_grid(plotlist = my_plots)
title <- ggdraw() + 
    draw_label(
      paste(my_species, "at 2050 and 100 for multiple climate change scenarios, model diff"),
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) # +
    # theme(
    #   # add margin on the left of the drawing canvas,
    #   # so title is aligned with left edge of first plot
    #   plot.margin = margin(15, 0, 10, 7)
    # )
save_grid <- plot_grid(title, grid_vis, ncol = 1, axis='b',rel_heights = c(0.03,1))

with(dat, AUC(obs = occurrences, pred = maxnet_pred, main = "Maxnet"))

# replaced with project_models function above

# ## 2050 RCP 8-8.5 ----
# scenario <- '2050_RCP8-8.5'
# vars_2050_85 <- gsub("BO2", "BO22_RCP85_2050", vars_sel)
# vars_2050_85 <- gsub('bdmin', 'ss', vars_2050_85)
# vars_2050_85 <- subset(vars_2050_85, vars_2050_85 %in% future_biooracle$layer_code)
# 
# # need to crop to layer with smallest extent
# 
# layers_2050_85 <- load_layers(layercodes = vars_2050_85, datadir = "../outputs/sdmpredictors/future", rasterstack = F)
# layers_2050_85 <- set_ext(layers_2050_85)
# # layers_2050_85 <- raster::addLayer(layers_2050_85,crop(load_layers(vars_2050_85[2]), layers_2050_85[[1]]))
# # layers_2050_85 <- raster::stack(layers_2050_85, crop(load_layers(vars_2050_85[3]), layers_2050_85[[1]]))
# 
# layers_2050_85 <- crop(layers_2050_85, mod_region, mask = T)
# plot(layers_2050_85)
# 
# # I think this shouldn't be done
# # primed_layers_2050_85$BO22_RCP85_2050_tempmax_ss[] <- unlist(lapply(primed_layers_2050_85$BO22_RCP85_2050_tempmax_ss[], subtract_one))
# # plot(primed_layers_2050_85$BO22_RCP85_2050_tempmax_ss)
# # plot(layers_2050_85$BO22_RCP85_2050_tempmax_ss)
# 
# names(layers_2050_85) <- gsub(pattern = "BO22_RCP85_2050", replacement = "BO2", x = names(layers_2050_85))
# names(layers_2050_85) <- gsub(pattern = "chlomin_ss", replacement = "chlomin_bdmin", x = names(layers_2050_85))
# names(layers_2050_85) <- gsub(pattern = "curvelmean_ss", replacement = "curvelmean_bdmin", x = names(layers_2050_85))
# 
# vars_sel_proj <- subset(names(layers_2050_85), names(layers_2050_85) %in% vars_sel)
# # make the projected models
# # would like to reuse maxnet models from previously, but not all layers are available in the future
# # so we are trimming down and need new corresponding models
# # NOW LAYERS ARE CONSISTENT SO REUSE
# 
# # compute and map maxnet predictions:
# # maxnet_mod_proj_2050_85 <- maxnet(p = dat$presence, data = dat[, vars_sel_proj], maxnet.formula(p = dat$presence, data = dat[, vars_sel_proj]))
# # 
# # maxnet_mod_proj <- maxnet_mod_proj_2050_85
# maxnet_pred_proj <- rast(predict(raster::stack(layers_2050_85), maxnet_mod, type = "cloglog"))
# 
# maxnet_prime_pred_proj_2050_85 <- rast(predict(raster::stack(layers_2050_85), maxnet_primed_mod, type = "cloglog"))
# maxnet_prime_diff_proj_2050_85 <- maxnet_prime_pred_proj_2050_85 - maxnet_pred_proj # positive values should be areas with increased predicted dist
# 
# plot_deets <- paste0(scenario, file_deets)
# png(paste0('../outputs/',my_species, '/GH_plots/', plot_deets,'.png'), width = 860, height = 560, units = 'px')
# ggplot() +
#   geom_spatraster(data = maxnet_prime_diff_proj_2050_85) +
#   scale_fill_gradient2(
#     low = "green",
#     mid = 'white',
#     high = "darkred",
#     midpoint = 0,
#     space = "Lab",
#     na.value = "grey50",
#     guide = "colourbar",
#     aesthetics = "fill",
#   ) +
#   geom_sf(data=pres_centroids,
#           size = .1,
#           color = alpha('blue', .25),
#   ) +
#   geom_sf(data=countries_diff) +
#   ggtitle(paste(my_species, 'projected to 2050 RCP85 - diff')) +
#   annotate(geom= 'label', label = paste(vars_sel_proj, sep = '\n', collapse = '\n'), x = Inf, y = -Inf, hjust = 1, vjust = 0) +
#   theme(axis.title.x=element_blank(), 
#         axis.title.y=element_blank(),
#         axis.text.x=element_blank(), #remove x axis labels
#         axis.ticks.x=element_blank(), #remove x axis ticks
#         axis.text.y=element_blank(),  #remove y axis labels
#         axis.ticks.y=element_blank()  #remove y axis ticks
#   )
# dev.off()
# 
# 
# ## 2050 RCP 2-4.5 ----
# scenario <- '2050_RCP2-4.5'
# vars_2050_45 <- gsub("BO2", "BO22_RCP45_2050", vars_sel)
# vars_2050_45 <- gsub('bdmin', 'ss', vars_2050_45)
# vars_2050_45 <- subset(vars_2050_45, vars_2050_45 %in% future_biooracle$layer_code)
# 
# 
# layers_2050_45 <- load_layers(layercodes = vars_2050_45[1], datadir = "../outputs/sdmpredictors/future")
# layers_2050_45 <- raster::stack(layers_2050_45, crop(load_layers(vars_2050_45[2]), layers_2050_45[[1]]))
# layers_2050_45 <- raster::stack(layers_2050_45, crop(load_layers(vars_2050_45[3]), layers_2050_45[[1]]))
# plot(layers_2050_45[[2]])
# 
# layers_2050_45 <- raster::stack(layers_2050_45)
# layers_2050_45 <- rast(layers_2050_45)
# layers_2050_45 <- crop(layers_2050_45, mod_region, mask = T)
# primed_layers_2050_45 <- layers_2050_45
# 
# primed_layers_2050_45$BO22_RCP45_2050_tempmax_ss[] <- unlist(lapply(primed_layers_2050_45$BO22_RCP45_2050_tempmax_ss[], subtract_one))
# plot(primed_layers_2050_45$BO22_RCP45_2050_tempmax_ss)
# plot(layers_2050_45$BO22_RCP45_2050_tempmax_ss)
# 
# names(layers_2050_45) <- gsub(pattern = "BO22_RCP45_2050", replacement = "BO2", x = names(layers_2050_45))
# names(layers_2050_45) <- gsub(pattern = "chlomin_ss", replacement = "chlomin_bdmin", x = names(layers_2050_45))
# names(primed_layers_2050_45) <- names(layers_2050_45)
# vars_sel_proj <- subset(names(layers_2050_45), names(layers_2050_45) %in% vars_sel)
# # make the projected models
# # reusing maxnet_mod and maxnet_primed_mod from before
# 
# # compute and map maxnet predictions:
# maxnet_mod_proj_2050_45 <- maxnet(p = dat$presence, data = dat[, vars_sel_proj], maxnet.formula(p = dat$presence, data = dat[, vars_sel_proj]))
# maxnet_mod_proj <- maxnet_mod_proj_2050_45
# maxnet_pred_proj_2050_45 <- rast(predict(raster::stack(layers_2050_45), maxnet_mod_proj, type = "cloglog"))
# 
# maxnet_prime_pred_proj_2050_45 <- rast(predict(raster::stack(primed_layers_2050_45), maxnet_mod_proj, type = "cloglog"))
# maxnet_prime_diff_proj_2050_45 <- maxnet_prime_pred_proj_2050_45 - maxnet_pred_proj_2050_45 # positive values should be areas with increased predicted dist
# 
# plot_deets <- paste0(threshold,'C_',scenario, '_', my_species)
# png(paste0('../outputs/',my_species, '/GH_plots/', plot_deets,'.png'), width = 860, height = 560, units = 'px')
# ggplot() +
#   geom_spatraster(data = maxnet_prime_diff_proj_2050_45) +
#   scale_fill_gradient2(
#     low = "green",
#     mid = 'white',
#     high = "darkred",
#     midpoint = 0,
#     space = "Lab",
#     na.value = "grey50",
#     guide = "colourbar",
#     aesthetics = "fill",
#   ) +
#   geom_sf(data=pres_centroids,
#           size = .1,
#           color = alpha('blue', .25),
#   ) +
#   geom_sf(data=countries_diff) +
#   ggtitle(paste(my_species, 'projected to 2050 RCP45 - diff')) +
#   annotate(geom= 'label', label = paste(vars_sel_proj, sep = '\n', collapse = '\n'), x = Inf, y = -Inf, hjust = 1, vjust = 0) +
#   theme(axis.title.x=element_blank(), 
#         axis.title.y=element_blank(),
#         axis.text.x=element_blank(), #remove x axis labels
#         axis.ticks.x=element_blank(), #remove x axis ticks
#         axis.text.y=element_blank(),  #remove y axis labels
#         axis.ticks.y=element_blank()  #remove y axis ticks
#   )
# dev.off()
# # plot(maxnet_pred_proj_2050_45, main= 'projected to 2050 RCP45')
# # plot(maxnet_prime_pred_proj_2050_45, main= 'projected to 2050 RCP45 PRIMED')
# 
# ## 2100 RCP 2-4.5 ----
# scenario <- '2100_RCP2-4.5'
# vars_2100_45 <- gsub("BO2", "BO22_RCP45_2100", vars_sel)
# vars_2100_45 <- gsub('bdmin', 'ss', vars_2100_45)
# vars_2100_45 <- subset(vars_2100_45, vars_2100_45 %in% future_biooracle$layer_code)
# 
# 
# layers_2100_45 <- load_layers(layercodes = vars_2100_45[1], datadir = "../outputs/sdmpredictors/future")
# layers_2100_45 <- raster::stack(layers_2100_45, crop(load_layers(vars_2100_45[2]), layers_2100_45[[1]]))
# layers_2100_45 <- raster::stack(layers_2100_45, crop(load_layers(vars_2100_45[3]), layers_2100_45[[1]]))
# plot(layers_2100_45[[2]])
# 
# layers_2100_45 <- raster::stack(layers_2100_45)
# layers_2100_45 <- rast(layers_2100_45)
# layers_2100_45 <- crop(layers_2100_45, mod_region, mask = T)
# primed_layers_2100_45 <- layers_2100_45
# 
# primed_layers_2100_45$BO22_RCP45_2100_tempmax_ss[] <- unlist(lapply(primed_layers_2100_45$BO22_RCP45_2100_tempmax_ss[], subtract_one))
# plot(primed_layers_2100_45$BO22_RCP45_2100_tempmax_ss)
# plot(layers_2100_45$BO22_RCP45_2100_tempmax_ss)
# 
# names(layers_2100_45) <- gsub(pattern = "BO22_RCP45_2100", replacement = "BO2", x = names(layers_2100_45))
# names(layers_2100_45) <- gsub(pattern = "chlomin_ss", replacement = "chlomin_bdmin", x = names(layers_2100_45))
# names(primed_layers_2100_45) <- names(layers_2100_45)
# vars_sel_proj <- subset(names(layers_2100_45), names(layers_2100_45) %in% vars_sel)
# # make the projected models
# # reusing maxnet_mod and maxnet_primed_mod from before
# 
# # compute and map maxnet predictions:
# maxnet_mod_proj_2100_45 <- maxnet(p = dat$presence, data = dat[, vars_sel_proj], maxnet.formula(p = dat$presence, data = dat[, vars_sel_proj]))
# maxnet_mod_proj <- maxnet_mod_proj_2100_45
# maxnet_pred_proj_2100_45 <- rast(predict(raster::stack(layers_2100_45), maxnet_mod_proj, type = "cloglog"))
# 
# maxnet_prime_pred_proj_2100_45 <- rast(predict(raster::stack(primed_layers_2100_45), maxnet_mod_proj, type = "cloglog"))
# maxnet_prime_diff_proj_2100_45 <- maxnet_prime_pred_proj_2100_45 - maxnet_pred_proj_2100_45 # positive values should be areas with increased predicted dist
# plot_deets <- paste0(threshold,'C_',scenario, '_', my_species)
# png(paste0('../outputs/',my_species, '/GH_plots/', plot_deets,'.png'), width = 860, height = 560, units = 'px')
# ggplot() +
#   geom_spatraster(data = maxnet_prime_diff_proj_2100_45) +
#   scale_fill_gradient2(
#     low = "green",
#     mid = 'white',
#     high = "darkred",
#     midpoint = 0,
#     space = "Lab",
#     na.value = "grey50",
#     guide = "colourbar",
#     aesthetics = "fill",
#   ) +
#   geom_sf(data=pres_centroids,
#           size = .1,
#           color = alpha('blue', .25),
#   ) +
#   geom_sf(data=countries_diff) +
#   ggtitle(paste(my_species, 'projected to 2100 RCP45 - diff')) +
#   annotate(geom= 'label', label = paste(vars_sel_proj, sep = '\n', collapse = '\n'), x = Inf, y = -Inf, hjust = 1, vjust = 0) +
#   theme(axis.title.x=element_blank(), 
#         axis.title.y=element_blank(),
#         axis.text.x=element_blank(), #remove x axis labels
#         axis.ticks.x=element_blank(), #remove x axis ticks
#         axis.text.y=element_blank(),  #remove y axis labels
#         axis.ticks.y=element_blank()  #remove y axis ticks
#   )
# dev.off()
# # plot(maxnet_pred_proj_2100_45, main= 'projected to 2100 RCP45')
# # plot(maxnet_prime_pred_proj_2100_45, main= 'projected to 2100 RCP45 PRIMED')
# 
# ## 2100 RCP8-8.5 ----
# scenario <- '2100_RCP8-8.5'
# vars_2100_85 <- gsub("BO2", "BO22_RCP85_2100", vars_sel)
# vars_2100_85 <- gsub('bdmin', 'ss', vars_2100_85)
# vars_2100_85 <- subset(vars_2100_85, vars_2100_85 %in% future_biooracle$layer_code)
# layers_2100_85 <- load_layers(layercodes = vars_2100_85[1], datadir = "../outputs/sdmpredictors/future")
# layers_2100_85 <- raster::stack(layers_2100_85, crop(load_layers(vars_2100_85[2]), layers_2100_85[[1]]))
# layers_2100_85 <- raster::stack(layers_2100_85, crop(load_layers(vars_2100_85[3]), layers_2100_85[[1]]))
# plot(layers_2100_85[[2]])
# 
# layers_2100_85 <- raster::stack(layers_2100_85)
# layers_2100_85 <- rast(layers_2100_85)
# layers_2100_85 <- crop(layers_2100_85, mod_region, mask = T)
# primed_layers_2100_85 <- layers_2100_85
# 
# primed_layers_2100_85$BO22_RCP85_2100_tempmax_ss[] <- unlist(lapply(primed_layers_2100_85$BO22_RCP85_2100_tempmax_ss[], subtract_one))
# plot(primed_layers_2100_85$BO22_RCP85_2100_tempmax_ss)
# plot(layers_2100_85$BO22_RCP85_2100_tempmax_ss)
# 
# names(layers_2100_85) <- gsub(pattern = "BO22_RCP85_2100", replacement = "BO2", x = names(layers_2100_85))
# names(layers_2100_85) <- gsub(pattern = "chlomin_ss", replacement = "chlomin_bdmin", x = names(layers_2100_85))
# names(primed_layers_2100_85) <- names(layers_2100_85)
# vars_sel_proj <- subset(names(layers_2100_85), names(layers_2100_85) %in% vars_sel)
# 
# # make the projected models
# # reusing maxnet_mod and maxnet_primed_mod from before
# 
# # compute and map maxnet predictions:
# maxnet_mod_proj_2100_85 <- maxnet(p = dat$presence, data = dat[, vars_sel_proj], maxnet.formula(p = dat$presence, data = dat[, vars_sel_proj]))
# maxnet_mod_proj <- maxnet_mod_proj_2100_85
# maxnet_pred_proj_2100_85 <- rast(predict(raster::stack(layers_2100_85), maxnet_mod_proj, type = "cloglog"))
# 
# maxnet_prime_pred_proj_2100_85 <- rast(predict(raster::stack(primed_layers_2100_85), maxnet_mod_proj, type = "cloglog"))
# maxnet_prime_diff_proj_2100_85 <- maxnet_prime_pred_proj_2100_85 - maxnet_pred_proj_2100_85 # positive values should be areas with increased predicted dist
# 
# plot_deets <- paste0(threshold,'C_',scenario, '_', my_species)
# png(paste0('../outputs/',my_species, '/GH_plots/', plot_deets,'.png'), width = 860, height = 560, units = 'px')
# ggplot() +
#   geom_spatraster(data = maxnet_prime_diff_proj_2100_85) +
#   scale_fill_gradient2(
#     low = "green",
#     mid = 'white',
#     high = "darkred",
#     midpoint = 0,
#     space = "Lab",
#     na.value = "grey50",
#     guide = "colourbar",
#     aesthetics = "fill",
#   ) +
#   geom_sf(data=pres_centroids,
#           size = .1,
#           color = alpha('blue', .25),
#   ) +
#   geom_sf(data=countries_diff) +
#   ggtitle(paste(my_species, 'projected to 2100 RCP85 - diff')) +
#   annotate(geom= 'label', label = paste(vars_sel_proj, sep = '\n', collapse = '\n'), x = Inf, y = -Inf, hjust = 1, vjust = 0) +
#   theme(axis.title.x=element_blank(), 
#         axis.title.y=element_blank(),
#         axis.text.x=element_blank(), #remove x axis labels
#         axis.ticks.x=element_blank(), #remove x axis ticks
#         axis.text.y=element_blank(),  #remove y axis labels
#         axis.ticks.y=element_blank()  #remove y axis ticks
#   )
# dev.off()
