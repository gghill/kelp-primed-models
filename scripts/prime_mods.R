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

my_species <- 'Saccharina latissima'

countries <- world(path = "../outputs/countries")
countries_write <- st_as_sf(countries)
st_write(countries_write, "../outputs/countries/countries.shp")

# Data read in ----
load('dataBrownAlgaePruned.RData')
pruned <- finalDataBaseBPruned
names(pruned)
head(pruned$acceptedName)
pruned_species <- pruned[pruned$acceptedName==my_species,]

pruned_species <- subset(pruned_species, select = c("decimalLongitude", "decimalLatitude", "country", "occurrenceStatus", "individualCount", "year", "coordinateUncertaintyInMeters", "sourceType", "datasetID", "basisOfRecord"))
pruned_species$decimalLongitude <- as.numeric(as.character(pruned_species$decimalLongitude))
pruned_species$decimalLatitude <- as.numeric(as.character(pruned_species$decimalLatitude))
pruned_species$year <- as.numeric(as.character(pruned_species$year))
pruned_species$coordinateUncertaintyInMeters <- as.numeric(as.character(pruned_species$coordinateUncertaintyInMeters))

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
write.csv(occurrences, paste0(occurrences_folder, "/occurrences_Slat_pruned.csv"), row.names = FALSE)

occurrences <- read.csv(paste0(occurrences_folder, "/occurrences_Slat_pruned.csv"))

# Layers
layers_ver2 <- list_layers(datasets = "Bio-ORACLE", version = 2)
# convert species occurrences table to a spatial object (like when you import a delimited text file into a GIS, you need to specify the names of the columns that contain the spatial coordinates (point geometry) and the geographic projection / coordinate reference system):
occ_points <- terra::vect(occurrences, geom = c("decimalLongitude", "decimalLatitude"), crs = "epsg:4326")

primed_layers <- c("BO2_tempmax_ss", 	
                   "BO2_chlomin_bdmin", "BO2_salinitymin_ss", "BO2_phosphatemin_bdmin", "BO2_nitratemin_bdmin")
layers_choice <- subset(layers_ver2,layer_code %in% primed_layers, select = c("name", "layer_code", "version"))
primed_layers <- rast(load_layers(layers_choice$layer_code, datadir = "../outputs/sdmpredictors"))

# Mod region
points_buff <- aggregate(buffer(occ_points, width = 200000))
countries <- disagg(countries)  # disaggregate into separate polygons for mainlands and islands
include <- c("Albania", "Algeria", "Bosnia and Herzegovina", "Croatia", 
             "Cyprus", "Egypt", "France", "Greece", "Israel", "Italy", "Lebanon", 
             "Libya", "Malta", "Monaco", "Montenegro", "Morocco", "Slovenia", 
             "Spain", "Syria", "Tunisia", "Turkey", "Canada", 
             "Canada", "Germany", "Denmark", 
             "Denmark", "Spain", "France", "United Kingdom", "United Kingdom", 
             "United Kingdom", "Isle of Man", 
             "Ireland", "Ireland", "Iceland", "Jersey", "Netherlands", "Netherlands", "Norway", "Portugal", "Russia", "Svalbard and Jan Mayen", 
             "Sweden", "Sweden", "United States", "United States", "Greenland") # including mediterranean states

pres_buff <- aggregate(buffer(subset(countries, countries$NAME_0 %in% include), width = 250000))
download_window <- ext(occ_points)
mod_region <- pres_buff

# combine point buffer region and country buffer region (if needed)
mod_region <- combineGeoms(mod_region, points_buff)

download_window <- ext(mod_region)
primed_layers_cut <- crop(primed_layers, mod_region, mask = T)
plot(download_window, border = "darkred", lwd = 2)
plot(pres_buff, col = "lightgrey", border = NA, add=T)
plot(primed_layers_cut[[1]], add=T)
points(occurrences[ , c("decimalLongitude", "decimalLatitude")], pch = 20)
plot(countries, border = "darkgrey", add = TRUE)
plot(mod_region, border = "blue", lwd = 3, add = TRUE)

primed_layers_mod <- primed_layers_cut
threshold <- 10
subtract_one <- function(x) {
  ifelse(!is.na(x) & x >= threshold, as.numeric(x - 1), x)
}
# primed_layers_mod[] <- unlist(lapply(primed_layers_mod[], subtract_one))

primed_layers_mod$BO2_tempmax_ss[] <- unlist(lapply(primed_layers_mod$BO2_tempmax_ss[], subtract_one))

plot(primed_layers_mod)

# closely inspect your species data vs. the size of the variables' pixels:
plot(primed_layers_mod[[1]], xlim = c(15, 20), ylim = c(69, 70))
points(occurrences, cex = 0.5)
# looking to have multiple points per pixel, not vice versa (pixels per point grid)

layers_folder <- paste0("../outputs/", my_species, "/layers")
if (!file.exists(layers_folder)) dir.create(layers_folder)

file_deets <- '_SlatTempNitPhos_PrimePrune_'

# prime sel is cut, but not forced (unmodified)
gridded_data <- gridRecords(rst = primed_layers_cut, pres.coords = occurrences[ , c("decimalLongitude", "decimalLatitude")])
nrow(gridded_data)==sum(!is.na(values(primed_layers_cut[[1]]))) # should be the same number
write.csv(gridded_data, paste0(occurrences_folder, paste0("/occurrences", file_deets, "sel.csv")), row.names = FALSE) # unforced, cut
writeRaster(primed_layers_cut, filename = paste0(layers_folder,"/", file_deets, "layers_cut.tif"), overwrite = TRUE)


# mod is after forcing
gridded_data <- gridRecords(rst = primed_layers_mod, pres.coords = occurrences[ , c("decimalLongitude", "decimalLatitude")])
nrow(gridded_data)==sum(!is.na(values(primed_layers_mod[[1]]))) # should be the same number
write.csv(gridded_data, paste0(occurrences_folder, paste0("/occurrences", file_deets, "mod.csv")), row.names = FALSE) # forced based on priming, cut
writeRaster(primed_layers_mod, filename = paste0(layers_folder, "/", file_deets, "layers_forced.tif"), overwrite = TRUE)

# Modeling ----
occ_folder <- paste0("../outputs/", my_species, "/occurrences")
# dat <- read.csv(paste0(occ_folder, "/occurrences_mod.csv"))
dat <- read.csv(paste0(occ_folder, paste0("/occurrences",file_deets,"sel.csv")))
dat_prime <- read.csv(paste0(occ_folder, "/occurrences_SlatTempNitPhos_PrimePrune_mod.csv"))

# Validation and Projection ----

future_layers <- list_layers_future()
future_biooracle <- subset(future_layers, dataset_code == "Bio-ORACLE")
vars_2050_85 <- gsub("BO2", "BO22_RCP85_2050", vars_sel)
vars_2050_85 <- gsub('bdmin', 'ss', vars_2050_85)
vars_2050_85 <- subset(vars_2050_85, vars_2050_85 %in% future_biooracle$layer_code)
vars_2100_85 <- gsub("BO2", "BO22_RCP85_2100", vars_sel)
vars_2100_85 <- gsub('bdmin', 'ss', vars_2100_85)
vars_2100_85 <- subset(vars_2100_85, vars_2100_85 %in% future_biooracle$layer_code)

layers_2050_85 <- load_layers(layercodes = vars_2050_85, datadir = "../outputs/sdmpredictors/future")
plot(layers_2050_85[[2]])
layers_2050_85[[2]] <- crop(layers_2050_85[[2]], layers_2050_85[[1]])
layers_2050_85[[3]] <- crop(layers_2050_85[[3]], layers_2050_85[[1]])

layers_2050_85 <- raster::stack(layers_2050_85)
layers_2050_85 <- rast(layers_2050_85)
layers_2050_85 <- crop(layers_2050_85, mod_region, mask = T)
primed_layers_2050_85 <- layers_2050_85

primed_layers_2050_85$BO22_RCP85_2050_tempmax_ss[] <- unlist(lapply(primed_layers_2050_85$BO22_RCP85_2050_tempmax_ss[], subtract_one))
plot(primed_layers_2050_85$BO22_RCP85_2050_tempmax_ss)
plot(layers_2050_85$BO22_RCP85_2050_tempmax_ss)

names(layers_2050_85) <- gsub(pattern = "BO22_RCP85_2050", replacement = "BO2", x = names(layers_2050_85))
names(layers_2050_85) <- gsub(pattern = "chlomin_ss", replacement = "chlomin_bdmin", x = names(layers_2050_85))
names(primed_layers_2050_85) <- names(layers_2050_85)
vars_sel_proj <- subset(names(layers_2050_85), names(layers_2050_85) %in% vars_sel)
# make the projected models
# reusing maxnet_mod and maxnet_primed_mod from before

# compute and map maxnet predictions:
maxnet_mod_proj <- maxnet(p = dat$presence, data = dat[, vars_sel_proj], maxnet.formula(p = dat$presence, data = dat[, vars_sel_proj]))
maxnet_pred_proj <- rast(predict(raster::stack(layers_2050_85), maxnet_mod_proj, type = "cloglog"))

maxnet_prime_pred_proj <- rast(predict(raster::stack(primed_layers_2050_85), maxnet_mod_proj, type = "cloglog"))
maxnet_prime_diff_proj <- maxnet_prime_pred_proj - maxnet_pred_proj # positive values should be areas with increased predicted dist

countries_diff <- subset(countries, countries$NAME_0 %in% include)

ggplot() +
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
  ggtitle(paste(my_species, 'projected to 2050 RCP85 - diff')) +
  annotate(geom= 'label', label = paste(vars_sel_proj, sep = '\n', collapse = '\n'), x = Inf, y = -Inf, hjust = 1, vjust = 0) +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  )
plot(maxnet_pred_proj, main= 'projected to 2050 RCP85')
plot(maxnet_prime_pred_proj, main= 'projected to 2050 RCP85 PRIMED')
