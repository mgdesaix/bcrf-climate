## Environmental variable selection

## Libraries
library(tidyverse)
library(terra)
library(corrplot)

##############################################################
## Raster processing
# Go through this process for both the 1991-2020 and 1961-1990 data sets

# read in rasters
tif.files.1961_1990 <- list.files("~/Downloads/Normal_1991_2020/", full.names = T)
elev.file <- "~/Downloads/elevation.tif"

# project rasters and crop to relevant geographic region
aw_1961_1990.ras <- terra::rast(tif.files.1961_1990) %>%
  project("epsg:4326") %>%
  terra::crop(ext(-109,-104,36,42))
aw_elevation.ras <- terra::rast(elev.file) %>%
  project("epsg:4326") %>%
  terra::crop(ext(-109,-104,36,42))

# get BCRF site data
bcrf.sites <- read_csv("../data/rofi-sites.csv") %>%
  dplyr::filter(!Near_Town %in% c("Humphrey Basin", "Virginia Lakes", "White Mountains"))

## create data frame of environmental data for sites

bcrf.env.1961_1990 <- terra::extract(aw_1961_1990.ras, bcrf.sites[,c("Lon", "Lat")])
bcrf.env.1991_2020 <- terra::extract(aw_1991_2020.ras, bcrf.sites[,c("Lon", "Lat")])
bcrf.elev <- terra::extract(aw_elevation.ras, bcrf.sites[,c("Lon", "Lat")])

bcrf.env.1961_1990$Elev <- bcrf.elev$elevation
bcrf.env.1991_2020$Elev <- bcrf.elev$elevation
bcrf.env.1961_1990 <- bcrf.env.1961_1990[,!colnames(bcrf.env.1961_1990) %in% c("ID", "MAR")] # remove MAR because it is NA for a site in one of the data sets
bcrf.env.1991_2020 <- bcrf.env.1991_2020[,!colnames(bcrf.env.1991_2020) %in% c("ID", "MAR")]

##############################################################
## Climate variation across Brown-capped Rosy-Finch sites

pca1 <- prcomp(bcrf.env.1961_1990, center = T, scale. = T)

pca1_varimp <- pca1$sdev/sum(pca1$sdev)
pca1_df <- data.frame(pca1$x) %>%
  add_column(Sites = bcrf.sites$Near_Town)

ggplot(pca1_df) +
  geom_point(aes(x = PC1, y = PC2, color = Sites))



##############################################################
## Environment correlations

env1_cor <- cor(bcrf.env.1991_2020)
env1_cor[abs(env1_cor) > 0.75] <- -1

corrplot(env1_cor)


