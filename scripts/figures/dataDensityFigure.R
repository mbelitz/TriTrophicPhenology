library(tidyverse)
library(sf)
library(terra)

# data figure

# read in hex cells
hex_sf <- vect("data/hexGrids/hex_grid_crop.shp") %>% 
  st_as_sf() %>% 
  st_transform(hex_sf, crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

usca <- rnaturalearth::ne_countries(country = c("United States of America", "Canada"),
                                    returnclass = "sf", scale = 10) %>% 
  st_transform(crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

# greenup data 
# read in greenup phenometrics and gdd data -----------------------------------

gdd <- data.table::fread('data/gdd_calcs.txt')

greenup <- readRDS("data/pheno-data-2020-08-25.rds") %>% 
  filter(!is.na(gr_mn)) %>% 
  distinct(cell,year, .keep_all = TRUE)

arrive <- greenup %>% 
  filter(!is.na(arr_GAM_mean))

na <- rnaturalearth::ne_countries(country = c("United States of America", "Canada"), returnclass = "sf")

greenup <- greenup %>% 
  filter(cell_lat > 24 & cell_lat < 50) %>% 
  filter(cell_lng > -95)

# process data -------------------------------------------------------------
#join fledge and GDD
# select limited number of data columns
gu_df <- ungroup(greenup) %>%
  dplyr::left_join(gdd) %>% 
  dplyr::select(cell, year, cell_lat, spring.gdd, gr_mn)

# center gdd for each cell
gu_df <- gu_df %>% 
  group_by(cell) %>% 
  mutate(spring_gdd_center = scale(spring.gdd, scale = F)) %>% # center gdd by cell
  filter(year >= 2002 & year <= 2017) # restrict to same years

# center lat across all cells
gu_df <- ungroup(gu_df) %>% 
  mutate(lat_center = scale(cell_lat, scale = F)) # center latitude across cell


# only include cells with at least three years of data
yearsGrouping <- gu_df %>% 
  group_by(cell) %>% 
  summarise(nYears = length(unique(year)))
CellToKeep <- filter(yearsGrouping, nYears >= 3)
gu_df <- gu_df %>% 
  filter(cell %in% CellToKeep$cell)

gu_gb <- gu_df %>% 
  reframe(nYears = length(unique(year)), .by = cell)

gu_sf <- left_join(hex_sf, gu_gb) %>% 
  filter(!is.na(nYears))

plot_a <- ggplot() +
  geom_sf(gu_sf, mapping = aes(fill = nYears), alpha = 0.8) +
  geom_sf(usca, mapping = aes(), fill = NA) +
  scale_fill_gradient2(low = "#b5c2b7", mid = "#8c93a8", high = "#62466b", 
                       limits = c(1,16), midpoint = 8) +
  labs(fill = "Years") +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
 # ggtitle("(A) Green-up") +
  theme_void() +
  theme(legend.box = "horizontal", legend.position = c(0.9,0.3),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))

plot_a

#############lepidoptera data
# read in lep phenometrics and gdd data -----------------------------------
leps <- read.csv("data/adult_bfly_phenometrics_noCountCircles_withFull2020Data.csv") %>% 
  rename(cell = HEXcell)

# process data -------------------------------------------------------------
leps_gdd <- left_join(leps, gdd)

# select limited number of data columns
leps_df <- leps_gdd %>% 
  dplyr::select(q5, year, cell, code, uniqObsDays, spring.gdd)

# center gdd for each cell
leps_df <- leps_df %>% 
  group_by(cell, code) %>% 
  mutate(spring_gdd_center = scale(spring.gdd, scale = F)) %>% # center gdd by cell
  filter(year >= 2002 & year <= 2017) # restrict to same years# 

leps_df <- filter(leps_df, !is.na(q5))
# only include cells with at least three years of sampling
yearsGrouping <- leps_df %>% 
  group_by(cell) %>% 
  summarise(nYears = length(unique(year)))
cellsToKeep <- filter(yearsGrouping, nYears >= 3)
leps_df <- leps_df %>% 
  filter(cell %in% cellsToKeep$cell)
# add cell_lat to dataframe
arr <- readRDS("data/pheno-data-2020-08-25.rds") # for greenup values
lat_df <- dplyr::distinct(arr, cell, .keep_all = T) %>% 
  select(cell, cell_lat)
leps_df <- left_join(leps_df, lat_df, by = "cell")
leps_df <- mutate(leps_df, CellCode = paste(cell, code, sep = "_"))

leps_gb <- ungroup(leps_df) %>% 
  reframe(nYears = length(unique(year)),
          nCode = length(unique(code)), .by = cell)

leps_sf <- left_join(hex_sf, leps_gb) %>% 
  filter(!is.na(nYears))

leps_sf_ows <- st_centroid(leps_sf)


plot_b <- ggplot() +
  geom_sf(leps_sf, mapping = aes(fill = nYears), alpha = 0.8) +
  geom_sf(usca, mapping = aes(), fill = NA) +
  geom_sf(leps_sf_ows, mapping = aes(size = nCode), shape = 18) +
  scale_fill_gradient2(low = "#b5c2b7", mid = "#8c93a8", high = "#62466b", limits = c(1,16), midpoint = 8) +
  scale_size_continuous(range = c(1,4), breaks = c(1,2,3)) +
  labs(fill = "Years", size = "Overwintering\nstages") +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
 # ggtitle("(B) Butterfly Emergence") +
  theme_void() +
  theme(legend.box = "horizontal", legend.position = c(1.1,0.3),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))

plot_b

###################### Arrival ################################################
###############################################################################
###############################################################################
arrival <- readRDS("data/pheno-data-2020-08-25.rds") %>% 
  filter(!is.na(arr_GAM_mean))

bird_pc <- readRDS("data/bird_PC_vals.rds")
# process data -------------------------------------------------------------
#join fledge and GDD
arr_gdd <- dplyr::left_join(arrival, gdd) %>% 
  left_join(bird_pc)

# select limited number of data columns
arr_df <- ungroup(arr_gdd) %>%
  mutate(sci_name = species) %>% 
  dplyr::select(sci_name, arr_GAM_mean, year, cell, species, PC1, cell_lat, spring.gdd)

# center gdd for each staion-species combination
arr_df <- arr_df %>% 
  group_by(cell, species) %>% 
  mutate(spring_gdd_center = scale(spring.gdd, scale = F)) %>% # center gdd by cell
  filter(year >= 2002 & year <= 2017) # restrict to same years

arr_df <- mutate(arr_df, SpeciesCell = paste(species, cell, sep = "_"))
yearsGrouping <- arr_df %>% 
  group_by(SpeciesCell) %>% 
  summarise(nYears = length(unique(year)))
SpeciesCellToKeep <- filter(yearsGrouping, nYears >= 3)
arr_df <- arr_df %>% 
  filter(SpeciesCell %in% SpeciesCellToKeep$SpeciesCell)

#for now remove GDD see if that is a problem
arr_df <- filter(arr_df, !is.na(spring.gdd))
# remove na's with lat -- this needs to be taken care of -- why am i missing latitudes for some cells?
arr_df <- filter(arr_df, !is.na(cell_lat))
#center latitutde by species
arr_df <- arr_df %>% 
  group_by(species) %>% 
  mutate(lat_center = scale(cell_lat, scale = F)) # center latitude by species

arr_gb <- ungroup(arr_df) %>% 
  reframe(nYears = length(unique(year)),
          nCode = length(unique(species)), .by = cell)

arr_sf <- left_join(hex_sf, arr_gb) %>% 
  filter(!is.na(nYears))

arr_sf_spp <- st_centroid(arr_sf)

# arrival data plot
plot_c <- ggplot() +
  geom_sf(arr_sf, mapping = aes(fill = nYears), alpha = 0.8) +
  geom_sf(usca, mapping = aes(), fill = NA) +
  geom_sf(arr_sf_spp, mapping = aes(size = nCode)) +
  scale_fill_gradient2(low = "#b5c2b7", mid = "#8c93a8", high = "#62466b", 
                       limits = c(1,16), midpoint = 8) +
  scale_size_continuous(range = c(1,5), limits = c(0,55), breaks = c(10,20,30,40,50)) +
  labs(fill = "Years", size = "Species") +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
 # ggtitle("(C) Bird Arrival") +
  theme_void() +
  theme(legend.box = "horizontal", legend.position = c(0.9,0.3),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))

plot_c


### bird fledge
###################### Fledge #################################################
###############################################################################
###############################################################################
fledge <- readRDS("data/MAPS-fledge-dates-2022-02-22.rds") %>% 
  dplyr::mutate(species = stringr::str_replace(sci_name, pattern = " ",
                                               replacement = "_"))

maps <- read.csv("data/MAPS_station_locations.csv")
# process data -------------------------------------------------------------
#join fledge and GDD
fledge_gdd <- dplyr::left_join(fledge, gdd)
# select limited number of data columns
fledge_df <- ungroup(fledge_gdd) %>% 
  dplyr::select(sci_name, juv_meanday, year, cell, station, species, PC1, spring.gdd, lat, lng)
# center gdd for each staion-species combination
fledge_df <- fledge_df %>% 
  group_by(station, species) %>% 
  mutate(spring_gdd_center = scale(spring.gdd, scale = F)) %>% # center gdd by cell
  filter(year >= 2002 & year <= 2017) # restrict to same years

fledge_df <- filter(fledge_df, !is.na(juv_meanday))
# only include species/stations with at least four years of data
fledge_df <- mutate(fledge_df, SpeciesStation = paste(species, station, sep = "_"))
yearsGrouping <- fledge_df %>% 
  group_by(SpeciesStation) %>% 
  summarise(nYears = length(unique(year)))
SpeciesStationToKeep <- filter(yearsGrouping, nYears >= 3)
fledge_df <- fledge_df %>% 
  filter(SpeciesStation %in% SpeciesStationToKeep$SpeciesStation)

#remove stations without GDD
fledge_df <- filter(fledge_df, !is.na(spring.gdd)) %>% 
  filter(lat > 24 & lat < 49) %>% # now filter to study location
  filter(lng > -95)

#center latitude by species
fledge_df <- fledge_df %>% 
  group_by(species) %>% 
  mutate(lat_center = scale(lat, scale = F)) # center gdd by cell

fledge_gb <- ungroup(fledge_df) %>% 
  reframe(nYears = length(unique(year)),
          nCode = length(unique(species)), .by = c(station, lat, lng))

fledge_sf <- st_as_sf(fledge_gb, coords = c(x = "lng", y = "lat"), crs = "WGS84") %>% 
  st_transform(crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") %>% 
  filter(!is.na(nYears))

fledge_cells <- ungroup(fledge_df) %>% 
  reframe(nYears = length(unique(year)),
          nCode = length(unique(species)), .by = cell)

fledge_cells_sf <- st_join(hex_sf, fledge_sf) %>% 
  filter(!is.na(nYears))

# fledge data plot
plot_d <- ggplot() +
  geom_sf(usca, mapping = aes(), fill = NA) +
  geom_sf(fledge_cells_sf, mapping = aes(), fill = NA, color = "grey75") +
  geom_sf(fledge_sf, mapping = aes(size = nCode, color = nYears)) +
  scale_color_gradient2(low = "#b5c2b7", mid = "#8c93a8", high = "#62466b", limits = c(1,16), midpoint = 8) +
  scale_size_continuous(range = c(1,5), limits = c(0,55), breaks = c(10,20,30,40,50)) +
  labs(color = "Years", size = "Species") +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
#  ggtitle("(D) Bird Breeding") +
  theme_void() +
  theme(legend.box = "horizontal", legend.position = c(0.9,0.3),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))

plot_d

#### ggarrange plots together
library(ggpubr)
ga <- ggarrange(plot_b, plot_c, plot_d,#labels = c("(A)\nGreen-up", "(B)\nButterfly\nEmergence",
                #          "(C)\nBird\nArrival", "(D)\nBird\nBreeding"),
                label.y = 0.95, nrow = 3, ncol = 1)

ga

ggsave(filename = "FigOutputs/Manuscript/dataDensityFigure.png", plot = ga, 
       width = 5, height = 10, dpi = 450)

