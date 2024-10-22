library(tidyverse)
library(cmdstanr)
library(MCMCvis)
library(tidybayes)
library(sf)
library(terra)

# read in hex cells
hex_sf <- vect("data/hexGrids/hex_grid_crop.shp") %>% 
  st_as_sf() %>% 
  st_transform(hex_sf, crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

usca <- rnaturalearth::ne_countries(country = c("United States of America", "Canada"),
                                    returnclass = "sf", scale = 10) %>% 
  st_transform(crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

###################### Greenup ################################################
###############################################################################
###############################################################################
# read in arrival phenometrics and gdd data -----------------------------------

gdd <- data.table::fread('data/gdd_calcs.txt')

greenup <- readRDS("data/pheno-data-2020-08-25.rds") %>% 
  filter(!is.na(gr_mn)) %>% 
  distinct(cell,year, .keep_all = TRUE)

arrive <- greenup %>% 
  filter(!is.na(arr_GAM_mean))

# reduce data areas east of 95 degrees W and 
# reducedata to areas North of 24 degrees N (Youngflesh et al., 2021) & S of 50 degrees N
greenup <- greenup %>% 
  filter(cell_lat > 24 & cell_lat < 50) %>% 
  filter(cell_lng > -95)

# process data -------------------------------------------------------------
#join fledge and GDD
# select limited number of data columns
gu_df <- ungroup(greenup) %>%
  dplyr::left_join(gdd) %>% 
  dplyr::select(cell, year, cell_lat, spring.gdd, gr_mn)

#ggplot(gu_df, aes(x = spring.gdd)) +
#  geom_histogram() +
#  facet_wrap(cell~species)

# center gdd for each cell
gu_df <- gu_df %>% 
  filter(year >= 2002 & year <= 2017) %>%  # restrict to same years
  group_by(cell) %>% 
  mutate(spring_gdd_center = scale(spring.gdd, scale = F)) %>% # center gdd by cell
  ungroup()
  
# center lat across all cells
gu_df <- ungroup(gu_df) %>% 
  mutate(lat_center = scale(cell_lat, scale = F))

# only include cells with at least three years of data
yearsGrouping <- gu_df %>% 
  group_by(cell) %>% 
  summarise(nYears = length(unique(year)))
CellToKeep <- filter(yearsGrouping, nYears >= 3)
gu_df <- gu_df %>% 
  filter(cell %in% CellToKeep$cell)

gu_df <- gu_df %>% 
  mutate(cellFactor = cell)

gu_df$cellFactor <- factor(gu_df$cellFactor)
gu_df$cellFactor <- droplevels(gu_df$cellFactor)
gu_df$cellFactor <- as.integer(gu_df$cellFactor)

#order by cell to make sure lat, cell match
gu_df <- gu_df %>%
  arrange(cell, year)

# now read in greenup model data
gu_data <- readRDS("greenup-2024-07-17/greenup-data-2024-07-17.rds")
gu_fit <- readRDS("greenup-2024-07-17/greenup-fit-2024-07-17.rds")

beta_df <- data.frame(beta = MCMCvis::MCMCpstr(gu_fit, params = 'beta', 
                                               fun = mean),
                      cell = 1:gu_data$M)

beta_df2 <- distinct(gu_df, cell, cell_lat) %>% 
  bind_cols(beta_df) %>% 
  rename(cell = 1)

beta_sf <- left_join(hex_sf, beta_df2) %>% 
  filter(!is.na(cell_lat))

greenup_sensitivity_plot <- ggplot() +
  geom_sf(beta_sf, mapping = aes(fill = beta)) +
  geom_sf(usca, mapping = aes(), fill = NA) +
  scale_fill_binned(type = "viridis", option = "mako",
                    breaks = c(-0.20,-0.15,-0.1,-0.08,-0.06,-0.04,-0.03,-0.02, -0.01,0),
                    limits = c(-0.25,0)) +
  labs(fill = "Sensitivity to GDD") +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
  ggtitle("(A) Green-up") +
  theme_void()
greenup_sensitivity_plot
###################### Leps ################################################
###############################################################################
###############################################################################
# read in lep phenometrics and gdd data -----------------------------------
leps <- read.csv("data/adult_bfly_phenometrics_noCountCircles_withFull2020Data.csv") %>% 
  rename(cell = HEXcell)

# process data -------------------------------------------------------------
leps_gdd <- left_join(leps, gdd)

# select limited number of data columns
leps_df <- leps_gdd %>% 
  dplyr::select(q5, year, cell, code, uniqObsDays, spring.gdd)

#ggplot(leps_df, aes(x = spring.gdd)) +
#  geom_histogram() +
#  facet_wrap(cell~code)

# restrict to same years
leps_df <- leps_df %>% 
  filter(year >= 2002 & year <= 2017) # restrict to same years# 

# remove rows without mean emergence estimate
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

# center latitude by code
leps_df <- leps_df %>% 
  group_by(code) %>% 
  mutate(lat_center = scale(cell_lat, scale = F)) %>% #center gdd by code
  ungroup() %>%
  group_by(cell, code) %>% 
  mutate(spring_gdd_center = scale(spring.gdd, scale = F)) %>% # center gdd by cell
  ungroup()

# center obs overall
leps_df <- leps_df %>% 
  mutate(obs_center = scale(uniqObsDays, scale = F))

# add CellCodeFactor and codeFactor
leps_df <- leps_df %>% 
  mutate(CellCodeFactor = CellCode,
         codeFactor = code)

# munge data for future model building
leps_df$CellCodeFactor <- factor(leps_df$CellCodeFactor)
leps_df$CellCodeFactor <- droplevels(leps_df$CellCodeFactor)
leps_df$CellCodeFactor <- as.integer(leps_df$CellCodeFactor)

leps_df$codeFactor <- factor(leps_df$codeFactor)
leps_df$codeFactor <- droplevels(leps_df$codeFactor)
leps_df$codeFactor <- as.integer(leps_df$codeFactor)

#order by CellCode to make sure code and lat match properly
leps_df <- leps_df %>%
  arrange(CellCodeFactor, year)

# now read in lepEmergence model data
leps_data <- readRDS("lepsEmergence-2024-07-17/lepEmergence-data-2024-07-17.rds")
leps_fit <- readRDS("lepsEmergence-2024-07-17/lepsEmergence-fit-2024-07-17.rds")

beta_df_leps <- data.frame(beta = MCMCvis::MCMCpstr(leps_fit, params = 'beta', 
                                               fun = mean),
                      CellCode = 1:leps_data$M)

beta_df_leps2 <- distinct(leps_df, cell, cell_lat, code, CellCode) %>% 
  bind_cols(beta_df_leps) %>% 
  rename(CellCode = 4)

beta_leps_sf <- left_join(hex_sf, beta_df_leps2) %>% 
  filter(!is.na(cell_lat)) %>% 
  filter(code == "RP") # change to RP

leps_sensitivity_plot <- ggplot() +
  geom_sf(beta_leps_sf, mapping = aes(fill = beta)) +
  geom_sf(usca, mapping = aes(), fill = NA) +
#  scale_fill_viridis_c(option = "mako", limits = c(-0.12, 0), na.value = "black") +
  scale_fill_binned(type = "viridis", option = "mako",
                    breaks = c(-0.20,-0.15,-0.1,-0.08,-0.06,-0.04,-0.03,-0.02, -0.01,0),
                    limits = c(-0.25,0)) +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
  labs(fill = "Sensitivity to GDD") +
  ggtitle("(B) Butterfly Emergence (Pupae OWS)") +
  theme_void()
leps_sensitivity_plot

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

# center gdd for each cell-species combination
arr_df <- arr_df %>% 
  filter(year >= 2002 & year <= 2017) %>% # restrict to same years
  group_by(cell, species) %>% 
  mutate(spring_gdd_center = scale(spring.gdd, scale = F)) %>% # center gdd by cell
  ungroup()

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
#center latitude by species
arr_df <- arr_df %>% 
  group_by(species) %>% 
  mutate(lat_center = scale(cell_lat, scale = F)) %>% # center latitude by species
  ungroup()

#order by SpeciesCell to make sure lat, Sp match, and check to make sure PC1 matches with Sp
arr_df <- arr_df %>%
  arrange(SpeciesCell, year)

# now read in arrival model data
arr_data <- readRDS("arrive-2024-07-17/arrive-data-2024-07-17.rds")
arr_fit <- readRDS("arrive-2024-07-17/arrive-fit-2024-07-17.rds")

beta_df_arr <- data.frame(beta = MCMCvis::MCMCpstr(arr_fit, params = 'beta', 
                                                    fun = mean),
                           SpeciesCell = 1:arr_data$M)

beta_df_arr2 <- distinct(arr_df, cell, cell_lat, species, SpeciesCell) %>% 
  bind_cols(beta_df_arr) %>% 
  rename(SpeciesCell = 4)

beta_arr_sf <- left_join(hex_sf, beta_df_arr2) %>% 
  filter(!is.na(cell_lat)) %>% 
  filter(species == "Geothlypis_trichas") # filter to 'most average' species

arrival_sensitivity_plot <- ggplot() +
  geom_sf(beta_arr_sf, mapping = aes(fill = beta)) +
  geom_sf(usca, mapping = aes(), fill = NA) +
#  scale_fill_viridis_c(option = "mako") +
  scale_fill_binned(type = "viridis", option = "mako",
                    breaks = c(-0.20,-0.15,-0.1,-0.08,-0.06,-0.04,-0.03,-0.02, -0.01,0),
                    limits = c(-0.25,0)) +
  labs(fill = "Sensitivity to GDD") +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
  ggtitle("(C) Bird Arrival (Common Yellowthroat)") +
  theme_void()
arrival_sensitivity_plot

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

# now read in fledge model data
fledge_data <- readRDS("fledge-2024-07-17/fledge-data-2024-07-17.rds")
fledge_fit <- readRDS("fledge-2024-07-17/fledge-fit-2024-07-17.rds")

beta_df_fledge <- data.frame(beta = MCMCvis::MCMCpstr(fledge_fit, params = 'beta', 
                                                   fun = mean),
                          SpeciesStation = 1:fledge_data$M)

beta_df_fledge2 <- distinct(fledge_df, cell, lat_center, SpeciesStation) %>% 
  bind_cols(beta_df_fledge) %>% 
  rename(SpeciesCell = 4)

beta_fledge_sf <- left_join(hex_sf, beta_df_fledge2) %>% 
  filter(species == "Geothlypis_trichas") # filter to 'most average' arrival species

# for double beta estimates in the same cell, take average?

fledge_sensitivity_plot <- ggplot() +
  geom_sf(beta_fledge_sf, mapping = aes(fill = beta)) +
  geom_sf(usca, mapping = aes(), fill = NA) +
#  scale_fill_viridis_c(option = "mako") +
  scale_fill_binned(type = "viridis", option = "mako",
                    breaks = c(-0.20,-0.15,-0.1,-0.08,-0.06,-0.04,-0.03,-0.02, -0.01,0),
                    limits = c(-0.25,0)) +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
  labs(fill = "Sensitivity to GDD") +
  ggtitle("(D) Bird Breeding (Common Yellowthroat)") +
  theme_void()
fledge_sensitivity_plot

#### ggarrange plots together
library(ggpubr)
ga <- ggarrange(greenup_sensitivity_plot, leps_sensitivity_plot,
                arrival_sensitivity_plot, fledge_sensitivity_plot, 
                common.legend = TRUE, legend = "right")

ga

ggsave(filename = "FigOutputs/gddSensitivity_pupaeOWS.png", plot = ga, 
       width = 9, height = 8, dpi = 450)


pp <- ggplot() +
  geom_sf(beta_leps_sf, mapping = aes(fill = beta)) +
  geom_sf(usca, mapping = aes(), fill = NA) +
  #  scale_fill_viridis_c(option = "mako", limits = c(-0.12, 0), na.value = "black") +
  scale_fill_binned(type = "viridis", option = "mako",
                    breaks = c(-0.20,-0.15,-0.1,-0.08,-0.06,-0.04,-0.03,-0.02, -0.01,0),
                    limits = c(-0.25,0)) +
  coord_sf(xlim = c(-22260, 3003338),
           ylim = c(-1664985,1754536)) +
  labs(fill = "Sensitivity to GDD") +
  ggtitle("Butterfly Emergence (Pupae OWS)") +
  theme_void()
pp

ggsave(filename = "FigOutputs/pupaeOWS_gddSensitivityToLat.png", 
       plot = pp, dpi = 450)
