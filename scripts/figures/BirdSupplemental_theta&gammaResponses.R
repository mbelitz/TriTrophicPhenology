library(tidyverse)
library(cmdstanr)
library(MCMCvis)
###################### Arrival ################################################
###############################################################################
###############################################################################
gdd <- data.table::fread('data/gdd_calcs.txt')

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

#ggplot(arr_df, aes(x = spring.gdd)) +
#  geom_histogram() +
#  facet_wrap(cell~species)

# center gdd for each cell-species combination
arr_df <- arr_df %>% 
  filter(year >= 2002 & year <= 2017) %>% # restrict to same years
  group_by(cell, species) %>% 
  mutate(spring_gdd_center = scale(spring.gdd, scale = F)) %>% # center gdd by cell
  ungroup()


#ggplot(arr_df, aes(x = spring_gdd_center)) +
#  geom_histogram() +
#  facet_wrap(cell~speies)
# only include species/cells with at least three years of data
arr_df <- mutate(arr_df, SpeciesCell = paste(species, cell, sep = "_"))
yearsGrouping <- arr_df %>% 
  group_by(SpeciesCell) %>% 
  summarise(nYears = length(unique(year)))
SpeciesCellToKeep <- filter(yearsGrouping, nYears >= 3)
arr_df <- arr_df %>% 
  filter(SpeciesCell %in% SpeciesCellToKeep$SpeciesCell)

#center latitutde by species
arr_df <- arr_df %>% 
  group_by(species) %>% 
  mutate(lat_center = scale(cell_lat, scale = F)) %>% # center latitude by species
  ungroup()

arr_df <- arr_df %>% 
  mutate(SpeciesCellFactor = SpeciesCell,
         speciesFactor = species)

# munge data for future model building
arr_df$SpeciesCellFactor <- factor(arr_df$SpeciesCellFactor)
arr_df$SpeciesCellFactor <- droplevels(arr_df$SpeciesCellFactor)
arr_df$SpeciesCellFactor <- as.integer(arr_df$SpeciesCellFactor)

arr_df$speciesFactor <- factor(arr_df$speciesFactor)
arr_df$speciesFactor <- droplevels(arr_df$speciesFactor)
arr_df$speciesFactor <- as.integer(arr_df$speciesFactor)

#order by SpeciesCell to make sure lat, Sp match, and check to make sure PC1 matches with Sp
arr_df <- arr_df %>%
  arrange(SpeciesCellFactor, year)

#scale PC1 to make sure that kappa is sens at mean PC1
PC1_df <- distinct(arr_df, species, PC1) %>%
  mutate(PC1_sc = scale(PC1, scale = F))

# now read in arrival model data, gamma is effect of GDD on fledge at mean lat for each species
arr_data <- readRDS("arrive-2024-07-17/arrive-data-2024-07-17.rds")
arr_fit <- readRDS("arrive-2024-07-17/arrive-fit-2024-07-17.rds")

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

#ggplot(fledge_df, aes(x = spring.gdd)) +
#  geom_histogram() +
#  facet_wrap(station~species)

# center gdd for each staion-species combination
fledge_df <- fledge_df %>% 
  filter(year >= 2002 & year <= 2017) %>% # restrict to same years
  group_by(station, species) %>% 
  mutate(spring_gdd_center = scale(spring.gdd, scale = F)) %>% # center gdd by cell
  ungroup()

#ggplot(fledge_df, aes(x = spring_gdd_center)) +
#  geom_histogram() +
#  facet_wrap(station~speies)
# remove rows without mean emergence estimate
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
  mutate(lat_center = scale(lat, scale = F)) %>% # center gdd by cell
  ungroup()

#make species station and species factor
fledge_df <- fledge_df %>% 
  mutate(SpeciesStationFactor = SpeciesStation,
         speciesFactor = species)


#ggplot(fledge_df, aes(x = lat_center)) +
#  geom_histogram() +
#  facet_wrap(~species)
# munge data for future model building
fledge_df$SpeciesStationFactor <- factor(fledge_df$SpeciesStationFactor)
fledge_df$SpeciesStationFactor <- droplevels(fledge_df$SpeciesStationFactor)
fledge_df$SpeciesStationFactor <- as.integer(fledge_df$SpeciesStationFactor)

fledge_df$speciesFactor <- factor(fledge_df$speciesFactor)
fledge_df$speciesFactor <- droplevels(fledge_df$speciesFactor)
fledge_df$speciesFactor <- as.integer(fledge_df$speciesFactor)

#order by SpeciesStation to make sure lat, Sp match, and check to make sure PC1 matches with Sp
fledge_df <- fledge_df %>%
  arrange(SpeciesStation, year)

#scale PC1 to make sure that kappa is sens at mean PC1
PC1_fledge_df <- distinct(fledge_df, species, PC1) %>%
  mutate(PC1_sc = scale(PC1, scale = F))

fledge_df <- left_join(fledge_df, PC1_fledge_df)

# now read in fledge model data
fledge_data <- readRDS("fledge-2024-07-17/fledge-data-2024-07-17.rds")
fledge_fit <- readRDS("fledge-2024-07-17/fledge-fit-2024-07-17.rds")

# supplemental response of leps by ows
gamma_df_arrival <- data.frame(theta = MCMCvis::MCMCpstr(arr_fit, params = 'gamma', 
                                                         func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                               theta =MCMCvis::MCMCpstr(arr_fit, params = 'gamma', 
                                                        fun = mean)) %>% 
  mutate(Pheno.Phase = "Bird Arrival") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)

gamma_df_arr2 <- distinct(arr_df, species) %>% 
  bind_cols(gamma_df_arrival) %>% 
  rename(species = 1,
         speciesFactor = 7)  %>% 
  left_join(bird_pc)

gamma_df_fledge <- data.frame(kappa = MCMCvis::MCMCpstr(fledge_fit, params = 'gamma', 
                                                        func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                              kappa =MCMCvis::MCMCpstr(fledge_fit, params = 'gamma', 
                                                       fun = mean))%>% 
  mutate(Pheno.Phase = "Bird Breeding") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)

gamma_df_fledge2 <- distinct(fledge_df, species) %>% 
  bind_cols(gamma_df_fledge) %>% 
  rename(species = 1,
         speciesFactor = 7)  %>% 
  left_join(bird_pc)

# latitude effect by species
theta_df_arrival <- data.frame(theta = MCMCvis::MCMCpstr(arr_fit, params = 'theta', 
                                                         func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                               theta =MCMCvis::MCMCpstr(arr_fit, params = 'theta', 
                                                        fun = mean)) %>% 
  mutate(Pheno.Phase = "Bird Arrival") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)

theta_df_arr2 <- distinct(arr_df, species) %>% 
  bind_cols(theta_df_arrival) %>% 
  rename(species = 1,
         speciesFactor = 7)  %>% 
  left_join(bird_pc)


theta_df_fledge <- data.frame(kappa = MCMCvis::MCMCpstr(fledge_fit, params = 'theta', 
                                                        func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                              kappa =MCMCvis::MCMCpstr(fledge_fit, params = 'theta', 
                                                       fun = mean))%>% 
  mutate(Pheno.Phase = "Bird Breeding") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)

theta_df_fledge2 <- distinct(fledge_df, species) %>% 
  bind_cols(theta_df_fledge) %>% 
  rename(species = 1,
         speciesFactor = 7)  %>% 
  left_join(bird_pc)


## effect of GDD at median latitude

gamma_effect_comb <- rbind(gamma_df_arr2, gamma_df_fledge2) %>% 
  mutate(species = str_replace(species, "_", " "))

gddEffectPlot <- ggplot(gamma_effect_comb, aes(y = species)) +
  geom_linerange(aes(x = median, xmin = CI_025, xmax = CI_975)) +
  geom_linerange(aes(x = median, xmin = CI_25, xmax = CI_75), linewidth = 1.5) +
  geom_point(aes(x = median), shape = 3, size = 2.5) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  labs(x = "Effect of GDD at mean latitude", y = "") +
  theme_classic() +
  theme(axis.text.y = element_text(face = "italic")) +
  facet_wrap(~speciesFactor)
gddEffectPlot

ggsave(filename = "FigOutputs/gddEffectBySpecies.png", dpi = 400, width = 8, height = 8)


## fct reorder version
#gddEffectPlot <- ggplot(gamma_effect_comb, aes(y = reorder(species, median))) +
#  geom_linerange(aes(x = median, xmin = CI_025, xmax = CI_975)) +
#  geom_linerange(aes(x = median, xmin = CI_25, xmax = CI_75), linewidth = 1.5) +
#  geom_point(aes(x = median), shape = 3, size = 2.5) +
#  geom_vline(aes(xintercept = 0), linetype = "dashed") +
#  labs(x = "Effect of GDD at median latitude", y = "") +
#  theme_classic() +
#  facet_wrap(~speciesFactor)


# plot effect of latitude on GDD sensitivity
theta_effect_comb <- rbind(theta_df_arr2, theta_df_fledge2) %>% 
  mutate(species = str_replace(species, "_", " "))


LatEffectPlot <- 
  ggplot(theta_effect_comb) +
  geom_linerange(aes(x = median, xmin = CI_025, xmax = CI_975, y = species)) +
  geom_linerange(aes(x = median, xmin = CI_25, xmax = CI_75, y = species), linewidth = 1.5) +
  geom_point(aes(x = median, y = species), shape = 3, size = 2.5) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  labs(x = "Effect of Latitude on GDD sensitivity", y = "") +
  theme_classic() +
  theme(axis.text.y = element_text(face = "italic")) +
  facet_wrap(~speciesFactor)

LatEffectPlot

ggsave(filename = "FigOutputs/latEffectBySpecies.png", dpi = 400, width = 8, height = 8)


