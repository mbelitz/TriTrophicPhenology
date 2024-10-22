##############
# fledge ~ GDD
##############

# specify run date --------------------------------------------------------

run_date <- '2024-07-17'


# load packages -----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(MCMCvis)


# read in lep phenometrics and gdd data -----------------------------------

gdd <- data.table::fread('data/gdd_calcs.txt')

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
PC1_df <- distinct(fledge_df, species, PC1) %>%
  mutate(PC1_sc = scale(PC1, scale = F))


# data for model ----------------------------------------------------------

DATA <- list(y = as.vector(fledge_df$juv_meanday),
             gdd = as.vector(fledge_df$spring_gdd_center),
             SpSt = as.vector(fledge_df$SpeciesStationFactor),
             Sp = as.vector(distinct(fledge_df, SpeciesStationFactor, speciesFactor))$speciesFactor,
             lat = as.vector(distinct(fledge_df, 
                                      SpeciesStationFactor, lat_center))$lat_center[,1],
             PC1 = PC1_df$PC1_sc[,1],
             N = length(fledge_df$juv_meanday),
             M = length(unique(fledge_df$SpeciesStationFactor)),
             O = length(unique(fledge_df$speciesFactor)))


# fit model ---------------------------------------------------------------
options(mc.cores = parallel::detectCores())
# set_cmdstan_path('c:/Users/Mike/Documents/.cmdstan/cmdstan-2.33.1/')
# DELTA <- 0.92
# TREE_DEPTH <- 12
# STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 4000

#compile model - this will tell you if the model syntax is valid before running the model
mod <- cmdstanr::cmdstan_model('scripts/stan/fledge.stan',
                               force_recompile = TRUE)

#sample
fit <- mod$sample(
  data = DATA,
  chains = CHAINS,
  iter_sampling = ITER / 2,
  iter_warmup = ITER / 2,
  parallel_chains = CHAINS,
  refresh = 500)
# max_treedepth = TREE_DEPTH
# adapt_delta = DELTA
# step_size = STEP_SIZE


# save summary space ------------------------------------------------------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(fit, 
                  round = 4,
                  file_name = paste0('fledge-results-', run_date),
                  dir = paste0('results'),
                  # mkdir = paste0('bird-ml-phylo-vint-measured-', run_date),
                  mkdir = paste0('fledge-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('fledge-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('fledge-data-', run_date),
                  cp_file = c(paste0('scripts/stan/fledge.stan'), 
                              paste0('scripts/fledge.R')),
                  cp_file_names = c(paste0('fledge-', run_date, '.stan'),
                                    paste0('fledge-', run_date, '.R')))

fig_dir <- paste0('results/fledge-', run_date, '/')

# #to load model object in for manipulation post fitting
 fit <- readRDS("fledge-2024-07-17/fledge-fit-2024-07-17.rds")

# #for interactive model diagnostics
# library(shinystan)
# shinystan::launch_shinystan(fit)


# make figs ---------------------------------------------------------------

# gamma = effect of GDD on fledge at mean lat for each species
# mu_theta = cross-species effect of latitude on GDD effect
# kappa = effect of GDD on fledge at mean species PC1
# phi = effect of PC1 on species-level GDD effect

MCMCvis::MCMCsummary(fit, 
                     params = c('kappa',
                                'phi'),
                     round = 3, pg0 = TRUE)

MCMCvis::MCMCplot(fit, 
                  params = 'gamma')

MCMCvis::MCMCplot(fit, 
                  params = 'theta')


MCMCvis::MCMCplot(fit, 
                  params = 'kappa')



MCMCvis::MCMCplot(fit, 
                  params = 'mu_theta')


MCMCvis::MCMCplot(fit, 
                  params = 'phi')


# visualize results
# visualize results
beta_df <- data.frame(beta = MCMCvis::MCMCpstr(fit, params = 'beta', 
                                               fun = mean),
                      CellCode = 1:DATA$M)
alpha_df <- data.frame(alpha = MCMCvis::MCMCpstr(fit, params = 'alpha', 
                                                 fun = mean),
                       CellCode = 1:DATA$M)

gamma_df <- data.frame(gamma = MCMCvis::MCMCpstr(fit, params = 'gamma', 
                                                 fun = mean),
                       species = 1:DATA$O)

mu_gamma_df <- data.frame(gamma = MCMCvis::MCMCpstr(fit, params = 'mu_gamma', 
                                                    fun = mean),
                          species = 1:DATA$O)

theta_df <- data.frame(gamma = MCMCvis::MCMCpstr(fit, params = 'theta', 
                                                 fun = mean),
                       species = 1:DATA$O)

## test for phylogenetic signal in residuals
# call in phylogenetic bird tree
library(rtrees); library(ape); library(phytools)
bird_list <- read.csv("data/phylo_names_key-2021-05-19.csv")

sp_list_df <- sp_list_df(sp_list = bird_list$phylo_sci_name, taxon = "bird")

bird_tree <- get_tree(sp_list = sp_list_df,
                      taxon = "bird", 
                      scenario = "at_basal_node",
                      show_grafted = TRUE)

consensusTree <- ls.consensus(bird_tree)
consensusTree <- ladderize(consensusTree)

# call in model fit
## test for phylogenetic signal in gamma
# call in phylogenetic bird tree
library(rtrees); library(ape); library(phytools)
bird_list <- read.csv("data/phylo_names_key-2021-05-19.csv")

sp_list_df <- sp_list_df(sp_list = bird_list$phylo_sci_name, taxon = "bird")

bird_tree <- get_tree(sp_list = sp_list_df,
                      taxon = "bird", 
                      scenario = "at_basal_node",
                      show_grafted = TRUE)

consensusTree <- ls.consensus(bird_tree)
consensusTree <- ladderize(consensusTree)

# call in model fit
resids_df <- left_join(gamma_df, mu_gamma_df) %>% 
  mutate(resids = mu_gamma - gamma)

species <- fledge_df$species
sci_name <- str_replace(fledge_df$sci_name, "_", " ")

species_df <- data.frame(species, sci_name) %>% 
  distinct(species,sci_name)

resids_df <- left_join(gamma_df, mu_gamma_df) %>% 
  mutate(resids = mu_gamma - gamma) %>% 
  left_join(species_df) %>% 
  left_join(bird_list)

resids_df <-resids_df %>% 
  mutate(phylo_sci_name = str_replace(phylo_sci_name, " ", "_"))

fledgePhy <- keep.tip(phy = consensusTree, tip = resids_df$phylo_sci_name)

resids_df <- setNames(object = resids_df$resids, nm = resids_df$phylo_sci_name)

phylosig(fledgePhy, resids_df, method = "lambda", test = TRUE) # no evidence of phylogenetic signal
#Phylogenetic signal lambda : 7.57987e-05 
#logL(lambda) : 144.288 
#LR(lambda=0) : -0.00121655 
#P-value (based on LR test) : 1 