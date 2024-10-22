##############
# arr ~ GDD
##############

# specify run date --------------------------------------------------------

run_date <- '2024-07-17'


# load packages -----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(MCMCvis)


# read in lep phenometrics and gdd data -----------------------------------

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

#ggplot(arr_df, aes(x = lat_center)) +
#  geom_histogram() +
#  facet_wrap(~species)
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


# data for model ----------------------------------------------------------

DATA <- list(y = as.vector(arr_df$arr_GAM_mean),
             gdd = as.vector(arr_df$spring_gdd_center),
             SpSt = as.vector(arr_df$SpeciesCellFactor),
             Sp = as.vector(distinct(arr_df, SpeciesCellFactor, speciesFactor))$speciesFactor,
             lat = as.vector(distinct(arr_df, 
                                      SpeciesCellFactor, lat_center))$lat_center[,1],
             PC1 = PC1_df$PC1_sc[,1],
             N = length(arr_df$arr_GAM_mean),
             M = length(unique(arr_df$SpeciesCellFactor)),
             O = length(unique(arr_df$speciesFactor)))


# fit model ---------------------------------------------------------------

options(mc.cores = parallel::detectCores())
set_cmdstan_path('c:/Users/Mike/Documents/.cmdstan/cmdstan-2.33.1/')
# DELTA <- 0.92
# TREE_DEPTH <- 12
# STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 4000

#compile model - this will tell you if the model syntax is valid before running the model
mod <- cmdstanr::cmdstan_model('scripts/stan/arrive.stan',
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
                  file_name = paste0('arrival-results-', run_date),
                  dir = paste0('results'),
                  # mkdir = paste0('bird-ml-phylo-vint-measured-', run_date),
                  mkdir = paste0('arrive-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('arrive-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('arrive-data-', run_date),
                  cp_file = c(paste0('scripts/stan/arrive.stan'), 
                              paste0('scripts/arrive.R')),
                  cp_file_names = c(paste0('arrive-', run_date, '.stan'),
                                    paste0('arrive-', run_date, '.R')))

fig_dir <- paste0('results/arrive-', run_date, '/')

# #to load model object in for manipulation post fitting
fit <- readRDS("arrive-2024-07-17/arrive-fit-2024-07-17.rds")

# #for interactive model diagnostics
#library(shinystan)
#shinystan::launch_shinystan(fit)
MCMCvis::MCMCplot(fit, 
                  params = 'beta')

MCMCvis::MCMCplot(fit, 
                  params = 'gamma')

MCMCvis::MCMCplot(fit, 
                  params = 'mu_gamma')


MCMCvis::MCMCplot(fit, 
                  params = 'theta')


MCMCvis::MCMCplot(fit, 
                  params = 'mu_theta')


MCMCvis::MCMCplot(fit, 
                  params = 'kappa')


MCMCvis::MCMCplot(fit, 
                  params = 'phi')


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

species <- arr_df$species
sci_name <- str_replace(arr_df$sci_name, "_", " ")

species_df <- data.frame(species, sci_name) %>% 
  distinct(species,sci_name)

resids_df <- left_join(gamma_df, mu_gamma_df) %>% 
  mutate(resids = mu_gamma - gamma) %>% 
  left_join(species_df) %>% 
  left_join(bird_list)

resids_df <-resids_df %>% 
  mutate(phylo_sci_name = str_replace(phylo_sci_name, " ", "_"))

arrivePhy <- keep.tip(phy = consensusTree, tip = resids_df$phylo_sci_name)

resids_df <- setNames(object = resids_df$resids, nm = resids_df$phylo_sci_name)


phylosig(arrivePhy, resids_df, method = "lambda", test = TRUE) # no evidence of phylogenetic signal
#Phylogenetic signal lambda : 7.53262e-05 
#logL(lambda) : 203.844 
#LR(lambda=0) : -0.000878779 
#P-value (based on LR test) : 1 