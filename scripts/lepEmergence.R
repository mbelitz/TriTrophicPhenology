##############
# emergence ~ GDD
##############

# specify run date --------------------------------------------------------

run_date <- '2024-07-17'


# load packages -----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(MCMCvis)


# read in lep phenometrics and gdd data -----------------------------------

gdd <- data.table::fread('data/gdd_calcs.txt')

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

# center gdd for each cell
leps_df <- leps_df %>% 
  filter(year >= 2002 & year <= 2017) # restrict to same years# 


#ggplot(leps_df, aes(x = spring_gdd_center)) +
#  geom_histogram() +
#  facet_wrap(cell~code)
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


# data for model ----------------------------------------------------------
DATA <- list(y = as.vector(leps_df$q5),
             gdd = as.vector(leps_df$spring_gdd_center),
             obs = as.vector(leps_df$obs_center), 
             CellCode = as.vector(leps_df$CellCodeFactor),
             code = as.vector(distinct(leps_df, 
                                       CellCodeFactor,codeFactor))$codeFactor,
             lat = as.vector(distinct(leps_df, 
                                      CellCodeFactor, lat_center))$lat_center[,1],
             N = length(leps_df$q5),
             M = length(unique(leps_df$CellCodeFactor)),
             O = length(unique(leps_df$codeFactor)))

# fit model ---------------------------------------------------------------

options(mc.cores = parallel::detectCores())
#set_cmdstan_path('c:/Users/Mike/Documents/.cmdstan/cmdstan-2.33.1/')
DELTA <- 0.99999
# TREE_DEPTH <- 12
STEP_SIZE <- 0.01
CHAINS <- 4
ITER <- 8000

#compile model - this will tell you if the model syntax is valid before running the model
mod <- cmdstanr::cmdstan_model('scripts/stan/leps.stan',
                               force_recompile = TRUE)
#sample
fit <- mod$sample(
  data = DATA,
  chains = CHAINS,
  iter_sampling = ITER / 2,
  iter_warmup = ITER / 2,
  parallel_chains = CHAINS,
  refresh = 500,
  adapt_delta = DELTA,
  step_size = STEP_SIZE)
# max_treedepth = TREE_DEPTH
# adapt_delta = DELTA
# step_size = STEP_SIZE


# save summary space ------------------------------------------------------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(fit, 
                  round = 4,
                  file_name = paste0('lepsEmergence-results-', run_date),
                  dir = paste0('results'),
                  # mkdir = paste0('bird-ml-phylo-vint-measured-', run_date),
                  mkdir = paste0('lepsEmergence-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('lepsEmergence-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('lepEmergence-data-', run_date),
                  cp_file = c(paste0('scripts/stan/leps.stan'), 
                              paste0('scripts/lepEmergence.R')),
                  cp_file_names = c(paste0('lepEmergence-', run_date, '.stan'),
                                    paste0('lepEmergence-', run_date, '.R')))

fig_dir <- paste0('results/lepEmergence-', run_date, '/')

# #to load model object in for manipulation post fitting
#fit <- readRDS("output/")

# #for interactive model diagnostics
#library(shinystan)
#lepEmergenece <- DATA$y
#shinystan::launch_shinystan(fit)
MCMCvis::MCMCplot(fit, 
                  params = 'beta')

MCMCvis::MCMCplot(fit, 
                  params = 'betaObs')

MCMCvis::MCMCplot(fit, 
                  params = 'gamma')

MCMCvis::MCMCplot(fit, 
                  params = 'mu_gamma')

MCMCvis::MCMCplot(fit, 
                  params = 'theta')