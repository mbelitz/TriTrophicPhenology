library(tidyverse)
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


gamma_df_arr <- data.frame(gamma = MCMCvis::MCMCpstr(arr_fit, params = 'gamma', 
                                                     func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                           gamma =MCMCvis::MCMCpstr(arr_fit, params = 'gamma', 
                                                    fun = median),
                          species = 1:arr_data$O) %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)
  
gamma_df_arr2 <- distinct(arr_df, species) %>% 
  bind_cols(gamma_df_arr) %>% 
  rename(species = 1,
         speciesFactor = 7)  %>% 
  left_join(PC1_df)

arr_df <- left_join(arr_df, PC1_df)

#extract posteriors, 
# kappa is effect of GDD on fledge at mean species PC1
# phi is effect of PC1 on species-level GDD effect
post_ch_arr <- MCMCvis::MCMCchains(arr_fit, params = c('kappa', 'phi'))

xsim_arr <- seq(min(arr_df$PC1_sc), max(arr_df$PC1_sc), by = 0.01)
ysim_arr <- matrix(NA, nrow = NROW(post_ch_arr), ncol = length(xsim_arr))
for (i in 1:NROW(post_ch_arr))
{
  #i <- 1
  ysim_arr[i,] <- post_ch_arr[i,1] + post_ch_arr[i,2] * xsim_arr
}

fit_plot_arr <- data.frame(xsim_arr = xsim_arr,
                       CI_025 = apply(ysim_arr, 2, function(x) quantile(x, probs = 0.025)),
                       CI_25 = apply(ysim_arr, 2, function(x) quantile(x, probs = 0.25)),
                       CI_75 = apply(ysim_arr, 2, function(x) quantile(x, probs = 0.75)),
                       CI_975 = apply(ysim_arr, 2, function(x) quantile(x, probs = 0.975)),
                       median = apply(ysim_arr, 2, function(x) quantile(x, probs = 0.5))
)

## arrival sensitivity to GDD
arrivalPlot <- ggplot() +
  geom_ribbon(data = fit_plot_arr,
              aes(x = xsim_arr, ymin = CI_025, ymax = CI_975),
              fill = 'grey', alpha = 0.7) +
  geom_point(gamma_df_arr2, mapping = aes(x = PC1_sc, y = median),
             shape = 3, size = 2.5, alpha = 0.9) +
  geom_errorbar(gamma_df_arr2, mapping = aes(x = PC1_sc, ymin = CI_025, ymax = CI_975), alpha = 0.9) +
  geom_errorbar(gamma_df_arr2, mapping = aes(x = PC1_sc, ymin = CI_25, ymax = CI_75), linewidth = 1.25, alpha = 0.9) +
  geom_smooth(fit_plot_arr, mapping = aes(x = xsim_arr, y = median), method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(x = "", y = "Sensitivity to GDD") +
  scale_x_continuous(limits = c(-3,3), breaks = -3:3) +
  theme_classic() +
  theme(text = element_text(size = 20),
      #  axis.ticks = element_blank(),
        axis.text.x = element_blank())
arrivalPlot


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

gamma_df_fledge <- data.frame(gamma = MCMCvis::MCMCpstr(fledge_fit, params = 'gamma', 
                                        func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
              gamma =MCMCvis::MCMCpstr(fledge_fit, params = 'gamma', 
                                       fun = median),
              species = 1:fledge_data$O) %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)

gamma_df_fledge2 <- distinct(fledge_df, species) %>% 
  bind_cols(gamma_df_fledge) %>% 
  rename(species = 1,
         speciesFactor = 7) %>% 
  left_join(PC1_fledge_df)

#extract posteriors
post_ch <- MCMCvis::MCMCchains(fledge_fit, params = c('kappa', 'phi'))

xsim <- seq(min(fledge_df$PC1_sc), max(fledge_df$PC1_sc), by = 0.01)
ysim <- matrix(NA, nrow = NROW(post_ch), ncol = length(xsim))
for (i in 1:NROW(post_ch))
{
  #i <- 1
  ysim[i,] <- post_ch[i,1] + post_ch[i,2] * xsim
}

fit_plot <- data.frame(xsim = xsim,
                      CI_025 = apply(ysim, 2, function(x) quantile(x, probs = 0.025)),
                      CI_25 = apply(ysim, 2, function(x) quantile(x, probs = 0.25)),
                      CI_75 = apply(ysim, 2, function(x) quantile(x, probs = 0.75)),
                      CI_975 = apply(ysim, 2, function(x) quantile(x, probs = 0.975)),
                      median = apply(ysim, 2, function(x) quantile(x, probs = 0.5))
                      )

# extract posterior values to report in results
phi_df_arrival <- data.frame(phi = MCMCvis::MCMCpstr(arr_fit, params = 'phi',
                                                         func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                               phi =MCMCvis::MCMCpstr(arr_fit, params = 'phi',
                                                        fun = median)) %>% 
  mutate(Pheno.Phase = "Bird Arrival") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)

# extract posterior values to report in results
phi_df_arrival_fledge <- data.frame(phi = MCMCvis::MCMCpstr(fledge_fit, params = 'phi',
                                                     func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                             phi =MCMCvis::MCMCpstr(fledge_fit, params = 'phi',
                                                    fun = median)) %>% 
  mutate(Pheno.Phase = "Bird Fledge") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)


# fledge Sensitivity to GDD plot across PC1
fledgePlot <- ggplot() +
  geom_ribbon(data = fit_plot,
              aes(x = xsim, ymin = CI_025, ymax = CI_975),
              fill = 'grey', alpha = 0.7) +
  geom_point(gamma_df_fledge2, mapping = aes(x = PC1_sc, y = median),
             shape = 3, size = 2.5, alpha = 0.9) +
  geom_errorbar(gamma_df_fledge2, mapping = aes(x = PC1_sc, ymin = CI_025, ymax = CI_975), alpha = 0.9) +
  geom_errorbar(gamma_df_fledge2, mapping = aes(x = PC1_sc, ymin = CI_25, ymax = CI_75), linewidth = 1.25, alpha = 0.9) +
  geom_smooth(fit_plot, mapping = aes(x = xsim, y = median), method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_x_continuous(limits = c(-3,3), breaks = -3:3) +
  labs(x = "Trait Score", y = "Sensitivity to GDD") +
  theme_classic() +
  theme(text = element_text(size = 20))

fledgePlot 


##PC1 Legends
a_legend <- ggplot()+ 
  geom_segment(aes(x = 1, y = -0.05, xend = 2, yend = -0.05),
                   arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = 1, y = -0.0575, xend = 2, yend = -0.0575),
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = 1, y = -0.065, xend = 2, yend = -0.065),
               arrow = arrow(length = unit(0.5, "cm"))) +
  
  geom_segment(aes(x = -1, y = -0.05, xend = -2, yend = -0.05),
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = -1, y = -0.0575, xend = -2, yend = -0.0575),
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = -1, y = -0.065, xend = -2, yend = -0.065),
               arrow = arrow(length = unit(0.5, "cm"))) +
  
  annotate("text", x = 0, y = -0.05, label = "Migration pace", size = 17/.pt) +
  annotate("text", x = 0, y = -0.0575, label = "Arrival date", size = 17/.pt) +
  annotate("text", x = 0, y = -0.065, label = "Winter latitude", size = 17/.pt) +
  
  annotate("text", x = 2.5, y = -0.05, label = "Fast", size = 17/.pt) +
  annotate("text", x = 2.5, y = -0.0575, label = "Late", size = 17/.pt) +
  annotate("text", x = 2.5, y = -0.065, label = "South", size = 17/.pt) +
  
  annotate("text", x = -2.5, y = -0.05, label = "Slow", size = 17/.pt) +
  annotate("text", x = -2.5, y = -0.0575, label = "Early", size = 17/.pt) +
  annotate("text", x = -2.5, y = -0.065, label = "North", size = 17/.pt) +
  
  scale_y_continuous(limits = c(-0.07, -0.045)) +
  theme_void() 

b_legend <- ggplot()+ 
  geom_segment(aes(x = 1, y = -0.05, xend = 2, yend = -0.05),
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = 1, y = -0.065, xend = 2, yend = -0.065),
               arrow = arrow(length = unit(0.5, "cm"))) +
  
  geom_segment(aes(x = -1, y = -0.05, xend = -2, yend = -0.05),
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = -1, y = -0.065, xend = -2, yend = -0.065),
               arrow = arrow(length = unit(0.5, "cm"))) +

  
  annotate("text", x = 0, y = -0.05, label = "Migration distance", size = 17/.pt) +
  annotate("text", x = 0, y = -0.065, label = "Breeding timing", size = 17/.pt) +
  
  annotate("text", x = 2.5, y = -0.05, label = "Long", size = 17/.pt) +
  annotate("text", x = 2.5, y = -0.065, label = "Late", size = 17/.pt) +
  
  annotate("text", x = -2.5, y = -0.05, label = "Short", size = 17/.pt) +
  annotate("text", x = -2.5, y = -0.065, label = "Early", size = 17/.pt) +
  
  scale_y_continuous(limits = c(-0.07, -0.045)) +
  theme_void()


# ggarrange
#library(ggpubr)

#ga <- ggarrange(arrivalPlot, fledgePlot,
#                a_legend, b_legend,
#                labels = c(" (A) Bird Arrival",
#                           "(B) Bird Breeding"),
#                ncol = 2, nrow = 2, label.x = 0, label.y = 1, 
#                widths = c(1,1,0.5,0.5),
#                heights = c(1,1,0.25, 0.25)
#                )
#ga


library(cowplot)
cp <- plot_grid(arrivalPlot, fledgePlot, a_legend, align = "v", 
          nrow = 3, ncol = 1,
          rel_heights = c(0.7,0.7,0.15),
          labels = c(" (A) Bird Arrival", "(B) Bird Breeding"),
          label_x = 0.1)      

ggsave(filename = "FigOutputs/Manuscript/PC1_SensitivityToGDD_vertical.png", 
       plot = cp, width = 6.8, height = 7.5, 
       dpi = 600)
