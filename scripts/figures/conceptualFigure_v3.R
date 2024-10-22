library(dplyr)
library(ggplot2)

# create dataframe for plotting predictive conceptual fig

dfa <- data.frame(
  median = c(-0.6,-0.5,-0.4,-0.4),
  low_ci = c(-0.7,-0.6,-0.5,-0.5),
  high_ci = c(-0.5,-0.4,-0.3,-0.3),
  pp = c("Vegetation Green-up",
         "Butterfly Emergence",
         "Bird Arrival",
         "Bird Breeding"),
  Scenario = rep("Scenario A",4)
)

dfb <- data.frame(
  median = c(-0.6,-0.5,-0.2,-0.2),
  low_ci = c(-0.7,-0.6,-0.3,-0.3),
  high_ci = c(-0.5,-0.4,-0.1,-0.1),
  pp = c("Vegetation Green-up",
         "Butterfly Emergence",
         "Bird Arrival",
         "Bird Breeding"),
  Scenario = rep("Scenario B",4)
)

## combine dataframes together into super dataframe
sdf <- rbind(dfa, dfb)
sdf$pp <- factor(sdf$pp, levels = rev(c("Vegetation Green-up",
                                        "Butterfly Emergence",
                                        "Bird Arrival",
                                        "Bird Breeding")))


pp <- ggplot(filter(sdf, Scenario == "Scenario B")) +
  geom_linerange(aes(x = median, xmin = low_ci, xmax = high_ci, y = pp)) +
  geom_point(aes(x = median, y = pp), size = 1.5) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  scale_x_continuous(limits = c(-0.7,0.25),
                     breaks = c(-0.5,-0.25,0,0.25),
                     labels = c("Earlier phenology", "", "No effect", "")) +
  labs(x = "Effect of GDD on phenology event", y = "") +
  theme_classic() +
  theme(axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
                                                       ends = "both")))

ggsave(filename = "FigOutputs/Manuscript/ConceptualFigure_v3.png", dpi = 400, 
       width = 4, height = 3)
