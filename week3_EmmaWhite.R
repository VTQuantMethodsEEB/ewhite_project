##################### Packages #####################
library(dplyr)
library(ggplot2)
library(egg)
######## Norris 2024 dataset
norris_dat <- read.csv("Norris2024_master_data.csv")
### Filter out individuals that died before or during experiments
norris_dat_filtered <- norris_dat %>% dplyr::filter(!(egg_id %in% c(32, 65, 107, 144, 162, 171, 225, 234,
                                                                    258, 297, 298, 307, 345, 359, 361, 365, 370,
                                                                    390, 402, 412, 416, 420, 429)))
# Remove the comments and notes columns
norris_dat_filtered_new <- norris_dat_filtered[-c(12,13)]
# Remove NAs from the entire dataset
norris_dat_filtered_new <- na.omit(norris_dat_filtered_new)

################# ggplots ##########################
View(norris_dat_filtered_new)
# I want to display trait value for each individual measured from each treatment and group by parental island source.
# Create boxplot with data points overlayed 
# Max Sprint Speed
bp1 <- ggplot(norris_dat_filtered_new, aes(x=treatment, y=hatchling_max_speed, color = parental_island)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width=0.2), size=2, alpha=0.6) +
  theme_classic() +
  xlab(NULL) +
  ylab("Hatchling Max Speed") +
  labs(color="Parental Island") +
  theme(legend.position = "none")

# Create boxplot with overlayed datapoints
# Log Hatchling SVL
bp2 <- ggplot(norris_dat_filtered_new, aes(x=treatment, y=log(hatch_svl), color = parental_island)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width=0.2), size=2, alpha=0.6) +
  theme_classic() +
  xlab(NULL) +
  ylab("Log Hatchling SVL") +
  labs(color="Parental Island") +
  theme(legend.position = "none")

# Log Sprint Trial SVL
bp3 <- ggplot(norris_dat_filtered_new, aes(x=treatment, y=log(sprint_svl), color = parental_island)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width=0.2), size=2, alpha=0.6) +
  theme_classic() +
  xlab(NULL) +
  ylab("Log Sprint Speed Trial SVL") +
  labs(color="Parental Island") +
  theme(legend.position = "none")

# Log Egg Mass
bp4 <- ggplot(norris_dat_filtered_new, aes(x=treatment, y= log(egg_mass), color = parental_island)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width=0.2), size=2, alpha=0.6) +
  theme_classic() +
  xlab(NULL) +
  ylab("Log Egg Mass") +
  labs(color="Parental Island")

# Log Hatchling Mass
bp5 <- ggplot(norris_dat_filtered_new, aes(x=treatment, y=log(hatch_mass), color = parental_island)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width=0.2), size=2, alpha=0.6) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Log Hatchling Mass") +
  labs(color="Parental Island") +
  theme(legend.position = "none")

# Log Water Loss (change in body mass)
bp6 <- ggplot(norris_dat_filtered_new, aes(x=treatment, y=log(water_loss_g), color = parental_island)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width=0.2), size=2, alpha=0.6) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Log Water Loss") +
  labs(color="Parental Island") +
  theme(legend.position = "none")

# I want to show all the boxplots in a grid format (6x6)
ggarrange(bp1, bp2, bp3, bp4, bp5, bp6)









