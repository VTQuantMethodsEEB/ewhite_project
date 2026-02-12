##################### Packages #####################
library(dplyr)
library(ggplot2)
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
# Create boxplot with data points overlayed 
# Max Sprint Speed
ggplot(norris_dat_filtered_new, aes(x=treatment, y=hatchling_max_speed, color = parental_island)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width=0.2), size=2, alpha=0.6) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Hatchling Max Speed") +
  labs(color="Parental Island")

# Create boxplot with overlayed datapoints
# Hatchling SVL
ggplot(norris_dat_filtered_new, aes(x=treatment, y=hatch_svl, color = parental_island)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width=0.2), size=2, alpha=0.6) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Hatch SVL") +
  labs(color="Parental Island")

# Sprint Trial SVL
ggplot(norris_dat_filtered_new, aes(x=treatment, y=sprint_svl, color = parental_island)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width=0.2), size=2, alpha=0.6) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Sprint Speed SVL") +
  labs(color="Parental Island")

# Egg Mass
ggplot(norris_dat_filtered_new, aes(x=treatment, y=egg_mass, color = parental_island)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width=0.2), size=2, alpha=0.6) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Egg Mass") +
  labs(color="Parental Island")

# Hatchling Mass
ggplot(norris_dat_filtered_new, aes(x=treatment, y=hatch_mass, color = parental_island)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width=0.2), size=2, alpha=0.6) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Hatchling Mass") +
  labs(color="Parental Island")

# Water Loss (change in body mass)
ggplot(norris_dat_filtered_new, aes(x=treatment, y=log(water_loss_g), color = parental_island)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width=0.2), size=2, alpha=0.6) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Water Loss") +
  labs(color="Parental Island")











