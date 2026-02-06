##################### Packages #####################
library(dplyr)
######## Norris 2024 dataset
norris_dat <- read.csv("Norris2024_master_data.csv")
############### Filtering data #####################
#### Examine data to remove (i.e., look at comments)
View(norris_dat)
### Filter out individuals that died before or during experiments
class(norris_dat$egg_id) # integer
norris_dat_filtered <- norris_dat %>% dplyr::filter(!(egg_id %in% c(32, 65, 107, 144, 162, 171, 225, 234,
                                                                  258, 297, 298, 307, 345, 359, 361, 365, 370,
                                                                  390, 402, 412, 416, 420, 429)))
View(norris_dat_filtered)
# Remove the comments and notes columns
norris_dat_filtered_new <- norris_dat_filtered[-c(12,13)]
View(norris_dat_filtered_new)











