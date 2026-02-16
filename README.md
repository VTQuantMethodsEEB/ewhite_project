# ewhite_project

Measuring GxE with reaction norms from studies of squamates.

The data is from the study, Norris et al. 2024, where the researchers investigated the developmental plasticity of phenotypic traits for two populations of brown anole lizards (Anolis sagrei) found on different islands (open and shade). Eggs from either of the two populations were exposed to two different incubation treatments which varied in temperature, moisture and substrate to match the two unique island environments. The dataset includes measurements of each hatchling’s body mass, egg mass, body size (SVL), maximum sprint speed, and water loss. The dataset also includes egg ID, cage number (corresponds to family litter), parental island (either open or shade), and incubation treatment (either open or shade). There are about ~400 hatchlings used in this study. I would like to use the reaction norm and genetic data obtained from this study (and possibly others) to measure genotype-by-environment interaction through a function in R. I want to investigate if there are differences in GxE amongst trait types, sex, and populations.

##WEEK 2## - Created new Rscript Norris2024_GxE. Loaded in data called Norris2024_master_data. Filtered data and removed NAs. Practiced with dplyr functions.

##WEEK 3## 
Rscript - week3_EmmaWhite.R
Data - Norris2024_master_data.csv.
Created boxplots comparing each traits across treatments and grouped by parental environment. 