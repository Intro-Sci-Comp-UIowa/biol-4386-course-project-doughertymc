## Vestibular Schwannoma Xenograft Radiation Data Analysis R Script
## Author: Mark Dougherty
## Date started: 2023-04-01
##  - note - ongoing work after initial creation date on 2023-03-26

# If not yet installed, run the following to install the tidyverse (required for this script):
# install.packages("tidyverse")

# Load tidyverse
library(tidyverse)

# Import raw data file from CSV to tibble using read_csv
vs_xeno_metabolomics_raw <- read_csv("C:/Users/mark1/Dropbox/BIOL_4386/Project_Folder/Formatted_Data/230125 Swnma NewXeno wGroupID GCLC sorted.csv")

## NEW CONSOLIDATED CODE 4.1.23 for data wrangling - calculates fold changes as desired
vs_xeno_curated <- vs_xeno_metabolomics_raw %>% select(!where(is_logical)) %>% rename(primary_tumor = "Corresponding Primary/Xenograft", dose = "RT Dose (Gy) - Xenografts")
metabolite_names_xeno <- colnames(vs_xeno_curated)[-(1:5)]
vs_xeno_fc <- vs_xeno_curated %>% select(primary_tumor, dose, all_of(metabolite_names_xeno)) %>% group_by(primary_tumor) %>% mutate(across(all_of(metabolite_names_xeno), ~ .x / mean(.x[dose==0], na.rm=TRUE)))

##
##-------------------------- Data wrangling progress marker as of HW3 submission ---------------------------
##-------------------------- Beyond here is work in progress ---------------------------
## next up: need to perform normality tests for all metabolites (Shapiro-Wilk test), then list which metabolites are normally distributed and which are lognormal. For the lognormal metabolites, they must then be log-transformed prior to further statistical analysis. Once this is odne, basic stats for correlation with radiation treatment dose (per metabolite), plus Two-Way ANOVA with Holm-Sidak test, then lastly produce array of box plots for most noteworthy metabolites

### THE FOLLOWING IS AN EXAMPLE PLOT BUT IS NOT YET THE FINAL GOAL PRODUCT.
## EXAMPLE BOXPLOT
select(vs_xeno_metabolomics_raw, !where(is_logical)) %>%
  ggplot(data = ., mapping = aes(x = `RT Dose (Gy) - Xenografts`, y = `2-Hydroxybutyrate`)) + 
  geom_boxplot(aes(group=`RT Dose (Gy) - Xenografts`)) +
  geom_point(color = "blue") + 
  labs(title = "2-Hydroxybutyrate", x= "RT Dose (Gy)", y = "Raw Value")


# Fold change boxplot for 2-Hydroxybutyrate
ggplot(vs_xeno_fc, aes(x = `dose`, y = `2-Hydroxybutyrate`)) +
  geom_boxplot(aes(group=`dose`)) +
  geom_point(aes(color = `primary_tumor`)) + 
  labs(title = "2-Hydroxybutyrate", x = "RT Dose (Gy)", y = "Fold Change")

# Fold change boxplot for 2-Hydroxybutyrate - with color labels
ggplot(vs_xeno_fc, aes(x = `dose`, y = `2-Hydroxybutyrate`)) +
  geom_boxplot(aes(group=`dose`)) +
  geom_point(aes(color = `primary_tumor`)) + 
  labs(title = "2-Hydroxybutyrate", x = "RT Dose (Gy)", y = "Fold Change")

# Fold change boxplot for 2-Hydroxybutyrate - with color labels and facet wrapping
ggplot(vs_xeno_fc, aes(x = `dose`, y = `2-Hydroxybutyrate`)) +
  geom_boxplot(aes(group=`dose`)) +
  geom_point(aes(color = `primary_tumor`)) + 
  facet_wrap(vars(`primary_tumor`)) +
  labs(title = "2-Hydroxybutyrate", x = "RT Dose (Gy)", y = "Fold Change")
