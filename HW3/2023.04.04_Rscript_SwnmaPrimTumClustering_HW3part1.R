## Vestibular Schwannoma Primary Tumor Clustering Data Analysis
## Author: Mark Dougherty
## Date started: 2023-04-04
##  - note - ongoing work

# Source of AutoPipe: https://github.com/falafel19/AutoPipe
# Reference: Masalha et al. Metabolic alterations in meningioma reflect the clinical course. BMC Cancer 21, 211 (2021). https://rdcu.be/c5yzG

# If not yet installed, run the following to install the tidyverse and devtools (required for this script):
# install.packages("tidyverse")
# install.packages("devtools")

# Load tidyverse & devtools
library(tidyverse)
library(devtools)

# Import raw data file from CSV to tibble using read_csv
vs_primary_metabolomics_raw <- read_csv("C:/Users/mark1/Dropbox/BIOL_4386/Project_Folder/Formatted_Data/230117_Swnma_PrimaryTumorsGCLC_09.2022_Batch.csv")
# Remove missing metabolite columns and save as curated tibble. Note that columns 1-6 are metadata.
vs_primary_curated <- vs_primary_metabolomics_raw %>% select(!where(is_logical))
metabolite_names_primary <- colnames(vs_primary_curated)[-(1:6)]
view(vs_primary_curated)

# If AutoPipe is not yet installed, use the following to install from GitHub using devtools:
install_github("falafel19/AutoPipe")
##NOTE TO SELF: AS of 4.5.23 AutoPipe not yet working - gives errors when attempt to install

## AutoPipe code (from https://github.com/falafel19/AutoPipe/README.md)
#Load data with Gene ENTREZ in rownames and samples in colnames
#data(y)
#dim(data)
#Optional: Read in clinical Infos with samples in rownames
#UnSuperClassifier(data,clinical_data=NULL,thr=2)

