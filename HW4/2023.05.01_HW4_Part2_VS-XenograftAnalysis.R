library(tidyverse)


# Import raw data file from CSV to tibble using read_csv
vs_xeno_metabolomics_raw <- read_csv("C:/Users/mark1/Dropbox/BIOL_4386/Project_Folder/Formatted_Data/230125 Swnma NewXeno wGroupID GCLC sorted.csv")

## NEW CONSOLIDATED CODE 4.1.23 for data wrangling - calculates fold changes as desired
vs_xeno_curated <- vs_xeno_metabolomics_raw %>% select(!where(is_logical)) %>% rename(primary_tumor = "Corresponding Primary/Xenograft", dose = "RT Dose (Gy) - Xenografts")
metabolite_names_xeno <- colnames(vs_xeno_curated)[-(1:5)]
vs_xeno_fc <- vs_xeno_curated %>% select(primary_tumor, dose, all_of(metabolite_names_xeno)) %>% group_by(primary_tumor) %>% 
  mutate(across(all_of(metabolite_names_xeno), ~ .x / mean(.x[dose==0], na.rm=TRUE)))


################# Outlier Detection - NOW WORKING AS OF 5.9.23 PM
##If not already installed, run: install.packages("outliers")
library(outliers)

# Run Grubbs' test for outliers on each column name in vs_xeno_fc that appears in the vector list metabolite_names_xeno
grubbs_results <- map_dfr(metabolite_names_xeno, ~{
  metabolite_col <- .x
  test_result_high <- grubbs.test(vs_xeno_fc[[metabolite_col]], opposite=FALSE, type=10)
  test_result_low <- grubbs.test(vs_xeno_fc[[metabolite_col]], opposite=TRUE, type=10)
  list(metabolite_name = metabolite_col,
       high_value = max(vs_xeno_fc[[metabolite_col]], na.rm=TRUE),
       p_value_high = test_result_high$p.value,
       low_value = min(vs_xeno_fc[[metabolite_col]], na.rm=TRUE),
       p_value_low = test_result_low$p.value)
})
# Filter for only the metabolites with p<0.01 on either high or low Grubbs test ('high' tests largest value, 'low' tests smallest value)
outlier_df <- grubbs_results %>% filter(., p_value_high<=0.01 | p_value_low<=0.01)
outlier_list <- outlier_df$metabolite_name
## NOTE THAT NONE OF THE LOW VALUES WERE SIGNIFICANT; 57 HIGH VALUES WERE SIGNIFICANT OUTLIERS PER GRUBBS TEST (alpha<0.01)
# Loop over each metabolite in outlier_df and replace the high value with NA in vs_xeno_fc
vs_xeno_fc_outliers_removed <- vs_xeno_fc
for (metabolite_name in outlier_df$metabolite_name) {
  vs_xeno_fc_outliers_removed <- vs_xeno_fc_outliers_removed %>%
    mutate(!!sym(metabolite_name) := if_else(!!sym(metabolite_name) == outlier_df$high_value[outlier_df$metabolite_name == metabolite_name], NA_real_, !!sym(metabolite_name)))
}

# APPENDIX (part 2.2): Print Grubbs' test results (grubbs_results) and the VS xenograft tibble with those outliers removed (vs_xeno_fc_outliers_removed)
view(grubbs_results)
view(vs_xeno_fc_outliers_removed)




######################
***********
  NOTE AS OF 5.9.23 AT 1830 ALL OF THE SUBSEQUENT ANALYSIS DOES NOT YET ACCOUNT FOR THE REMOVAL OF OUTLIERS
######################



################# Test for Normality with Shapiro-Wilk test  (NOTE: THIS DOES NOT YET ACCOUNT FOR VALUES REMOVED BY OUTLIER TEST)
# Initialize an empty tibble to store the p-values from the Shapiro-Wilk test
shapiro_pvalues <- tibble(metabolite = character(),
                          p_value = double())

# Loop over the outcome variables and perform the Shapiro-Wilk test
for (i in 1:length(metabolite_names_xeno)) {
  # Extract the outcome variable
  outcome_var <- metabolite_names_xeno[i]
  # Perform the Shapiro-Wilk test
  shapiro_test <- shapiro.test(vs_xeno_fc[[outcome_var]])
  # Store the variable name and p-value in the tibble
  shapiro_pvalues <- shapiro_pvalues %>% 
    add_row(metabolite = outcome_var, p_value = shapiro_test$p.value)
}
# Sort the tibble by p-values and filter for p-values > 0.05 (non-normally distributed) and < 0.05 (normally distributed)
# non-normal list:
shapiro_pvalues_nonnormal <- shapiro_pvalues %>% arrange(p_value) %>% filter(p_value > 0.05)
# normally distributed:
shapiro_pvalues_normal <- shapiro_pvalues %>% arrange(p_value) %>% filter(p_value <= 0.05)
# Create tibble of only normally distributed metabolites
vs_xeno_fc_normal <- vs_xeno_fc %>% group_by(primary_tumor, dose) %>% select(all_of(shapiro_pvalues_normal$metabolite))


################ Transformation of LogNormal metabolites
# Create new tibble with only the columns to log transform
vs_xeno_fc_nonnormal <- vs_xeno_fc %>% group_by(primary_tumor, dose) %>% select(all_of(shapiro_pvalues_nonnormal$metabolite))

# Log-transform the columns
vs_xeno_fc_log_transformed <- vs_xeno_fc_nonnormal %>% mutate_if(is.numeric, ~ ifelse(. > 0, log(.), NA))

# IF DESIRED: Rename the columns to indicate they are log-transformed
# colnames(vs_xeno_fc_log_transformed) <- paste0(colnames(vs_xeno_fc_log_transformed), "_logtransformed")
# IF DESIRED CAN RE-BIND THESE VALUES TO THE NORMALLY DISTRIBUTED VALUES USING CODE SIMILAR TO THE FOLLOWING:
## Combine the log-transformed columns with the rest of the original tibble
# my_tibble_transformed <- bind_cols(my_tibble %>% select(-all_of(var_names_to_log)), my_tibble_log_transformed)
## View the resulting tibble
# my_tibble_transformed

## APPENDIX (part 2.4): View results of normal and lognormal data (Appendix)
view(vs_xeno_fc_normal)
view(vs_xeno_fc_log_transformed)


################STATISTICAL TESTS (note: does not yet account for outliers)
######### Correlation with radiation
# Create an empty list to store the correlations
my_correlations <- list()
# Loop over the output variables and compute the correlations
for (i in 3:ncol(vs_xeno_fc)) {
  output_var <- names(vs_xeno_fc)[i]
  cor_test <- cor.test(vs_xeno_fc$dose, vs_xeno_fc[[output_var]], method = "pearson")
  my_correlations[[output_var]] <- cor_test$estimate
}
# Combine the correlations into a data frame
cor_df <- data.frame(output_var = names(my_correlations),
                     correlation = unlist(my_correlations))
cor_df_filtered_sorted <- cor_df %>% arrange(desc(correlation)) %>% filter(correlation > 0.2)

#### REPEAT THE ABOVE FOR ONLY THE NORMALLY DISTRIBUTED METABOLITES
my_correlations <- list()
for (i in 3:ncol(vs_xeno_fc_normal)) {
  output_var <- names(vs_xeno_fc_normal)[i]
  cor_test <- cor.test(vs_xeno_fc_normal$dose, vs_xeno_fc_normal[[output_var]], method = "pearson")
  my_correlations[[output_var]] <- cor_test$estimate
}
# Combine the correlations into a data frame
cor_df_normal <- data.frame(output_var = names(my_correlations),
                     correlation = unlist(my_correlations))
cor_df_normal_filtered_sorted <- cor_df_normal %>% arrange(desc(correlation)) %>% filter(correlation > 0.2)

#### REPEAT THE ABOVE FOR ONLY LOG-TRANSFORMED METABOLITES
my_correlations <- list()
for (i in 3:ncol(vs_xeno_fc_log_transformed)) {
  output_var <- names(vs_xeno_fc_log_transformed)[i]
  cor_test <- cor.test(vs_xeno_fc_log_transformed$dose, vs_xeno_fc_log_transformed[[output_var]], method = "pearson")
  my_correlations[[output_var]] <- cor_test$estimate
}
# Combine the correlations into a data frame
cor_df_log_normal <- data.frame(output_var = names(my_correlations),
                            correlation = unlist(my_correlations))
cor_df_log_normal_filtered_sorted <- cor_df_log_normal %>% arrange(desc(correlation)) %>% filter(correlation > 0.2)

#Create data frame with both normal and log transformed correlations that are >0.2 by pearson test and sorted from high to low
cor_df_all_filtered_sorted <- rbind(cor_df_normal_filtered_sorted, cor_df_log_normal_filtered_sorted) %>% arrange(desc(correlation))


# APPENDINX: view all the metabolites with correlation coeff > 0.2, sorted greatest to least:
view(cor_df_all_filtered_sorted)
# Print top 10 metabolites by correlation values:
cor_df_all_filtered_sorted[1:10,1]




################ Two-way ANOVA with Holm-Sidak test *****************INCOMPLETE - AS OF 5.9.23 PM, STILL LACKS POST-HOC ANALYSIS
######ALSO NOTE: DOES NOT ACCOUNT FOR OUTLIERS YET (as of 5.9.23 PM)
library(broom)
# Normally distributed metabolites:
## First pivot longer to reformat the output columns for anova function
normal_data_long <- vs_xeno_fc_normal %>%
  pivot_longer(cols = 3:ncol(.), names_to = "metabolite", values_to = "value")
# Then run the two-way ANOVA by the first two columns  across all metabolites
####NOTE: primary_tumor + dose does not include the interaction term between primary_tumor and dose (my_anova). If desire the interaction, change to primary_tumor * dose inside aov function
normal_anova <- normal_data_long %>%
  group_by(metabolite) %>%
  do(tidy(aov(value ~ primary_tumor + dose, data = .)))
## Then filter and sort for just the significant metabolites by radiation dose:
normal_anova_dose_significant <- normal_anova %>% 
  filter(term == "dose" & p.value <= 0.05) %>%
  arrange(p.value) %>%
  mutate(normality = 'normal')

# Log-transformed metabolites:
## First pivot longer to reformat the output columns for anova function
logtransform_data_long <- vs_xeno_fc_log_transformed %>%
  pivot_longer(cols = 3:ncol(.), names_to = "metabolite", values_to = "value")
# Then run the two-way ANOVA by the first two columns  across all metabolites
####NOTE: primary_tumor + dose does not include the interaction term between primary_tumor and dose (my_anova). If desire the interaction, change to primary_tumor * dose inside aov function
logtransform_anova <- logtransform_data_long %>%
  group_by(metabolite) %>%
  do(tidy(aov(value ~ primary_tumor + dose, data = .)))
## Then filter and sort for just the significant metabolites by radiation dose:
logtransform_anova_dose_significant <- logtransform_anova %>% 
  filter(term == "dose" & p.value <= 0.05) %>%
  arrange(p.value) %>% 
  mutate(normality = 'lognormal')

all_anova_dose_significant <- rbind(normal_anova_dose_significant, logtransform_anova_dose_significant) %>% 
  arrange(p.value)
view(all_anova_dose_significant)

significant_metabolites <- all_anova_dose_significant$metabolite
significant_metabolites_normal <- normal_anova_dose_significant$metabolite
significant_metabolites_logtransform <- logtransform_anova_dose_significant$metabolite


### APPENDIX: View list of significant metabolites based on Two-Way ANOVA (no post-hoc test yet)
view(significant_metabolites)


******STILL NEED TO ADD/FIGURE OUT POST HOC ANALYSIS (holm sidak)

******************

spermidine_anova <- anova(lm(vs_xeno_fc$Spermidine ~ dose + primary_tumor, data = vs_xeno_fc))
## DISCOVERY AS OF 5.3.23 AT 1637: spermidine_anova shows spermidine P values significant, but then when do posthoc pairwise tests with Holm correction 
### there is actually no pairwise group that is significant. In contrast, if there is no correction (p.adjust.method="none") then the 0-10 and 0-20 are 
### significant. Thus, the anova is only significant prior to correcting for multiple comparisons.
##### CONCLUSION: need to figure out how to systematically test this of all the variables that were significant on plain two-way ANOVA 
spermidine_posthoc_dose <- pairwise.t.test(vs_xeno_fc$Spermidine, vs_xeno_fc$dose, p.adjust.method = "holm", pool.sd = FALSE)
spermidine_posthoc_primary_tumor <- pairwise.t.test(vs_xeno_fc$Spermidine, vs_xeno_fc$primary_tumor, p.adjust.method = "holm", pool.sd = FALSE)


spermidine_posthoc_dose <- pairwise.t.test(vs_xeno_fc$Spermidine, vs_xeno_fc$dose, p.adjust.method = "none", pool.sd = FALSE)
















************
################ Graphs of metabolites that are significantly correlated with radiation dose (limit to ~top 20 candidates -- select those that are significant based on ANOVA) 
  ###########################INCOMPLETE AS OF 5.9.23 - ONLY HAVE A FEW GRAPHS DONE
# beta-Hydroxybutyrate fold change boxplot with color labels
ggplot(vs_xeno_fc, aes(x = `dose`, y = `beta-Hydroxybutyrate (3-Hydroxybutyrate)`)) +
  geom_boxplot(aes(group=`dose`)) +
  geom_point(aes(color = `primary_tumor`)) + 
  labs(title = "beta-Hydroxybutyrate", x = "RT Dose (Gy)", y = "Fold Change")

# Lysine fold change boxplot with color labels
ggplot(vs_xeno_fc, aes(x = `dose`, y = `Lysine`)) +
  geom_boxplot(aes(group=`dose`)) +
  geom_point(aes(color = `primary_tumor`)) + 
  labs(title = "Lysine", x = "RT Dose (Gy)", y = "Fold Change")

# Fold change boxplot for Spermidine - with color labels
ggplot(vs_xeno_fc, aes(x = `dose`, y = `Spermidine`)) +
  geom_boxplot(aes(group=`dose`)) +
  geom_point(aes(color = `primary_tumor`)) + 
  labs(title = "Spermidine", x = "RT Dose (Gy)", y = "Fold Change")



####################NOT YET SURE HOW TO AUTOMATE THE PRODUCTION OF MULTIPLE GRAPHS (perhaps not worth it)
#### Create a loop to make Fold Change Boxplots for all metabolites that are significant based on normal_anova_dose_significant and logtransform_anova_dose_significant
significant_metabolites_normal <- normal_anova_dose_significant$metabolite
significant_metabolites_logtransform <- logtransform_anova_dose_significant$metabolite


***THIS FOR LOOP DOESNT WORK YET: - perhaps should try using 'across' instead??
for (i in 1:length(significant_metabolites_normal)) {
    ggplot(vs_xeno_fc_normal, aes(x = `dose`, y = names(vs_xeno_fc_normal) == significant_metabolites_normal[i])) +
    geom_boxplot(aes(group=`dose`)) +
    geom_point(aes(color = `primary_tumor`)) + 
    labs(title = significant_metabolites_normal[i], x = "RT Dose (Gy)", y = "Fold Change")
}

ggplot(vs_xeno_fc_normal, aes(x = `dose`, y = names(vs_xeno_fc_normal) == significant_metabolites_normal[3])) +
  geom_boxplot(aes(group=`dose`)) +
  geom_point(aes(color = `primary_tumor`)) + 
  labs(title = significant_metabolites_normal[3], x = "RT Dose (Gy)", y = "Fold Change")

*************


### THE FOLLOWING ARE EXAMPLE PLOTS BUT NOT YET THE FINAL GOAL PRODUCT.
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

