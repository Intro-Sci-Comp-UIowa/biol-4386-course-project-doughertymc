Folder “R Analysis” ("C:\Users\mark1\Dropbox\Research\Schwannoma & Meningioma Basic Research\Metabolomics Data\R Analysis") Created 1.25.2023 as first stab at organizing metabolomics data for analysis using R. I anticipate this will be especially helpful for Xenograft Radiation analysis that requires multiple steps of normalization (thus MetaboAnalyst cannot easily handle).
Specifically I want to take raw metabolomics data from Excel files
- Convert values to fold changes by dividing by the mean of the control samples derived from the same original primary tumor
- Then analyze each metabolite for outlier values with Grubb’s test) and produce list of those outliers with SampleID, metabolite, Mouse Group/Patient of origin, and Radiation Dose specified, creating a documented list of which values were removed AND a clean list without those outliers
- Then Shapiro-Wilk normality test, document results, and LogTransform values of metabolites that are LogNormal, producing documented list of normality for each metabolite and then new dataset with the appropriate transformations as needed
- Then analysis (list subject to change):
o ANOVA or T-test
o PCA
o K-means clustering or other unsupervised clustering
o Heatmap
o Correlational analysis
- Then produce graphs/figures
o Bar charts with individual values and error bars

NOTES FROM YOUTUBE VIDEO TUTORIAL
Packages: ggplot2, rio (R input output), rmarkdown; pacman (package manager)
