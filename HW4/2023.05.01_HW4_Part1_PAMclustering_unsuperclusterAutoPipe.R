# Install AutoPipe from GitHub using devtools (devtools must be loaded):
# First must install dependencies if not already installed
library(BiocManager)
BiocManager::install(c('siggenes', 'annotate', 'fgsea', 'org.Hs.eg.db', 'ConsensusClusterPlus', 'clusterProfiler'))
install_github("falafel19/AutoPipe")
library(AutoPipe)


## *********************NOTE TO SELF 4.28.23 10:57AM-CURRENTLY THE IMPORT INCLUDES SPINAL SCHWANNOMAS TOO - WILL WANT TO RUN WITH AND WITHOUT THOSE INCLUDED

# Import raw data file from CSV to tibble using read_csv
vs_primary_metabolomics_raw <- read_csv("C:/Users/mark1/Dropbox/BIOL_4386/Project_Folder/Formatted_Data/230117_Swnma_PrimaryTumorsGCLC_09.2022_Batch.csv")
# Remove missing metabolite columns and save as curated tibble. Note that columns 1-6 are metadata.
vs_primary_curated <- vs_primary_metabolomics_raw %>% select(!where(is_logical))
metabolite_names_primary <- colnames(vs_primary_curated)[-(1:6)]
# Note that must remove sample S15 to run clustering because S15 is missing all LCMS data (~half of metabolites), AND must remove S46 because it gets classified as a cluster on its own so the heatmap function throws an error with it included. Notably on prior MetaboAnalyst analysis, had found that S46 is an outlier on the PCA when including all metabolites due to crazy outlier value for Malonate. May be worthwhile to remove Malonate and re-test with S46 included, although not critical because only adds one more sample.
vs_primary_curated <- vs_primary_curated %>% filter(Sample_Label != "S15") %>% filter(Sample_Label != "S46")
view(vs_primary_curated)

vs_clinical_data <- vs_primary_curated[c(1,3,4,5,6)]
vs_clinical_data_df <- vs_clinical_data %>% column_to_rownames(var = "Sample_Label") %>% as.data.frame(.)

#Remove metabolite columns with missing data, and remove metadata columns 2:6
vs_primary_no_missing <- vs_primary_curated[-(2:6)] %>% select(where(~all(!is.na(.))))
view(vs_primary_no_missing)

#Must convert from tibble to dataframe for AutoPipe::TopPAM to work; this also converts the column Sample_Label to the actual row labels for the dataframe
vs_primary_df <- vs_primary_no_missing %>% column_to_rownames(var = "Sample_Label") %>% as.data.frame(.)
class(vs_primary_df)
vs_transposed <- t(vs_primary_df)
# Run AutoPipe's TopPAM feature to calcluate optimal number of clusters using PAM clustering. NOTE this requires metabolites in ROWS and samples in COLUMNS
res <- AutoPipe::TopPAM(vs_transposed, max_clusters = 15, TOP=139, B=100, clusterboot=FALSE)
# TopPAM result -> 2 groups are best for PAM clustering, but one of them is just sample S46

me_TOP <- res[[1]]
dim(me_TOP)
number_of_k <- res[[3]]
File_genes <- AutoPipe::Groups_Sup(me_TOP, me = vs_transposed, number_of_k,TRw=-1)
groups_men=File_genes[[2]]

AutoPipe::Supervised_Cluster_Heatmap(groups_men = groups_men, gene_matrix=File_genes[[1]], TOP_Cluster= 15, method="PAMR", show_sil=TRUE, show_clin=T, genes_to_print = 4, print_genes=T, TOP=50, GSE=F, plot_mean_sil=T, stats_clust=res[[2]], threshold=1, samples_data = vs_clinical_data_df)

