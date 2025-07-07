library(Seurat)
library(Matrix)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ------------------------------ Read blood, paclitaxel plus atezolizumab -------------------------------------
# List of blood samples treated with paclitaxel plus atezolizumab
folder_names <- c("Post_P001_b","Post_P005_b","Post_P012_b","Post_P014_b","Post_P017_b","Post_P019_b",
                   "Pre_P001_b","Pre_P002_b","Pre_P004_b","Pre_P005_b","Pre_P007_b","Pre_P010_b",
                   "Pre_P012_b","Pre_P014_b","Pre_P016_b","Pre_P017_b","Pre_P019_b",
                   "Prog_P002_b","Prog_P004_b","Prog_P007_b","Prog_P010_b","Prog_P016_b")  

# ------------------------------ Read blood, paclitaxel -------------------------------------
# The following block (commented out) can be used to read different sets of blood samples
# folder_names <- c("Post_P011_b","Post_P018_b","Post_P020_b","Post_P022_b", "Post_P023_b","Post_P024_b","Post_P025_b",
#                   "Pre_P008_b","Pre_P011_b","Pre_P013_b","Pre_P018_b","Pre_P020_b","Pre_P022_b","Pre_P023_b","Pre_P024_b","Pre_P025_b",
#                   "Prog_P008_b","Prog_P013_b")  

# ------------------------------ Read tumor, paclitaxel plus atezolizumab -------------------------------------
# This is another set of tumor samples (commented out)
# folder_names <- c("Post_P002_t","Pre_P004_t","Post_P005_t","Post_P012_t","Post_P016_t","Post_P017_t","Post_P019_t",
#                   "Pre_P002_t","Pre_P005_t","Pre_P007_t","Pre_P012_t","Pre_P016_t","Pre_P017_t","Pre_P019_t",
#                   "Prog_P007_t")  

# ------------------------------ Read tumor, paclitaxel -------------------------------------
# In this case, tumor samples treated with paclitaxel are being used.
folder_names <- c("Post_P003_t","Post_P018_t","Post_P020_t","Post_P022_t","Post_P023_t","Post_P025_t",
                   "Pre_P013_t","Pre_P018_t","Pre_P020_t","Pre_P022_t","Pre_P023_t","Pre_P025_t",
                   "Prog_P013_t")  

# ------------------------------------------------------------------------------------
seurat_list <- list()  # Initialize an empty list to store Seurat objects

# Loop through each sample folder and process the data
for (i in seq_along(folder_names)) {
  folder_name <- folder_names[i]
  sample_path <- file.path("D:/Thesis/Proposal/Analysis/WGCNA/Blood-Patient/Blood-Patient-Paclitaxel&Atizo-Analysis", folder_name)
  
  # Read the data in 10X Genomics format
  data <- Read10X(data.dir = sample_path)
  seurat_object <- CreateSeuratObject(counts = data, min.features = 100, project = folder_name)
  
  # Add metadata for identification and conditions
  seurat_object$orig.ident <- folder_name
  seurat_object$sample <- tolower(strsplit(folder_name, "_")[[1]][1])  # Set condition based on sample name
  
  # Store the Seurat object in the list
  seurat_list[[folder_name]] <- seurat_object
}

# --------- Merge the Seurat objects from different samples into one large Seurat object
merged_samples <- merge(
  x = seurat_list[[1]],  # Start with the first Seurat object
  y = seurat_list[-1],   # Merge the rest of the Seurat objects from the list
  add.cell.ids = names(seurat_list),  # Use folder names as unique cell IDs
  project = "Blood_samples"  # Set the project name for the merged object
)

# Prepare a list of samples for trait assignment (post-treatment and progressive samples)
samples <- folder_names
post_samples <- samples[grepl("Post|Prog", samples)]  # Identify post-treatment and progressive samples
trait <- ifelse(grepl("Prog", post_samples), "2", "1")  # Assign trait value "2" for progressive, "1" for post-treatment

# Remove unnecessary variables to save memory
rm(seurat_list)
rm(seurat_object)
rm(data)

# --------- Pseudo bulk transformation: Aggregate gene counts across all cells in each sample
raw_counts <- GetAssayData(merged_samples, layer = "counts")  # Extract raw count data
head(raw_counts)  # Preview the data
dim(raw_counts)   # Check the dimensions (genes x cells)

# Count the number of cells in each sample
cell_counts <- table(merged_samples$orig.ident)
cell_counts

# Aggregate the gene counts by sample
pseudo_bulk <- t(rowsum(t(as.matrix(raw_counts)), group = merged_samples$orig.ident))
head(pseudo_bulk)  # Preview aggregated counts
dim(pseudo_bulk)   # New dimensions: genes x samples

# Normalize the counts by dividing by the number of cells in each sample
mean_counts <- sweep(pseudo_bulk, 2, cell_counts, FUN = "/")
head(mean_counts)
dim(mean_counts)

# Update pseudo_bulk with normalized counts
pseudo_bulk <- mean_counts
rm(mean_counts)
rm(raw_counts)
rm(combined)

# --------- Split samples based on pre-treatment and post-treatment for constructing two WGCNA networks
sample_names <- colnames(pseudo_bulk)
sample_group <- ifelse(grepl("^Pre", sample_names), "pre", "post")  # Categorize into "pre" or "post"
table(sample_group)

# Get the pre-treatment and post-treatment samples
pre_samples <- colnames(pseudo_bulk)[sample_group == "pre"]
post_samples <- colnames(pseudo_bulk)[sample_group == "post"]

# Create expression matrices for pre and post-treatment data
datExpr_pre <- as.data.frame(t(pseudo_bulk[, pre_samples]))
datExpr_post <- as.data.frame(t(pseudo_bulk[, post_samples]))

# ---------- Variance filtering: Remove genes with zero variance across both pre and post groups
var_pre  <- apply(datExpr_pre, 2, var)
var_post <- apply(datExpr_post, 2, var)

# Identify and remove genes with zero variance in both groups
common_zero_var_genes <- names(var_pre)[var_pre == 0 & var_post == 0]
datExpr_pre  <- datExpr_pre[, !colnames(datExpr_pre) %in% common_zero_var_genes]  # Remove from pre-treatment data
View(datExpr_pre)  # View the filtered data
datExpr_post <- datExpr_post[, !colnames(datExpr_post) %in% common_zero_var_genes]  # Same filtering for post-treatment data

# ---------------------------
# Read gene information for further analysis
data1 <- read.csv("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Blood-Patient-Paclitaxel&Atizo-Analysis\\TOP10-Clusters-final-Blood-Pacli-plus Azi.csv")

# Create an output directory if it doesn't exist
output_dir <- "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Blood-Patient-Paclitaxel&Atizo-Analysis\\Results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Extract gene information (gene names)
gene_info <- data1[, c("gene")]
head(gene_info)

# Select genes that are present in both the gene list and the pre-treatment data
genes_to_use <- intersect(gene_info, colnames(datExpr_pre))
datExpr_pre1 <- datExpr_pre[, genes_to_use]

# Get the number of modules in the pre-treatment network
num_modules <- ncol(MEs_pr)
head(MEs_pr)  # Preview the eigengenes
module_names <- colnames(MEs_pr)  # Get module names
module_names

# Select common genes
genes_to_use <- intersect(gene_info, colnames(datExpr_pre))  # Select genes common to both data
datExpr_pre1 <- datExpr_pre[, genes_to_use]

# Eigen Gene Calculation: Match samples between eigengenes and selected genes
common_samples <- intersect(rownames(MEs_pr), rownames(datExpr_pre1))
MEs_mat <- MEs_pr[common_samples, ]  # Subset the eigengene matrix
expr_mat <- datExpr_pre1[common_samples, ]  # Subset the gene expression matrix

# Check the number of common samples
length(common_samples)
head(rownames(MEs_pr))  # Preview eigengene sample names
head(rownames(datExpr_pre1))  # Preview expression matrix sample names

# Module Membership (kME): Correlate gene expression with module eigengenes
module_membership <- cor(expr_mat, MEs_mat, use = "pairwise.complete.obs")

# Save the module membership data to a CSV file
write.csv(module_membership, file = "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Blood-Patient-Paclitaxel&Atizo-Analysis\\Results\\ModuleMembership_NewData.csv")

# Calculate the weighted correlation matrix: Multiply module membership by module importance
weighted_cor_matrix <- module_membership * module_importance  

# Save the weighted correlation matrix to a CSV file
write.csv(weighted_cor_matrix, file = "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Blood-Patient-Paclitaxel&Atizo-Analysis\\Results\\Weighted_ModuleCorrelation_Pre_vs_Post.csv")

# Find the maximum membership for each gene and save to CSV
max_membership_per_gene <- apply(weighted_cor_matrix, 1, max)

# Create a data frame with genes and their maximum module membership
max_membership_df <- data.frame(
  gene = rownames(weighted_cor_matrix),
  Max_Membership = max_membership_per_gene
)

# Save the results to a CSV file
write.csv(max_membership_df, file = "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Blood-Patient-Paclitaxel&Atizo-Analysis\\Results\\Max_Module_Membership_Per_Gene.csv", row.names = FALSE)
