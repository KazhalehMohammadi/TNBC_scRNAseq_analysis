library(Seurat)
library(Matrix)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ------------------------------ Read blood, paclitaxel plus atezolizumab -------------------------------------
# The folder names below refer to samples taken from patients who received treatment with paclitaxel plus atezolizumab
folder_names <- c("Post_P001_b","Post_P005_b","Post_P0012_b","Post_P0014_b","Post_P0017_b","Post_P0019_b",
                  "Pre_P001_b","Pre_P002_b","Pre_P004_b","Pre_P005_b","Pre_P007_b","Pre_P0010_b",
                  "Pre_P0012_b","Pre_P0014_b","Pre_P0016_b","Pre_P0017_b","Pre_P0019_b",
                  "Prog_P002_b","Prog_P004_b","Prog_P007_b","Prog_P0010_b","Prog_P0016_b")  

# ------------------------------ Read blood, paclitaxel -------------------------------------
# Alternative blood sample sets for paclitaxel treatment (this section is commented out)
# folder_names <- c("Post_P0011_b","Post_P0018_b","Post_P0020_b","Post_P0022_b", "Post_P0023_b","Post_P0024_b","Post_P0025_b",
#                   "Pre_P008_b","Pre_P0011_b","Pre_P0013_b","Pre_P0018_b","Pre_P0020_b","Pre_P0022_b","Pre_P0023_b","Pre_P0024_b","Pre_P0025_b",
#                   "Prog_P008_b","Prog_P0013_b")  

# ------------------------------ Read tumor, paclitaxel plus atezolizumab -------------------------------------
# Tumor sample sets for paclitaxel plus atezolizumab treatment (this section is commented out)
# folder_names <- c("Post_P002_t","Pre_P004_t","Post_P005_t","Post_P012_t","Post_P016_t","Post_P017_t","Post_P019_t",
#                   "Pre_P002_t","Pre_P005_t","Pre_P007_t","Pre_P012_t","Pre_P016_t","Pre_P017_t","Pre_P019_t",
#                   "Prog_P007_t")  

# ------------------------------ Read tumor, paclitaxel -------------------------------------
folder_names <- c("Post_P003_t","Post_P018_t","Post_P020_t","Post_P022_t","Post_P023_t","Post_P025_t",
                  "Pre_P013_t","Pre_P018_t","Pre_P020_t","Pre_P022_t","Pre_P023_t","Pre_P025_t",
                  "Prog_P013_t")  
# ------------------------------------------------------------------------------------
seurat_list <- list()
for (i in seq_along(folder_names)) {
  folder_name <- folder_names[i]
  sample_path <- file.path("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient\\Tumor-Patient-Paclitaxel-Analysis", folder_name)
  
  # Read data using the 10X Genomics format
  data <- Read10X(data.dir = sample_path)
  seurat_object <- CreateSeuratObject(counts = data, min.features = 100, project = folder_name)
  seurat_object$orig.ident <- folder_name  # Add original identifier for each sample
  seurat_object$sample <- tolower(strsplit(folder_name, "_")[[1]][1])  # Create a "sample" attribute from the folder name
  seurat_list[[folder_name]] <- seurat_object  # Store each Seurat object in the list
}

# --------- Merge all Seurat objects into one single object
merged_samples <- merge(
  x = seurat_list[[1]],  # Start with the first Seurat object
  y = seurat_list[-1],   # Merge the rest of the Seurat objects in the list
  add.cell.ids = names(seurat_list),   # Add the folder names as unique cell IDs
  project = "Tumor_samples"  # Set project name
)

# Prepare a list of samples, with "Post" or "Prog" samples categorized separately for later use
samples <- folder_names
post_samples <- samples[grepl("Post|Prog", samples)]
trait <- ifelse(grepl("Prog", post_samples), "2", "1")  # Assign a "2" for progressive samples, and "1" for post samples

# Clean up to save memory
rm(seurat_list)
rm(seurat_object)
rm(data)

# --------- Pseudo bulk transformation: summing counts by sample (Pre/Post) and normalizing
raw_counts <- GetAssayData(merged_samples, layer = "counts")  # Extract count data from the Seurat object
head(raw_counts)
dim(raw_counts)            # Check dimensions (genes x cells)
cell_counts <- table(merged_samples$orig.ident)  # Count cells per sample
cell_counts
pseudo_bulk <- t(rowsum(t(as.matrix(raw_counts)), group = merged_samples$orig.ident))  # Sum by group (sample)
head(pseudo_bulk)
dim(pseudo_bulk)   # Check dimensions (genes x sample groups)
mean_counts <- sweep(pseudo_bulk, 2, cell_counts, FUN = "/")  # Normalize by the number of cells in each group
head(mean_counts)
dim(mean_counts)
pseudo_bulk <- mean_counts  # Update pseudo_bulk with the normalized counts
rm(mean_counts)
rm(raw_counts)
rm(combined)

# --------- Split samples into pre-treatment and post-treatment for constructing two WGCNA networks
sample_names <- colnames(pseudo_bulk)

sample_group <- ifelse(grepl("^Pre", sample_names), "pre", "post")  # Categorize samples into "pre" and "post"
table(sample_group)

pre_samples <- colnames(pseudo_bulk)[sample_group == "pre"]  # Select pre-treatment samples
post_samples <- colnames(pseudo_bulk)[sample_group == "post"]  # Select post-treatment samples

datExpr_pre <- as.data.frame(t(pseudo_bulk[, pre_samples]))  # Create expression data frame for pre-treatment
datExpr_post <- as.data.frame(t(pseudo_bulk[, post_samples]))  # Create expression data frame for post-treatment

# ---------- Variance filtering: Remove genes with zero variance across both pre and post groups
var_pre <- apply(datExpr_pre, 2, var)
var_post <- apply(datExpr_post, 2, var)
common_zero_var_genes <- names(var_pre)[var_pre == 0 & var_post == 0]  # Identify genes with zero variance in both groups
datExpr_pre <- datExpr_pre[, !colnames(datExpr_pre) %in% common_zero_var_genes]  # Remove these genes
View(datExpr_pre)  # Check the filtered data
datExpr_post <- datExpr_post[, !colnames(datExpr_post) %in% common_zero_var_genes]  # Apply the same filtering to post-treatment data

# ---------------------------
data1 <- read.csv("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient\\Tumor-Patient-Paclitaxel-Analysis\\TOP10-Clusters-final-Tumor-Pacli.csv")  # Read gene information
# data2 <- read.csv("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\TOP10-Clusters-final-Blood-Pacli.csv")  # Uncomment if needed
output_dir <- "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient\\Tumor-Patient-Paclitaxel-Analysis\\Results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)  # Create the output directory if it doesn't exist
gene_info <- data1[, c("gene")]  # Extract gene names from the file
head(gene_info)
genes_to_use <- intersect(gene_info, colnames(datExpr_pre))  # Find genes that are present in both the gene list and the pre-treatment data
datExpr_pre1 <- datExpr_pre[, genes_to_use]  # Subset pre-treatment data to only these genes

# Get module information (Eigengenes)
num_modules <- ncol(MEs_pr)  # Number of modules in the pre-treatment network
head(MEs_pr)  # Display a preview of module eigengenes
module_names <- colnames(MEs_pr)  # List of module names
module_names

# Common Gene Selection
genes_to_use <- intersect(gene_info, colnames(datExpr_pre))  # Select genes that are common to the gene list and data
datExpr_pre1 <- datExpr_pre[, genes_to_use]  # Subset data to only the selected genes

# Eigengenes Genes: Find common samples between eigengenes and the selected genes
common_samples <- intersect(rownames(MEs_pr), rownames(datExpr_pre1))
MEs_mat <- MEs_pr[common_samples, ]  # Subset eigengenes
expr_mat <- datExpr_pre1[common_samples, ]  # Subset gene expression matrix

# Calculate module membership (kME) by correlating expression matrix with eigengene matrix
module_membership <- cor(expr_mat, MEs_mat, use = "pairwise.complete.obs")  # Correlation between genes and modules

# Output module membership to CSV
write.csv(module_membership, file = "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient\\Tumor-Patient-Paclitaxel-Analysis\\Results\\ModuleMembership_NewData.csv")

# Calculate weighted correlation matrix by multiplying module membership with module importance
weighted_cor_matrix <- module_membership * module_importance  # Element-wise multiplication
weighted_cor_matrix  # Display the weighted correlation matrix
write.csv(weighted_cor_matrix, file = "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient\\Tumor-Patient-Paclitaxel-Analysis\\Results\\Weighted_ModuleCorrelation_Pre_vs_Post.csv.csv")

# Find the maximum module membership for each gene
max_membership_per_gene <- apply(weighted_cor_matrix, 1, max)
max_membership_df <- data.frame(
  gene = rownames(weighted_cor_matrix),
  Max_Membership = max_membership_per_gene  # Maximum membership value for each gene
)

write.csv(max_membership_df, file = "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient\\Tumor-Patient-Paclitaxel-Analysis\\Results\\Max_Module_Membership_Per_Gene.csv", row.names = FALSE)
