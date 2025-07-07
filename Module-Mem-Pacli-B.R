library(Seurat)
library(Matrix)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ------------------------------ Read blood, paclitaxel plus atezolizumab -------------------------------------
# List of sample names for blood samples treated with paclitaxel plus atezolizumab
folder_names <- c("Post_P011_b","Post_P018_b","Post_P020_b","Post_P022_b", "Post_P023_b","Post_P024_b","Post_P025_b",
                   "Pre_P008_b","Pre_P011_b","Pre_P013_b","Pre_P018_b","Pre_P020_b","Pre_P022_b","Pre_P023_b","Pre_P024_b","Pre_P025_b",
                   "Prog_P008_b","Prog_P013_b")  

# ------------------------------ Read tumor, paclitaxel plus atezolizumab -------------------------------------
# The following sample set (commented) refers to tumor samples treated with paclitaxel plus atezolizumab.
# folder_names <- c("Post_P002_t","Pre_P004_t","Post_P005_t","Post_P012_t","Post_P016_t","Post_P017_t","Post_P019_t",
#                   "Pre_P002_t","Pre_P005_t","Pre_P007_t","Pre_P012_t","Pre_P016_t","Pre_P017_t","Pre_P019_t",
#                   "Prog_P007_t")  

# ------------------------------ Read tumor, paclitaxel -------------------------------------
# A different set of tumor samples treated with paclitaxel (this is being used in this case)
folder_names <- c("Post_P003_t","Post_P018_t","Post_P020_t","Post_P022_t","Post_P023_t","Post_P025_t",
                   "Pre_P013_t","Pre_P018_t","Pre_P020_t","Pre_P022_t","Pre_P023_t","Pre_P025_t",
                   "Prog_P013_t")  

# ------------------------------------------------------------------------------------
seurat_list <- list()
for (i in seq_along(folder_names)) {
  folder_name <- folder_names[i]
  sample_path <- file.path("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Blood-Patient-Paclitaxel-Analysis", folder_name)
  
  # Reading the data using the 10X Genomics format
  data <- Read10X(data.dir = sample_path)
  seurat_object <- CreateSeuratObject(counts = data, min.features = 100, project = folder_name)
  
  # Add metadata to identify the sample origin and its condition (pre, post, or progressive)
  seurat_object$orig.ident <- folder_name
  seurat_object$sample <- tolower(strsplit(folder_name, "_")[[1]][1])
  
  # Store Seurat object in a list
  seurat_list[[folder_name]] <- seurat_object
}

# --------- Merge the Seurat objects from different samples
merged_samples <- merge(
  x = seurat_list[[1]],  # Start with the first Seurat object
  y = seurat_list[-1],   # Merge the rest of the Seurat objects in the list
  add.cell.ids = names(seurat_list),  # Use the folder names as unique cell IDs
  project = "Tumor_samples"  # Assign a project name to the merged object
)

# Create a list of sample identifiers, and assign trait values ("1" for post and "2" for progressive samples)
samples <- folder_names
post_samples <- samples[grepl("Post|Prog", samples)]
trait <- ifelse(grepl("Prog", post_samples), "2", "1")

# Remove unnecessary variables to save memory
rm(seurat_list)
rm(seurat_object)
rm(data)

# --------- Pseudo bulk transformation: Aggregating gene counts across cells by sample
raw_counts <- GetAssayData(merged_samples, layer = "counts")  # Extract the raw count data
head(raw_counts)  # Preview the data
dim(raw_counts)   # Display the dimensions (genes x cells)

# Count the number of cells per sample
cell_counts <- table(merged_samples$orig.ident)
cell_counts

# Aggregate the counts across all cells in each sample (pre/post-treatment groups)
pseudo_bulk <- t(rowsum(t(as.matrix(raw_counts)), group = merged_samples$orig.ident))
head(pseudo_bulk)
dim(pseudo_bulk)  # Check the new dimensions (genes x samples)

# Normalize by the number of cells in each sample
mean_counts <- sweep(pseudo_bulk, 2, cell_counts, FUN = "/")
head(mean_counts)
dim(mean_counts)

# Update pseudo_bulk with the normalized counts
pseudo_bulk <- mean_counts
rm(mean_counts)
rm(raw_counts)
rm(combined)

# --------- Split samples into pre-treatment and post-treatment groups for constructing two WGCNA networks
sample_names <- colnames(pseudo_bulk)
sample_group <- ifelse(grepl("^Pre", sample_names), "pre", "post")  # Categorize samples as "pre" or "post"
table(sample_group)

pre_samples <- colnames(pseudo_bulk)[sample_group == "pre"]
post_samples <- colnames(pseudo_bulk)[sample_group == "post"]

# Create expression matrices for pre and post treatment
datExpr_pre <- as.data.frame(t(pseudo_bulk[, pre_samples]))
datExpr_post <- as.data.frame(t(pseudo_bulk[, post_samples]))

# ---------- Variance filtering: Remove genes with zero variance across both pre and post groups
var_pre  <- apply(datExpr_pre, 2, var)
var_post <- apply(datExpr_post, 2, var)

# Identify genes that have zero variance in both pre and post groups
common_zero_var_genes <- names(var_pre)[var_pre == 0 & var_post == 0]
datExpr_pre  <- datExpr_pre[, !colnames(datExpr_pre) %in% common_zero_var_genes]  # Remove those genes from pre-treatment data
View(datExpr_pre)  # Check the filtered data
datExpr_post <- datExpr_post[, !colnames(datExpr_post) %in% common_zero_var_genes]  # Apply the same filtering to post-treatment data

# ---------------------------
data1 <- read.csv("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Blood-Patient-Paclitaxel-Analysis\\TOP10-Clusters-final-Blood-Pacli.csv")  # Read gene information from CSV
# data2 <- read.csv("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\TOP10-Clusters-final-Blood-Pacli.csv")  # Uncomment if necessary

# Create an output directory if it doesn't already exist
output_dir <- "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Blood-Patient-Paclitaxel-Analysis\\Results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Extract gene information (gene names)
gene_info <- data1[, c("gene")]
head(gene_info)

# Select genes that are present in both the gene information list and the pre-treatment data
genes_to_use <- intersect(gene_info, colnames(datExpr_pre))
datExpr_pre1 <- datExpr_pre[, genes_to_use]

# Get the number of modules in the pre-treatment network
num_modules <- ncol(MEs_pr)
head(MEs_pr)
module_names <- colnames(MEs_pr)  # Get the names of the modules
module_names

# Common Gene Selection
genes_to_use <- intersect(gene_info, colnames(datExpr_pre))  # Select genes that are common in both the data and gene info
datExpr_pre1 <- datExpr_pre[, genes_to_use]  # Subset the data to include only those genes

# Eigen Gene Calculation: Matching samples between eigengenes and selected genes
common_samples <- intersect(rownames(MEs_pr), rownames(datExpr_pre1))
MEs_mat <- MEs_pr[common_samples, ]  # Subset the eigengene matrix
expr_mat <- datExpr_pre1[common_samples, ]  # Subset the gene expression matrix

length(common_samples)  # Check the number of common samples
head(rownames(MEs_pr))  # Preview row names (samples) of eigengenes
head(rownames(datExpr_pre1))  # Preview row names (samples) of the expression matrix

# Module membership (kME): Correlation between gene expression and module eigengenes
module_membership <- cor(expr_mat, MEs_mat, use = "pairwise.complete.obs")

# Save module membership to a CSV file
write.csv(module_membership, file = "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Blood-Patient-Paclitaxel-Analysis\\Results\\ModuleMembership_NewData.csv")

# Weighted correlation matrix: Multiply module membership by module importance
weighted_cor_matrix <- module_membership * module_importance  

# Save the weighted correlation matrix to a CSV file
write.csv(weighted_cor_matrix, file = "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Blood-Patient-Paclitaxel-Analysis\\Results\\Weighted_ModuleCorrelation_Pre_vs_Post.csv")

# Find the maximum membership for each gene and save to CSV
max_membership_per_gene <- apply(weighted_cor_matrix, 1, max)

max_membership_df <- data.frame(
  gene = rownames(weighted_cor_matrix),
  Max_Membership = max_membership_per_gene
)

# Save the maximum module membership per gene to a CSV file
write.csv(max_membership_df, file = "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Blood-Patient-Paclitaxel-Analysis\\Results\\Max_Module_Membership_Per_Gene.csv", row.names = FALSE)
