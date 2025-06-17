library(Seurat)
library(Matrix)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
# #------------------------------read blood, paclitaxel plus atezolizumab-------------------------------------
# folder_names <- c("Post_P001_b","Post_P005_b","Post_P0012_b","Post_P0014_b","Post_P0017_b","Post_P0019_b",
#                   "Pre_P001_b","Pre_P002_b","Pre_P004_b","Pre_P005_b","Pre_P007_b","Pre_P0010_b",
#                   "Pre_P0012_b","Pre_P0014_b","Pre_P0016_b","Pre_P0017_b","Pre_P0019_b",
#                   "Prog_P002_b","Prog_P004_b","Prog_P007_b","Prog_P0010_b","Prog_P0016_b")  
#------------------------------read blood, paclitaxel-------------------------------------
# folder_names <- c("Post_P0011_b","Post_P0018_b","Post_P0020_b","Post_P0022_b", "Post_P0023_b","Post_P0024_b","Post_P0025_b",
#                   "Pre_P008_b","Pre_P0011_b","Pre_P0013_b","Pre_P0018_b","Pre_P0020_b","Pre_P0022_b","Pre_P0023_b","Pre_P0024_b","Pre_P0025_b",
#                   "Prog_P008_b","Prog_P0013_b")  
#------------------------------read tumor, paclitaxel plus atezolizumab-------------------------------------
# folder_names <- c("Post_P002_t","Pre_P004_t","Post_P005_t","Post_P012_t","Post_P016_t","Post_P017_t","Post_P019_t",
#                   "Pre_P002_t","Pre_P005_t","Pre_P007_t","Pre_P012_t","Pre_P016_t","Pre_P017_t","Pre_P019_t",
#                   "Prog_P007_t")  
# #------------------------------read tumor, paclitaxel-------------------------------------
 folder_names <- c("Post_P003_t","Post_P018_t","Post_P020_t","Post_P022_t","Post_P023_t","Post_P025_t",
                   "Pre_P013_t","Pre_P018_t","Pre_P020_t","Pre_P022_t","Pre_P023_t","Pre_P025_t",
                   "Prog_P013_t")  
#------------------------------------------------------------------------------------
seurat_list <- list()
for (i in seq_along(folder_names)) {
  folder_name <- folder_names[i]
  sample_path <- file.path("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient\\Tumor-Patient-Paclitaxel-Analysis", folder_name)
  
  data <- Read10X(data.dir = sample_path)
  seurat_object <- CreateSeuratObject(counts = data, min.features = 100, project = folder_name)
  seurat_object$orig.ident <- folder_name
  seurat_object$sample <- tolower(strsplit(folder_name, "_")[[1]][1])
  seurat_list[[folder_name]] <- seurat_object
}
#---------merge
merged_samples <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),   # Take the folder names
  project = "Tumor_samples"
)

samples <- folder_names
post_samples <- samples[grepl("Post|Prog", samples)]
trait <- ifelse(grepl("Prog", post_samples), "2", "1")

rm(seurat_list)
rm(seurat_object)
rm(data)
#---------pseudo bulk
raw_counts <- GetAssayData(merged_samples, layer = "counts")
head(raw_counts)
dim(raw_counts)            # gene number* cell number
cell_counts <- table(merged_samples$orig.ident) #number of cell
cell_counts
pseudo_bulk <- t(rowsum(t(as.matrix(raw_counts)), group = merged_samples$orig.ident))
head(pseudo_bulk)
dim(pseudo_bulk)   # gene number* Sample number
mean_counts <- sweep(pseudo_bulk, 2, cell_counts, FUN = "/")
head(mean_counts)
dim(mean_counts)
pseudo_bulk<-mean_counts
rm(mean_counts)
rm(raw_counts)
rm(combined)
#----------split samples based on pre and post for construct two wgcna network
sample_names <- colnames(pseudo_bulk)

sample_group <- ifelse(grepl("^Pre", sample_names), "pre", "post")
table(sample_group)

pre_samples <-  colnames(pseudo_bulk)[sample_group == "pre"]
post_samples <- colnames(pseudo_bulk)[sample_group == "post"]

datExpr_pre <- as.data.frame(t(pseudo_bulk[, pre_samples]))
datExpr_post <- as.data.frame(t(pseudo_bulk[, post_samples]))
#----------
var_pre  <- apply(datExpr_pre, 2, var)
var_post <- apply(datExpr_post, 2, var)
common_zero_var_genes <- names(var_pre)[var_pre == 0 & var_post == 0]
datExpr_pre  <- datExpr_pre[, !colnames(datExpr_pre) %in% common_zero_var_genes]
View(datExpr_pre)
datExpr_post <- datExpr_post[, !colnames(datExpr_post) %in% common_zero_var_genes]
#---------------------------
data1 <- read.csv("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient\\Tumor-Patient-Paclitaxel-Analysis\\TOP10-Clusters-final-Tumor-Pacli.csv")
#data2 <- read.csv("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\TOP10-Clusters-final-Blood-Pacli.csv")
output_dir <- "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient\\Tumor-Patient-Paclitaxel-Analysis\\Results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)  # ایجاد مسیر در صورت عدم وجود
gene_info <- data1[, c("gene")]
head(gene_info)
genes_to_use <- intersect(gene_info, colnames(datExpr_pre))
datExpr_pre1 <- datExpr_pre[, genes_to_use]

num_modules <- ncol(MEs_pr)
head(MEs_pr)
module_names <- colnames(MEs_pr)
module_names

# Common Gene Selection
genes_to_use <- intersect(gene_info, colnames(datExpr_pre))
datExpr_pre1 <- datExpr_pre[, genes_to_use]

# Eigengenes Genes
common_samples <- intersect(rownames(MEs_pr), rownames(datExpr_pre1))
MEs_mat <- MEs_pr[common_samples, ]
expr_mat <- datExpr_pre1[common_samples, ]

# module membership (kME)
module_membership <- cor(expr_mat, MEs_mat, use = "pairwise.complete.obs")

# Output
write.csv(module_membership, file = "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient\\Tumor-Patient-Paclitaxel-Analysis\\Results\\ModuleMembership_NewData.csv")
weighted_cor_matrix <- module_membership * module_importance  # ضرب سطری به صورت برداری در R
weighted_cor_matrix
write.csv(weighted_cor_matrix, file = "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient\\Tumor-Patient-Paclitaxel-Analysis\\Results\\Weighted_ModuleCorrelation_Pre_vs_Post.csv")

max_membership_per_gene <- apply(weighted_cor_matrix, 1, max)
max_membership_df <- data.frame(
  gene = rownames(weighted_cor_matrix),
  Max_Membership = max_membership_per_gene
)

write.csv(max_membership_df, file = "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient\\Tumor-Patient-Paclitaxel-Analysis\\Results\\Max_Module_Membership_Per_Gene.csv", row.names = FALSE)


