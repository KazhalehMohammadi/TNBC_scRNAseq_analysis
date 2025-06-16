library(Seurat)
library(Matrix)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
# #------------------------------read blood, paclitaxel plus atezolizumab-------------------------------------
# folder_names <- c("Post_P001_b","Post_P005_b","Post_P012_b","Post_P014_b","Post_P017_b","Post_P019_b",
#                   "Pre_P001_b","Pre_P002_b","Pre_P004_b","Pre_P005_b","Pre_P007_b","Pre_P010_b",
#                   "Pre_P012_b","Pre_P014_b","Pre_P016_b","Pre_P017_b","Pre_P019_b",
#                   "Prog_P002_b","Prog_P004_b","Prog_P007_b","Prog_P010_b","Prog_P016_b")  
# #------------------------------read blood, paclitaxel-------------------------------------
# folder_names <- c("Post_P011_b","Post_P018_b","Post_P020_b","Post_P022_b", "Post_P023_b","Post_P024_b","Post_P025_b",
#                   "Pre_P008_b","Pre_P011_b","Pre_P013_b","Pre_P018_b","Pre_P020_b","Pre_P022_b","Pre_P023_b","Pre_P024_b","Pre_P025_b",
#                   "Prog_P008_b","Prog_P013_b")  
# #------------------------------read tumor, paclitaxel plus atezolizumab-------------------------------------
# folder_names <- c("Post_P002_t","Post_P005_t","Post_P012_t","Post_P016_t","Post_P017_t","Post_P019_t",
#                   "Pre_P002_t","Pre_P004_t","Pre_P005_t","Pre_P007_t","Pre_P012_t","Pre_P016_t","Pre_P017_t","Pre_P019_t",
#                   "Prog_P007_t")  
#------------------------------read tumor, paclitaxel-------------------------------------
folder_names <- c("Post_GSM7073191_t","Post_P020_t","Post_P022_t","Post_P025_t",
                  "Pre_P013_t","Pre_P020_t","Pre_P022_t","Pre_P025_t",
                  "Prog_P013_t","Post_P003_t","Pre_P018_t","Post_P018_t","Pre_P023_t","Post_P023_t")  
#------------------------------------------------------------------------------------
seurat_list <- list()
for (i in seq_along(folder_names)) {
  folder_name <- folder_names[i]
  sample_path <- file.path("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient-Paclitaxel-Analysis", folder_name)
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
  project = "Blood_samples"
)

samples <- folder_names
post_samples <- samples[grepl("Post|Prog", samples)]
trait <- ifelse(grepl("Prog", post_samples), "2", "1")

rm(seurat_list)
rm(seurat_object)
rm(data)
#---------pseudo bulk
raw_counts <- GetAssayData(merged_samples, layer = "counts")
dim(raw_counts)            # gene number* cell number
cell_counts <- table(merged_samples$orig.ident) #number of cell 
pseudo_bulk <- t(rowsum(t(as.matrix(raw_counts)), group = merged_samples$orig.ident))
dim(pseudo_bulk)
mean_counts <- sweep(pseudo_bulk, 2, cell_counts, FUN = "/")
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
datExpr_post <- datExpr_post[, !colnames(datExpr_post) %in% common_zero_var_genes]
#------------------------------Soft-threshold Plot

par(mfrow = c(1,2))  

# برای pre
plot(sft_pre$fitIndices[,1], -sign(sft_pre$fitIndices[,3])*sft_pre$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R²", 
     type="n", main = "Scale independence (pre)")
text(sft_pre$fitIndices[,1], -sign(sft_pre$fitIndices[,3])*sft_pre$fitIndices[,2],
     labels=powers, cex=0.7, col="red")
abline(h=0.9, col="red")  # خط مرجع برای R² = 0.9

# برای pre
plot(sft_pre$fitIndices[,1], sft_pre$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = "Mean connectivity (pre)")
text(sft_pre$fitIndices[,1], sft_pre$fitIndices[,5], labels=powers, cex=0.7, col="red")

#-------------------------------------------------------------------------------
powers = c(1:20, seq(22, 100, by=2))
sft_pre  = pickSoftThreshold(datExpr_pre, powerVector = powers, verbose = 5)
sft_post = pickSoftThreshold(datExpr_post, powerVector = powers, verbose = 5)
softPower_pre <- 24  
softPower_post <- 22 
# ============================
# --- Pre
adjacency_pre <- adjacency(datExpr_pre, power = softPower_pre)
TOM_pre <- TOMsimilarity(adjacency_pre)
dissTOM_pre <- 1 - TOM_pre
geneTree_pre <- hclust(as.dist(dissTOM_pre), method = "average")
# --- Post
adjacency_post <- adjacency(datExpr_post, power = softPower_post)
TOM_post <- TOMsimilarity(adjacency_post)
dissTOM_post <- 1 - TOM_post
geneTree_post <- hclust(as.dist(dissTOM_post), method = "average")
# ============================
minModuleSize <- 30
dynamicMods_pre <- cutreeDynamic(dendro = geneTree_pre,distM = dissTOM_pre,deepSplit = 2,
                                 pamRespectsDendro = FALSE,minClusterSize = minModuleSize)
moduleColors_pre <- labels2colors(dynamicMods_pre)
png("dendrogram_pre.png", width = 800, height = 600)
plotDendroAndColors(geneTree_pre, moduleColors_pre,"Module colors (pre)", dendroLabels = FALSE, 
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()

dynamicMods_post <- cutreeDynamic(dendro = geneTree_post,distM = dissTOM_post,deepSplit = 2,
                                  pamRespectsDendro = FALSE,minClusterSize = minModuleSize)
moduleColors_post <- labels2colors(dynamicMods_post)
png("dendrogram_post.png", width = 800, height = 600)
plotDendroAndColors(geneTree_post, moduleColors_post,"Module colors (Post)",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
# ============================eigen gene
MEList_pre <- moduleEigengenes(datExpr_pre, colors = moduleColors_pre)
MEs_pr <- MEList_pre$eigengenes

MEList_post <- moduleEigengenes(datExpr_post, colors = moduleColors_post)
MEs_post <- MEList_post$eigengenes
# ============================
num_modules_pre <- ncol(MEs_pr)  # Number of modules in the Pre network
num_modules_post <- ncol(MEs_post)  # Number of modules in the Post network
num_samples <- nrow(MEs_pr)  # Number of samples
module_importance <- numeric(num_modules_pre)
for (i in 1:num_modules_pre ) {
  current_importance <- 0
  for (j in 1:num_modules_post ) {
    cor_ME_pre_post <- cor(MEs_pr[, i], MEs_post[, j], use = "pairwise.complete.obs")
    cor_ME_post_trait <- cor(MEs_post[, j], trait, use = "pairwise.complete.obs")
    current_importance <- current_importance + cor_ME_pre_post * cor_ME_post_trait
  }
  module_importance[i] <- current_importance
}

# Display the results
module_importance

threshold <- 10  
important_modules <- which(abs(module_importance) >  threshold)
important_modules
#=============================
my_data <- read.csv("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient-Paclitaxel-Analysis\\TOP10-Clusters-final.csv")
gene_info <- my_data[, c("gene", "Pre_avg_log2FC")]
num_modules <- ncol(MEs_pr)  
num_genes <- nrow(gene_info)  
correlation_results <- data.frame(
  gene = gene_info$gene,  
  Pre_expression = gene_info$Pre_avg_log2FC
)

for (i in 1:num_modules) {
  module_cor <- apply(gene_info[, "Pre_avg_log2FC", drop = FALSE], 1, function(gene_row) {
    cor(MEs_pr[, i], gene_row, use = "pairwise.complete.obs")
  })
  correlation_results[paste("Module_", i, "_correlation", sep = "")] <- module_cor
}
for (i in 1:num_modules) {
  correlation_results[paste("Module_", i, "_importance", sep = "")] <- module_importance[i]
}
head(correlation_results)
write.csv(correlation_results, "D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Tumor-Patient-Paclitaxel-Analysis\\Results\\correlation_results.csv", row.names = FALSE)



