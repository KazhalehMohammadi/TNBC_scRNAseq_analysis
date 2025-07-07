library(Seurat)
library(Matrix)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ------------------------------ Read blood, paclitaxel plus atezolizumab -------------------------------------
folder_names <- c("Post_P001_b","Post_P005_b","Post_P012_b","Post_P014_b","Post_P017_b","Post_P019_b",
                  "Pre_P001_b","Pre_P002_b","Pre_P004_b","Pre_P005_b","Pre_P007_b","Pre_P010_b",
                  "Pre_P012_b","Pre_P014_b","Pre_P016_b","Pre_P017_b","Pre_P019_b",
                  "Prog_P002_b","Prog_P004_b","Prog_P007_b","Prog_P010_b","Prog_P016_b","Post_P011_b","Post_P018_b",
                  "Post_P020_b","Post_P022_b","Post_P023_b","Post_P024_b","Post_P025_b",                   
                  "Pre_P008_b","Pre_P011_b","Pre_P013_b","Pre_P018_b","Pre_P020_b","Pre_P022_b",
                  "Pre_P023_b","Pre_P024_b","Pre_P025_b",
                  "Prog_P008_b","Prog_P013_b")  

# ------------------------------ Read blood, paclitaxel -------------------------------------
# folder_names <- c("Post_P0011_b","Post_P0018_b","Post_P0020_b","Post_P0022_b", "Post_P0023_b","Post_P0024_b","Post_P0025_b",
#                   "Pre_P008_b","Pre_P0011_b","Pre_P0013_b","Pre_P0018_b","Pre_P0020_b","Pre_P0022_b","Pre_P0023_b","Pre_P0024_b","Pre_P0025_b",
#                   "Prog_P008_b","Prog_P0013_b")  

# ------------------------------ Read tumor, paclitaxel plus atezolizumab -------------------------------------
# folder_names <- c("Post_P002_t","Post_P005_t","Post_P0012_t","Post_P0016_t","Post_P0017_t","Post_P0019_t",
#                   "Pre_P002_t","Pre_P004_t","Pre_P005_t","Pre_P007_t","Pre_P0012_t","Pre_P0016_t","Pre_P0017_t","Pre_P0019_t",
#                   "Prog_P007_t")  

# ------------------------------ Read tumor, paclitaxel -------------------------------------
# folder_names <- c("Post_P003_t","Post_P0018_t","Post_P0020_t","Post_P0022_t","Post_P0023_t","Post_P0024_t","Post_P0025_t",
#                   "Pre_P0013_t","Pre_P0018_t","Pre_P0020_t","Pre_P0022_t","Pre_P0023_t","Pre_P0025_t",
#                   "Prog_P0013_t")  

# ------------------------------------------------------------------------------------
seurat_list <- list()
for (i in seq_along(folder_names)) {
  folder_name <- folder_names[i]
  sample_path <- file.path("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Important_modules_Blood", folder_name)
  data <- Read10X(data.dir = sample_path)
  seurat_object <- CreateSeuratObject(counts = data, min.features = 100, project = folder_name)
  seurat_object$orig.ident <- folder_name
  seurat_object$sample <- tolower(strsplit(folder_name, "_")[[1]][1])
  seurat_list[[folder_name]] <- seurat_object
}

# --------- Merge the samples
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

# --------- Pseudo bulk transformation
raw_counts <- GetAssayData(merged_samples, layer = "counts")
dim(raw_counts)            # gene number * cell number
cell_counts <- table(merged_samples$orig.ident) # number of cells 
pseudo_bulk <- t(rowsum(t(as.matrix(raw_counts)), group = merged_samples$orig.ident))
dim(pseudo_bulk)
mean_counts <- sweep(pseudo_bulk, 2, cell_counts, FUN = "/")
dim(mean_counts)
pseudo_bulk <- mean_counts
rm(mean_counts)
rm(raw_counts)
rm(combined)

# --------- Split samples based on pre and post for constructing two WGCNA networks
sample_names <- colnames(pseudo_bulk)

sample_group <- ifelse(grepl("^Pre", sample_names), "pre", "post")
table(sample_group)

pre_samples <-  colnames(pseudo_bulk)[sample_group == "pre"]
post_samples <- colnames(pseudo_bulk)[sample_group == "post"]

datExpr_pre <- as.data.frame(t(pseudo_bulk[, pre_samples]))
datExpr_post <- as.data.frame(t(pseudo_bulk[, post_samples]))

# ---------- Variance filtering
var_pre  <- apply(datExpr_pre, 2, var)
var_post <- apply(datExpr_post, 2, var)
common_zero_var_genes <- names(var_pre)[var_pre == 0 & var_post == 0]
datExpr_pre  <- datExpr_pre[, !colnames(datExpr_pre) %in% common_zero_var_genes]
datExpr_post <- datExpr_post[, !colnames(datExpr_post) %in% common_zero_var_genes]

#-------------------------------------------------------------------------------
powers = c(1:20, seq(22, 100, by=2))
sft_pre  = pickSoftThreshold(datExpr_pre, powerVector = powers, verbose = 5)
sft_post = pickSoftThreshold(datExpr_post, powerVector = powers, verbose = 5)

#---------------------------------------------- Soft-threshold Plot
# Install and load ggplot2 if needed
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# Function to create data frame and find appropriate power
prepare_sft_plot <- function(sft) {
  df <- data.frame(
    Power = sft$fitIndices[,1],
    SFT_R2 = -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
    MeanConnectivity = sft$fitIndices[,5]
  )
  
  # Find the first power that has an R² value greater than 0.9
  selected_power <- df$Power[which(df$SFT_R2 > 0.9)[1]]
  
  list(df = df, selected_power = selected_power)
}

# Prepare data for pre and post
pre_result <- prepare_sft_plot(sft_pre)
post_result <- prepare_sft_plot(sft_post)

# Plot for pre
ggplot(pre_result$df, aes(x=Power)) +
  geom_line(aes(y=SFT_R2), color="blue", size=1) +
  geom_point(aes(y=SFT_R2), color="blue") +
  geom_line(aes(y=MeanConnectivity/10), color="green", size=1) +
  geom_point(aes(y=MeanConnectivity/10), color="green") +
  geom_hline(yintercept=0.9, color="red", linetype="dashed") +
  geom_point(data=subset(pre_result$df, Power == pre_result$selected_power),
             aes(x=Power, y=SFT_R2), color="red", size=4, shape=21, fill="red") +
  scale_y_continuous(
    name = "Scale Free Topology Model Fit (R²)",
    sec.axis = sec_axis(~.*10, name="Mean Connectivity")
  ) +
  ggtitle(paste("Soft-thresholding power selection (Pre)\nSelected power =", pre_result$selected_power)) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
ggsave("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Important_modules_Blood\\Important_modules-Results\\pre_soft_thresholding_plot.png", width=8, height=6, dpi=300)

# Plot for post
ggplot(post_result$df, aes(x=Power)) +
  geom_line(aes(y=SFT_R2), color="blue", size=1) +
  geom_point(aes(y=SFT_R2), color="blue") +
  geom_line(aes(y=MeanConnectivity/10), color="green", size=1) +
  geom_point(aes(y=MeanConnectivity/10), color="green") +
  geom_hline(yintercept=0.9, color="red", linetype="dashed") +
  geom_point(data=subset(post_result$df, Power == post_result$selected_power),
             aes(x=Power, y=SFT_R2), color="red", size=4, shape=21, fill="red") +
  scale_y_continuous(
    name = "Scale Free Topology Model Fit (R²)",
    sec.axis = sec_axis(~.*10, name="Mean Connectivity")
  ) +
  ggtitle(paste("Soft-thresholding power selection (Post)\nSelected power =", post_result$selected_power)) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Important_modules_Blood\\Important_modules-Results\\post_soft_thresholding_plot.png", width=8, height=6, dpi=300)

#-------------------- SOFT POWER
softPower_pre <- 36  
softPower_post <- 44 

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
dynamicMods_pre <- cutreeDynamic(dendro = geneTree_pre, distM = dissTOM_pre, deepSplit = 2,
                                 pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
moduleColors_pre <- labels2colors(dynamicMods_pre)
png("dendrogram_pre.png", width = 800, height = 600)
plotDendroAndColors(geneTree_pre, moduleColors_pre, "Module colors (pre)", dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)
ggsave("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Important_modules_Blood\\Important_modules-Results\\dendrogram_pre.png", width=8, height=6, dpi=300)
dev.off()

dynamicMods_post <- cutreeDynamic(dendro = geneTree_post, distM = dissTOM_post, deepSplit = 2,
                                  pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
moduleColors_post <- labels2colors(dynamicMods_post)
png("dendrogram_post.png", width = 800, height = 600)
plotDendroAndColors(geneTree_post, moduleColors_post, "Module colors (Post)", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)
ggsave("D:\\Thesis\\Proposal\\Analysis\\WGCNA\\Blood-Patient\\Important_modules_Blood\\Important_modules-Results\\dendrogram_post.png", width=8, height=6, dpi=300)
dev.off()

# ============================ Eigen gene calculation
MEList_pre <- moduleEigengenes(datExpr_pre, colors = moduleColors_pre)
MEs_pr <- MEList_pre$eigengenes

MEList_post <- moduleEigengenes(datExpr_post, colors = moduleColors_post)
MEs_post <- MEList_post$eigengenes

# ============================
num_modules_pre <- ncol(MEs_pr)  # Number of modules in the Pre network
num_modules_post <- ncol(MEs_post)  # Number of modules in the Post network
num_samples <- nrow(MEs_pr)  # Number of samples
module_importance <- numeric(num_modules_pre)
for (i in 1:num_modules_pre) {
  current_importance <- 0
  for (j in 1:num_modules_post) {
    cor_ME_pre_post <- cor(MEs_pr[, i], MEs_post[, j], use = "pairwise.complete.obs")
    cor_ME_post_trait <- cor(MEs_post[, j], trait, use = "pairwise.complete.obs")
    current_importance <- current_importance + cor_ME_pre_post * cor_ME_post_trait
  }
  module_importance[i] <- current_importance
}

# Display the results
module_importance
