# Pipeline by Ruth De Paula
## Proteomics data preprocessing (for ROSMAP) ##
# Data was downloaded from Synapse.
# Very important to use dataframe preprocessed using TAMPOR, otherwise the data will be very biased by batch effects.

setwd("C:/Users/ruthb/Documents/BCM/PhD project/Results")
options(digits=22)  ## necessary to ensure R is going to read big key numbers

# Load required libraries
library(openxlsx)
library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(plotly)
library(tidyr)
library(dplyr)
library(lmtest)
library(ggplot2) # graphs
library(ggrepel) # graphs

### Proteomics preprocessing ----
# Not regressed (TAMPOR) - proteomics round 1 (n=400 samples)
data <- fread("Proteomics/ROSMAP data/Round1-2.Unregressed_batch-corrected_TMT_reporter_abundance.csv")

metadata_r1 <- read.table("Proteomics/ROSMAP data/round1_ROSMAP_sample_metadata_Johnson_2020.txt", sep="\t", header=T)
names(metadata_r1)[names(metadata_r1) == 'ROSMAP_projid'] <- 'projid'
sample_key <- read.table("AMP-AD_rosmap_id_key.csv", sep=",", header=T)[,1:2]
metadata_r1 <- merge(sample_key, metadata_r1, by="projid")

metadata_r1$search_name <- metadata_r1$TMT_Batch.Channel
colnames(data) <- gsub("^b0", "b", colnames(data))

# Create a named vector mapping TMT_Batch.Channel to individualID
channel_to_id <- setNames(metadata_r1$individualID, metadata_r1$TMT_Batch.Channel)

# Replace colnames in `data` that match any TMT_Batch.Channel
colnames(data) <- sapply(colnames(data), function(col) {
  if (col %in% names(channel_to_id)) {
    channel_to_id[[col]]
  } else {
    col
  }
})

# Reorganize
names <- data$V1
data <- as.data.frame(t(data[,-1]))
colnames(data) <- names


# Still need to sum all peptides that belong to the same protein and impute.

# Step 1: Impute missing values with the column-wise minimum
data_imputed <- data  # copy structure
# data_imputed <- data_protein  # for option 1C
for (j in seq_len(ncol(data_imputed))) {
  col_min <- min(data_imputed[[j]], na.rm = TRUE)
  data_imputed[[j]][is.na(data_imputed[[j]])] <- col_min
}

# Step 2: Remove peptides (columns) with >30% missing *before imputation*
# Use original data to calculate missingness
missing_fraction <- colMeans(is.na(data))
data_filtered <- data_imputed[, missing_fraction <= 0.3]

# Step 3: Collapse peptides into proteins by summing columns with same protein name

# Extract protein names from colnames (e.g., "VAMP1|P23763" → "VAMP1")
protein_names <- sub("\\|.*", "", colnames(data_filtered))

# Remove peptide columns with protein name "0"
valid_cols <- protein_names != "0"
data_filtered <- data_filtered[, valid_cols]
protein_names <- protein_names[valid_cols]

# Use aggregate approach to sum peptides for the same protein (column-wise grouping)
# Transpose → group rows → sum → transpose back
data_transposed <- t(data_filtered)  # now: rows = peptides, cols = individuals
data_grouped <- rowsum(data_transposed, group = protein_names)
data_protein <- as.data.frame(t(data_grouped))  # now: rows = individuals, cols = proteins


# Step 4: Log2 transform (no mean-centering)
data_log2 <- as.data.frame(log2(data_protein))

path <- "C:/Users/ruthb/Documents/BCM/PhD project/Results/Proteomics/ROSMAP data/"
saveRDS(data_protein, paste0(path,"Round1-2.Unregressed_batch-corrected_NA_removed.rds"))  # 7596 total peptides imput / 5084 no NA
saveRDS(data_filtered, paste0(path,"Round1-2.Unregressed_batch-corrected_imputed_proteins.rds"))  # 7458 total proteins imput
saveRDS(data_log2, paste0(path,"Round1-2.Unregressed_batch-corrected_imputed_proteins_log2.rds"))  # 7596 peptides imput / 5084 no NA / 7458 total proteins imput


# Quick PCA plot
metadata <- metadata_r1[,c("individualID","Batch")]
data_log2$individualID <- rownames(data_log2)
data_log2 <- merge(metadata, data_log2, by="individualID")

# Run PCA (center = TRUE ensures proper PCA, do not scale if you want to preserve intensity scale)
pca_res <- prcomp(as.matrix(data_log2[,-c(1:2)]), center = TRUE, scale. = TRUE)

# Plot the first two PCs using base R
group_vector <- data_log2$Batch

plot(pca_res$x[, 1], pca_res$x[, 2],
     col = as.factor(group_vector),
     pch = 19,
     xlab = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "%)"),
     ylab = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "%)"),
     main = "PCA Colored by Batch")
# legend("topright", legend = levels(as.factor(group_vector)), col = 1:length(unique(group_vector)), pch = 19)


# Now, 3d PCA plot
pca_df <- as.data.frame(pca_res$x[, 1:3])  # PC1–PC3
pca_df$Sample <- data_log2$individualID
pca_df$Batch <- data_log2$Batch

# 3D interactive plot
fig <- plot_ly(
  data = pca_df,
  x = ~PC1, y = ~PC2, z = ~PC3,
  # color = ~Batch,
  # colors = "Set2",
  type = 'scatter3d',
  mode = 'markers+text',
  text = ~Sample,
  textposition = 'top center',
  marker = list(size = 5, color = 'steelblue')
) %>%
  layout(
    title = "3D PCA of Protein Abundances",
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )

fig

# Individual R5211056 seems to be an obvious outlier, and maybe R2343042
