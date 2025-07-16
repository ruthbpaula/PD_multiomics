# Pipeline by Ruth De Paula
## Metabolomics and proteomics differential "expression" ##
# This pipeline should work for both ROSMAP brains and BCM plasma samples (both use Metabolon's HD4 panel).
# This can be similarly adapted to proteomics data analysis, since the data is also based on mass-spectrometry.

setwd("C:/Users/ruthb/Documents/BCM/PhD project/Results/Metabolomics/BCM_DATA/CAND (RAPPID and DNA bank) Plasma/New_CAND_400+")
options(digits=22)  ## necessary to ensure R is going to read big key numbers

############## LIBRARIES ##############
library(aod)
library(lmtest)
library(car)
library(tidyr)
library(ggpubr)
library(dplyr)
library(openxlsx)
library(ggplot2) # graphs
library(hrbrthemes) # graphs
library(viridis) # graphs
library(ComplexHeatmap) # graphs
library(limma) # graphs
library(pdist) # graphs
library("RColorBrewer") # graphs
library(circlize) # graphs
library(seriation) # graphs
library(ggfortify) # graphs
library(ggrepel) # graphs
#######################################

## Load data ----
metab_ab_norm_filtered <- read.xlsx("CAND_400+_plasma_Metabolon_batch_norm_imput_common.xlsx")  ## New complete dataset
rownames(metab_ab_norm_filtered) <- metab_ab_norm_filtered$PARENT_SAMPLE_NAME
metab_ab_norm_filtered <- metab_ab_norm_filtered[,-1]

missingness <- read.table("CAND_400+_count_blanks.txt", header=T)
missing_metab <- missingness[missingness$n_blanks>135,]$CHEM_ID  # 30%+ missing in merged table
metab_ab_norm_filtered <- metab_ab_norm_filtered[,!(colnames(metab_ab_norm_filtered) %in% missing_metab)] # 1171 metabolites

# Data is not yet log2 transformed. Now, do log2 transformation and calculate z-scores
metab_ab_norm_filtered <- as.data.frame(scale(log2(metab_ab_norm_filtered)))  # scale by column
# metab_ab_norm_filtered <- as.data.frame(log2(metab_ab_norm_filtered))

# Rescuing metadata
metadata_key <- read.xlsx("C:/Users/ruthb/Documents/BCM/PhD project/Results/BCM_Patient_Metadata/CAND_400+/sample_correspondence_CAND.xlsx")
colnames(metadata_key)[2] <- "PLASMA_ID"

metadata <- read.xlsx("C:/Users/ruthb/Documents/BCM/PhD project/Results/BCM_Patient_Metadata/CAND_400+/updated_CAND_450_metadata_July_2024.xlsx")
colnames(metadata)[11] <- "GROUP_NAME"
metadata$GROUP_NAME <- sub("NCI", "Control", metadata$GROUP_NAME)
colnames(metadata)[6] <- "PLASMA_ID"
colnames(metadata)[7] <- "AGE_PLASMA"
colnames(metadata)[8] <- "GENDER"

metadata_2 <- metadata[,c("PLASMA_ID","GROUP_NAME","AGE_PLASMA","GENDER")]
metadata_2 <- merge(metadata_key, metadata_2, by="PLASMA_ID")

metab_ab_norm_filtered$PARENT_SAMPLE_NAME <- rownames(metab_ab_norm_filtered)
metab_ab_norm_filtered <- merge(metadata_2, metab_ab_norm_filtered, by="PARENT_SAMPLE_NAME")

metab_ab_norm_filtered$GENDER <- as.factor(metab_ab_norm_filtered$GENDER)
metab_ab_norm_filtered$GROUP_NAME <- factor(metab_ab_norm_filtered$GROUP_NAME) # levels = c("Control", "PD") -> Control as ref

### **OPTIONAL - For calculating and saving residuals (USE ALL SUBJECTS IN THIS STEP) ----
# This may be important in the future to generate boxplots and heatmaps.
residualized_z_scores <- data.frame(matrix(nrow=nrow(metab_ab_norm_filtered),ncol=ncol(metab_ab_norm_filtered)-5))
colnames(residualized_z_scores) <- colnames(metab_ab_norm_filtered)[6:ncol(metab_ab_norm_filtered)]
rownames(residualized_z_scores) <- metab_ab_norm_filtered$PARENT_SAMPLE_NAME

for (i in 6:ncol(metab_ab_norm_filtered)) {
  result0 <- scale(as.data.frame(lm(metab_ab_norm_filtered[,i] ~ AGE_PLASMA + GENDER, data=metab_ab_norm_filtered)$residuals))
  residualized_z_scores[,i-5] <- result0[1:nrow(result0),]
}

residualized_z_scores <- cbind(metab_ab_norm_filtered[,1:5], residualized_z_scores)

saveRDS(residualized_z_scores, "residualized_z_scores_CAND_400+_all_chemID_age_and_sex_updated_meta_NOT_ZSCORED.rds")


## Differential expression itself ----
linear_LSD_classif <- data.frame(matrix(nrow=ncol(metab_ab_norm_filtered),ncol=2))
linear_LSD_classif[,1] <- colnames(metab_ab_norm_filtered)
colnames(linear_LSD_classif) <- c("CHEM_ID","lrt_test")

# OPTIONAL - Exclude subjects that are not of interest
metab_ab_norm_filtered <- metab_ab_norm_filtered[!(metab_ab_norm_filtered$GROUP_NAME %in% c("MCI","AD","OD")),]

metab_ab_norm_filtered$GROUP_NAME <- as.factor(metab_ab_norm_filtered$GROUP_NAME)
metab_ab_norm_filtered$GENDER <- as.factor(metab_ab_norm_filtered$GENDER)

for (i in 6:ncol(metab_ab_norm_filtered)) {
  model_full <- lm(metab_ab_norm_filtered[,i] ~ GROUP_NAME + AGE_PLASMA + GENDER, data=metab_ab_norm_filtered)
  model_reduced <- lm(metab_ab_norm_filtered[,i] ~ AGE_PLASMA + GENDER, data=metab_ab_norm_filtered)
  
  # LRT test
  lrt_res <- lrtest(model_full, model_reduced)  ## If this p-value is bigger than 0.05, it means the full and reduced model can be 
  ## considered the same and then diagnosis does not differentiate metabolite levels
  linear_LSD_classif[i,2] <- lrt_res[2,5]
}

# FDR correction of LRT test
linear_LSD_classif <- as.data.frame(linear_LSD_classif[6:nrow(linear_LSD_classif),])
linear_LSD_classif <- cbind(linear_LSD_classif, p.adjust(linear_LSD_classif$lrt_test, method = "BH"))
colnames(linear_LSD_classif) <- c("CHEM_ID", "lrt_test", "padj_lrt")
rownames(linear_LSD_classif) <- linear_LSD_classif$CHEM_ID

# After running LRT, calculate mean metabolite level per metabolite, for control and for case separately. 
# Use that to evaluate whether a certain metabolite is upregulated or downregulated in the dataset.
control_data <- metab_ab_norm_filtered[metab_ab_norm_filtered$GROUP_NAME %in% "Control",]
PD_data <- metab_ab_norm_filtered[metab_ab_norm_filtered$GROUP_NAME %in% "PD",]

linear_LSD_classif$PD_mean <- NA
linear_LSD_classif$ctrl_mean <- NA
linear_LSD_classif$PD_dysreg <- NA

for (i in 1:nrow(linear_LSD_classif)) {
  metab_name <- rownames(linear_LSD_classif)[i]
  linear_LSD_classif[i,"PD_mean"] <- mean(PD_data[,metab_name])
  linear_LSD_classif[i,"ctrl_mean"] <- mean(control_data[,metab_name])
  
  if(linear_LSD_classif[i,"PD_mean"] > linear_LSD_classif[i,"ctrl_mean"]){
    linear_LSD_classif[i,"PD_dysreg"] <- "up"
  }
  
  if(linear_LSD_classif[i,"PD_mean"] < linear_LSD_classif[i,"ctrl_mean"]){
    linear_LSD_classif[i,"PD_dysreg"] <- "down"
  }
}


#### Calculating log2FC ----
# First, load imputed but not log2/z-scored dataset
metab_not_transformed <- read.xlsx("CAND_400+_plasma_Metabolon_batch_norm_imput_common.xlsx")
rownames(metab_not_transformed) <- metab_not_transformed$PARENT_SAMPLE_NAME
# OPTIONAL - Exclude metab with a lot of missingness
missingness <- read.table("CAND_400+_count_blanks.txt", header=T)
missing_metab <- missingness[missingness$n_blanks>135,]$CHEM_ID  # 30%+ missing in merged table
metab_not_transformed <- metab_not_transformed[,!(colnames(metab_not_transformed) %in% missing_metab)] # 1171 metabolites

# Rescue case/control information
path_classif_meta <- metab_ab_norm_filtered[,c("PARENT_SAMPLE_NAME","GROUP_NAME")]
metab_not_transformed <- merge(path_classif_meta, metab_not_transformed, by="PARENT_SAMPLE_NAME")

# Then, prepare data and calculate log2FC
control_data <- metab_not_transformed[metab_not_transformed$GROUP_NAME %in% "Control",]
PD_data <- metab_not_transformed[metab_not_transformed$GROUP_NAME %in% "PD",]

linear_LSD_classif$PD_mean_no_transf <- NA
linear_LSD_classif$ctrl_mean_no_transf <- NA
linear_LSD_classif$log2FC <- NA

for (i in 1:nrow(linear_LSD_classif)) {
  metab_name <- rownames(linear_LSD_classif)[i]
  linear_LSD_classif[i,"PD_mean_no_transf"] <- mean(PD_data[,metab_name])
  linear_LSD_classif[i,"ctrl_mean_no_transf"] <- mean(control_data[,metab_name])
  
  linear_LSD_classif[i,"log2FC"] <- log2(linear_LSD_classif[i,"PD_mean_no_transf"]/linear_LSD_classif[i,"ctrl_mean_no_transf"])
}

# write.table(linear_LSD_classif, "linear_LSD_classif_CAND_400+_PD_vs_ctrl_updated_meta.txt", sep="\t", col.names = T, row.names = F)

# The above-mentioned method to calculate log2FC does not take covariates into consideration, but the resulting log2FC look better when plotting Volcano plots.
# If you need both the p-values and the log2FC values to be corrected for covariates, you need to extract the beta from the linear regression model.
# colnames(metab_ab_norm_filtered)[-c(1:5)] <- lapply(colnames(metab_ab_norm_filtered)[-c(1:5)], gsub, pattern='^', replacement="X")
linear_LSD_classif$beta_log2FC <- NA
linear_LSD_classif$linear_p <- NA

for (i in 1:nrow(linear_LSD_classif)) {
  metab_name <- rownames(linear_LSD_classif)[i]
  
  model <- lm(metab_ab_norm_filtered[,metab_name] ~ GROUP_NAME + AGE_PLASMA + GENDER, data=metab_ab_norm_filtered)
  linear_LSD_classif[i,"beta_log2FC"] <- coef(model)[2]
  linear_LSD_classif[i,"linear_p"] <- summary(model)$coefficients[2,4]
}

linear_LSD_classif$padj_linear <- p.adjust(linear_LSD_classif$linear_p, method = "BH")

write.table(linear_LSD_classif, "linear_LSD_classif_CAND_400+_PD_vs_ctrl_updated_meta_w_beta_log2FC_NOT_ZSCORED.txt", sep="\t", col.names = T, row.names = F)


### Volcano plot highlighting key metabolites ----
# Load differential expression results
# linear_LSD_classif <- read.xlsx("multiple_regression_models/linear_LSD_classif_PD_vs_NO_PD-AD_Metabolon_normalized_data.xlsx", rowNames = T)
rownames(linear_LSD_classif) <- lapply(rownames(linear_LSD_classif), gsub, pattern='^', replacement="X")
rownames(linear_LSD_classif) <- lapply(rownames(linear_LSD_classif), gsub, pattern=' ', replacement="")
linear_LSD_classif$metabolite_name <- rownames(linear_LSD_classif)

# Define list of metabolites of interest
metab.list <- c("X1342","X100004634","X100006129","X100006361","X999912410","X100006360","X117","X999918779","X100022484")  # 6 most significant

# OPTIONAL - Subset result
# full_result <- linear_LSD_classif
linear_LSD_classif <- full_result
linear_LSD_classif <- linear_LSD_classif[linear_LSD_classif$padj_lrt > 0.0000000005,]

# OPTIONAL - Plot p-values instead of padj
colnames(linear_LSD_classif)[3] <- "padj_temp"
colnames(linear_LSD_classif)[2] <- "padj_lrt"

## Plot (red dots are upregulated DEPs; blue dots are downregulated DEPs; grey dots are not DEPs)
options(digits=2)
ggplot(linear_LSD_classif, aes(y= -log10(linear_LSD_classif$padj_lrt), x=linear_LSD_classif$log2FC,
                               color = ifelse(linear_LSD_classif$padj_lrt < 0.05 & linear_LSD_classif$log2FC > 0, "red",
                                              ifelse(linear_LSD_classif$padj_lrt < 0.05 & linear_LSD_classif$log2FC < 0, "blue",
                                                     ifelse(linear_LSD_classif$padj_lrt > 0.05, "gray", NA))))) +
# ggplot(linear_LSD_classif, aes(y= -log10(linear_LSD_classif$padj_lrt), x=linear_LSD_classif$beta_log2FC,
#                                color = ifelse(linear_LSD_classif$padj_lrt < 0.05 & linear_LSD_classif$beta_log2FC > 0, "red",
#                                               ifelse(linear_LSD_classif$padj_lrt < 0.05 & linear_LSD_classif$beta_log2FC < 0, "blue",
#                                                      ifelse(linear_LSD_classif$padj_lrt > 0.05, "gray", NA))))) +
  geom_point(size = 1.2) +
  scale_color_manual(name = NULL, values = c('blue','gray','red'), labels = c('LFC>0','LFC<0','NS')) +
  theme_minimal() +
  theme(legend.position = 'none', plot.title = element_text(size = 20), text = element_text(size = 14)) +
  labs(y='-log10(p)', x = 'log2FC', title = paste0('n metab p<0.05= ', length(which(linear_LSD_classif$padj_lrt<0.05)),' / ', nrow(linear_LSD_classif))) +
  geom_label_repel(aes(label=ifelse(rownames(linear_LSD_classif) %in% metab.list, rownames(linear_LSD_classif),'')), size = 2.6, color = 'black', 
                   segment.size =0.5 , point.padding = 0.8, min.segment.length = 1, segment.color = 'grey70', box.padding = 0.3, max.overlaps = 10000) +
  geom_point(aes(stroke = ifelse(rownames(linear_LSD_classif) %in% metab.list, 2, 0))) +
  geom_vline(xintercept = 0, size = 0.6, color = 'black', linetype = 'dashed') #+
#geom_hline(yintercept = -log10(0.05), size = 0.6, color = 'black', linetype = 'dashed')
options(digits=22)

# Save displayed plot
ggsave("BCM_volcano_plot_LRT_age-collection_sex.pdf")
