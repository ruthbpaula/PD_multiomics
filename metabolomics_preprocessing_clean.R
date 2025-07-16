# Pipeline by Ruth De Paula
## Metabolomics data preprocessing ##
# This pipeline should work for both ROSMAP brains and BCM plasma samples (both use Metabolon's HD4 panel).
# Since ROSMAP data is older (before 2020) and imputed data was not provided in Synapse portal, extra steps are needed:
# - Removal of "NIST" standards, and calculation of imputed values/imputation
# CAND data is newer (2022-2023), and imputed table is already provided by Metabolon (imputation as minimum).
# In that case, we will use their imputed table and simply remove the metabolites that initially had too many missing values (>30%).

setwd("C:/Users/ruthb/Documents/BCM/PhD project/Results/Metabolomics/ROSMAP_DATA/MetabolonHD4_brain_syn26007830_syn26427298")

options(digits=22)  ## necessary to ensure R is going to read big key numbers

############## LIBRARIES ##############
library(openxlsx) # excel things
library(ggplot2)
#######################################

## Basic pipeline for metabolomics data preprocessing:
# 0. Load data and (only for ROSMAP) remove NIST samples
# 1. Plot data distribution before any preprocessing
# 2. Select metabolites (present in ~70% of the patients)
# 3. (only for ROSMAP) Imputation of missing values
# 4. Data transformation. Approach:
# --- log2Fold + mean centering/z-score

#_____________________________________________________________________________________________________________________
### 0. Load data and remove "NIST" columns (these are standard plasma samples that will not be used right now) ----
metabolon_data <- as.data.frame(read.table("ROSMAP_Metabolon_HD4_Brain514_assay_data_syn25985690.csv", header = T, sep=","))
metabolon_data <- metabolon_data[!(metabolon_data$individualID %in% c("NIST01","NIST02","NIST03","NIST04","NIST05","NIST06","NIST07")),]

# ______________________________________________________________________________________
### 1. Plot data distribution before preprocessing ----
metabolon_data_chem_2 <- metabolon_data

# Pick a single metabolite as example
hist(as.numeric(metabolon_data_chem_2$X35),  ## non-normal distribution
     xlab='Metabolite abundance across subjects')

# ______________________________________________________________________________________
### 2. Select metabolites (expressed in ~70% of individuals) ----
# 70% threshold is = 360 individuals (70% of 514)

# Count number of cells that are not NA per column/metabolite
metab_noNA_counts <- data.frame(matrix(nrow=1,ncol=ncol(metabolon_data_chem_2)))
colnames(metab_noNA_counts) <- colnames(metabolon_data_chem_2)

for (i in 9:ncol(metabolon_data_chem_2)) {
  metab_noNA_counts[1,i] <- sum(!is.na(metabolon_data_chem_2[,i]))
}

metab_noNA_counts <- as.data.frame(t(metab_noNA_counts[,9:ncol(metabolon_data_chem_2)]))

# Select metabolites that are present in 360 subjects or more
metab_noNA_counts_70 <- as.data.frame(metab_noNA_counts[metab_noNA_counts$V1 > 359,])  ## 685 of 1055 metabolites remained (~65%)
metab_noNA_counts_70 <- metab_noNA_counts_70[order(-metab_noNA_counts_70[,1]), , drop = FALSE]
metab_noNA_counts_sort <- metab_noNA_counts[order(-metab_noNA_counts[,1]), , drop = FALSE]
metab_noNA_counts_sort$label <- rownames(metab_noNA_counts_sort)
metab_noNA_counts_70$label <- metab_noNA_counts_sort[1:685,2]

string <- as.character(metab_noNA_counts_70$label)
metabolon_data_chem_3 <- metabolon_data_chem_2[, string, drop = FALSE]
rownames(metabolon_data_chem_3) <- metabolon_data_chem_2$individualID

# ______________________________________________________________________________________
### 3. Imputation of missing values ----

# metabolon_data_chem_3 <- read.table("metabolon_data_chem_3_chem_id.txt", sep="\t", header=T)
# colnames(metabolon_data_chem_3) <- gsub(".*\\.", "", colnames(metabolon_data_chem_3))

## Conservative approach 1: impute with minimum value
metab_impute_2 <- data.frame(matrix(nrow=1,ncol=ncol(metabolon_data_chem_3)))
colnames(metab_impute_2) <- colnames(metabolon_data_chem_3)

for (i in 1:ncol(metabolon_data_chem_3)) {
  metab_impute_2[1,i] <- min(as.numeric(metabolon_data_chem_3[,i]), na.rm = TRUE)
}

# Replace NA with imputation value
metabolon_data_chem_4 <- metabolon_data_chem_3

for (i in 1:ncol(metabolon_data_chem_4)) {
  subset_df <- as.data.frame(metabolon_data_chem_4[,i])
  subset_df[is.na(subset_df)] <- metab_impute_2[1,i]
  metabolon_data_chem_4[,i] <- subset_df
}

# saveRDS(metabolon_data_chem_4, "metabolon_data_chem_4_imput_as_min_chem_id.rds")


# Data distribution - after imputation
# Pick a single metabolite as example
hist(as.numeric(metabolon_data_chem_4$X35),  ## non-normal distribution
     xlab='Metabolite abundance across subjects')

# ______________________________________________________________________________________
### 4. Data transformation ----
## Conservative approach: log2Fold + z-score (for mean centering, use scale = FALSE)
metabolon_data_chem_normalized <- log2(metabolon_data_chem_4)
metabolon_data_chem_normalized <- as.data.frame(scale(metabolon_data_chem_normalized))  ## z-scores

# Data distribution - after transformation
# Pick a single metabolite as example
hist(as.numeric(metabolon_data_chem_normalized$X35),  ## non-normal distribution
     xlab='Metabolite abundance across subjects')

## Test data for normality
shapiro.test(metabolon_data_chem_normalized$X35)

saveRDS(metabolon_data_chem_normalized, "metabolon_data_chem_normalized_2_imput_as_min_chem_id.rds")

