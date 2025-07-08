# Pipeline by Ruth De Paula
## Additive and recessive models for variant effect analysis ##

setwd("/mnt/belinda_local/ruth/data/omics_integration/CAND_ADPD_450")
options(digits=22)  ## necessary to ensure R is going to read big key numbers

############## LIBRARIES ##############
library(openxlsx)
library(ggplot2)
library(biomaRt)
library(ropls)
library(ggbeeswarm)  # For geom_quasirandom()
library(data.table)
#######################################

### Additive model statistics - For SPTSSB paper ----
# First, load filtered variant spreadsheet for gene of interest
filtered_vars_hom <- read.csv("SPTSSB_coding_variants_with_sample_names_annot_homo.txt", sep="\t")
filtered_vars_het <- read.csv("SPTSSB_coding_variants_with_sample_names_annot_hetero.txt", sep="\t")
filtered_vars_hom <- read.csv("ARSA_coding_variants_with_sample_names_annot_homo.txt", sep="\t") # need to generate
filtered_vars_het <- read.csv("ARSA_coding_variants_with_sample_names_annot_hetero.txt", sep="\t") # need to generate

# SPTSSB variants
filtered_vars_hom <- filtered_vars_hom[filtered_vars_hom$POS %in% 161359842 & filtered_vars_hom$ALT %in% "G",]  # rs1450522
filtered_vars_het <- filtered_vars_het[filtered_vars_het$POS %in% 161359842 & filtered_vars_het$ALT %in% "G",]  # rs1450522

# ARSA variants
# filtered_vars_hom <- filtered_vars_hom[filtered_vars_hom$POS %in% c(51063477,51064416) & filtered_vars_hom$ALT %in% "C",]  # rs6151429 and rs2071421, hg19
# filtered_vars_het <- filtered_vars_het[filtered_vars_het$POS %in% c(51063477,51064416) & filtered_vars_het$ALT %in% "C",]  # rs6151429 and rs2071421, hg19
filtered_vars_hom <- filtered_vars_hom[filtered_vars_hom$POS %in% c(50625988) & filtered_vars_hom$ALT %in% "C",]  # rs2071421, hg38 - the other one is not counted in this file
filtered_vars_het <- filtered_vars_het[filtered_vars_het$POS %in% c(50625988) & filtered_vars_het$ALT %in% "C",]  # rs2071421, hg38 - the other one is not counted in this file

# Load id key
id_key <- read.xlsx("sample_correspondence_pVCF_Lichtarge_lab.xlsx")

# Get lists of carriers - 1 variant in the gene / merged variants
hom_carriers <- unlist(strsplit(filtered_vars_hom[,"Samples"], "/"))
het_carriers <- unlist(strsplit(filtered_vars_het[,"Samples"], "/"))
all_carriers <- c(hom_carriers, het_carriers)

# Create dataframe with subjects and carrier status
selected_samples <- id_key[id_key$Cohort %in% c("CAND-Methodist","CAND"),][,c(1,2)]
variant_df <- as.data.frame(matrix(nrow=nrow(selected_samples), ncol=3))
variant_df[,c(1:2)] <- selected_samples[,c(1:2)]
colnames(variant_df) <- c("wgs_id","DNA_id","carrier_status")

for (i in 1:nrow(variant_df)) {
  
  # if (variant_df[i,"wgs_id"] %in% hom_carriers) {variant_df[i,"carrier_status"] <- 2} # additive analysis
  # if (variant_df[i,"wgs_id"] %in% het_carriers) {variant_df[i,"carrier_status"] <- 1} # additive analysis
  if (variant_df[i,"wgs_id"] %in% hom_carriers) {variant_df[i,"carrier_status"] <- 1} # recessive analysis
  if (variant_df[i,"wgs_id"] %in% het_carriers) {variant_df[i,"carrier_status"] <- 0} # recessive analysis
  # if (variant_df[i,"wgs_id"] %in% all_carriers) {variant_df[i,"carrier_status"] <- 1} # dominant analysis
  
}
variant_df[is.na(variant_df$carrier_status),]$carrier_status <- 0

# Rescue metadata
metadata <- read.xlsx("/mnt/belinda_local/ruth/data/metabolomics/New_CAND_400+/updated_CAND_450_metadata_July_2024.xlsx")[,c(3,5,6,7,9,10)]
colnames(metadata) <- c("DNA_id","plasma_id","age_collection","sex","age_onset","diagnosis")
metadata$sex <- gsub(0, "female", metadata$sex)
metadata$sex <- gsub(1, "male", metadata$sex)
variant_df <- merge(metadata, variant_df, by="DNA_id")

# Rescue metabolites
metab_levels <- readRDS("/mnt/belinda_local/ruth/data/metabolomics/New_CAND_400+/metab_ab_norm_filtered_CAND_400+_2.rds")
# metab_levels <- readRDS("/mnt/belinda_local/ruth/data/metabolomics/New_CAND_400+/metab_ab_norm_filtered_CAND_400+_NO_ZSCORE.rds")
metab_levels$PARENT_SAMPLE_NAME <- rownames(metab_levels)

sample_corresp <- read.xlsx("/mnt/belinda_local/ruth/data/metabolomics/New_CAND_400+/sample_correspondence_CAND.xlsx")
colnames(sample_corresp)[2] <- "plasma_id"
metab_levels <- merge(sample_corresp, metab_levels, by="PARENT_SAMPLE_NAME")

variant_metab_df <- merge(variant_df, metab_levels, by="plasma_id")  # just 147 PD and 148 controls remain
colnames(variant_metab_df)[-c(1:9)] <- paste0("X", colnames(variant_metab_df)[-c(1:9)])


### Prepare for analysis ----
# variant_metab_sub <- variant_metab_df[variant_metab_df$diagnosis %in% c("PD","NCI"),] # For PD vs Control  # 2 PD and 2 Ctrl w/ metabolomics missing (bc don't have genetic data)
variant_metab_sub <- variant_metab_df[variant_metab_df$diagnosis %in% c("PD"),] # For SPTSSB carrier status

stats_df <- as.data.frame(matrix(nrow=ncol(variant_metab_sub[,-c(1:9)])*4, ncol=6))
# rownames(stats_df) <- rep(colnames(variant_metab_sub[,-c(1:9)]), 4)
colnames(stats_df) <- c("CHEM_ID","rn","Estimate","Std. Error","t value","Pr(>|t|)")
stats_df$CHEM_ID <- rep(colnames(variant_metab_sub[,-c(1:9)]), each=4)

for (i in colnames(variant_metab_sub[,-c(1:9)])) {
  
  # Linear model
  model <- lm(variant_metab_sub[,i] ~ carrier_status + age_collection + sex, data=variant_metab_sub)
  # model <- lm(variant_metab_sub[,i] ~ diagnosis + age_collection + sex, data=variant_metab_sub)
  summary_result <- summary(model)
  result <- as.data.table(summary_result$coefficients, keep.rownames=T)
  
  stats_df[stats_df$CHEM_ID %in% i,-1] <- result
  
}


# Rescue metabolite information
data_dictionary <- read.xlsx("/mnt/belinda_local/ruth/data/metabolomics/New_CAND_400+/CAND_400+_plasma_Metabolon_chemical_annotation_common.xlsx")
data_dictionary$CHEM_ID <- paste0("X", data_dictionary$CHEM_ID)

stats_df_2 <- merge(data_dictionary[,c("CHEM_ID","SUPER_PATHWAY","SUB_PATHWAY","CHEMICAL_NAME")], stats_df, by="CHEM_ID")

write.table(stats_df_2, "recessive_model_SPTSSB_rs1450522_PD_only.txt", row.names = F, col.names = T, sep="\t")


sphingo_in_PD_and_ttest <- c("X100009038","X100008957","X100009025","X100015789","X100015751","X100009030")

stats_df_3 <- stats_df_2[stats_df_2$CHEM_ID %in% sphingo_in_PD_and_ttest,]
