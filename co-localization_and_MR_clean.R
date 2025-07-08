options(digits=22)  ## necessary to ensure R is going to read big key numbers

### 1. Co-localization with GWAS-eQTLs and GWAS-mQTLs ----
setwd("/mnt/belinda_local/ruth/data/omics_integration/colocalization_analysis")

# GitHub file: https://github.com/gan-orlab/GALC/blob/main/coloc (Senkevich et al., 2023)
library(coloc)          # Original coloc package, which is used to run coloc.abf() and sensitivity()
library(colochelpR)     # Helper package with wrapper functions
library(data.table)     # Loaded for fread(), which permits fast loading of files with many rows
library(ggpubr)         # For formatting of plots
library(here)           # For file path construction
library(purrr)
library(dplyr)
library(dbplyr)
library(tibble)
library(locuscomparer)
library(stringr)
library(ggplot2)
library(ggrepel)

#### Preparing datasets for coloc ----
# Will use hg19 files, to be consistent w/ SMR MR analysis later (SMR already provides formatted QTLs in hg19).

# Coloc GWAS-QTL file mandatory columns should be (first half comes from GWAS file, and second half comes from QTL file):
# "b_MA.gwas","varbeta.gwas","p_value.gwas","b_MA.eqtl","varbeta.eqtl","MAF.eqtl","p_value.eqtl","snp (rsid)"

# If you don't have "b_MA", "MAF", and "varbeta" columns, calculate them this way:
# SPTSSB_eqtl$b_MA <- SPTSSB_eqtl$Effect
# SPTSSB_eqtl[SPTSSB_eqtl$Freq1 > 0.50, ]$b_MA <- SPTSSB_eqtl[SPTSSB_eqtl$Freq1 > 0.50, ]$b_MA *-1
# 
# SPTSSB_eqtl$MAF <- SPTSSB_eqtl$Freq1
# SPTSSB_eqtl[SPTSSB_eqtl$Freq1 > 0.50, ]$MAF <- 1 - SPTSSB_eqtl[SPTSSB_eqtl$Freq1 > 0.50, ]$MAF
# 
# SPTSSB_eqtl$varbeta = SPTSSB_eqtl$StdErr*SPTSSB_eqtl$StdErr

# BESD eQTL files were downloaded from here: https://yanglab.westlake.edu.cn/software/smr/#DataResource
# eQTLs were extracted from SMR BESD files, using the following command line:
# smr-1.3.1-linux-x86_64/smr --beqtl-summary BrainMeta_cis_eqtl_summary/BrainMeta_cis_eQTL_chr3 --query 1 --snp-chr 3 --from-snp-kb 160577 --to-snp-kb 161678 --out SPTSSB_locus_cis_eQTL_BrainMeta_query_ALL_GENES

# nerve <- read.table("nerve_tibial_hg19_df_for_SPTSSB_coloc_ALL_GENES.txt", header=T, sep="\t") # nerve tibial GTEx SPTSSB eQTL
# hypothalamus <- read.table("hypothalamus_hg19_df_for_SPTSSB_coloc_ALL_GENES.txt", header=T, sep="\t") # hypothalamus GTEx SPTSSB eQTL
# cerebellum <- read.table("cerebellum_hg19_df_for_SPTSSB_coloc_ALL_GENES.txt", header=T, sep="\t") # cerebellum GTEx SPTSSB eQTL
# blood <- read.table("whole_blood_hg19_df_for_SPTSSB_coloc_ALL_GENES.txt", header=T, sep="\t") # whole blood GTEx SPTSSB eQTL
# accumbens <- read.table("nucleus_accumbens_hg19_df_for_SPTSSB_coloc_ALL_GENES.txt", header=T, sep="\t") # nucleus accumbens GTEx SPTSSB eQTL
# cortex <- read.table("brain_cortex_hg19_df_for_SPTSSB_coloc_ALL_GENES.txt", header=T, sep="\t") # brain cortex BrainMeta SPTSSB eQTL

blood <- read.table("whole_blood_hg19_df_for_SPTSSB_coloc_ALL_GENES_case-control.txt", header=T, sep="\t") # whole blood GTEx SPTSSB eQTL
cortex <- read.table("brain_cortex_hg19_df_for_SPTSSB_coloc_ALL_GENES_case-control.txt", header=T, sep="\t") # brain cortex BrainMeta SPTSSB eQTL
cortex <- cortex[cortex$Probe %in% "ENSG00000196542.8",]

# For single-cell, no genes were detected for SPTSSB locus in endothelial_cells
excitatory <- read.table("excitatory_neurons_hg19_df_for_SPTSSB_coloc_ALL_GENES_case-control.txt", header=T, sep="\t") # brain DLPFC ROSMAP snRNA-seq eQTL
# microglia <- read.table("microglia_hg19_df_for_SPTSSB_coloc_ALL_GENES_case-control.txt", header=T, sep="\t") # brain DLPFC ROSMAP snRNA-seq eQTL; unfortunately, no SPTSSB eQTLs detected
oligo <- read.table("oligodendrocytes_hg19_df_for_SPTSSB_coloc_ALL_GENES_case-control.txt", header=T, sep="\t") # brain DLPFC ROSMAP snRNA-seq eQTL
astro <- read.table("astrocytes_hg19_df_for_SPTSSB_coloc_ALL_GENES_case-control.txt", header=T, sep="\t") # brain DLPFC ROSMAP snRNA-seq eQTL
inhibitory <- read.table("inhibitory_neurons_hg19_df_for_SPTSSB_coloc_ALL_GENES_case-control.txt", header=T, sep="\t") # brain DLPFC ROSMAP snRNA-seq eQTL
# OPC <- read.table("OPC_hg19_df_for_SPTSSB_coloc_ALL_GENES_case-control.txt", header=T, sep="\t") # brain DLPFC ROSMAP snRNA-seq eQTL; unfortunately, no SPTSSB eQTLs detected

# mQTL summary stats were downloaded from here: https://omicscience.org/apps/mgwas/mgwas.table.php
# heptanoate <- read.table("heptanoate_hg19_df_for_SPTSSB_coloc.txt", header=T, sep="\t") # plasma Surendran mQTL
# oleate <- read.table("oleate-vaccenate_hg19_df_for_SPTSSB_coloc.txt", header=T, sep="\t") # plasma Surendran mQTL
# myristoleate <- read.table("myristoleate_hg19_df_for_SPTSSB_coloc.txt", header=T, sep="\t") # plasma Surendran mQTL
# sphinganine <- read.table("sphinganine_hg19_df_for_SPTSSB_coloc.txt", header=T, sep="\t") # plasma Surendran mQTL

heptanoate <- read.table("heptanoate_hg19_df_for_SPTSSB_coloc_case-control.txt", header=T, sep="\t") # plasma Surendran mQTL - allele 2 is effect allele
oleate <- read.table("oleate-vaccenate_hg19_df_for_SPTSSB_coloc_case-control.txt", header=T, sep="\t") # plasma Surendran mQTL - allele 2 is effect allele
myristoleate <- read.table("myristoleate_hg19_df_for_SPTSSB_coloc_case-control.txt", header=T, sep="\t") # plasma Surendran mQTL - allele 2 is effect allele
sphinganine <- read.table("sphinganine_hg19_df_for_SPTSSB_coloc_case-control.txt", header=T, sep="\t") # plasma Surendran mQTL - allele 2 is effect allele
CERd20_1_24_0 <- read.table("CERd20-1_24-0_hg19_df_for_SPTSSB_coloc_case-control.txt", header=T, sep="\t") # plasma Cadby mQTL
CERd20_1_24_1 <- read.table("CERd20-1_24-1_hg19_df_for_SPTSSB_coloc_case-control.txt", header=T, sep="\t") # plasma Cadby mQTL
dhCERd18_0_18_0 <- read.table("dhCERd18-0_18-0_hg19_df_for_SPTSSB_coloc_case-control.txt", header=T, sep="\t") # plasma Cadby mQTL
SM43_1 <- read.table("SM43-1_hg19_df_for_SPTSSB_coloc_case-control.txt", header=T, sep="\t") # plasma Cadby mQTL
DE22_6 <- read.table("DE22-6_hg19_df_for_SPTSSB_coloc_case-control.txt", header=T, sep="\t") # plasma Cadby mQTL

# Select QTL list of interest
df <- CERd20_1_24_1
unique(df$Gene)  # see which genes of interest have eQTLs
QTL <- "CERd20_1_24_1_mQTL_SPTSSB_PD-GWAS"

# Select gene eQTLs of interest
df_sub <- df[df$Gene == "SPTSSB", ]
df_sub <- df[df$Gene == "PPM1L", ]
df_sub <- df[df$Gene == "B3GALNT1", ]
df_sub <- df[df$Gene == "NMD3", ]
# df_sub <- df[df$Gene == "OTOL1", ]

# Remove zeroed/duplicated SNPs
df_sub <- df_sub[!(df_sub$snp %in% df_sub$snp[duplicated(df_sub$snp)]), ]


#### Run coloc ----
## GWAS vs. QTL
res <- coloc.abf(dataset1=list(beta=df_sub$b_MA.gwas, varbeta=df_sub$varbeta.gwas, type="cc", snp=df_sub$snp, pvalue=df_sub$p_value),
                 dataset2=list(beta=df_sub$b_MA.eqtl, varbeta=df_sub$varbeta.eqtl, N=nrow(df_sub), MAF=as.numeric(df_sub$MAF.eqtl), type="quant", snp=df_sub$snp, pvalue=df_sub$P.value))

## eQTL vs. mQTL
# Re-format input
eQTL <- excitatory[,c(3,15:ncol(excitatory))] # snRNA-seq
# eQTL <- cortex[,c(2:3,15:ncol(cortex))] # bulk RNA-seq
mQTL <- heptanoate[,c(3,15:ncol(heptanoate))]
df <- merge(eQTL, mQTL, by="SNP_ID")
df_sub <- df[df$Gene.x == "SPTSSB", ]
df_sub <- df_sub[!(df_sub$snp %in% df_sub$snp[duplicated(df_sub$snp)]), ]
# Run
res <- coloc.abf(dataset1=list(beta=df_sub$b_MA.eqtl.x, varbeta=df_sub$varbeta.eqtl.x, N=nrow(df_sub), MAF=as.numeric(df_sub$MAF.eqtl.x), type="quant", snp=df_sub$snp, pvalue=df_sub$P.value.x),
                 dataset2=list(beta=df_sub$b_MA.eqtl.y, varbeta=df_sub$varbeta.eqtl.y, N=nrow(df_sub), MAF=as.numeric(df_sub$MAF.eqtl.y), type="quant", snp=df_sub$snp, pvalue=df_sub$P.value.y))


results <- as.data.frame(t(res$summary))
results

all_res_output_name <- paste0("Cadby_", QTL, "_case-control.tsv")
# all_res_output_name <- paste0("GTEx-Surendran_", QTL, ".tsv")
write.table(results, all_res_output_name, sep = "\t")#, quote = F, row.names = F)

#### Stats per SNP instead of per region ----
# In colocalization analysis, the most important or significant SNPs are typically those 
# that have the highest posterior probabilities of being shared causal variants across traits. 
# Specifically, these SNPs are those with high posterior probabilities of colocalization 
# (i.e., shared causal variants), particularly from models like PP.H4 (common causal variant).

o <- order(res$results$SNP.PP.H4,decreasing=TRUE)
cs <- cumsum(res$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]
top_snps_index <- res$results[o,][1:w,]$snp
top_snps_stats <- res$results[o,][1:w,][,c("snp","SNP.PP.H4")] # High PP.H4.abf values (e.g., >0.5 or higher) suggest strong evidence for shared causality.
# top_snps_index <- res$results[o,][1:20,]$snp
# top_snps_stats <- res$results[o,][1:20,][,c("snp","SNP.PP.H4")] # High PP.H4.abf values (e.g., >0.5 or higher) suggest strong evidence for shared causality.

top_snps <- merge(top_snps_stats, df_sub, by="snp")

# save a table with SNPs
top_snps_output_name <- paste0("Cadby_", QTL, "_case-control_top_snps.tsv")
# top_snps_output_name <- paste0("GTEx-Surendran_", QTL, "_top_snps.tsv")
write.table(top_snps, top_snps_output_name, sep = "\t", quote = F)#, row.names = F)

#### Plot LocusCompare plots ----
# select top SNP
top_var <- top_snps[which.max(top_snps[[2]]), ]$snp
top_var

# Optional - instead of using top var as LD reference, use variant of choice
top_var <- "rs1450522" # SPTSSB

# prepare datasets for locuscompare
## GWAS vs. QTL
df_sub_eqtl <- df_sub %>%
  dplyr::select(snp, P.value) %>%
  dplyr::mutate(rsid = as.character(snp)) %>%
  dplyr::rename(pval = P.value)
df_sub_gwas <- df_sub %>%
  dplyr::select(snp, p_value) %>%
  dplyr::mutate(rsid = as.character(snp)) %>%
  dplyr::rename(pval = p_value)

# plot significant colocalization
plot <- locuscompare(in_fn1=df_sub_gwas,
                     in_fn2=df_sub_eqtl,
                     title1="PD GWAS",
                     title2="Cer(d20:1/24:1) mQTL",
                     # title2="QTL",
                     snp=top_var,
                     genome="hg19",
                     population = "EUR")
plot


## eQTL vs. mQTL
df_sub_mqtl <- df_sub %>%
  dplyr::select(snp, P.value.y) %>%
  dplyr::mutate(rsid = as.character(snp)) %>%
  dplyr::rename(pval = P.value.y)
df_sub_eqtl <- df_sub %>%
  dplyr::select(snp, P.value.x) %>%
  dplyr::mutate(rsid = as.character(snp)) %>%
  dplyr::rename(pval = P.value.x)

# plot significant colocalization
plot <- locuscompare(in_fn1=df_sub_eqtl,
                     in_fn2=df_sub_mqtl,
                     title1="Brain cortex eQTL",
                     title2="Plasma heptanoate mQTL",
                     # title2="QTL",
                     snp=top_var,
                     genome="hg19",
                     population = "EUR")
plot


# save plot
ggsave(filename = paste0("SPTSSB_", QTL, "_locuscompare_hg19_rs1450522_ref_case-control.pdf"),
       plot = plot + theme(legend.position = "none") +
         geom_text_repel(max.overlaps = 1000000), # Increase the limit
       device = "pdf",
       width = 12,
       height = 6)


#________________________________________________________________________________
### 2. SMR/HEIDI Mendelian Randomization ----
setwd("/mnt/belinda_local/ruth/data/omics_integration/MR_analysis")

# Outcome: Clinical PD status; 
# Instruments: GWAS significant SNPs; 
# Exposure: mQTLs or eQTLs

# MR pipeline (from SMR):
# Github of source code: https://github.com/NIH-CARD/NDD_SMR/blob/main/scripts/SMR_v8.ipynb
# Program download: https://yanglab.westlake.edu.cn/software/smr/#Download
# Instructions: https://yanglab.westlake.edu.cn/software/smr/#SMR&HEIDIanalysis
# Data visualization tool: https://yanglab.westlake.edu.cn/smr-portal/viewer (very similar to LocusZoom)


#### Formatting GWAS summary stats ----
# Basically, GWAS file mandatory columns should be:
# "rsid","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value","N_datasets"

# gwas_all <- read.table("../colocalization_analysis/SPTSSB_gwas_locus_for_coloc-MR_hg19.txt", header=T)
# gwas_sig <- read.table("../colocalization_analysis/SPTSSB_gwas_locus_for_coloc-MR_0_00000005_hg19.txt", header=T)
gwas_all <- read.table("../colocalization_analysis/SPTSSB_gwas_locus_for_coloc-MR_hg19_case-control.txt", header=T)
gwas_sig <- read.table("../colocalization_analysis/SPTSSB_gwas_locus_for_coloc-MR_0_00000005_hg19_case-control.txt", header=T)

mygwas <- gwas_sig

# Optional - Replacing empty rsid fields (those are not SNVs but are still short mutations)
mygwas[mygwas$SNP_ID %in% "chr3:160958484:A:AAT",]$rsid <- "rs34371992"  # effect allele:other allele
mygwas[mygwas$SNP_ID %in% "chr3:161093790:C:CA",]$rsid <- "rs5853931"  # effect allele:other allele
mygwas[mygwas$SNP_ID %in% "chr3:161096040:T:TTTG",]$rsid <- "rs149705199"  # effect allele:other allele
mygwas[mygwas$SNP_ID %in% "chr3:160600563:G:A",]$rsid <- "rs935497"  # effect allele:other allele


mygwas <- mygwas[,c("rsid","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value","N_datasets")]
colnames(mygwas) <- c("SNP","A1","A2","freq","b","se","p","N")

#### Formatting QTL summary stats (only necessary for Surendran mQTLs and ROSMAP single-cell eQTLs) ----
smr_path <- "/mnt/belinda_70t/ruth_data/QTLs/for_SMR_coloc_graphs/smr-1.3.1-linux-x86_64/smr"

SPTSSB_eqtl_2 <- read.table("/mnt/belinda_local/ruth/data/omics_integration/colocalization_analysis/Surendran_mGWAS_SPTSSB_data/SPTSSB_PD_GWAS_var_+-_500kb_heptanoate_only.txt", sep = '\t', header = T, stringsAsFactors = TRUE)
# SPTSSB_eqtl_2 <- read.table("/mnt/belinda_local/ruth/data/omics_integration/colocalization_analysis/Surendran_mGWAS_SPTSSB_data/extract_SPTSSB_160577630_161577630_M52285.txt", sep = '\t', header = T, stringsAsFactors = TRUE) # oleate-vaccenate
# SPTSSB_eqtl_2 <- read.table("/mnt/belinda_70t/ruth_data/QTLs/celltype-eqtl-sumstats.Exc.SPTSSB_locus_genes.tsv", header = T)  # ROSMAP snRNA-seq for excitatory neurons

# For single-cell eQTLs
SPTSSB_eqtl_2 <- SPTSSB_eqtl_2[,c(3:12)]
colnames(SPTSSB_eqtl_2) <- c("gene","rsid","CHR","pos.hg38","Allele2","Allele1","Freq1","Effect","StdErr","P.value")
SPTSSB_eqtl_1  <- read.table("/mnt/belinda_local/ruth/data/omics_integration/colocalization_analysis/Surendran_mGWAS_SPTSSB_data/sptssb_heptanoate_500kb+-coordinates_hg38_conversion_w_rsids.txt", sep = '\t', header = T, stringsAsFactors = TRUE)[,c(2,4)]
colnames(SPTSSB_eqtl_1)[1] <- "Pos"
SPTSSB_eqtl_2 <- merge(SPTSSB_eqtl_1, SPTSSB_eqtl_2, by="pos.hg38")


colnames(SPTSSB_eqtl_2)[2] <- "Pos"

# For mQTLs other than heptanoate - rescue rsids
rsids <- read.table("/mnt/belinda_local/ruth/data/omics_integration/colocalization_analysis/Surendran_mGWAS_SPTSSB_data/SPTSSB_PD_GWAS_var_+-_500kb_heptanoate_only.txt", sep = '\t', header = T, stringsAsFactors = TRUE)[,c("Pos","rsid")]
SPTSSB_eqtl_2 <- merge(rsids, SPTSSB_eqtl_2, by="Pos")


SPTSSB_eqtl_2$Allele1 <- toupper(SPTSSB_eqtl_2$Allele1)
SPTSSB_eqtl_2$Allele2 <- toupper(SPTSSB_eqtl_2$Allele2)

SPTSSB_eqtl <- SPTSSB_eqtl_2

# Only for mQTLs
SPTSSB_eqtl$gene <- "ENSG00000196542"  # SPTSSB


SPTSSB_eqtl$FDR <- p.adjust(SPTSSB_eqtl$P.value, method = "BH")
SPTSSB_eqtl$t_stat <- SPTSSB_eqtl$Effect / SPTSSB_eqtl$StdErr

SPTSSB_eqtl <- SPTSSB_eqtl[,c("rsid","gene","Effect","t_stat","P.value","FDR")]
colnames(SPTSSB_eqtl) <- c("SNP","gene","beta","t.stat","p.value","FDR")

# Remove columns with either empty or duplicated rsids
SPTSSB_eqtl <- SPTSSB_eqtl[!(SPTSSB_eqtl$SNP %in% ""),]
SPTSSB_eqtl <- SPTSSB_eqtl[!(SPTSSB_eqtl$SNP %in% SPTSSB_eqtl$SNP[duplicated(SPTSSB_eqtl$SNP)]), ] # if there are no duplicated values, result will be empty

# For Surendran mQTLs only - Need to flip beta signal (because they flipped effect and reference alleles)
SPTSSB_eqtl$beta <- -(SPTSSB_eqtl$beta)

write.table(SPTSSB_eqtl, "heptanoate_SPTSSB_for_SMR_flipped.txt", row.names=F, col.names=T, quote=F, sep="\t")

# For single-cell eQTLs only - Need to save eQTLs per gene
temp <- SPTSSB_eqtl[SPTSSB_eqtl$gene %in% "ENSG00000196542",]
write.table(temp, "excitatory_neurons_SPTSSB_for_SMR.txt", row.names=F, col.names=T, quote=F, sep="\t")


##### Creating .besd files (only necessary for Surendran mQTLs/single-cell eQTLs) ----
# Instructions: https://yanglab.westlake.edu.cn/software/smr/#DataManagement

## Make a BESD file from Matrix QTL output

# First, ----for mQTLs----, we need to properly name the genes where the SNVs fall onto
# gene_rescue <- unique(read.table("VEP_SPTSSB_locus_SNV-gene_rescue_Surendran.txt")[,c(1,6:7)])
# gene_rescue <- gene_rescue[grep("^ENSG", gene_rescue$V7), ]  # only over half of SNPs get rescued...

PPM1L_coord <- c(160473995, 160788817) # start, end
B3GALNT1_coord <- c(160801670, 160823160)
NMD3_coord <- c(160939098, 160969795)
SPTSSB_coord <- c(161062579, 161089871)
OTOL1_coord <- c(161214595, 161221730)

# Gene names and coordinates as a list
genes <- list(
  PPM1L = list(coord = PPM1L_coord, gene_id = "ENSG00000163590"),
  B3GALNT1 = list(coord = B3GALNT1_coord, gene_id = "ENSG00000169255"),
  NMD3 = list(coord = NMD3_coord, gene_id = "ENSG00000169251"),
  SPTSSB = list(coord = SPTSSB_coord, gene_id = "ENSG00000196542"),
  OTOL1 = list(coord = OTOL1_coord, gene_id = "ENSG00000182447")
)

# Load mQTL data
mQTL <- read.table("heptanoate_SPTSSB_for_SMR_flipped.txt", header = TRUE)
mQTL$gene <- NA  # Initialize gene column

# Assign genes to SNPs
for (i in 1:nrow(mQTL)) {
  rsid <- mQTL[i, "SNP"]
  pos <- SPTSSB_eqtl_2[SPTSSB_eqtl_2$rsid %in% rsid,]$Pos
  
  # Calculate distances to all genes
  distances <- sapply(genes, function(gene) {
    if (pos >= gene$coord[1] && pos <= gene$coord[2]) {
      return(0)  # SNP falls within the gene
    } else if (pos < gene$coord[1]) {
      return(gene$coord[1] - pos)  # Distance to gene start
    } else {
      return(pos - gene$coord[2])  # Distance to gene end
    }
  })
  
  # Assign the nearest gene
  nearest_gene <- names(distances)[which.min(distances)]
  mQTL[i, "gene"] <- genes[[nearest_gene]]$gene_id
}

write.table(mQTL, "heptanoate_SPTSSB_for_SMR_flipped.txt", row.names = F, col.names = T, sep="\t", quote = F)

# Create BESD itself
cmd <- paste0(smr_path, " --eqtl-summary heptanoate_SPTSSB_for_SMR_flipped.txt --matrix-eqtl-format --make-besd --out heptanoate_SPTSSB_for_SMR_flipped ")
# cmd <- paste0(smr_path, " --eqtl-summary excitatory_neurons_SPTSSB_for_SMR.txt --matrix-eqtl-format --make-besd --out excitatory_neurons_SPTSSB_for_SMR ")

cmd
system(cmd)

# After make BESD file command, I realized the.esi and .epi files were not properly formatted.
# Need to re-format it, so it looks like this:
# # .esi:
# chromosome  SNP genetic_distance(can be 0) position  effect_allele other_allele  freq
# 1    rs1001  0   744055  A   G   0.23
# 1    rs1002  0   765522  C   G   0.06
# 1    rs1003  0   995669  T   C   0.11
# ......
# # .epi (only used for graphs):
# chromosome  probeID(exon or transc) genetic_distance(can be 0)  physical_position geneID  gene_strand
# 1    probe1001   0   924243  Gene01  +
# 1    probe1002   0   939564  Gene02  -
# 1    probe1003   0   1130681 Gene03  -
# ......
# Gene range list: https://www.cog-genomics.org/static/bin/plink/glist-hg38

esi_file <- read.table("heptanoate_SPTSSB_for_SMR_flipped.esi")

positions <- SPTSSB_eqtl_2[,c("rsid","Pos","Allele1","Allele2","Freq1")]
# positions <- positions[!(positions$rsid %in% ""),]
positions <- positions[!(positions$rsid %in% positions$rsid[duplicated(positions$rsid)]), ]
positions <- unique(positions)

esi_file$V1 <- 3
esi_file$V4 <- positions$Pos

# When alleles are NOT flipped
esi_file$V5 <- positions$Allele1
esi_file$V6 <- positions$Allele2
esi_file$V7 <- positions$Freq1

# When alleles are flipped
esi_file$V5 <- positions$Allele2
esi_file$V6 <- positions$Allele1
esi_file$V7 <- 1-positions$Freq1

write.table(esi_file, "heptanoate_SPTSSB_for_SMR_flipped.esi", row.names = F, col.names = F, sep = "\t", quote = F)

# Rescue info for genes of interest
epi_file <- read.table("heptanoate_SPTSSB_for_SMR_flipped.epi")

# Rescue only for the genes with eQTL within the .besd file
epi_file[1,] <- c(3, "ENSG00000163590", 0, 160631406, "PPM1L", "+")  # coordinate in the middle of the gene
epi_file[2,] <- c(3, "ENSG00000169255", 0, 160812415, "B3GALNT1", "-")
epi_file[3,] <- c(3, "ENSG00000169251", 0, 160954447, "NMD3", "+")
epi_file[4,] <- c(3, "ENSG00000196542", 0, 161076225, "SPTSSB", "-")
epi_file[5,] <- c(3, "ENSG00000182447", 0, 161218162, "OTOL1", "+")
write.table(epi_file, "heptanoate_SPTSSB_for_SMR_flipped.epi", row.names = F, col.names = F, sep = "\t", quote = F)


#### Running MR itself (SMR and HEIDI test) ----
# Instructions: https://yanglab.westlake.edu.cn/software/smr/#SMR&HEIDIanalysis

binaries_path <- "/mnt/belinda_local/ruth/home/plink_binaries/g1000_eur_hg19"
# gwas_path <- "./gwas_SPTSSB_for_SMR_hg19.ma"
gwas_path <- "./gwas_SPTSSB_for_SMR_hg19_case-control.ma"
heptanoate_path <- "./heptanoate_SPTSSB_for_SMR_flipped"
# oleate_path <- "./oleate-vaccenate_SPTSSB_for_SMR"
nerve_path <- "/mnt/belinda_70t/ruth_data/QTLs/for_GTEx_coloc_graphs/Nerve_Tibial/Nerve_Tibial"
hypothalamus_path <- "/mnt/belinda_70t/ruth_data/QTLs/for_GTEx_coloc_graphs/Brain_Hypothalamus/Brain_Hypothalamus"
cortex_path <- "/mnt/belinda_70t/ruth_data/QTLs/for_SMR_coloc_graphs/BrainMeta_cis_eqtl_summary/BrainMeta_cis_eQTL_chr3"
blood_path <- "/mnt/belinda_70t/ruth_data/QTLs/for_GTEx_coloc_graphs/Whole_Blood/Whole_Blood"
# excitatory_path <- "./excitatory_neurons_SPTSSB_for_SMR_all_genes"
excitatory_path <- "./excitatory_neurons_SPTSSB_for_SMR"

# Multi-SNP-based SMR test
# Below shows an option to combine the information from all the SNPs in a region that pass a p-value threshold 
# (the default value is 5.0e-8 which can be modified by the flag --peqtl-smr)
# Note that SMR always considers GWAS phenotype as outcome, so it is unclear whether it can test reverse causation...
# But for the version with 2 molecular traits, it is indeed possible to test reverse causation by flipping the order the exposure and outcome
# summary stats are presented in the command line.
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", heptanoate_path, " --gwas-summary ", gwas_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_PD_GWAS_plasma_heptanoate_FLIPPED_mQTLs_multiSNP_0.05_case-control ")
# cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", oleate_path, " --gwas-summary ", gwas_path, " oleate-vaccenate_SPTSSB_for_SMR --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_PD_GWAS_plasma_oleate-vaccenate_mQTLs_multiSNP_0.05 ")

cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", nerve_path, " --gwas-summary ", gwas_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_PD_GWAS_nerve_tibial_SPTSSB_eQTLs_multiSNP_0.05 ") # genome build needs to match!
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", hypothalamus_path, " --gwas-summary ", gwas_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_PD_GWAS_hypothalamus_SPTSSB_eQTLs_multiSNP_0.05 ") # genome build needs to match!
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", cortex_path, " --gwas-summary ", gwas_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_PD_GWAS_cortex_SPTSSB_eQTLs_multiSNP_0.05_case-control ") # genome build needs to match!
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", blood_path, " --gwas-summary ", gwas_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_PD_GWAS_blood_SPTSSB_eQTLs_multiSNP_0.05 ") # genome build needs to match!

cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", excitatory_path, " --gwas-summary ", gwas_path, " excitatory_neurons_SPTSSB_for_SMR_all_genes --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_PD_GWAS_excitatory_neurons_SPTSSB_all_genes_eQTLs_multiSNP_0.05_case-control ")

cmd
system(cmd)

# --bfile (plink binaries) can be replaced by --bld, which reads LD information from a binary file in BLD format
# You may add "--maf 0.01" to filter SNPs by MAF (in the reference sample) - OPTIONAL
# You may run with multiple threads: "--thread-num 10" - OPTIONAL

# SMR analysis of two molecular traits
# Here we provide an option to test the pleiotropic association between two molecular traits using summary data.
# Exposure should come first, outcome should come second

# eQTL -> mQTL
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", nerve_path, " --beqtl-summary ", heptanoate_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_nerve_to_heptanoate_multiSNP_0.05")
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", hypothalamus_path, " --beqtl-summary ", heptanoate_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_hypothalamus_to_heptanoate_multiSNP_0.05")
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", cortex_path, " --beqtl-summary ", heptanoate_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_cortex_to_heptanoate_multiSNP_0.05")
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", blood_path, " --beqtl-summary ", heptanoate_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_blood_to_heptanoate_multiSNP_0.05")
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", excitatory_path, " --beqtl-summary ", heptanoate_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_excitatory_to_heptanoate_multiSNP_0.05_SPTSSB_only")

# cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", nerve_path, " --beqtl-summary ", oleate_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_nerve_to_oleate-vaccenate_multiSNP_0.05")
# cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", hypothalamus_path, " --beqtl-summary ", oleate_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_hypothalamus_to_oleate-vaccenate_multiSNP_0.05")
# cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", cortex_path, " --beqtl-summary ", oleate_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_cortex_to_oleate-vaccenate_multiSNP_0.05")
# cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", blood_path, " --beqtl-summary ", oleate_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_blood_to_oleate-vaccenate_multiSNP_0.05")

cmd
system(cmd)

# mQTL -> eQTL (reverse causation)
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", heptanoate_path, " --beqtl-summary ", nerve_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_heptanoate_to_nerve_multiSNP_0.05")
# cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", oleate_path, " --beqtl-summary ", nerve_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_oleate_to_nerve_multiSNP_0.05")
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", heptanoate_path, " --beqtl-summary ", hypothalamus_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_heptanoate_to_hypothalamus_multiSNP_0.05")
# cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", oleate_path, " --beqtl-summary ", hypothalamus_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_oleate_to_hypothalamus_multiSNP_0.05")
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", heptanoate_path, " --beqtl-summary ", cortex_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_heptanoate_to_cortex_multiSNP_0.05")
# cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", oleate_path, " --beqtl-summary ", cortex_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_oleate_to_cortex_multiSNP_0.05")
cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", heptanoate_path, " --beqtl-summary ", blood_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_heptanoate_to_blood_multiSNP_0.05")
# cmd <- paste0(smr_path, " --bfile ", binaries_path, " --beqtl-summary ", oleate_path, " --beqtl-summary ", blood_path, " --peqtl-smr 0.05 --peqtl-heidi 0.05 --smr-multi --out smr_oleate_to_blood_multiSNP_0.05")

cmd
system(cmd)

# Data visualization tool: https://yanglab.westlake.edu.cn/smr-portal/viewer
# Formatting files for plotting: https://yanglab.westlake.edu.cn/software/smr/#SMRlocusplot19
