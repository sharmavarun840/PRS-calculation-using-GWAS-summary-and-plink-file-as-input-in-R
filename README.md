# PRS-calculation-using-GWAS-summary-and-plink-file-as-input-in-R and using principle components startifaction as confounding factor for PRS 
The calculation of polygenic risk scores from GWAS summary statistics and genotype data in PLINK format typically involves several steps, including the following:

Load the summary statistics file into R and extract the necessary columns, such as the SNP ID, effect size, and p-value.

# Load necessary packages
library(data.table)

# Read in the GWAS summary statistics file
gwas_summary <- fread("gwas_summary_statistics.txt")

# Extract the necessary columns from the summary statistics file
snp_id <- gwas_summary$SNP
beta <- gwas_summary$BETA
pval <- gwas_summary$PVAL


Load the genotype data in PLINK format into R using the plink() function from the plink package and data visulisation.

# Load the plink package
library(plink)

# Load the genotype data in PLINK format
geno_data <- plink("genotype_data")

Use the snp_id vector from the summary statistics file to extract the corresponding genotypes from the genotype data. This can be done using the snp_data() function from the plink package.

# Extract the genotypes for the SNPs in the snp_id vector
snp_genotypes <- snp_data(geno_data, id = snp_id)


Calculate the polygenic risk score for each individual by summing the product of the effect size and genotype for each SNP.

# Calculate the polygenic risk score for each individual
prs <- rowSums(snp_genotypes * beta)


Optionally, adjust the polygenic risk scores for population stratification using a method such as principal component analysis.

# Load the necessary packages
library(prcomp)

# Perform principal component analysis on the genotype data
pca_result <- prcomp(snp_genotypes)

# Adjust the polygenic risk scores for population stratification
prs_adj <- prs - (pca_result$x[, 1:10] %*% solve(pca_result$rotation[1:10, 1:10])) %*% pca_result$rotation[1:10, ]

Please note that this is just an example of how to calculate polygenic risk scores in R, and there may be other approaches or considerations to take into account depending on the specific data and analysis goals.

