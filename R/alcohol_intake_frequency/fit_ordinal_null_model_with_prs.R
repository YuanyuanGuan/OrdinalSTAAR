##########################################################
# fit ordinal null model with regenie binary multitrait
# Yuanyuan Guan
# 19/06/2025
##########################################################

# sinfo
# srun -p maths-gpu -w c90 --mem=300G --cpus-per-task=8 --time=24:00:00 --pty bash
# conda activate ukbb
# R

library(data.table)
library(Matrix)
library(stringr)
library(STAAR)
library(GENESIS)
library(GMMAT)
library(STAARpipeline)
library(SeqVarTools)
library(tidyr)
library(dplyr)
library(ordinal)
library(ordinalSTAAR)

# step1: generate full data ------------------------------------------------

load("/datapool/maths/UKBB/Data/Phenotype/main_survey/bd_676623.Rdata")

# PC1-PC20
pcs <- fread("/datapool/maths/UKBB/Data/PC_GRM/sGRM_and_PCs/output.pca.score")
colnames(pcs)[1:2] <- c("userId", "userID2")
colnames(pcs)[3:22] <- paste0("PC", 1:20)

# Phenotype
pheno <- bd[, c(c(0, 769) + 1)]
colnames(pheno) <- c("userId", "alcohol_intake_frequency")

# Covariates
covariate <- bd[, c(c(0, 22, 10318) + 1)]
colnames(covariate)[1:3] <- c("userId", "sex", "age")

# merge PCs, phenotype and covariates
t1 <- merge(pheno, covariate, by = "userId")
fullDat <- merge(t1, pcs, by = "userId")

# data filter
gId <- get(load("/datapool/maths/UKBB/Data/Phenotype/sample_id/UKB_500K_WGS_sampleid.RData"))
fullDat <- fullDat[
  (fullDat$userId %in% gId) &
    !is.na(fullDat$alcohol_intake_frequency) &
    !is.na(fullDat$PC1) & !is.na(fullDat$age) & !is.na(fullDat$sex),
]

# data transformation
fullDat$age2 <- (fullDat$age)^2
fullDat <- fullDat[fullDat$alcohol_intake_frequency != "Prefer not to answer", ]
fullDat$alcohol_intake_frequency <- ordered(
  factor(
    fullDat$alcohol_intake_frequency,
    levels = c("Never", "Special occasions only",
               "One to three times a month", "Once or twice a week",
               "Three or four times a week", "Daily or almost daily")
  )
)
fullDat$alcohol_ordinal <- as.integer(fullDat$alcohol_intake_frequency)

ncat <- max(fullDat$alcohol_ordinal, na.rm=TRUE) # 6
for(i in 1:(ncat - 1)) {
  alcohol_bin <- paste0("alcohol_bin_", i)
  fullDat[[alcohol_bin]] <- ifelse(fullDat$alcohol_ordinal <= i, 1, 0)
}

# data standardization
cols_to_scale <- c("age", "age2", paste0("PC", 1:10))
fullDat[, (cols_to_scale) := lapply(.SD, scale), .SDcols = cols_to_scale]

save(fullDat, file="/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/WGS_binary_alcohol_intake_frequency_fullDat.20250601.Rdata")
# load("/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/WGS_binary_alcohol_intake_frequency_fullDat.20250601.Rdata")
# load("D:/desktop/multiclass/Nullmodel/alcohol_intake_frequency/LOCO/WGS_alcohol_intake_frequency_fullDat.20250619.Rdata")

# step2: generate regenie input files -------------------------------------

# generate phenotype and covariates
pheno_cols <- c("userId", "userId", paste0("alcohol_bin_", 1:5))
pheno_file <- fullDat[, ..pheno_cols]
setnames(pheno_file, old = 1:2, new = c("FID", "IID"))
fwrite(pheno_file,
       file = "/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/phenotype_alcohol_intake_frequency.20250601.txt",
       sep = "\t",
       row.names = FALSE)

covar_cols <- c("userId", "userId", "sex", "age", "age2", paste0("PC", 1:10))
covar_file <- fullDat[, ..covar_cols]
setnames(covar_file, old = 1:2, new = c("FID", "IID"))
fwrite(covar_file,
       file = "/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/covariates.20250517.txt",
       sep = "\t",
       row.names = FALSE)

# merge phenotype and covariates
pheno_covar <- merge(pheno_file, covar_file, by = c("FID", "IID"))
pheno_covar$FID_IID <- paste0(pheno_covar$FID, "_", pheno_covar$IID)

# step3: regenie ----------------------------------------------------------

# conda activate regenie_env
# regenie \
# --step 1 \
# --bed /datapool/maths/UKBB/UKB_genotypes_QCed_ZHL/chrall \
# --phenoFile /datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/phenotype_alcohol_intake_frequency.20250601.txt \
# --covarFile /datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/covariates.20250517.txt \
# --covarColList sex,age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
# --catCovarList sex \
# --bt \
# --bsize 1000 \
# --out /datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination \
# --threads 8 \
# --print-prs

# Step 4: Read and Reshape WIDE format PRS files, then Merge ----------------------------------------------------------

output_path <- "/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/binary_regenie_combination/"

prs_reshaped_list <- list()
ncat <- max(fullDat$alcohol_ordinal, na.rm=TRUE) # 6
ncat <- 6

# Loop through each binary phenotype to read and reshape its PRS file
for (i in 1:(ncat - 1)) {
  # --- 1. Construct file name ---
  prs_file_name <- paste0(output_path, "binary_regenie_combination_", i, ".prs")

  if (file.exists(prs_file_name)) {

    # --- 2. Read the WIDE format PRS file ---
    prs_wide_data <- fread(prs_file_name, header = TRUE, data.table = FALSE)
    colnames(prs_wide_data)[1] <- "chr"

    prs_long_data <- prs_wide_data %>%
      pivot_longer(
        cols = -chr,
        names_to = "FID_IID",
        values_to = paste0("prs_bin_", i)
      ) %>%
      select(FID_IID, !!paste0("prs_bin_", i))

    prs_reshaped_list[[i]] <- prs_long_data

    message(paste("Successfully read and reshaped PRS for alcohol_bin_", i))

  } else {
    warning(paste("PRS file not found:", prs_file_name))
  }
}

# --- Step 5: Merge all reshaped PRS results ----------------------------------------------------------
fullDat <- fullDat %>%
  select(
    userId, alcohol_intake_frequency, sex, age, age2,
    paste0("PC", 1:10)
  )
fullDat$FID_IID <- paste(fullDat$userId, fullDat$userId, sep = "_")

merged_prs_data <- Reduce(
  function(x, y) full_join(x, y, by = "FID_IID"),
  prs_reshaped_list
)

fullDat <- fullDat %>%
  left_join(merged_prs_data, by = "FID_IID") %>%
  select(-FID_IID)

all_pred_cols <- paste0("prs_bin_", 1:(ncat - 1))
all_pred <- fullDat %>% select(all_of(all_pred_cols))

# Handle NAs
na_count <- sum(is.na(all_pred))
if(na_count > 0) {
  warning(paste(na_count, "NA values found in the merged PRS columns. Removing affected rows."))
  fullDat <- fullDat[complete.cases(all_pred), ]
  all_pred <- all_pred[complete.cases(all_pred), ]
}


# Step 6: PCA -------------------------------------------------------------
if (any(!is.finite(as.matrix(all_pred)))) {
  warning("Non-finite values found in PRS columns. Removing affected rows.")
  finite_rows <- rowSums(is.finite(as.matrix(all_pred))) == ncol(all_pred)
  fullDat <- fullDat[finite_rows, ]
  all_pred <- all_pred[finite_rows, ]
}

pca_result <- prcomp(all_pred, center = TRUE, scale. = TRUE)

keep_num <- 5
pca_components <- as.data.frame(pca_result$x[, 1:keep_num])
colnames(pca_components) <- paste0("prs_pc", 1:keep_num)
# head(pca_components)

data_for_null_model <- bind_cols(fullDat, pca_components)
save(data_for_null_model,
     file = "/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/latent_pheno/WGS_alcohol_intake_frequency_fullDat_with_PRS_PCs.20250910.Rdata")


# --- Step 7: Fit Final Ordinal Null Model -----------------------------------------------------------
load("D:\\desktop\\Proj\\Multiclass\\Nullmodel\\alcohol_intake_frequency\\PRS\\WGS_alcohol_intake_frequency_fullDat_with_PRS_PCs.20250619.Rdata")
load("/mnt/project/UKB_500K_WGS_staarpipeline/Multiclass/alcohol_intake_frequency/WGS_alcohol_intake_frequency_fullDat_with_PRS_PCs.20250619.Rdata")

base_covar_cols <- c("sex", "age", "age2", paste0("PC", 1:10))
prs_pc_cols <- paste0("prs_pc", 1:5)

outcome_column <- "alcohol_intake_frequency"
id_column <- "userId"

obj.STAAR.UKB.alcohol_intake_frequency <- NullModel_Ordinal(
  phenofile = data_for_null_model,
  outcomeCol = outcome_column,
  sampleCol = id_column,
  covCol = base_covar_cols,
  PRSCol = prs_pc_cols,
  verbose = TRUE
)

save(
  obj.STAAR.UKB.alcohol_intake_frequency,
  file = "D:/desktop/multiclass/Nullmodel/alcohol_intake_frequency/PRS/obj.STAAR.UKB.alcohol_intake_frequency.20250909.Rdata"
)

message("--- Analysis Complete ---")
