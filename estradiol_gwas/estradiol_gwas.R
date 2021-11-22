library(tobitGwas)
library(data.table)
library(magrittr)
library(stringr)

genotypes <- read_genotype_directory("genotypes")

covariates_f <- data.table::fread("~/hormone_gwas/females/covariates.txt")

# Create dummy variables for categorical data
covars_f_transformed <- covariates_f %>%
  cbind(., model.matrix( ~ HRT - 1, .)[,-1]) %>%
  cbind(., model.matrix( ~ OC - 1, .)[,-1]) %>%
  cbind(., model.matrix( ~ SMOKE - 1, .)[,-1])

# Delete categorical variables
covars_f_transformed[, HRT := NULL]
covars_f_transformed[, OC := NULL]
covars_f_transformed[, SMOKE := NULL]
covars_f_transformed[, HYSTER := as.logical(HYSTER)]

# Load phenotype and merge with covariates
pheno <- data.table::fread("pheno.txt")

# Keep individuals with complete data, only
common_ids_f <- genotypes$IID %>%
  intersect(covars_f_transformed$IID) %>%
  intersect(phenoIID)
genotypes_f <- genotypes[IID %in% common_ids_f]
covars_f_transformed <- covars_f_transformed[IID %in% common_ids_f]
pheno_f <- pheno[IID %in% common_ids_f]

# Ensure ordering
setkey(genotypes_f, IID)
setkey(covars_f_transformed, IID)
setkey(pheno_f, IID)

# The normalized detection limit is the highest normalized estradiol value that was undetectable
detection_limit <- pheno_f[oest == 175, max(oest.norm)]

# Only use variables for regression that make sense
not_predictors <- c("IID", "FID")
predictor_cols_f <- !(colnames(covars_f_transformed) %in% not_predictors)
predictors_f <- covars_f_transformed[, ..predictor_cols_f] %>% as.matrix()

# Regression happens here

result <- regress_all(
  geno = genotypes_f[, -c(1,2)],
  pheno = pheno_f$oest.norm,
  covars = predictors_f,
  detection_limit = detection_limit,
  hde_test = FALSE
)

# Write results to disk
write.table(result, "out/results_f.txt", sep="\t", quote=F, row.names=F)

# Do the same stuff for men
covariates_m <- data.table::fread("~/hormone_gwas/males/covariates.txt")
covars_m_transformed <- cbind(covariates_m, model.matrix( ~ SMOKE  - 1, covariates_m)[,-1])
covars_m_transformed[, SMOKE := NULL]

common_ids_m <- pheno$IID %>%
  intersect(covariates_m$IID) %>%
  intersect(genotypes$IID)

genotypes_m <- genotypes[IID %in% common_ids_m]
covars_m_transformed <- covars_m_transformed[IID %in% common_ids_m]
pheno_m <- pheno[IID %in% common_ids_m]

#  The normalized detection limit is different for men b.c. the distribution is different
detection_limit_m <- pheno_m[oest == 175, max(oest.norm)]

predictor_cols_m <- !(colnames(data_m) %in% not_predictors)
predictors_m <- data_m[, ..predictor_cols_m] %>% as.matrix()

# Regression

result_m <- regress_all(
  geno = genotypes_m[, -c(1,2)],
  pheno = pheno_m$oest.norm,
  covars = predictors_m,
  detection_limit = detection_limit,
  hde_test = FALSE
)

write.table(result_m, "out/results_m.txt", sep="\t", quote=F, row.names=F)
