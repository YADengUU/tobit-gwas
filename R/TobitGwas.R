#!Rscript

library(magrittr)
library(data.table)
library(parallel)

regress <- function(geno, pheno, covars, lower) {
  VGAM::vglm( pheno ~ geno + covars, tobit(Lower=lower))
}

genotypes <- list.files(path="genotypes", pattern="raw", full.names=TRUE) %>%
  lapply(fread) %>%
  lapply(extract, , c(1,2,7)) %>%
  Reduce(f=merge)

covariates_f <- fread("~/hormone_gwas/females/covariates.txt")

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
pheno <- fread("pheno.txt")
data_f <- merge(covars_f_transformed, pheno)

# Keep individuals with complete data, only
common_ids_f <- intersect(genotypes$IID, data_f$IID)
genotypes_f <- genotypes[IID %in% common_ids_f]
data_f <- data_f[IID %in% common_ids_f]

# The normalized detection limit is the highest normalized estradiol value that was undetectable
detection_limit <- data_f[oest == 175, max(oest.norm)]

# Only use variables for regression that make sense
not_predictors <- c("IID", "FID", "oest", "oest.norm")
predictor_cols_f <- !(colnames(data_f) %in% not_predictors)
predictors_f <- data_f[, ..predictor_cols_f] %>% as.matrix()

# Regression happens here
models <- genotypes_f[, -c(1,2)] %>% mclapply(regress, data_f$oest.norm, predictors_f, detection_limit, mc.cores=8)
summaries <- mclapply(models, VGAM::summaryvglm, HDEtest=FALSE, mc.cores=8)
result <- lapply(summaries, coef) %>% lapply(extract, "geno", ) %>% do.call("rbind", .) %>% as.data.table(keep.rownames="SNP")

# Make the results pretty
snp_split <- str_split_fixed(result$SNP, "_", 2)
result[, ID := snp_split[,1]]
result[, A1 := snp_split[,2]]

# Write results to disk
write.table(result, "out/results_f.txt", sep="\t", quote=F, row.names=F)

# Do the same stuff for men
covariates_m <- fread("~/hormone_gwas/males/covariates.txt")
covars_m_transformed <- cbind(covariates_m, model.matrix( ~ SMOKE  - 1, covariates_m)[,-1])
covars_m_transformed[, SMOKE := NULL]

common_ids_m <- intersect(pheno$IID, covariates_m$IID) %>% intersect(genotypes$IID)
data_m <- merge(covars_m_transformed, pheno)[IID %in% common_ids]

#  The normalized detection limit is different for men b.c. the distribution is different
detection_limit_m <- data_m[oest == 175, max(oest.norm)]
genotypes_m <- genotypes[IID %in% common_ids_m]

predictor_cols_m <- !(colnames(data_m) %in% not_predictors)
predictors_m <- data_m[, ..predictor_cols_m] %>% as.matrix()

# Regression
models_m <- mclapply(genotypes_m[, -c(1,2)], regress, data_m$oest.norm, predictors_m, detection_limit_m, mc.cores=8)
summaries_m <- mclapply(models_m, summary, mc.cores=8)
result_m <- lapply(summaries_m, coef) %>% lapply(extract, "geno", ) %>% do.call("rbind", .) %>% as.data.table(keep.rownames="SNP")

# Prettification
snp_split_m <- stringr::str_split_fixed(result_m$SNP, "_", 2)
result_m[, ID := snp_split_m[,1]]
result_m[, A1 := snp_split_m[,2]]
write.table(result_m, "out/results_m.txt", sep="\t", quote=F, row.names=F)
