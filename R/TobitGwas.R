regress_snp <- function(geno, pheno, covars, detection_limit, hde_test = FALSE) {
  coefficients <-
    VGAM::vglm(
      pheno ~ geno + covars,
      tobit(Lower=detection_limit)) %>%
    VGAM::summaryvglm(model, HDEtest = hde_test) %>%
    coef()
  coefficients["geno", ]
}

regress_all <- function(geno, pheno, covars, detection_limit, hde_test = FALSE) {
  result <- geno %>%
    parallel::mclapply(
      FUN = regress,
      pheno = data_f$oest.norm,
      covars = covars,
      detection_limit = detection_limit,
      hde_test = hde_test,
      mc.cores = 8) %>%
    do.call("rbind", .) %>%
    as.data.table(keep.rownames="SNP")

  snp_split <- stringr::str_split_fixed(result$SNP, "_", 2)
  result[, ID := snp_split[,1]]
  result[, A1 := snp_split[,2]]
  result[, SNP := NULL]
  data.table::setcolorder(result, c("ID", "A1"))
  result
}


uninteresting_columns <- c("MAT", "PAT", "SEX", "PHENOTYPE")
extract_genotype_columns <- function(genotype_table){
  interesting_columns <- !(colnames(genotype_table) %in% uninteresting_columns)
  genotype_table[, ..interesting_columns]
}

collect_genotypes <- function(files) {
  files %>%
    lapply(data.table::fread) %>%
    lapply(extract_genotype_columns) %>%
    Reduce(f=merge)
}
