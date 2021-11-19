#' Run Tobit Regression on a single SNP
#'
#' Runs a Tobit-I censored regression on the provided SNP.
#'
#' @param geno A numeric vector of dosages of one SNP for all individuals
#' @param pheno A numeric vector of phenotypes for all individuals.
#' Measurements below detection limit must be set to a value less than or equal detection limit.
#' @param covars N x M covariate matrix (rows are individuals)
#' @param detection_limit The detection limit of the phenotype
#' @param hde_test If TRUE, a test for the Hauck-Donner effect will be performed.
#' The result is currently not returned
#'
#' @return A named vector containing estimate, standard error, z and p values of the regression.
#' @export
regress_snp <- function(geno, pheno, covars, detection_limit, hde_test = FALSE) {
  coefficients <-
    VGAM::vglm(
      pheno ~ geno + covars,
      tobit(Lower=detection_limit)) %>%
    VGAM::summaryvglm(model, HDEtest = hde_test) %>%
    coef()
  coefficients["geno", ]
}

#' Perform a Tobit regression on all SNPs
#'
#' Performs a censored-regression association study on all provided SNPs.
#' This is parallelized by default using [parallel::mclapply].
#'
#' @param geno Genotype matrix of the cohort. Rows are individuals
#' @param pheno A numeric vector of phenotypes for all individuals.
#' Measurements below detection limit must be set to a value less than or equal detection limit.
#' @param covars N x M covariate matrix (rows are individuals)
#' @param detection_limit The detection limit of the phenotype
#' @param hde_test If TRUE, a test for the Hauck-Donner effect will be performed.
#' The result is currently not returned
#' @param mc.cores The number of cores to use for parallelization, see [parallel::mclapply].
#'
#' @return A [data.table::data.table] containing regression results for all SNPs
#' @export
regress_all <- function(geno, pheno, covars, detection_limit, hde_test = FALSE, mc.cores = getOption("mc.cores", 2L)) {
  result <- geno %>%
    parallel::mclapply(
      FUN = regress,
      pheno = data_f$oest.norm,
      covars = covars,
      detection_limit = detection_limit,
      hde_test = hde_test,
      mc.cores = mc.cores) %>%
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


#' Extracts columns from raw genotype files that we actually need for regression
#'
#' @param genotype_table The genotype file parsed as a [data.table::data.table]
#'
#' @return A [data.table::data.table] with only IDs and dosages remaining
extract_genotype_columns <- function(genotype_table){
  interesting_columns <- !(colnames(genotype_table) %in% uninteresting_columns)
  genotype_table[, ..interesting_columns]
}


#' Implements an outer join using [base::merge].
#'
#' We need this because [base::Reduce] can't pass additional arguments to its function.
#'
#' @return An outer join of a and b
outer_join <- function(a,b){
  merge(a,b, all=TRUE)
}

#' Reads raw genotype files from PLINK and combines them into a single genotype matrix.
#'
#' @param files A vector of filenames to read
#'
#' @return A [data.table::data.table] of dosages for individuals (rows) and SNPs (columns).
#' Individuals that are missing from some files will have their dosages set to NA for the SNPs in these files.
#'
#' @export
read_genotypes <- function(files) {
  files %>%
    lapply(data.table::fread) %>%
    lapply(extract_genotype_columns) %>%
    Reduce(f=outer_join())
}

#' Reads a directory of raw genotype files.
#'
#' @param directory The directory to read.
#'
#' @return A [data.table::data.table] of dosages for individuals (rows) and SNPs (columns).
#' Individuals that are missing from some files will have their dosages set to NA for the SNPs in these files.
#'
#' @export
read_genotype_directory <- function(directory) {
  list.files(path=directory, pattern="raw", full.names = T) %>% read_genotypes()
}

#' Rank Transformation for Tobit Regression
#'
#' Rank transforms values that do not fit a normal distribution to fit a
#' normal distribution and calculates the detection limit of the transformed data.
#'
#' @param values The raw values to transform
#' @param detection_limit The detection limit of the method used to measure the values
#' @param na_is_below_dt Indicates whether all NAs measurements below the detection limit.
#'
#' @details
#' If `detection_limit` is not given, the minimum observed value is assumed
#' to be the detection limit
#'
#' If `na_is_below_dt` is TRUE, all NAs are assumed to be measurements below
#' detection and will be set to the detection limit before transformation.
#' If FALSE, NAs are not modified and it is assumed that all measurements
#' below detection limit have been set to a value less than or equal to the
#' detection limit.
#'
#' @return A named list with two elements.
#' `values` is a vector of transformed measurements.
#' `detection_limit` is the transformed detection limit.
#' @export
rank_transform <- function(values, detection_limit = min(values, na.rm = TRUE), na_is_below_dt=FALSE){
  if (na_is_below_dt){
    values[is.na(values)] <- detection_limit
  }
  defined <- values[!is.na(values)]
  tmp.1 <- (order(defined)-0.5)/length(defined)
  tmp.2 <- rep(NA,length(values))
  tmp.2[!is.na(values)] <- tmp.1
  transformed <- qnorm(tmp.2,0,1)
  transformed_dt <- transformed[values <= detection_limit] %>%
    max(na.rm = TRUE)
  list(values = transformed, detection_limit = transformed_dt)
}
