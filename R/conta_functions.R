# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' @useDynLib conta
#' @importFrom Rcpp sourceCpp
#' @import data.table parallel ggplot2
NULL

#' Return nucleotide bases
#' @param context_mode whether to get bases from 3-base context
#'
#' @export
get_bases <- function(context_mode = FALSE) {
  if (!context_mode)
    return(c("A", "T", "G", "C"))
  else {
    bases <- expand.grid(get_bases(), get_bases(), get_bases())
    return(paste(bases$Var1, bases$Var2, bases$Var3, sep = ""))
  }
}

#' Test whether the two contexts are compatible with a one base variance in the
#' middle base
#' @param context1 context of the first case
#' @param context2 context of the second case
#'
#' @export
is_compatible_context <- function(context1, context2) {
  if (substring(context1, 1, 1) == substring(context2, 1, 1) &&
      substring(context1, 3, 3) == substring(context2, 3, 3))
    return(TRUE)
  else
    return(FALSE)
}

#' Return initial test range
#'
#' @export
get_initial_range <- function() {
  return(c(1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3,
           5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1))
}

#' Return an empty results table to output when quitting early
#' @param sample sample name
#' @importFrom utils packageVersion
#' @export
empty_result <- function(sample) {
  vals <- list(conta_version = as.character(packageVersion("conta")),
               sample = sample)
  na_names <- c("conta_call", "cf", "sum_log_lr", "avg_log_lr", "snps", "depth",
                "pos_lr_all", "pos_lr_x", "pos_lr_chr_cv", "y_count",
                "y_norm_count", "y_frac", "pregnancy",
                "excluded_regions", "error_rate",
                "T>A", "G>A", "C>A", "A>T", "G>T", "C>T",
                "A>G", "T>G", "C>G", "A>C", "T>C", "G>C")
  na_vals <- rep(NA_character_, length(na_names))
  names(na_vals) <- na_names
  empty_result <- do.call(data.table, c(vals, na_vals))
  return(empty_result)
}

#' Fail if certain conditions are not met
#'
#' @param dat data.table with SNP counts and genotypes
#'
#' @export
fail_test <- function(dat) {
  if ( is.null(dat) || nrow(dat) == 0 || is.null(dat$gt) ||
       nrow(dat[gt == "0/0", ]) == 0 ||
       nrow(dat[gt == "1/1", ]) == 0 ||
       nrow(dat[gt == "0/1", ]) == 0) {
    msg <- "Either no SNPs passed filters or not all genotypes were called."
    if (interactive()) {
      stop(msg)
    } else {
      message(msg)
      quit(save = "no", status = 2, runLast = FALSE)
    }
  }
}

#' Return expected contamination fraction
#'
#' Average ratio of heterozygotes to homzoygotes for contaminant SNPs between
#' 100 individual samples was determined to be 80%. For these SNPs, only one
#' allele (half the fraction) contributes to the contamination. For the
#' remaining 20% of SNPs which are homozygotes, both alleles contribute to the
#' contamination. This 80% number is an approximation and may be changed based
#' on prior information.
#'
#' @param cf numeric contamination fraction
#' @param het_rate ratio of heterozygote to homozygote contaminated alleles
#' @return expected contamination fraction
#'
#' @export
get_exp_cf <- function(cf, het_rate = 0.8) {
  return(het_rate * (cf / 2) + (1 - het_rate) * cf)
}

#' Simulate contamination
#'
#' It works by adding in each allele proportional to specified
#' contamination fraction and maf. So it actually simulates a
#' sample randomly ignoring linkage disequilibrium.
#'
#' @param wgs data.frame containing counts and metrics per SNP
#' @param cf numeric contamination fraction
#'
#' @return data.frame containing counts and metrics per SNP
#'
#' @importFrom stats rpois runif
#' @export
sim_conta <- function(wgs, cf) {

  if (cf == 0) return(wgs)

  # Calculate allele depths for each chromosome
  allele1 <- ifelse(runif(nrow(wgs)) > wgs$maf, wgs$major, wgs$minor)
  allele2 <- ifelse(runif(nrow(wgs)) > wgs$maf, wgs$major, wgs$minor)

  # Add in contamination
  bases <- get_bases()
  for (b in bases) {
    locs1 <- which(allele1 == b)
    allele1_probs <- wgs[locs1, "depth"] * cf / 2
    allele1_draws <- rpois(allele1_probs[, .N], allele1_probs[, depth])
    wgs[locs1, b] <- wgs[locs1, b, with = FALSE][[1]] + allele1_draws

    locs2 <- which(allele2 == b)
    allele2_probs <- wgs[locs2, "depth"] * cf / 2
    allele2_draws <- rpois(allele2_probs[, .N], allele2_probs[, depth])
    wgs[locs2, b] <- wgs[locs2, b, with = FALSE][[1]] + allele2_draws
  }

  # Recalculate depth, ratio and counts
  wgs <- ratio_and_counts(wgs)

  return(wgs)
}

#' Simulate LOH and/or contamination
#'
#' This simulation returns a simple data structure with metrics required to
#' test LOH functions. It works by simulating 4 chromosomes (2 for host, and
#' 2 for contaminant). It also simulates the ratio of the host chromosomes which
#' is analogous to loss of heterozygosity.
#'
#' @param n number of SNPs to simulate
#' @param min_maf minimum maf to simulate (and maximum is 1 - maf)
#' @param dp_min min epth to simulate
#' @param dp_max max depth to simulate
#' @param er_min minimum error rate to simulate
#' @param er_max maximum error rate to simulate
#' @param delta LOH deviation from heterozygosity
#' @param alpha contamination level
#' @param seed numeric set seed for simulation
#'
#' @return data.frame containing counts and metrics per SNP
#'
#' @export
simulate_loh_conta <- function(n, min_maf, dp_min, dp_max, er_min, er_max,
                               delta, alpha, seed = 1359) {
  set.seed(seed)

  # Add maf check
  if (min_maf > 0.5)
    min_maf <- 0.5

  # Simulate maf for each SNP (at range 25 to 75%)
  maf <- runif(n, min = min_maf, max = 1 - min_maf)

  # Simulate error rate for each SNP randomly at specified range
  er <- runif(n, min = er_min, max = er_max)

  # Simulate depth for each SNP randomly around specified target
  dp <- rpois(n, runif(n, dp_min, dp_max))

  # Simulate alleles for host
  a1 <- ifelse(runif(n) > maf, 1, 0)
  a2 <- ifelse(runif(n) > maf, 1, 0)

  # Simulate alleles for contaminant
  a3 <- ifelse(runif(n) > maf, 1, 0)
  a4 <- ifelse(runif(n) > maf, 1, 0)

  # Simulate depth for host
  dp1 <- rpois(n, dp * (1 - alpha) * (0.5 - delta / 2))
  dp2 <- rpois(n, dp * (1 - alpha) * (0.5 + delta / 2))

  # Simulate depth for contaminant
  dp3 <- rpois(n, dp * alpha / 2)
  dp4 <- rpois(n, dp * alpha / 2)

  # Draw minor allele depths for host
  ad1 <- rpois(n, dp1 * abs(a1 - er))
  ad1_b <- rpois(n, dp1 * (1 - abs(a1 - er)))
  ad2 <- rpois(n, dp2 * abs(a2 - er))
  ad2_b <- rpois(n, dp2 * (1 - abs(a2 - er)))

  # Draw minor allele depths for contaminant
  ad3 <- rpois(n, dp3 * abs(a3 - er))
  ad3_b <- rpois(n, dp3 * (1 - abs(a3 - er)))
  ad4 <- rpois(n, dp4 * abs(a4 - er))
  ad4_b <- rpois(n, dp4 * (1 - abs(a4 - er)))

  # Sum the ads
  ad <- ad1 + ad2 + ad3 + ad4
  dp <- ad + ad1_b + ad2_b + ad3_b + ad4_b

  return(data.table(depth = dp, minor_count = ad, maf = maf, er = er))
}

#' Load conta files
#'
#' For a given record, mt and snps, read vfn or gt file, and return it.
#'
#' @param conta_loc location of conta vfn or gt file
#' @param snps targeted dbsnp data.table
#' @return data.table with conta genotypes
#'
#' @export
load_conta_file <- function(conta_loc, snps = NULL) {

  if ((!is.null(snps)) & !is.na(conta_loc) &
      (file.exists(conta_loc) ||
       (requireNamespace("grails3r") && grails3r::s3_file_exists(conta_loc)) ||
       as.logical(do.call(aws.s3::head_object,
                          c(file, aws.signature::locate_credentials()))))) {

    # Read targeted genotypes and rsid from vfn file format
    conta_loc_dt <- read_data_table(conta_loc, sep = ",", showProgress = FALSE)

    # Add rsid if missing.
    if (!is.null(conta_loc_dt))
      conta_loc_dt[, rsid := snps$ID]

  } else if (!is.na(conta_loc) &
             (file.exists(conta_loc) ||
              (requireNamespace("grails3r") &&
               grails3r::s3_file_exists(conta_loc)) ||
              as.logical(do.call(aws.s3::head_object,
                                 c(file, locate_credentials()))))) {

    # Read conta genotypes and rsid from gt file format
    conta_loc_dt <- read_data_table(conta_loc, showProgress = FALSE)

    # Vf is the variant allele frequency, needs to be set of wgs files
    conta_loc_dt[, vf := vr / dp]

  }
  else {

    conta_loc_dt <- NULL

  }

  return(conta_loc_dt)
}

#' Calculate concordance between two samples' genotypes
#'
#' TODO (edamato): support gt files from other pipelines
#'
#' @param dt1 data.table of genotype set 1
#' @param dt2 data.table of genotype set 2
#'
#' @export
genotype_concordance <- function(dt1, dt2) {

  # Merge the two data tables to put each SNP in a single row
  m1 <- merge(dt1, dt2, by = "rsid")

  # Calculate genotype concordance which is the fraction of matching genotypes
  return(m1[, mean(gt.x == gt.y, na.rm = TRUE)])
}

#' Get folders under a specified s3 directory
#'
#' Returns the set of folders under a given s3 path. get_bucket_df normally
#' returns all folders and files recursively, this method parses get_bucket_df
#' output to return only folders directly under the specified dir.
#'
#' @param s3_path to display the contents
#' @export
get_s3_folders <- function(s3_path) {

  # Stop if this is not an s3 path
  if (!startsWith(s3_path, "s3://"))
    stop("get_s3_ls() requires a valid s3 path")

  # Retireve bucket text
  bucket_text <- strsplit(s3_path, "/")[[1]][3]

  # Retrieve prefix text
  pre_text <- paste(strsplit(s3_path, "/")[[1]][c(-1, -2, -3)],
                    sep = "/", collapse = "/")

  # Get all files (and folders) under specified bucket/prefix
  bc_all_files <- aws.s3::get_bucket_df(bucket = bucket_text,
                                        prefix = pre_text, max = Inf)

  # Number of folders under prefix
  pl <- length(strsplit(pre_text, "/")[[1]])

  # Get all paths directly under prefix (omit files and folder in deeper levels)
  paths <- unique(sapply(strsplit(bc_all_files$Key, "/"),
                  function(x) {
                    if (length(x) > (pl + 1)) x[[pl + 1]]
                    }))

  # Add a slash in the end before returning to denote a folder
  return(paste(unlist(paths), "/", sep = ""))
}

#' Set numeric equivalent of chromosomes for sorting purposes
#'
#' @section TODO: Handle chromosomes outside X, Y and M/MT
#'
#' @param dat data.table to add numeric chromosome
#' @export
set_numeric_chrs <- function(dat) {

  suppressMessages(require(dplyr))

  datt <- dat %>%
    mutate(
      chrom_int = substring(chrom, 4)) %>%
    mutate(
      chrom_int = case_when(
        chrom_int == "X" ~ 23,
        chrom_int == "Y" ~ 24,
        chrom_int == "M" ~ 25,
        chrom_int == "MT" ~ 25,
        TRUE ~ as.numeric(chrom_int))
    )
  return(data.table(datt))
}
