#' @useDynLib conta
#' @importFrom Rcpp sourceCpp
#' @import data.table parallel ggplot2 grails3r
NULL

#' Return nucleotide bases
#'
#' @export
get_bases <- function() {
  return(c("A", "T", "G", "C"))
}

#' Return initial test range
#'
#' @export
get_initial_range <- function() {
  return(c(1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3,
           5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1))
}

#' Return an empty results table to output when quitting early
#'
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
    msg <- "Either no SNPs passed filters or no genotypes were called."
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
#' Average ratio of heterozygotes for informative (contaminant) SNPs between
#' 100 individual samples was determined to be 80%. For these SNPs, only one
#' allele (half the fraction) contributes to the contamination. For the
#' remaining 20% of SNPs which are homozygotes, both alleles contribute to the
#' contamination.
#'
#' @param cf numeric contamination fraction
#' @return expected contamination fraction
#'
#' @export
get_exp_cf <- function(cf) {
  return(0.8 * cf / 2 + 0.2 * cf)
}

#' Simulate contamination
#'
#' It works by adding in each allele proportional to specified
#' contamination fraction and maf. So it actually simulates a
#' sample randomly ignoring linkage disequilibrium.
#'
#' @param wgs data.frame containing counts and metrics per SNP
#' @param cf numeric contamination fraction
#' @param seed numeric set seed for simulation
#'
#' @return data.frame containing counts and metrics per SNP
#'
#' @export
sim_conta <- function(wgs, cf, seed = 1359) {

  if (cf == 0) return(wgs)

  set.seed(seed)

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

#' Simulate loh and/or contamination
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
#' @param delta loh deviation from heterozygosity
#' @param alpha contamination level
#' @param seed numeric set seed for simulation
#'
#' @return data.frame containing counts and metrics per SNP
#'
#' @export
simulate_loh_conta <- function(n, min_maf, dp_min, dp_max, er_min, er_max,
                               delta, alpha, seed = 1359) {

  set.seed(seed)

  # Simulate maf for each SNP (at range 25 to 75%)
  maf <- runif(n, min = min_maf, max = 1 - min_maf)

  # Simulate error rate for each SNP randomly at specified range
  er <- runif(n, min = er_min, max = er_max)

  # Simulate depth for each SNP randomly around specified target
  dp <- rpois(n, runif(n, dp_min, dp_max))

  # Simulate alleles for each chromosome
  a1 <- ifelse(runif(n) > maf, 1, 0)
  a2 <- ifelse(runif(n) > maf, 1, 0)

  # Simulate contamination in case necessary
  a3 <- ifelse(runif(n) > maf, 1, 0)
  a4 <- ifelse(runif(n) > maf, 1, 0)

  # Draw allelic depths based on parameters
  # Allele 1 and Allele 2 contribute 50% unless there is LOH, in that case
  # Allele 1 contributes a proportion of specified by delta (loh coefficient)
  ad1 <- rpois(n, a1 * dp * (0.5 - delta) * (1 - alpha) * abs(a1 - er))
  ad2 <- rpois(n, a2 * dp * (0.5 + delta) * (1 - alpha) * abs(a2 - er))

  # Allele 3 and 4 contributes only if there is contamination (alpha)
  ad3 <- rpois(n, a3 * dp * alpha * abs(a3 - er))
  ad4 <- rpois(n, a4 * dp * alpha * abs(a4 - er))

  # Sum the ads
  ad <- ad1 + ad2 + ad3 + ad4

  # Set contamination probs based on maf
  cp <- ifelse(a1 == 1 & a2 == 1, 1 - maf ^ 2,
               ifelse(a1 == 0 & a2 == 0, 1 - (1 - maf) ^ 2,
                      1 - 2 * maf * (1 - maf)))

  return(data.table(depth = dp, minor_count = ad, maf = maf, er = er, cp = cp))
}
