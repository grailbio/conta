#' @useDynLib conta
#' @importFrom Rcpp sourceCpp
#' @import data.table parallel ggplot2 grails3r
NULL

#' Clean given object and perform garbage collection
#'
#' @export
clean <- function(object) {
  rm(object)
  gc()
}

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
