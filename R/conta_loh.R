# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' Log likelihood ratio for loss of heterozygosity hypothesis.
#'
#' At a given level of delta which is the allelic imbalance or LOH level as an
#' absolute difference from 0.5, calculate the likelihood of observing the
#' data given a specific ad and depth.
#'
#' @param ad alternative allele depth for a given SNP
#' @param depth total depth for a given SNP
#' @param maf minor allele frequencies for host sample
#' @param er error rate for a given SNP
#' @param delta loss of heterozygosity level to be tested
#' @param blackswan blackswan term sets a limit on probability for each event
#' @return log likelihood ratio for LOH vs no LOH
#'
#' @export
llr_loh <- function(ad, depth, maf, er, delta, blackswan) {
  l_min <- blackswan / depth
  l_het <- l_min + (1 - l_min) * l_loh(ad, depth, maf, maf, 0, 0, er)
  l_loh <- l_min + (1 - l_min) * l_loh(ad, depth, maf, maf, delta, 0, er)
  return(log(l_loh / l_het))
}

#' Avg of log likelihood ratios for loss of heterozygosity.
#'
#' Avg log likelihood ratio across a given set of SNPs
#'
#' @param dat data table containing required fields
#' @param delta loss of heterozygosity level to be tested
#' @param blackswan blackswan term sets a limit on probability for each event
#' @return sum of log likelihood ratios for LOH vs no LOH
#'
#' @export
avg_llr_loh <- function(dat, delta, blackswan) {
  weighted.mean(llr_loh(dat$minor_count, dat$depth, dat$maf, dat$er,
                       delta, blackswan), dat$depth, na.rm = T)
}


#' Log likelihood ratio of contamination
#'
#' At a given level of cf (contamination fraction), calculate
#' the likelihood of observing the data given ad and depth.
#'
#' @param ad alternative allele depth for a given SNP
#' @param depth total depth for a given SNP
#' @param maf1 minor allele frequencies for host sample
#' @param maf2 minor allele frequencies for candidate contaminant
#' @param er error rate for a given SNP
#' @param alpha contamination fraction to be tested
#' @param cp contamination probability calculated from minor allele frequency
#' @param blackswan blackswan term sets a limit on probability for each event
#' @return log likelihood ratio for contamination vs no contamination
#'
#' @export
llr_cont <- function(ad, depth, maf1, maf2, er, alpha, delta, blackswan) {
  l_min <- blackswan / depth
  l_het <- l_min + (1 - l_min) * l_loh(ad, depth, maf1, maf2, delta, 0, er)
  l_cont <- l_min + (1 - l_min) * l_loh(ad, depth, maf1, maf2, delta, alpha, er)
  return(log(l_cont / l_het))
}


#' Avg log likelihood ratios for contamination
#'
#' Average log likelihood ratio across a given set of SNPs
#'
#' @param dat data table containing required fields
#' @param alpha contamination level level to be tested
#' @param blackswan blackswan term sets a limit on probability for each event
#' @return sum of log likelihood ratios for LOH vs no LOH
#'
#' @export
avg_llr_cont <- function(dat, alpha, blackswan, delta = 0) {
  weighted.mean(llr_cont(dat$minor_count, dat$depth, dat$maf, dat$maf,
                        dat$er, alpha, delta, blackswan),
                dat$depth, na.rm = T)
}

#' Likelihood of observing the data with LOH
#'
#' Calculates the posterior probability of observing a given alternative allele
#' depth and total depth given a contamination level (alpha) and a loss of
#' heterozygosity deviation level (delta). Prior probabilities on the host and
#' contaminate genotypes are calculated using the given minor allele frequencies
#' for the host (maf1) and contaminate (maf2).
#'
#' @param ad alternative allele depth for a given SNP
#' @param dp total depth for a given SNP
#' @param maf1 minor allele frequencies for host sample
#' @param maf2 minor allele frequencies for candidate contaminant
#' @param d LOH level as deviation from 0.5
#' @param a alpha contamination level
#' @param er error rate
#'
#' @importFrom stats dbinom
#' @importFrom stats dnbinom
#' @export
l_loh <- function(ad, dp, maf1, maf2, d = 0, a = 0, er = 3e-4) {
  pg2(0, 0, maf1, maf2) * dbinom(ad, dp, er) +
    pg2(0, 1, maf1, maf2) * dbinom(ad, dp, er + a / 2) +
    pg2(0, 2, maf1, maf2) * dbinom(ad, dp, er + a) +
    pg2(1, 0, maf1, maf2) * (0.5 * dnbinom(ad, dp - ad, 0.5 - d / 2 - a / 2) +
                             0.5 * dnbinom(ad, dp - ad, 0.5 + d / 2 - a / 2)) +
    pg2(1, 1, maf1, maf2) * (0.25 * dnbinom(ad, dp - ad, 0.5 - d / 2 - a / 2) +
                             0.25 * dnbinom(ad, dp - ad, 0.5 + d / 2 - a / 2) +
                             0.25 * dnbinom(ad, dp - ad, 0.5 - d / 2 + a / 2) +
                             0.25 * dnbinom(ad, dp - ad, 0.5 + d / 2 + a / 2)) +
    pg2(1, 2, maf1, maf2) * (0.5 * dnbinom(ad, dp - ad, 0.5 - d / 2 + a / 2) +
                             0.5 * dnbinom(ad, dp - ad, 0.5 + d / 2 + a / 2)) +
    pg2(2, 0, maf1, maf2) * dbinom(ad, dp, 1 - er - a) +
    pg2(2, 1, maf1, maf2) * dbinom(ad, dp, 1 - er - a / 2) +
    pg2(2, 2, maf1, maf2) * dbinom(ad, dp, 1 - er)
}

#' Prior probability of genotype pair
#'
#' @param gt1 Genotype of first sample's SNP
#'            0, 1, or 2 corresponding to hom reference, het, hom alternative
#' @param gt2 Genotype of second sample's SNP
#'            0, 1, or 2 corresponding to hom reference, het, hom alternative
#' @param maf1 minor allele frequency first sample's SNP
#' @param maf2 minor allele frequency of second sample's SNP
#' @export
pg2 <- function(gt1, gt2, maf1, maf2) {
  pg(gt1, maf1) * pg(gt2, maf2)
}

#' Prior probability of genotype
#'
#' @param gt 0, 1, or 2 corresponding to hom reference, het, hom alternative
#' @param maf minor allele frequency
#' @export
pg <- function(gt, maf) {
  ifelse(gt == rep(0, length(maf)), (1 - maf) * (1 - maf),
          ifelse(gt == rep(1, length(maf)), 2 * (1 - maf) * maf,
                  ifelse(gt == rep(2, length(maf)), maf * maf, NA)))
}

#' Calculate per bin (large chromosomal segments) stats
#'
#' @param dat data.table containing counts and metrics per SNP
#' @param save_dir save directory
#' @param sample sample name
#' @param min_lr likelihood ratio threshold to call LOH
#' @param blackswan blackswan term sets a limit on probability for each event
#' @param min_loh to worry about affecting contamination results,
#'     if LOH is less than this fraction do not call
#' @param min_snps minimum number of snps required to evaluate a region for LOH
#' @param max_snps max snps to evaluate in a region, if greater number of
#'     SNPs are present, subsample to this many SNPs.
#' @param seed random seed
#'
#' @return data.table containing stats per bin
#'
#' @importFrom utils write.table
#' @export
get_per_bin_loh <- function(dat, save_dir, sample, min_lr, blackswan,
                            min_loh = 0.15, min_snps = 20, max_snps = 200,
                            min_auto_loh = 0.4, seed = 1359) {

  results <- data.table(chrom = character(), chunk = numeric(),
                        snps = numeric(),
                        loh_mle = numeric(), loh_val = numeric(),
                        cont_mle = character(), cont_val = numeric(),
                        loh = logical())

  # combs contains combination of chrom and chunk
  combs <- expand.grid(unique(dat$chrom), unique(dat$chunk))
  combs <- combs[order(combs[, 1]), ]

  for (i in 1:nrow(combs)) {

    # Get the data in this specific chunk
    dat_chunk <- dat[chrom == combs[i, 1] & chunk == combs[i, 2], ]

    # Set dat as a subset of dat if it exceeds a pre-determined size
    if (nrow(dat_chunk) > max_snps) {
      set.seed(seed)
      dat_chunk <- dat_chunk[sort(sample(1:nrow(dat_chunk), max_snps)), ]
    }

    # Optimize LOH starting on a grid (if there are enough SNPs)
    loh_opt <- optimize_after_grid(dat_chunk, "avg_llr_loh", blackswan)

    # If loh_opt is not null, its objective (avg. likelihood ratio) is above
    # a specified threshold, its level is above a minimum LOH level, then
    # calculate contamination mle for the same region. Do not call a region as
    # LOH otherwise, also do not call if there aren't enough SNPs to access.
    if (!is.na(loh_opt$objective)
        && dat_chunk[, .N] >= min_snps
        && loh_opt$objective >= min_lr
        && loh_opt$maximum >= min_loh) {

      cont_opt <- optimize_after_grid(dat_chunk, "avg_llr_cont", blackswan)

      dfres <- data.table(chrom = combs[i, 1], chunk = combs[i, 2],
                          snps = dat_chunk[, .N],
                          loh_mle = loh_opt$maximum,
                          loh_val = loh_opt$objective,
                          cont_mle = cont_opt$maximum,
                          cont_val = cont_opt$objective,
                          loh = (loh_opt$maximum >= min_auto_loh ||
                                   (loh_opt$objective > cont_opt$objective &
                                      loh_opt$objective > min_lr)))
    } else {

      dfres <- data.table(chrom = combs[i, 1], chunk = combs[i, 2],
                          snps = dat_chunk[, .N],
                          loh_mle = loh_opt$maximum,
                          loh_val = loh_opt$objective,
                          cont_mle = NA, cont_val = NA,
                          loh = FALSE)
    }

    results <- rbind(results, dfres)
  }

  write.table(results,
              paste(save_dir, "/", sample, ".loh_regions.tsv", sep = ""),
              sep = "\t", quote = F, col.names = T, row.names = F)

  return(results)
}


#' optimize further after grid search
#` It first performs grid search, and find-search further using optimize()
#'
#' @param dat_chunk data.table containing counts and metrics per SNP
#' @param LR_func Log-likelihood ratio function
#' @param blackswan blackswan term sets a limit on probability for each event
#' @return maxima and its objective
#'
#' @importFrom stats optimize
#' @export
optimize_after_grid <- function(dat_chunk, LR_func, blackswan) {

  # Get initial range
  cf <- get_initial_range()

  grid_lr <- lapply(cf, function(x) match.fun(LR_func)(dat_chunk, x,
                                                       blackswan = blackswan))

  # Cf that gives the max result from the initial grid search
  vmax <- which.max(grid_lr)

  # If vmax is all NAs it will return a logical vector of length zero, in that
  # case set it as the 2nd element to optimize in the minimum range.
  if (length(vmax) == 0)
    return(list(maximum = NA, objective = NA))

  # If max cf is the last or the first value, shift to prevent out of bounds
  vmax <- ifelse(vmax == length(cf), vmax - 1, ifelse(vmax == 1, 2, vmax))

  opt_val <- optimize(match.fun(LR_func), lower = cf[vmax - 1],
                      upper = cf[vmax + 1], tol = 0.001, maximum = T,
                      dat = dat_chunk, blackswan = blackswan)

  return(opt_val)
}


#' Remove regions with high LOH
#'
#' @param dat data.table SNP data
#' @param bin_stats data.table LOH stats for each bin
#' @return data.table data with SNPs on LOH regions removed
#'
#' @export
exclude_high_loh_regions <- function(dat, bin_stats) {
  remove_bin <- bin_stats[loh == TRUE]
  if (nrow(remove_bin) > 0) {
    remove_bin_id <- paste(remove_bin$chrom, "-", remove_bin$chunk)
    dat_ids <- paste(dat$chrom, "-", dat$chunk)
    dat <- dat[!(dat_ids %in% remove_bin_id), ]
  }
  return(dat)
}
