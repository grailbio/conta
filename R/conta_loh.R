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
#' @param maf minor allele frequency for a given SNP
#' @param er error rate for a given SNP
#' @param delta loss of heterozygosity level to be tested
#' @param blackswan blackswan term sets a limit on probability for each event
#' @return log likelihood ratio for LOH vs no LOH
#'
#' @export
lr_loh <- function(ad, depth, maf, er, delta, blackswan) {
  min_lh <- blackswan / depth
  phet <- min_lh + (1 - min_lh) * p_loh(ad, depth, maf, 0, 0, er)
  ploh <- min_lh + (1 - min_lh) * p_loh(ad, depth, maf, delta, 0, er)
  return(log(ploh / phet))
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
avg_lr_loh <- function(dat, delta, blackswan) {
  mean(lr_loh(dat$minor_count, dat$depth, dat$maf, dat$er,
              delta, blackswan), na.rm = T)
}


#' Log likelihood ratio of contamination
#'
#' At a given level of cf (contamination fraction), calculate
#' the likelihood of observing the data given ad and depth.
#'
#' @param ad alternative allele depth for a given SNP
#' @param depth total depth for a given SNP
#' @param maf minor allele frequency for a given SNP
#' @param er error rate for a given SNP
#' @param alpha contamination fraction to be tested
#' @param cp contamination probability calculated from minor allele frequency
#' @param blackswan blackswan term sets a limit on probability for each event
#' @return log likelihood ratio for contamination vs no contamination
#'
#' @export
lr_cont <- function(ad, depth, maf, er, alpha, cp, blackswan) {
  min_lh <- blackswan / depth
  phet <- min_lh + (1 - min_lh) * p_loh(ad, depth, maf, 0, 0, er)
  pcont <- min_lh + (1 - min_lh) * p_cont(ad, depth, maf, er, alpha, cp)
  return(log(pcont / phet))
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
avg_lr_cont <- function(dat, alpha, blackswan) {
  mean(lr_cont(dat$minor_count, dat$depth, dat$maf,
               dat$er, alpha, dat$cp, blackswan), na.rm = T)
}

#' Probability of observing the data with LOH
#'
#' This function tests the probability of observing a given number of alt
#' reads at a SNP position for a given delta (LOH) and alpha (contamination)
#' When delta = 0 and alpha = 0, we are essentially testing the probability
#' of observing the data from a normal site (without LOH or contamination).
#'
#' The probabilities are conditioned on the genotype. Note the genotype priors
#' add up to 1. Genotype priors can be one of 3 categories:
#' 1) (1-maf) * (1-maf) is the probability of observing a ref/ref
#' 2) 2 * (1-maf) * maf is the probability of observing a ref/alt, and
#' 3) maf * maf is the probability of observing an alt/alt
#'
#' Error rate er is a small number, and since the application here
#' is LOH detection, it should matter only for low coverage (<20) sites where
#' heterozygosity can be achieved from errors on a homozygote site.
#'
#'
#' @param ad alternative allele depth for a given SNP
#' @param depth total depth for a given SNP
#' @param maf minor allele frequency for a given SNP
#' @param delta loh level as deviation from 0.5
#' @param alpha contamination level
#' @param er error rate
#'
#' @importFrom stats dbinom
#' @export
p_loh <- function(ad, depth, maf, delta = 0, alpha = 0, er = 3e-4) {
  (1 - maf) * (1 - maf) * dbinom(ad, depth, er + alpha) +
    2 * (1 - maf) * maf * (0.5 * dbinom(ad, depth, 0.5 - delta)
                           + 0.5 * dbinom(ad, depth, 0.5 + delta)) +
    maf * maf * dbinom(ad, depth, 1 - er - alpha)
}

#' Probability of observing the data with contamination.
#'
#' First term signifies the probability of observing the data from a SNP that
#' has the same genotypes as the host, whereas the second part signifies the
#' probability of observing the data from an actual contaminated site.
#' When there is a contamination, it will either manifest as a homozygote allele
#' coming into an heterozygote site, or as an heterozygote or opposite
#' homozygote allele coming into the host homozygote site.
#'
#' @param ad alternative allele depth for a given SNP
#' @param depth total depth for a given SNP
#' @param maf minor allele frequency for a given SNP
#' @param er error rate for a given SNP
#' @param alpha contamination fraction
#' @param cp contamination probability

#' @export
p_cont <- function(ad, depth, maf, er, alpha, cp) {
  (1 - cp) * p_loh(ad, depth, maf, 0, 0, er) +
    cp * p_loh(ad, depth, maf, alpha / 2, get_exp_cf(alpha), er)
}

#' Calculate per bin (large chromosomal segments) stats
#'
#' @param dat data.table containing counts and metrics per SNP
#' @param save_dir save directory
#' @param sample sample name
#' @param cutoff threshold to call loh
#' @param blackswan blackswan term sets a limit on probability for each event
#' @param min_loh to worry about, if loh is less than this fraction, it is ok
#' @param min_snps minimum number of snps required to evaluate a region for loh
#' @return data.table containing stats per bin
#'
#' @importFrom utils write.table
#' @export
get_per_bin_loh <- function(dat, save_dir, sample, cutoff, blackswan,
                            min_loh = 0.05, min_snps = 10) {

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

    # Optimize loh starting on a grid (if there are enough SNPs)
    loh_opt <- optimize_after_grid(dat_chunk, "avg_lr_loh", blackswan)

    # If loh_opt is not null, its objective (avg. likelihood ratio) is above
    # a specified threshold, its level is above a minimum loh level, then
    # calculate contamination mle for the same region. Do not call a region as
    # LOH otherwise, also do not call if there aren't enough SNPs to access.
    if (!is.na(loh_opt$objective)
        && dat_chunk[, .N] >= min_snps
        && loh_opt$objective >= cutoff
        && loh_opt$maximum >= min_loh) {

      cont_opt <- optimize_after_grid(dat_chunk, "avg_lr_cont", blackswan)

      dfres <- data.table(chrom = combs[i, 1], chunk = combs[i, 2],
                          snps = dat_chunk[, .N],
                          loh_mle = loh_opt$maximum,
                          loh_val = loh_opt$objective,
                          cont_mle = cont_opt$maximum,
                          cont_val = cont_opt$objective,
                          loh = loh_opt$objective > cont_opt$objective &
                            loh_opt$objective > cutoff)
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
