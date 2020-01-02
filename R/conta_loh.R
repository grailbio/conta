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
#' @param hm het mean frequency from baseline
#' @param M het alpha + beta parameters for beta binomial alt formulation
#' @param er error rate for a given SNP
#' @param delta loss of heterozygosity level to be tested
#' @param blackswan blackswan term sets a limit on probability for each event
#' @return log likelihood ratio for LOH vs no LOH
#'
#' @export
llr_loh <- function(ad, depth, maf, hm, M, er, delta, blackswan) {
  l_min <- blackswan / depth
  l_het <- l_min + (1 - l_min) * l_loh(ad, depth, maf, maf, hm, M, 0, 0, er)
  l_loh <- l_min + (1 - l_min) * l_loh(ad, depth, maf, maf, hm, M, delta, 0, er)
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
  weighted.mean(llr_loh_all(dat, delta, blackswan), log(dat$depth),
                na.rm = TRUE)
}

#' Log likelihood ratios for loss of heterozygosity for a set of snps
#'
#' Log likelihood ratio across a given set of SNPs
#'
#' @param dat data table containing required fields
#' @param delta loss of heterozygosity level to be tested
#' @param blackswan blackswan term sets a limit on probability for each event
#' @return sum of log likelihood ratios for LOH vs no LOH
#'
#' @export
llr_loh_all <- function(dat, delta, blackswan) {
  llr_loh(dat$minor_count, dat$depth, dat$maf, dat$het_mean, dat$het_sum,
          dat$er, delta, blackswan)
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
#' @param hm het mean frequency from baseline
#' @param M het alpha + beta parameters for beta binomial alt formulation
#' @param er error rate for a given SNP
#' @param alpha contamination fraction to be tested
#' @param delta LoH delta to be tested
#' @param blackswan blackswan term sets a limit on probability for each event
#' @return log likelihood ratio for contamination vs no contamination
#'
#' @export
llr_cont <- function(ad, depth, maf1, maf2, hm, M, er, alpha, delta, blackswan) {
  l_min <- blackswan / depth
  l_het <- l_min + (1 - l_min) * l_loh(ad, depth, maf1, maf2, hm, M, delta, 0, er)
  l_cont <- l_min + (1 - l_min) * l_loh(ad, depth, maf1, maf2, hm, M, delta, alpha, er)
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
  weighted.mean(llr_cont_all(dat, alpha, blackswan), log(dat$depth),
                na.rm = TRUE)
}

#' Log likelihood ratios for contamination for a set of SNPs
#'
#' Average log likelihood ratio across a given set of SNPs
#'
#' @param dat data table containing required fields
#' @param alpha contamination level level to be tested
#' @param blackswan blackswan term sets a limit on probability for each event
#' @return sum of log likelihood ratios for LOH vs no LOH
#'
#' @export
llr_cont_all <- function(dat, alpha, blackswan, delta = 0) {
  llr_cont(dat$minor_count, dat$depth, dat$maf, dat$maf, dat$het_mean,
           dat$het_sum, dat$er, alpha, delta, blackswan)
}

#' Likelihood of observing the data with LOH
#'
#' Calculates the posterior probability of observing a given alternative allele
#' depth and total depth given a contamination level (alpha) and a loss of
#' heterozygosity deviation level (delta). Prior probabilities on the host and
#' contaminate genotypes are calculated using the given minor allele frequencies
#' for the host (maf1) and contaminate (maf2).
#'
#' Description of expectation of alt allele:
#' Host hom ref, cont hom ref: just the error rate
#' Host hom ref, cont het: error rate times expectation of alt
#' Host hom ref, cont hom alt: error rate times twice the expectation of alt
#' Host het, cont hom ref: case 1) d - a/2: delta pushes the expected
#' het frequency up, while contamination pulls it down since contaminant does
#' not have the alt allele
#' Host het, cont hom ref: case 2) -d -a/2: same thing delta is in the other
#' direction. Since we don't know which chromosome is affected by LoH, we test
#' probability of both delta and -delta, and multiply each by half
#' Host het, cont het: just d or -d, because host and contaminant have the
#' same alleles.
#' Host het, cont hom alt: Opposite of case 1 and 2, where contamination is
#' actually pulling the allele frequency up because contaminant is hom alt.
#' Host hom alt cases are mirror images of host hom ref with corresponding cont
#'
#' @param ad alternative allele depth for a given SNP
#' @param dp total depth for a given SNP
#' @param maf1 minor allele frequencies for host sample
#' @param maf2 minor allele frequencies for candidate contaminant
#' @param hm mean heterozygote frequency from baseline
#' @param M het alpha + beta parameters for beta binomial alt formulation
#' @param d LOH level as deviation from 0.5
#' @param a alpha contamination level
#' @param er error rate
#' @TODO try an alternate formulation of l_loh modeling the total alt and
#'       ref counts from host and contaminant
#' @importFrom stats dbinom
#' @export
l_loh <- function(ad, dp, maf1, maf2, hm, M = 100, d = 0, a = 0, er = 3e-4) {
  pg2(0, 0, maf1, maf2) * dbinom(ad, dp, er) +
    pg2(0, 1, maf1, maf2) * dbinom(ad, dp, er + hm * a) +
    pg2(0, 2, maf1, maf2) * dbinom(ad, dp, er + 2 * hm * a) +
    pg2(1, 0, maf1, maf2) * (0.5 * dbetabinom(ad, dp, nloh(hm, d - a/2), M) +
                           0.5 * dbetabinom(ad, dp, nloh(hm, -d - a/2), M)) +
    pg2(1, 1, maf1, maf2) * (0.5 * dbetabinom(ad, dp, nloh(hm, d), M) +
                           0.5 * dbetabinom(ad, dp, nloh(hm, -d), M)) +
    pg2(1, 2, maf1, maf2) * (0.5 * dbetabinom(ad, dp, nloh(hm, d + a/2), M) +
                           0.5 * dbetabinom(ad, dp, nloh(hm, -d + a/2), M)) +
    pg2(2, 0, maf1, maf2) * dbinom(ad, dp, 1 - er - 2 * (1 - hm) * a) +
    pg2(2, 1, maf1, maf2) * dbinom(ad, dp, 1 - er - (1 - hm) * a) +
    pg2(2, 2, maf1, maf2) * dbinom(ad, dp, 1 - er)
}

#' Normalize expected allele frequency by LoH or contamination
#'
#' This function calculates the expected allele frequency for a SNP based on
#' additional contamination or LoH. hm is the original expected minor allele
#' fraction. If it is say 0.2, we expect 0.2 fraction of alleles to be from
#' minor allele and 0.8 fraction from reference. This expectation is further
#' modified by delta value which can be due to LoH or contamination. 0.5 + d
#' is the expected modifier for allele1, and 0.5 - d is the modifier for
#' allele2. If we are testing d = 0.2, with hm = 0.2, we expect 0.2 * 0.7 = 0.14
#' fraction of alleles to come from minor allele and 0.8 * 0.3 = 0.24 from
#' reference allele. When normalized this corresponds to 0.14 / (0.14 + 0.24) =
#' ~0.37. This is the expected minor allele frequency. For LoH likelihood, we
#' would also test delta = -0.2 at the same time since SNPs are not phased.
#' In that case, we get 0.2 * 0.3 = 0.06 for minor allele and 0.8 * 0.7 = 0.56
#' for reference, which when normalized is 0.06 / (0.06 + 0.56) = 0.1 minor
#' allele frequency.
#'
#'
#' @param hm mean heterozygote frequency baseline for a specific SNP
#' @param d delta from level of heterozygote where 0.5 is no deviation case.
#'     If there is no deviation (d = 0), then function will return a value
#'     equal to hm. Otherwise it will increase or decrease hm based on the
#'     level of difference of this value from 0.5 considering the original hm.
#' @param d_min d should not be less than this number
#' @param d_max d should not be more than this number
#'
#' @export
nloh <- function(hm, d, d_min = -0.49, d_max = 0.49) {

  # d should be less than 0.49 and greater than 0.01
  d <- pmax(d_min, pmin(d_max, d))

  return(hm * (0.5 + d) / (hm * (0.5 + d) + (1 - hm) * (0.5 - d)))
}

#' Calculate probability density for alternative parameterization of beta
#' binomial distribution
#'
#' @param k total successes binomial
#' @param n total trials binomial
#' @param mu expected mean from beta prior (alpha / (alpha + beta))
#' @param M sum from beta prior (alpha + beta)
#'
#' @export
dbetabinom <- function(k, n, mu, M) {
  betaprob <- choose(n, k) * beta(k + mu * M, n - k + M * (1 - mu)) /
    beta(mu * M, M * (1 - mu))
  return(betaprob)
}

#' Calculate log probability density for alternative parameterization of beta
#' binomial distribution
#'
#' @param k total successes binomial
#' @param n total trials binomial
#' @param mu mean from beta (alpha / (alpha + beta))
#' @param M sum from beta (alpha + beta)
#' @param mu_min mu should not be less than this number
#' @param mu_max mu should not be more than this number
#'
#' @export
ldbetabinom <- function(k, n, mu, M, mu_min = 0.01, mu_max = 0.99) {

  mu <- pmax(mu_min, pmin(mu_max, mu))
  lbetaprob <- lchoose(n, k) + lbeta(k + mu * M, n - k + M * (1 - mu)) -
    lbeta(mu * M, M * (1 - mu))
  return(lbetaprob)
}

#' Return initial test range for LoH
#'
#' @export
get_initial_loh_range <- function() {
  return(c(0.1, 0.2, 0.3, 0.4, 0.45, 0.49))
}

#' Skewed sigmoid with center around provided mean
#'
#' @param x input value
#' @param m the value when x = 0
#'
#' @export
ssig <- function(x, m = 0.5) {
  1 / (1 + exp(-x) * (1 - m) / m)
}

#' Biased sigmoid function to convert a given number to a proportion
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
#' @param conta_cf contamination fraction to be tested
#' @param min_loh to worry about affecting contamination results,
#'     if LOH is less than this fraction do not call
#' @param min_snps minimum number of snps required to evaluate a region for LOH
#' @param max_snps max snps to evaluate in a region, if greater number of
#'     SNPs are present, subsample to this many SNPs.
#' @param het_mean_min minimum het allele frequency to call
#' @param seed random seed
#'
#' @return data.table containing stats per bin
#'
#' @importFrom utils write.table
#' @export
get_per_bin_loh <- function(dat, save_dir, sample, min_lr = 0.01, blackswan, conta_cf,
                            min_loh = 0.2, min_snps = 20, max_snps = 200,
                            het_mean_min = 0.25, seed = 1359) {

  results <- data.table(chrom = character(), chunk = numeric(),
                        snps = numeric(),
                        loh_mle = numeric(), loh_val = numeric(),
                        cont_mle = character(), cont_val = numeric(),
                        loh = logical())

  # combs contains combination of chrom and chunk
  combs <- expand.grid(unique(dat$chrom), unique(dat$chunk))
  combs <- combs %>%
    arrange(Var1, Var2)

  for (i in 1:nrow(combs)) {

    # Get the data in this specific chunk
    dat_chunk <- dat %>%
      dplyr::filter(chrom == combs[i, 1],
                    chunk == combs[i, 2])

    # Remove SNPs without baseline or having means below or above a threshold
    dat_chunk <- as.data.table(dat_chunk) %>%
      dplyr::filter(!is.na(het_mean) &
                      het_mean > het_mean_min &
                      het_mean < (1 - het_mean_min))

    # Set dat as a subset of dat if it exceeds a pre-determined size
    if (nrow(dat_chunk) > max_snps) {
      set.seed(seed)
      dat_chunk <- dat_chunk[sort(sample(1:nrow(dat_chunk), max_snps)), ]
    }

    # Optimize LOH starting on a grid (if there are enough SNPs)
    loh_opt <- optimize_after_grid(dat_chunk, "avg_llr_loh", blackswan,
                                   cf = get_initial_loh_range())

    # If loh_opt is not null, its objective (avg. likelihood ratio) is above
    # a specified threshold, its level is above a minimum LOH level, then
    # calculate contamination mle for the same region. Do not call a region as
    # LOH otherwise, also do not call if there aren't enough SNPs to access.
    if (!is.na(loh_opt$objective) &&
        nrow(dat_chunk) >= min_snps &&
        loh_opt$objective >= min_lr &&
        abs(loh_opt$maximum) >= min_loh) {

      #cont_opt <- optimize_after_grid(dat_chunk, "avg_llr_cont", blackswan)
      conta_opt <- avg_llr_cont(dat_chunk, alpha = conta_cf, blackswan)

      dfres <- data.table(chrom = combs[i, 1], chunk = combs[i, 2],
                          snps = nrow(dat_chunk),
                          loh_mle = loh_opt$maximum,
                          loh_val = loh_opt$objective,
                          cont_mle = conta_cf,
                          cont_val = conta_opt,
                          loh = (loh_opt$objective > conta_opt &
                                      loh_opt$objective > min_lr))
    } else {

      dfres <- data.table(chrom = combs[i, 1], chunk = combs[i, 2],
                          snps = nrow(dat_chunk),
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
#' @param cf contamination fraction to be tested
#' @return maxima and its objective
#'
#' @importFrom stats optimize
#' @export
optimize_after_grid <- function(dat_chunk, LR_func, blackswan, cf) {

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
