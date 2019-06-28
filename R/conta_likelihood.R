# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' Calculate log likelihood ratio
#'
#' Calculate ll-ratio between cf+noise vs noise.
#' Binomial model for each probability requires average error mu.
#'
#' @param cp numeric contamination prior based on maf and genotype
#' @param depth numeric number of reads
#' @param cf numeric tested contamination fraction
#' @param mu numeric average error rate for this type of SNP
#' @param ad numeric number of alternative alleles called
#' @param blackswan numeric limit likelihood of individual hypotheses
#' @param limit_lh numeric limit maximum for minimum lh
#'
#' @return numeric log likelihood ratio
#'
#' @importFrom stats dbinom
#' @export
log_lr <- function(cp, depth, cf, mu, ad, blackswan = 0.01, limit_lh = 0.01) {

  # Probability of het contamination
  binom_cont <- dbinom(ad, size = depth, prob = (mu + cf), log = FALSE)

  # Probability of noise
  binom_noise <- dbinom(ad, size = depth, prob = mu, log = FALSE)

  # Minimum likelihood for each hypothesis to set a limit
  min_lh <- pmin(limit_lh, blackswan / depth ^ 2)

  # Likelihood of contamination (lower limit at min_lh)
  lh <- (1 - min_lh) * (cp * binom_cont + (1 - cp) * binom_noise) + min_lh

  # Likelihood of no contamination (lower limit at min_lh)
  lh0 <- (1 - min_lh) * binom_noise + min_lh

  # Log likelihood ratio
  return(log(lh / lh0))
}

#' Calculate likelihood ratio for each SNP
#'
#' @param cf numeric contamination fraction to be tested
#' @param dat data.table containing counts and metrics per SNP, hets filtered
#' @param blackswan blackswan term for maximum likelihood
#' @param outlier_frac fraction of outliers to remove
#' @param het_rate ratio of heterozygote to homozygote alt contaminant SNPs
#'
#' @return numeric avg. log-likelihood ratio for the given cf for each SNP
#'
#' @importFrom stats quantile weighted.mean sd
#' @export
get_lr <- function(cf, dat, blackswan = 0.05, outlier_frac = 0.005,
                   het_rate = 0.8) {

  # Remove het genotype and calculate log likelihood ratio for each SNP
  dat <- dat[gt != "0/1"]
  dat$lr <- log_lr(dat$cp, dat$depth, get_exp_cf(cf, het_rate = het_rate),
                   dat$er, dat$vr, blackswan)

  # Remove outliers
  dat <- dat[lr < quantile(lr, 1 - outlier_frac) &
               lr > quantile(lr, outlier_frac)]

  return(dat)
}

#' Conta log likelihood ratio
#'
#' @param cf numeric contamination fraction to be tested
#' @param dat data.table containing counts and metrics per SNP, hets filtered
#' @param save_dir directory to save results
#' @param filename_prefix file name prefix for printing purposes
#' @param loh whether to plot in LOH mode
#' @param blackswan blackswan term for maximum likelihood
#' @param outlier_frac fraction of outliers to remove
#' @param het_rate ratio of heterozygote to homozygote alt contaminant SNPs
#'
#' @return numeric avg. log-likelihood ratio for the given cf or a more detailed
#'   result if the save_dir and filename_prefix are specified.
#'
#' @importFrom stats quantile weighted.mean sd
#' @export
conta_llr <- function(cf, dat, save_dir = NA, filename_prefix = NA,
                      loh = FALSE, blackswan = 0.05, outlier_frac = 0.005,
                      het_rate = 0.8) {

  # Likelihood for each SNP to be contaminated at cf level
  dat <- get_lr(cf, dat, blackswan, outlier_frac, het_rate)

  # Calculate average likelihoods of contamination per SNP per depth
  if (is.na(save_dir) | is.na(filename_prefix)) {
    return(max(weighted.mean(dat[, lr], dat[, depth], na.rm = TRUE), 0))
  } else {

    sum_log_lr <- sum(dat[, lr], na.rm = TRUE)
    avg_log_lr <- weighted.mean(dat[, lr], dat[, depth], na.rm = TRUE)
    pos_lr_all <- mean(dat[, lr] > 0, na.rm = TRUE)

    # Generate stats per chr
    per_chr <- dat[, .(.N, depth = mean(depth, na.rm = TRUE),
                       avg_lr = weighted.mean(lr, depth, na.rm = TRUE),
                       pos_lr = sum(lr > 0, na.rm = TRUE),
                       pos_ratio = mean(lr > 0, na.rm = TRUE),
                       vfn = mean(vfn, na.rm = TRUE),
                       vr = mean(vr, na.rm = TRUE)),
                   by = chrom]

    # Generate likelihood ratio plots
    if (loh) {
      plot_lr(save_dir, filename_prefix, dat, per_chr,
              ext_chr_table = "per_chr.loh.tsv",
              ext_loh_table = "per_bin.loh.tsv",
              ext_loh_plot = "bin.lr.loh.png")
    } else {
      plot_lr(save_dir, filename_prefix, dat, per_chr)
    }

    # Pos lr ratio on X chr tells us whether the contaminant is male or female
    pos_lr_x <- per_chr[chrom == "X" | chrom == "chrX", pos_ratio]
    pos_cv <- as.numeric(per_chr[, .(sd(pos_ratio, na.rm = TRUE) /
                            mean(pos_ratio, na.rm = TRUE))])

    return(data.table( cf = round(cf, 5),
                       sum_log_lr = max(sum_log_lr, 0),
                       avg_log_lr = max(avg_log_lr, 0),
                       hom_snps = dat[, .N],
                       depth = round(dat[, mean(depth)], 0),
                       pos_lr_all = pos_lr_all,
                       pos_lr_x = ifelse(length(pos_lr_x) > 0, pos_lr_x, NA),
                       pos_lr_chr_cv = pos_cv,
                       stringsAsFactors = FALSE ))
  }
}

#' Optimize a set of likelihood metrics across SNPs using optimize function.
#'
#'
#' @param dat data.table containing counts and metrics per SNP, hets filtered
#' @param lr_th numeric likelihood threshold to make a call
#' @param save_dir character location to save the results
#' @param filename_prefix character file name prefix
#' @param loh whether to plot in LOH mode
#' @param blackswan blackswan term for maximum likelihood
#' @param min_cf minimum contamination fraction to call
#' @param cf_correction cf correction which is subtracted from calculated cf
#' @param outlier_frac fraction of outliers to remove
#'
#' @return data.table of likelihoods and related metrics
#'
#' @importFrom stats optimize
#' @export
optimize_likelihood <- function(dat, lr_th, save_dir = NA, filename_prefix = NA,
                                loh = FALSE, blackswan, min_cf, cf_correction,
                                outlier_frac) {

  # Do an initial grid search
  cf <- get_initial_range()
  grid_lr <- parallel::mclapply(cf, function(x)
    conta_llr(x, dat, blackswan = blackswan, outlier_frac = outlier_frac))

  # Cf that gives the max result from the initial grid search
  vmax <- which.max(grid_lr)

  # If max cf is the last or the first value, shift to prevent out of bounds
  vmax <- ifelse(vmax == length(cf), vmax - 1, ifelse(vmax == 1, 2, vmax))

  # Optimize avg. likelihood around the cf that gave best result
  opt_val <- optimize(conta_llr, lower = cf[vmax - 1], upper = cf[vmax + 1],
                      dat = dat, blackswan = blackswan,
                      outlier_frac = outlier_frac,
                      tol = 1e-6, maximum = TRUE)

  # Get optimized result
  result <- conta_llr(opt_val$maximum, dat, save_dir, filename_prefix,
                     loh = loh, blackswan = blackswan,
                     outlier_frac = outlier_frac)

  # Correct cf by an empiricially calculated factor (depends on characteristics
  # of the samples such as depth and error rate variance). A cf_correction
  # factor should be calculated for accurate cf readouts at levels of
  # contamination close to the limit of detection, but it is not required.
  result$cf <- max(result$cf - cf_correction, 0)

  # Make a call based on likelihood threshold
  conta_call <- (result$avg_log_lr >= lr_th) & (result$cf >= min_cf)

  # Plot the likelihood distribution and the maximum likelihood (in red cross)
  # Do not plot if this was not a post-loh run
  if (conta_call & loh)
    plot_max_likelihood(cf, grid_lr, opt_val, save_dir, filename_prefix)

  return(cbind(conta_call, result))
}

#' Calculate likelihood ratio based on known source genotypes
#'
#' @param gt1 genotype/variant call table of the host
#' @param gt2 genotype/variant call table of the contaminant
#' @param cf numeric contamination fraction
#' @param blackswan blackswan term for maximum likelihood
#' @param outlier_frac fraction of outliers to remove
#' @param match_prob prior probability of contamination
#' @param detailed_results return likelihood ratios per SNP
#'
#' @return likelihood ratio from the source
#'
#' @importFrom stats weighted.mean
#' @export
get_source_llr <- function(gt1, gt2, cf, blackswan = 1,
                           outlier_frac = 0.005, match_prob = 0.9,
                           detailed_results = FALSE) {

  # Merge host with source candidate
  m1 <- merge(gt1, gt2, by = "rsid")

  # Prepare data table for conta
  dat <- data.table(rsid = m1$rsid,
                    cp = ifelse(m1$gt.x != m1$gt.y, match_prob, 1 - match_prob),
                    gt = m1$gt.x,
                    depth = m1$dp.x,
                    er = m1$er.x,
                    vr = m1$vr.x)

  if (detailed_results) {
    return(get_lr(cf = cf, dat = dat, blackswan = blackswan,
                  outlier_frac = outlier_frac))
  } else {
    return(conta_llr(cf = cf, dat = dat, blackswan = blackswan,
                     outlier_frac = outlier_frac))
  }
}

#' Calculate genotype concordance between two samples
#'
#' @param gt1 genotype/variant call table of the host
#' @param gt2 genotype/variant call table of the contaminant
#' @return numeric fraction of matching genotypes between two samples
#'
#' @export
get_genotype_concordance <- function(gt1, gt2) {

  # Merge host with source candidate
  m1 <- merge(gt1, gt2, by = "rsid")

  return(mean(m1$gt.x == m1$gt.y, na.rm = TRUE))
}
