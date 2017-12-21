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
#'
#' @return numeric log likelihood ratio
#'
#' @export
log_lr <- function(cp, depth, cf, mu, ad, blackswan = 0.01) {

  # Probability of het contamination
  binom_cont <- dbinom(ad, size = depth, prob = (mu + cf), log = FALSE)

  # Probability of noise
  binom_noise <- dbinom(ad, size = depth, prob = mu, log = FALSE)

  # Minimum likelihood for each hypothesis to set a limit
  min_lh <- blackswan / depth

  # Likelihood of contamination
  lh <- (1 - min_lh) * (cp * binom_cont + (1 - cp) * binom_noise) + min_lh

  # Likelihood of no contamination
  lh0 <- (1 - min_lh) * binom_noise + min_lh

  # Log likelihood ratio
  return(log(lh / lh0))
}

#' Calculate a set of likelihood metrics across SNPs
#'
#' @param cf numeric contamination fraction to be tested
#' @param dat data.table containing counts and metrics per SNP, hets filtered
#' @param EE data.frame error model
#' @param save_dir directory to save results
#' @param sample name of the sample for printing purposes
#' @param loh whether to plot in loh mode
#' @param blackswan blackswan term for maximum likelihood
#' @param out_frac fraction of outliers to remove
#' @param min_maf minimum allele frequency used to filter SNPs
#'
#' @return numeric avg. log-likelihood ratio for the given cf or a more detailed
#'   result if the save_dir and sample name are specified.
#'
#' @export
conta_lr <- function(cf, dat, EE, save_dir = NA, sample = NA,
                     loh = FALSE, blackswan, out_frac = 0.01, min_maf = 0.25) {

  # Likelihood for each SNP to be contaminated at this level
  lr <- log_lr(dat$cp, dat$depth, get_exp_cf(cf), dat$er, dat$vr, blackswan)

  # Calculate average likelihoods of contamination per SNP per depth
  if (is.na(save_dir) | is.na(sample)) {

    # Remove outliers
    cp <- dat[lr < quantile(lr, 1 - out_frac) & lr > quantile(lr, out_frac), cp]
    lr <- lr[lr < quantile(lr, 1 - out_frac) & lr > quantile(lr, out_frac)]
    return(max(weighted.mean(lr, cp, na.rm = TRUE), 0))
  } else {

    # Likelihood for each SNP to be contaminated at this level
    dat$lr <- lr

    # Remove outliers
    dat <- dat[lr < quantile(lr, 1 - out_frac) & lr > quantile(lr, out_frac)]

    sum_log_lr <- sum(dat[, lr], na.rm = TRUE)
    avg_log_lr <- weighted.mean(dat[, lr], dat[, cp], na.rm = TRUE)
    pos_lr_all <- mean(dat[, lr] > 0, na.rm = TRUE)

    # Generate stats per chr
    per_chr <- dat[, .(.N, depth = mean(depth, na.rm = TRUE),
                       avg_lr = weighted.mean(lr, cp, na.rm = TRUE),
                       pos_lr = sum(lr > 0, na.rm = TRUE),
                       pos_ratio = mean(lr > 0, na.rm = TRUE),
                       vfn = mean(vfn, na.rm = TRUE),
                       vr = mean(vr, na.rm = TRUE)),
                   by = chrom]

    # Generate likelihood ratio plots
    if (loh) {
      plot_lr(save_dir, sample, dat, per_chr,
              ext_chr_table = "per_chr.loh.tsv",
              ext_loh_table = "per_bin.loh.tsv",
              ext_loh_plot = "bin.lr.loh.png",
              min_maf = min_maf)
    } else {
      plot_lr(save_dir, sample, dat, per_chr, min_maf = min_maf)
    }

    # Pos lr ratio on X chr tells us whether the contaminant is male or female
    pos_lr_x <- per_chr[chrom == "chrX", pos_ratio]
    pos_cv <- as.numeric(per_chr[, .(sd(pos_ratio, na.rm = TRUE) /
                            mean(pos_ratio, na.rm = TRUE))])

    return(data.table( cf = round(cf, 5),
                       sum_log_lr = max(sum_log_lr, 0),
                       avg_log_lr = max(avg_log_lr, 0),
                       snps = dat[, .N],
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
#' @param EE data.frame error model
#' @param lr_th numeric likelihood threshold to make a call
#' @param save_dir character location to save the results
#' @param sample character sample name
#' @param loh whether to plot in LOH mode
#' @param blackswan blackswan term for maximum likelihood
#' @param min_cf minimum contamination fraction to call
#' @param cf_correction cf correction which is subtracted from calculated cf
#' @param outlier_frac fraction of outlier SNPs (based on likelihood) to remove
#'
#' @return data.table of likelihoods and related metrics
#'
#' @export
optimize_likelihood <- function(dat, EE, lr_th, save_dir = NA, sample = NA,
                                loh = FALSE, blackswan, min_cf, cf_correction,
                                outlier_frac) {

  # Do an initial grid search
  cf <- get_initial_range()
  grid_lr <- mclapply(cf, function(x) conta_lr(x, dat, EE,
                                               blackswan = blackswan))

  # Cf that gives the max result from the initial grid search
  vmax <- which.max(grid_lr)

  # If max cf is the last or the first value, shift to prevent out of bounds
  vmax <- ifelse(vmax == length(cf), vmax - 1, ifelse(vmax == 1, 2, vmax))

  # Optimize avg. likelihood around the cf that gave best result
  opt_val <- optimize(conta_lr, lower = cf[vmax - 1], upper = cf[vmax + 1],
                     dat = dat, EE = EE, blackswan = blackswan,
                     tol = 1e-6, maximum = TRUE)

  # Get optimized result
  result <- conta_lr(opt_val$maximum, dat, EE, save_dir, sample,
                     loh = loh, blackswan = blackswan, out_frac = outlier_frac)

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
    plot_max_likelihood(cf, grid_lr, opt_val, save_dir, sample)

  return(cbind(conta_call, result))
}

#' Calculate likelihood ratio based on known source genotypes
#'
#' @param gt1 genotype/variant call table of the host
#' @param gt2 genotype/variant call table of the contaminant
#' @param cf numeric contamination fraction
#' @param blackswan blackswan term for maximum likelihood
#'
#' @return likelihood ratio from the source
#'
#' @export
get_source_lr <- function(gt1, gt2, cf, blackswan) {

  # Merge host with source candidate
  m1 <- merge(gt1, gt2, by = "rsid")

  # Set up contamination probabilities with some arbitrary high and low
  cp <- ifelse(m1$gt.x != m1$gt.y, 0.9, 0.1)

  return(weighted.mean(log_lr(cp, m1$dp.x, get_exp_cf(cf),
                              m1$er.x, m1$vr.x, blackswan), cp, na.rm = TRUE))
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
