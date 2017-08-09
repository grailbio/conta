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
log_lr <- function(cp, depth, cf, mu, ad, blackswan = 10 ^ -2) {

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
#' @param save.dir directory to save results
#' @param sample name of the sample for printing purposes
#' 
#' @return numeric avg. log-likelihood ratio for the given cf or a more detailed
#'   result if the save.dir and sample name are specified.
#'
#' @export
conta_lr <- function(cf, dat, EE, save.dir = NA, sample = NA) {

  # Likelihood for each SNP to be contaminated at this level
  lr <- log_lr(dat$cp, dat$depth, get_exp_cf(cf), dat$er, dat$vr)

  # Calculate average likelihoods of contamination per SNP per depth
  if (is.na(save.dir) | is.na(sample)) {
    return(max(mean(lr, na.rm = TRUE), 0))
  } else {

    # Likelihood for each SNP to be contaminated at this level
    dat$lr <- lr
    sum_log_lr <- sum(dat[, lr], na.rm = TRUE)
    avg_log_lr <- mean(dat[, lr], na.rm = TRUE)
    pos_lr_all <- mean(dat[, lr] > 0, na.rm = TRUE)

    # Generate stats per chr
    per_chr <- dat[, .(.N, depth = mean(depth, na.rm = TRUE),
                       avg_lr = mean(lr, na.rm = TRUE),
                       pos_lr = sum(lr > 0, na.rm = TRUE),
                       pos_ratio = mean(lr > 0, na.rm = TRUE),
                       avg_lr_ratio = mean(lr, na.rm = TRUE) / avg_log_lr,
                       vfn = mean(vfn, na.rm = TRUE),
                       vr = mean(vr, na.rm = TRUE)),
                   by = chrom]

    plot_lr(save.dir, sample, dat, per_chr)

    # Pos lr ratio on X chr tells us whether the contaminant is male or female
    pos_lr_X <- per_chr[chrom == "chrX", pos_ratio]
    pos_cv <- as.numeric(per_chr[, .(sd(pos_ratio, na.rm = TRUE) /
                            mean(pos_ratio, na.rm = TRUE))])
    return(data.table( cf = round(cf, 5),
                       sum_log_lr = max(sum_log_lr, 0),
                       avg_log_lr = max(avg_log_lr, 0),
                       snps = dat[, .N],
                       depth = round(dat[, mean(depth)], 0),
                       pos_lr_all = pos_lr_all,
                       pos_lr_X = ifelse(length(pos_lr_X) > 0, pos_lr_X, NA),
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
#' @param save.dir character location to save the results
#' @param sample character sample name
#' 
#' @return data.table of likelihoods and related metrics
#'
#' @export
optimize_likelihood <- function(dat, EE, lr_th, save.dir = NA, sample = NA) {

  # Do an initial grid
  cf <- get_initial_range()
  grid_lr <- mclapply(cf, function(x) conta_lr(x, dat, EE))

  # Cf that gives the max result from the initial grid search
  vmax <- which.max(grid_lr)

  # If max cf is the last or the first value, shift to prevent out of bounds
  vmax <- ifelse(vmax == length(cf), vmax - 1, ifelse(vmax == 1, 2, vmax))

  # Optimize avg. likelihood around the cf that gave best result
  opt_val <- optimize(conta_lr, lower = cf[vmax - 1], upper = cf[vmax + 1],
                     dat = dat, EE = EE, tol = 1e-6, maximum = TRUE)
  
  # Get optimized result
  result <- conta_lr(opt_val$maximum, dat, EE, save.dir, sample)
  
  # Make a call based on likelihood threshold
  conta_call <- result$avg_log_lr >= lr_th
  
  # Plot the likelihood distribution and the maximum likelihood (in red cross)
  if (conta_call)
    plot_max_likelihood(cf, grid_lr, opt_val, save.dir, sample)

  return(cbind(conta_call, result))
}

#' Calculate likelihood ratio based on known source genotypes
#'
#' @param gt1 genotype/variant call table of the host
#' @param gt2 genotype/variant call table of the contaminant
#' @param cf numeric contamination fraction
#' @return likelihood ratio from the source
#'
#' @export
get_source_lr <- function(gt1, gt2, cf) {

  # Merge host with source candidate
  m1 <- merge(gt1, gt2, by = "rsid")

  # Set up contamination probabilities with some arbitrary high and low
  cp <- ifelse(m1$gt.x != m1$gt.y, 0.9, 0.1)

  return(mean(log_lr(cp, m1$dp.x, get_exp_cf(cf), m1$er.x, m1$vr.x),
              na.rm = TRUE))
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
