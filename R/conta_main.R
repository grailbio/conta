# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' Run conta.
#'
#' Reads a given counts file, processes it and calculates
#' contamination likelihood for a set of ranges.
#'
#' @param tsv_file input tsv file
#' @param metrics_file input tsv metrics file in long format
#' @param sample experiment name (basename for outputs)
#' @param save_dir output folder
#' @param lr_th min avg. likelihood ratio per SNP to make a call, this
#'     number is highly dependent on the data type used and should be optimized
#'     by the user based on a sensitivity - specificity trade-off.
#' @param sim_level if non-zero, a contaminant at this level will
#'     be added to the current sample. The contaminant will be randomly
#'     generated from minor allele frequencies of the SNPs in the population.
#' @param baseline baseline file for blacklist or noise model
#' @param min_depth minimum depth for a SNP to be considered
#' @param max_depth maximum depth for a SNP to be considered
#' @param loh_lr_cutoff minimum likelihood ratio to call a region as LOH
#' @param loh_delta_cutoff minimum delta (het deviation) to call a region as LOH
#' @param loh_min_snps minimum number of SNPs in a region to consider LOH
#' @param loh_max_snps maxmimum number of SNPs in a region to use for LOH, if
#'     there are more SNPs, they are subsampled to this number
#' @param cf_correction cf correction calculated from empirical data
#' @param min_maf minimum minor allele frequency to include a SNP
#' @param min_cf minimum contamination fraction to call
#' @param blackswan blackswan term for maximum likelihood estimation
#' @param outlier_frac fraction of outlier SNPs (based on depth) to remove
#' @param tsv_rev_file input tsv file for reverse strand reads, if this option
#'     is provided then first tsv_file is considered as positive strand reads
#' @param cores number of cores to be used for parallelization
#' @param subsample Either NA (use all SNPs) or number of SNPs to subsample to
#' @param context_mode whether to run with errors calculated in 3-base context
#' @param seed random seed
#'
#' @return none
#'
#' @importFrom utils packageVersion
#' @importFrom utils write.table
#' @export
conta_main <- function(tsv_file, sample, save_dir, metrics_file = "",
                       lr_th = 0.001, sim_level = 0, baseline = NA,
                       min_depth = 10, max_depth = 10000, loh_lr_cutoff = 0.01,
                       loh_delta_cutoff = 0.3, loh_min_snps = 20,
                       loh_max_snps = 1000, min_maf = 0.01, subsample = NA,
                       cf_correction = 0, min_cf = 0.0001, blackswan = 1,
                       outlier_frac = 0.002, tsv_rev_file = NA,
                       cores = 2, context_mode = FALSE,
                       chr_y_male_threshold = 0.0005, seed = 1359) {

  options("digits" = 8)
  options("mc.cores" = cores)
  set.seed(seed)

  # Create output folder if it doesn't exist
  dir.create(save_dir, showWarnings = FALSE)

  # Write out a default TSV result file in case we exit early from fail_test
  empty_result <- empty_result(sample)
  out_file <- file.path(save_dir, paste(sample, "conta.tsv", sep = "."))
  write.table(empty_result, file = out_file, sep = "\t", row.names = FALSE,
              quote = FALSE)

  # Read baseline if available
  if (!is.na(baseline))
    baseline <- read_data_table(baseline)

  # Prep snp counts
  dat <- read_and_prep(tsv_file, tsv_rev_file, baseline)

  # Set dat as a subset of dat if it exceeds a pre-determined size
  if (!is.na(subsample) & nrow(dat) > subsample) {
    dat <- dat[sort(sample(1:nrow(dat), subsample)), ]
  }

  # Simulate contamination if sim_level is non-0
  dat <- sim_conta(dat, sim_level)

  # Original depth by chr
  plot_depth_by_chr(dat, save_dir, sample, min_depth, ext_plot = "depth.png")

  # Add in more useful fields and filter based on depth and outlier fractions
  dat <- annotate_and_filter(dat, min_depth = min_depth, max_depth = max_depth,
                             out_frac = outlier_frac)

  # Fail if there is no data, or one of the genotypes is never observed
  fail_test(dat)

  # Depth by chr after filters
  plot_depth_by_chr(dat, save_dir, sample, min_depth,
                    ext_plot = "filtered.depth.png")

  # Calculate substitution rates per base and add them to SNP data table
  EE <- calculate_error_model(dat, save_dir, sample,
    context_mode = context_mode)
  dat <- add_error_rates(dat, EE, context_mode)

  # Remove low maf positions (they were kept for error model calculations)
  dat <- dat[ maf > min_maf & maf < (1 - min_maf), ]

  # Add experimental bayesian genotype estimation.
  dat <- bayesian_genotype(dat)

  # Fail if there is no data, or one of the genotypes is never observed
  fail_test(dat)

  # Obtain results and plot lr without LOH filter applied yet
  result <- optimize_likelihood(dat, lr_th, save_dir, sample,
                                loh = FALSE, blackswan, min_cf, cf_correction,
                                outlier_frac = outlier_frac)

  # Remove loss of heterozygosity regions
  bin_stats <- get_per_bin_loh(dat, save_dir, sample, min_lr = loh_lr_cutoff,
                               blackswan = blackswan,
                               min_loh = loh_delta_cutoff,
                               min_snps = loh_min_snps,
                               max_snps = loh_max_snps)
  dat_loh <- exclude_high_loh_regions(dat, bin_stats)

  # TODO: Re-calculate substitution rates after LOH exclusion
  # Note this requires SNPs removed from analyses after LOH removal and
  # re-calculation of error rates. Perhaps give up on visualization without
  # LOH filters (because that step would be costly with all the SNPs)
  #EE <- calculate_error_model(dat_loh, save_dir, sample,
  #  context_mode = context_mode)
  #dat <- add_error_rates(dat_loh, EE, context_mode)

  # Plot minor allele ratio plot (.vr) with LOH
  plot_minor_ratio(dat, dat_loh, save_dir, sample)

  # Fail if there is no data, or one of the genotypes is never observed
  fail_test(dat_loh)

  # Calculate likelihood, max likelihood, and whether to call it
  result <- optimize_likelihood(dat_loh, lr_th, save_dir,
                                sample, loh = TRUE, blackswan, min_cf,
                                cf_correction, outlier_frac = outlier_frac)

  # Calculate chr X and Y counts from metrics file if provided
  chr_y_stats <- chr_stats(metrics_file, "chrY")
  chr_x_stats <- chr_stats(metrics_file, "chrX")
  x_het_snps <- nrow(dat[gt == "0/1" & (chrom == "X" | chrom == "chrX"), ])

  # Make a sex/gender call for the host
  sex <- get_sex_call(chr_y_stats, chr_y_male_threshold, result)

  # Make a pregnancy call for the contaminant
  pregnancy <- get_pregnancy_call(result, sex,
                                  chr_y_stats, chr_y_male_threshold)

  # Final results table
  max_result <- data.table(conta_version = packageVersion("conta"),
                           sample = sample,
                           sex = sex,
                           format(result, digits = 5, trim = TRUE),
                           y_count = round(chr_y_stats$count, 4),
                           y_norm_count = round(chr_y_stats$normalized_count, 6),
                           y_fraction_covered = round(chr_y_stats$fraction_covered, 4),
                           x_count = round(chr_x_stats$count, 4),
                           x_norm_count = round(chr_x_stats$normalized_count, 6),
                           x_fraction_covered = round(chr_x_stats$fraction_covered, 4),
                           x_het_snps = x_het_snps,
                           pregnancy = pregnancy,
                           excluded_regions = bin_stats[loh == TRUE, .N],
                           error_rate = round(mean(EE$er), 7),
                           round(t(data.frame(EE$er,
                                 row.names = rownames(EE))), digits = 7))

  # Write conta results
  write.table(max_result, file = out_file, sep = "\t", row.names = FALSE,
              quote = FALSE)

  # Write genotype files for all SNPs
  write_gt_file(dat, max_result, blackswan, outlier_frac,
                file.path(save_dir, paste(sample, "gt.tsv", sep = ".")))

  # Write genotype files for LOH-excluded SNPs
  write_gt_file(dat_loh, max_result, blackswan, outlier_frac,
                file.path(save_dir, paste(sample, "gt.loh.tsv", sep = ".")))
}
