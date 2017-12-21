#' Run conta.
#'
#' Reads a given counts file, processes it and calculates
#' contamination likelihood for a set of ranges.
#'
#' @param tsv_file character input tsv file name
#' @param bin_file character input bin file name
#' @param cnv_file character input cnv file name
#' @param sample character name of experiment (basename of file usually)
#' @param save_dir character folder where the output will go to
#' @param lr_th numeric min avg. likelihood ratio per SNP to make a call, this
#'     number is highly dependent on the data type used and should be optimized
#'     by the user based on a sensitivity - specificity trade-off.
#' @param sim_level numeric if non-zero, a mix of this level will
#'     be added to the current sample. The contaminant will be randomly
#'     generated from minor allele frequencies of the SNPs in the population.
#' @param baseline data.table noise model
#' @param min_depth minimum depth for a SNP to be considered by conta
#' @param max_depth maximum depth for a SNP to be considered by conta
#' @param loh_cutoff minimum maximum likelihood to call a region as LOH
#' @param loh_delta_cutoff minimum delta (het deviation) to call a region as LOH
#' @param cf_correction cf correction calculated from empirical data
#' @param min_maf minimum minor allele frequency to include a SNP
#' @param min_cf minimum cf to call
#' @param blackswan blackswan term for maximum likelihood estimation
#' @param outlier_frac fraction of outlier SNPs (based on depth) to remove
#' @param cores number of cores to be used for parallelization
#'
#' @return none
#'
#' @export
conta_main <- function(tsv_file, sample, save_dir, bin_file = NA, cnv_file = NA,
                      lr_th = 0.05, sim_level = 0, baseline = NA, min_depth = 5,
                      max_depth = 10000, loh_cutoff = 0.01,
                      loh_delta_cutoff = 0.05, min_maf = 0.25,
                      cf_correction = 0, min_cf = 0.00025,
                      blackswan = 0.05, outlier_frac = 0.01, cores = 2) {

  options("digits" = 8)
  options("mc.cores" = cores)

  # Create output folder if it doesn't exist
  dir.create(save_dir, showWarnings = FALSE)

  # Prep snp counts
  dat <- read_and_prep(tsv_file)

  # Simulate contamination if sim_level is non-0
  dat <- sim_conta(dat, sim_level)

  # Original depth by chr
  plot_depth_by_chr(dat, save_dir, sample, ext_plot = "depth.png")

  # Add in more useful fields and filter based on depth, maf and other alleles
  dat <- annotate_and_filter(dat, min_depth = min_depth, max_depth = max_depth,
                             out_frac = outlier_frac)

  # Depth by chr after filters
  plot_depth_by_chr(dat, save_dir, sample, ext_plot = "filtered.depth.png")

  # Fail if there is no data, or one of the genotypes is never observed
  fail_test(dat)

  # Calculate substitution rates per base and add them to SNP data table
  EE <- calculate_error_model(dat, save_dir, sample)
  dat <- add_error_rates(dat, EE)

  # Remove low maf positions (they were kept for error model calculations)
  dat <- dat[ maf > min_maf & maf < (1 - min_maf), ]

  # Fail if there is no data, or one of the genotypes is never observed
  fail_test(dat)

  # Plot minor allele ratio
  plot_minor_ratio(dat, save_dir, sample)

  # Obtain results and plot lr without loh filter applied yet
  result <- optimize_likelihood(dat[gt != "0/1"], EE, lr_th, save_dir, sample,
                                loh = FALSE, blackswan, min_cf, cf_correction,
                                outlier_frac = outlier_frac)

  # Remove LOH regions
  bin_stats <- get_per_bin_loh(dat, save_dir, sample, loh_cutoff,
                               blackswan, min_loh = loh_delta_cutoff)
  dat_loh <- exclude_high_loh_regions(dat, bin_stats)

  # Fail if there is no data, or one of the genotypes is never observed
  fail_test(dat_loh)

  # Calculate likelihood, max likelihood, and whether to call it
  result <- optimize_likelihood(dat_loh[gt != "0/1"], EE, lr_th, save_dir,
                                sample, loh = TRUE, blackswan, min_cf,
                                cf_correction, outlier_frac = outlier_frac)

  # Read Y chr counts from bin file
  y_count <- count_ychr(bin_file)

  # Read z-scores from cna QC file
  final_stein <- get_final_stein(cnv_file)
  final_mapd <- get_final_mapd(cnv_file)

  # A male pregnancy is a female host with lower likelihood on X chr
  # and Y chr count above expected (expected for female is ~0.002). We do not
  # currently call a female pregnancy.
  female_host <- !is.na(y_count) & y_count < 0.2
  male_contamination <- result$conta_call & !is.na(y_count) & y_count >= 0.005
  pregnancy <- ifelse(!female_host | !male_contamination, NA,
                      ifelse(result$pos_lr_x <= result$pos_lr_all / 2,
                             TRUE, FALSE))

  # Results to be written
  max_result <- data.table(conta_version = packageVersion("conta"),
                           sample = sample,
                           format(result, digits = 5),
                           y_count = round(y_count, 4),
                           pregnancy = pregnancy,
                           excluded_regions = bin_stats[loh == TRUE, .N],
                           error_rate = round(mean(EE$er), 7),
                           round(t(data.frame(EE$er,
                                 row.names = rownames(EE))), digits = 7))

  # Write conta results to file
  out_file <- file.path(save_dir, paste(sample, "conta.tsv", sep = "."))
  write.table(max_result, file = out_file, sep = "\t", row.names = FALSE,
              quote = FALSE)

  # Write genotypes and possible contaminant reads
  gt_sum <- rbind(dat[, .(rsid, cp, dp = depth, er, gt, vr)])
  out_file_gt <- file.path(save_dir, paste(sample, "gt.tsv", sep = "."))
  write.table(format(gt_sum, digits = 6), file = out_file_gt, sep = "\t",
              col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Write genotypes and possible contaminant reads
  gt_sum_loh <- rbind(dat_loh[, .(rsid, cp, dp = depth, er, gt, vr)])
  out_file_gt_loh <- file.path(save_dir, paste(sample, "gt.loh.tsv", sep = "."))
  write.table(format(gt_sum_loh, digits = 6), file = out_file_gt_loh,
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}
