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
#' @param cores number of cores to be used for parallelization
#'
#' @return none
#'
#' @export
conta_main <- function(tsv_file, sample, save_dir, bin_file = NA, cnv_file = NA,
                      lr_th = 0.05, sim_level = 0, baseline = NA, min_depth = 5,
                      cores = 2) {

  options("digits" = 8)
  options("mc.cores" = cores)

  # Create output folder if it doesn't exist
  dir.create(save_dir, showWarnings = FALSE)

  # Prep snp counts
  dat <- read_and_prep(tsv_file)

  # Simulate contamination if sim_level is non-0
  dat <- sim_conta(dat, sim_level)

  # Add in more useful fields and filter based on depth, maf and other alleles
  dat <- annotate_and_filter(dat, min_depth = min_depth)

  if ( nrow(dat) == 0 ) stop(paste("No SNPs passed filters"))

  # Calculate substitution rates per base and add them to SNP data table
  EE <- calculate_error_model(dat, save_dir, sample)
  dat <- add_error_rates(dat, EE)

  # Remove low maf positions (they were kept for error model calculations)
  dat <- dat[ maf > 0.25, ]

  # Gather hets vs homs
  dat_hom <- dat[gt == "1/1" | gt == "0/0"]
  dat_het <- dat[gt == "0/1"]

  # Calculate likelihood, max likelihood, and whether to call it
  result <- optimize_likelihood(dat_hom, EE, lr_th, save_dir, sample)

  # Read Y chr counts from bin file
  y_count <- count_ychr(bin_file)

  # Read z-scores from cna QC file
  final_stein <- get_final_stein(cnv_file)
  final_mapd <- get_final_mapd(cnv_file)
  call_pass_stein <- ifelse(is.na(final_stein), NA,
                            ifelse(final_stein < 0.5 & result$conta_call,
                                   TRUE, FALSE))

  # Plot het distortion and gather mad of hd metric
  loh_metric <- plot_het_dist(dat_het, save_dir, sample)

  # A male pregnancy is a female host with lower likelihood on X chr
  # and Y chr count above expected (expected for female is ~0.002). We do not
  # currently call a female pregnancy.
  female_host <- !is.na(y_count) & y_count < 0.2
  male_contamination <- result$conta_call & !is.na(y_count) & y_count >= 0.005
  pregnancy <- ifelse(!female_host | !male_contamination, NA,
                      ifelse(result$pos_lr_X <= result$pos_lr_all / 2,
                             TRUE, FALSE))

  # Results to be written
  max_result <- data.table(sample = sample,
                           call_pass_stein = call_pass_stein,
                           format(result, digits = 5),
                           y_count = round(y_count, 4),
                           loh_metric = round(loh_metric, 3),
                           pregnancy = pregnancy,
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
}
