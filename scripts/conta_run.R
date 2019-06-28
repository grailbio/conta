#!/usr/bin/Rscript

# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

suppressMessages(library(conta))
suppressMessages(library(optparse))

# constants

# functions
main <- function() {

  start_time <- proc.time()

  option_list <- list(
    make_option(c("-i", "--tsv_file"), type = "character", default = NULL,
                help = "input pileup file with counts on SNP positions",
                metavar = "character"),
    make_option(c("-a", "--tsv_rev_file"), type = "character", default = NULL,
                help = "input pileup file with reverse strand counts (optional),
                when provided, first input file is assumed to contain
                forward (+) strand counts)", metavar = "character"),
    make_option(c("-b", "--biometrics_file"), type = "character", default = "",
                help = "bio-metrics file", metavar = "character"),
    make_option(c("-d", "--dbSNP_file"), type = "character", default = NULL,
                help = "input vcf file from dbSNP", metavar = "character"),
    make_option(c("-c", "--dbSNP_rev_file"), type = "character", default = NULL,
                help = "input vcf file from rev dbSNP", metavar = "character"),
    make_option(c("-s", "--sample"), type = "character", default = "test",
                help = "Sample ID / file name to prefix out", metavar = "character"),
    make_option(c("-l", "--lr_th"), type = "numeric", default = 0.001,
                help = "Likelihood ratio threshold", metavar = "numeric"),
    make_option(c("-m", "--sim_level"), type = "numeric", default = 0,
                help = "Add simulated cf. 0 means none.", metavar = "numeric"),
    make_option(c("-p", "--min_depth"), type = "numeric", default = 10,
                help = "Minimum depth for a SNP.", metavar = "numeric"),
    make_option(c("-g", "--max_depth"), type = "numeric", default = 10000,
                help = "Maximum depth for a SNP.", metavar = "numeric"),
    make_option(c("-z", "--min_cf"), type = "numeric", default = 0.0001,
                help = "Minimum cf to call.", metavar = "numeric"),
    make_option(c("-q", "--min_maf"), type = "numeric", default = 0.01,
                help = "Minimum maf for a SNP.", metavar = "numeric"),
    make_option(c("-o", "--save_dir"), type = "character", default = "out",
                help = "output folder", metavar = "character"),
    make_option(c("-r", "--remove_maf_file"), type = "logical", default = TRUE,
                help = "remove maf file", metavar = "logical"),
    make_option(c("-n", "--baseline"), type = "character", default = NA,
                help = "baseline file for blacklist or noise model",
                metavar = "character"),
    make_option(c("-f", "--loh_cutoff"), type = "numeric", default = 0.001,
                help = "loh likelihood ratio cutoff", metavar = "numeric"),
    make_option(c("-e", "--loh_delta_cutoff"), type = "numeric", default = 0.3,
                help = "loh delta (deviation) cutoff", metavar = "numeric"),
    make_option(c("", "--loh_auto_delta_cutoff"), type = "numeric",
                default = 0.4, help = "loh delta (deviation) auto cutoff",
                metavar = "numeric"),
    make_option(c("", "--loh_min_snps"), type = "numeric", default = 20,
                help = "minimum number of SNPs per chromosomal region",
                metavar = "numeric"),
    make_option(c("", "--loh_max_snps"), type = "numeric", default = 1000,
                help = "maximum number of SNPs per chromosomal region",
                metavar = "numeric"),
    make_option(c("-w", "--blackswan"), type = "numeric", default = 1,
                help = "black swan term for MLE", metavar = "numeric"),
    make_option(c("-v", "--outlier_frac"), type = "numeric", default = 0.002,
                help = "fraction of outliers to remove", metavar = "numeric"),
    make_option(c("-x", "--cf_correction"), type = "numeric", default = 0,
                help = "conta fraction correction", metavar = "numeric"),
    make_option(c("", "--subsample"), type = "numeric", default = NA,
                help = "subsample to number of SNPs", metavar = "numeric"),
    make_option(c("-u", "--cores"), type = "numeric", default = 2,
                help = "cpu cores", metavar = "numeric"),
    make_option(c("-y", "--non_dbSNP"), type = "logical", default = FALSE,
                help = "Include non-dbSNP for error rate", metavar = "logical"),
    make_option(c("", "--fasta"), type = "character", default = NA,
                help = "input genome fasta file for nucleotide context capture,
                        if provided, conta switches to analyzing errors
                        in tri-nucleotide contexts", metavar = "character"),
    make_option(c("", "--chr_y"), type = "numeric", default = 0.0005,
                help = "Chromosome Y normalized read threshold for male calls"),
    make_option(c("", "--seed"), type = "integer", default = 1234,
                help = "seed for simulations", metavar = "integer"))

  opt_parser <- OptionParser(option_list = option_list);
  opt <- parse_args(opt_parser);
  filename_prefix <- opt$sample
  sample_id       <- opt$sample

  if (is.null(opt$tsv_file) | is.null(opt$dbSNP_file)) {
    print_help(opt_parser)
    stop("At least a pileup file and dbSNP must be supplied.\n", call. = FALSE)
  }

  # Create output dir if it doesn't exist
  dir.create(opt$save_dir, showWarnings = FALSE)

  # Intersect counts file with dbSNP
  maf_file <- paste(opt$save_dir, "/", filename_prefix, ".maf.tsv", sep = "")

  if (file.exists(maf_file)) {
    # File is already interesected, skip it
    message(paste("Skipping intersection since ", maf_file, "already exists"))
  } else {
    message(paste("Intersecting", maf_file, "with", opt$dbSNP_file))
    conta::intersect_snps(opt$tsv_file, maf_file, opt$dbSNP_file, opt$non_dbSNP,
                          opt$fasta)
  }

  # Intersect counts with minus strand counts if provided
  if (!is.null(opt$tsv_rev_file) & !is.null(opt$dbSNP_rev_file)) {
    maf_file2 <- paste(opt$save_dir, "/", filename_prefix, ".maf.rev.tsv", sep = "")

    if (file.exists(maf_file2)) {
      # File is already interesected, skip it
      message(paste("Skipping intersection since ",
                    maf_file2, "already exists"))
    } else {
      message(paste("Intersecting", maf_file2, "with", opt$dbSNP_rev_file))
      conta::intersect_snps(opt$tsv_rev_file, maf_file2, opt$dbSNP_rev_file,
                            opt$non_dbSNP, opt$fasta)
    }
  } else {
    maf_file2 <- NA
  }

  # Run conta
  message(paste("Starting conta"))
  conta::conta_main(tsv_file = maf_file,
                    tsv_rev_file = maf_file2,
                    sample_id = sample_id,
                    save_dir = opt$save_dir,
                    metrics_file = opt$biometrics_file,
                    lr_th = opt$lr_th,
                    sim_level = opt$sim_level,
                    baseline = opt$baseline,
                    min_depth = opt$min_depth,
                    max_depth = opt$max_depth,
                    loh_lr_cutoff = opt$loh_cutoff,
                    loh_delta_cutoff = opt$loh_delta_cutoff,
                    loh_auto_delta_cutoff = opt$loh_auto_delta_cutoff,
                    loh_min_snps = opt$loh_min_snps,
                    loh_max_snps = opt$loh_max_snps,
                    min_maf = opt$min_maf,
                    subsample = opt$subsample,
                    cf_correction = opt$cf_correction,
                    min_cf = opt$min_cf,
                    blackswan = opt$blackswan,
                    outlier_frac = opt$outlier_frac,
                    cores = opt$cores,
                    context_mode = !is.na(opt$fasta),
                    chr_y_male_threshold = opt$chr_y,
                    seed = opt$seed)
  message(paste("Done."))

  # Remove maf file
  if (opt$remove_maf_file) {
      message(paste("Removing annotated pileup files.",
                     "Set --remove_maf_file FALSE to keep them."))
      file.remove(maf_file)
    if (!is.na(maf_file2))
      file.remove(maf_file2)
  }

  proc.time() - start_time
}

if (sys.nframe() == 0) {
  main()
}
