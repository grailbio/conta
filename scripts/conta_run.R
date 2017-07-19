#!/usr/bin/Rscript

suppressMessages(library(conta))
suppressMessages(library(optparse))

# constants

# functions
main <- function() {

  start_time <- proc.time()

  option_list <- list(
    make_option(c("-i", "--tsv_file"), type = "character", default = NULL,
                help = "input tsv file with SNP counts", metavar = "character"),
    make_option(c("-b", "--bin_file"), type = "character", default = NA,
                help = "bin tsv file from bedtools", metavar = "character"),
    make_option(c("-c", "--cnv_file"), type = "character", default = NA,
                help = "CNV QC file with CNV calls", metavar = "character"),
    make_option(c("-d", "--dbSNP_file"), type = "character", default = NULL,
                help = "input vcf file from dbSNP", metavar = "character"),
    make_option(c("-s", "--sample"), type = "character", default = "test",
                help = "sample name to prefix out", metavar = "character"),
    make_option(c("-l", "--lr_th"), type = "numeric", default = 0.05,
                help = "Likelihood ratio threshold", metavar = "numeric"),
    make_option(c("-m", "--sim_level"), type = "numeric", default = 0,
                help = "Add simulated cf. 0 means none.", metavar = "numeric"),
    make_option(c("-p", "--min_depth"), type = "numeric", default = 5,
                help = "Minimum depth for a SNP.", metavar = "numeric"),
    make_option(c("-o", "--save_dir"), type = "character", default = "out",
                help = "output folder", metavar = "character"),
    make_option(c("-r", "--remove_maf_file"), type = "logical", default = TRUE,
                help = "remove maf file", metavar = "logical"),
    make_option(c("-n", "--baseline"), type = "character", default = NA,
                help = "baseline file", metavar = "character"),
    make_option(c("-u", "--cores"), type = "numeric", default = 2,
                help = "cpu cores", metavar = "numeric"))

  opt_parser <- OptionParser(option_list = option_list);
  opt <- parse_args(opt_parser);

  if (is.null(opt$tsv_file) | is.null(opt$dbSNP_file)) {
    print_help(opt_parser)
    stop("At least input file and dbSNP must be supplied.\n", call. = FALSE)
  }

  # Create output dir if it doesn't exist
  dir.create(opt$save_dir, showWarnings = FALSE)

  # Intersect counts file with dbSNP
  maf_file <- paste(opt$save_dir, "/", opt$sample, ".maf.tsv", sep = "")

  if (file.exists(maf_file)) {
    # File is already interesected, skip it
    message(paste("Skipping intersection since ", maf_file, "already exists"))
  } else {
    conta::intersect_snps(opt$tsv_file, maf_file, opt$dbSNP_file, FALSE)
  }

  # Run conta
  conta::conta_main(tsv_file = maf_file, sample = opt$sample,
                    save_dir = opt$save_dir, bin_file = opt$bin_file,
                    cnv_file = opt$cnv_file, lr_th = opt$lr_th,
                    sim_level = opt$sim_level, baseline = opt$baseline,
                    min_depth = opt$min_depth, cores = opt$cores)

  # Remove maf file
  if (opt$remove_maf_file)
    file.remove(maf_file)

  proc.time() - start_time
}

if (sys.nframe() == 0) {
  main()
}
