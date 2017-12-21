#!/usr/bin/Rscript
suppressMessages(library(conta))
suppressMessages(library(data.table))
suppressMessages(library(optparse))

main <- function() {

  start_time <- proc.time()

  option_list <- list(
    make_option(c("-p", "--pair_file"), type = "character", default = NULL,
                help = "input match of cfDNA and gDNA", metavar = "character"),
    make_option(c("-o", "--output_file"), type = "character", default = NULL,
                help = "output concordance metrics", metavar = "character"),
    make_option(c("-s", "--shiny_file"), type = "character", default = NULL,
                help = "input locations and metrics", metavar = "character"),
    make_option(c("-d", "--dbsnp_file"), type = "character", default = NULL,
                help = "dbsnp for targeted analysis", metavar = "character"),
    make_option(c("-c", "--corrections"), type = "character", default = NULL,
                help = "swap corrections", metavar = "character"),
    make_option(c("-n", "--signatures"), type = "character", default = NULL,
                help = "signature SNPs", metavar = "character"),
    make_option(c("-r", "--randomize"), type = "character", default = FALSE,
                help = "random pairs for thresholding", metavar = "character")
  )

  opt_parser <- OptionParser(option_list = option_list);
  opt <- parse_args(opt_parser);

  if (is.null(opt$pair_file) | is.null(opt$output_file) |
      is.null(opt$shiny_file) | is.null(opt$dbsnp_file)) {
    print_help(opt_parser)
    stop("All four input and output files must be provided.\n", call. = FALSE)
  }

  # Run conta
  conta::swap(pairing_table_name = opt$pair_file,
              out_table_name = opt$output_file,
              shiny_loc = opt$shiny_file,
              dbsnp_targeted = opt$dbsnp_file,
              randomize = opt$randomize,
              corrections = opt$corrections,
              signatures = opt$signatures)

  proc.time() - start_time

}

if (sys.nframe() == 0) {
  main()
}
