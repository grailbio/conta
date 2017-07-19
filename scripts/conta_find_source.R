#!/usr/bin/Rscript

suppressMessages(library(conta))
suppressMessages(library(optparse))

# constants

# functions
main <- function() {

  option_list <- list(
  make_option(c("-b", "--base"), type = "character", default = NULL,
              help = "base folder with conta results", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "out file with source results", metavar = "character"),
  make_option(c("-u", "--cores"), type = "numeric", default = 8,
              help = "cpu cores", metavar = "numeric"))

  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser);

  if (is.null(opt$base) | is.null(opt$out)) {
    print_help(opt_parser)
    stop("Both base folder and output fil must be supplied.\n", call. = FALSE)
  }
  # Run conta
  conta::conta_source(opt$base, opt$out, opt$cores)
}

if (sys.nframe() == 0) {
  main()
}
