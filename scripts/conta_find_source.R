#!/usr/bin/Rscript

# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

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
  make_option(c("-s", "--subfolder"), type = "character", default = "",
              help = "specify if conta is in subfolder", metavar = "character"),
  make_option(c("-t", "--threshold"), type = "numeric", default = NA,
              help = "specify conta threshold", metavar = "numeric"),
  make_option(c("-w", "--blackswan"), type = "numeric", default = 1,
              help = "black swan term for MLE", metavar = "numeric"),
  make_option(c("-q", "--outlier_frac"), type = "numeric", default = 0.01,
              help = "fraction of outliers to remove", metavar = "numeric"),
  make_option(c("-m", "--source_threshold"), type = "numeric", default = 0.01,
              help = "specify source threshold", metavar = "numeric"),
  make_option(c("-u", "--cores"), type = "numeric", default = 8,
              help = "cpu cores", metavar = "numeric"))

  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser);

  if (is.null(opt$base) | is.null(opt$out)) {
    print_help(opt_parser)
    stop("Both base folder and output file must be supplied.\n", call. = FALSE)
  }
  # Run conta
  conta::conta_source(opt$base, opt$out,
                      subfolder = opt$subfolder,
                      threshold = opt$threshold,
                      blackswan = opt$blackswan,
                      outlier_frac = opt$outlier_frac,
                      source_threshold = opt$source_threshold,
                      cores = opt$cores)
}

if (sys.nframe() == 0) {
  main()
}
