#!/usr/bin/Rscript

# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

suppressMessages(library(conta))
suppressMessages(library(optparse))

# functions
main <- function() {

  start_time <- proc.time()

  option_list <- list(
    make_option(c("-i", "--gt-file-loc"), type = "character", default = NULL,
                help = "input pileup file with counts on SNP positions",
                metavar = "character", dest = "gt_file_input"),
    make_option(c("-l", "--gt-loh-file-loc"), type = "character", default = NULL,
                help = "input pileup file with counts on SNP positions",
                metavar = "character", dest = "gt_loh_file_input"),
    make_option(c("-s", "--sample"), type = "character", default = "test",
                help = "Sample ID / file name to prefix out", metavar = "character"),
    make_option(c("-o", "--out-dir"), type = "character", default = NULL,
                help = "Path to output directory",
                metavar = "character", dest = "out_dir")
    )

  opt_parser <- OptionParser(option_list = option_list);
  opt <- parse_args(opt_parser);

  if (is.null(opt$gt_file_input) | is.null(opt$out_dir)) {
    print_help(opt_parser)
    stop("conta/graphics: both input and output are required.\n", call. = FALSE)
  }

  # Create output dir if it doesn't exist
  dir.create(opt$out_dir, showWarnings = FALSE)

  # Run conta
  message(paste("Starting conta"))
  conta::conta_graphics(gt_file = opt$gt_file_input,
                        gt_loh_file = opt$gt_loh_file_input,
                        sample = opt$sample,
                        save_dir = opt$out_dir)
  message(paste("Done"))

  proc.time() - start_time
}

if (sys.nframe() == 0) {
  main()
}
