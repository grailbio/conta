#!/usr/bin/Rscript

# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

parse_args <- function() {
  option_list <- list(
    optparse::make_option(c("-i", "--gt-loh-file-loc"), type = "character", default = NULL,
                          help = "TSV of sample gt.loh.", metavar = "character",
                          dest = "gt_file_input"),
    optparse::make_option(c("-v", "--vcf"), type = "character", default = NULL,
                          help = "1000 genomes vcf file", metavar = "character",
                          dest = "vcf"),
    optparse::make_option(c("-o", "--out-dir"), type = "character", default = NULL,
                          help = "Path to output directory", metavar = "character",
                          dest = "out_dir"))

  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser);

  if (is.null(opt$gt_file_input) || is.null(opt$out_dir)) {
    stop("Both sample input and output dir must be supplied.\n")
  }

  return(opt)
}

main <- function() {
  opt <- parse_args()
  if (!dir.exists(opt$out_dir)) {
    dir.create(opt$out_dir)
  }

  # Run conta ethnicity workflow
  conta::conta_ethnicity(gt_loh_file_loc = opt$gt_file_input,
                         vcf = opt$vcf,
                         outdir = opt$out_dir)
}

if (sys.nframe() == 0) {
  main()
}
