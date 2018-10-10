#!/usr/bin/Rscript

# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

parse_args <- function() {
  option_list <- list(
  optparse::make_option(c("-i", "--simulation-input"), type = "character", default = NULL,
              help = "TSV of simulated contamination results", metavar = "character",
              dest = "simulation_input"),
  optparse::make_option(c("-o", "--out-dir"), type = "character", default = NULL,
              help = "Path to output directory", metavar = "character",
              dest = "out_dir"),
  optparse::make_option(c("-e", "--extreme-threshold"), type = "numeric", default = 1.5,
              help = paste("Samples with uncontaminated avg LLR above this ",
                           "threshold will be filtered. First filter applied.",
                           "Defaults to 1.5"),
              metavar = "numeric",
              dest = "extreme_threshold"),
  optparse::make_option(c("-q", "--quantile-threshold"), type = "numeric", default = 0.98,
              help = paste("Quantile threshold of samples to keep. Second filter",
                           "applied. Defaults to 0.98."),
              metavar = "numeric",
              dest = "quantile_threshold"),
  optparse::make_option(c("-m", "--mad-threshold"), type = "numeric", default = 3,
              help = paste("Samples with uncontaminated avg LLR more than",
                           "`mad-threshold` MADs away from the median will be",
                           "filtered. Final threshold applied. Defaults to 3."),
              metavar = "numeric",
              dest = "mad_threshold"))

  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser);

  if (is.null(opt$simulation_input) || is.null(opt$out_dir)) {
    stop("Both simulation input and output dir must be supplied.\n")
  }

  return(opt)
}

main <- function() {
  opt <- parse_args()
  if (!dir.exists(opt$out_dir)) {
    dir.create(opt$out_dir)
  }

  # Parse simulation results
  sim_data <- conta::read_sim_results(opt$simulation_input)
  conta::conta_threshold(sim_data,
                         opt$out_dir,
                         extreme_level = opt$extreme_threshold,
                         filter_quantile = opt$quantile_threshold,
                         mad_thresh = opt$mad_threshold)
}

if (sys.nframe() == 0) {
  main()
}
