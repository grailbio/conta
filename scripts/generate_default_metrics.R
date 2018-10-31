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
    make_option(c("-o", "--save_dir"), type = "character", default = "out",
                help = "output folder", metavar = "character"),
    make_option(c("-s", "--sample"), type = "character", default = "test",
                help = "sample name to prefix out", metavar = "character"))

  opt_parser <- OptionParser(option_list = option_list);
  opt <- parse_args(opt_parser);

  # Create output dir if it doesn't exist
  dir.create(opt$save_dir, showWarnings = FALSE)

  empty_result <- empty_result(opt$sample)
  out_file <- file.path(opt$save_dir, paste(opt$sample, "conta.tsv", sep = "."))
  utils::write.table(empty_result, file = out_file, sep = "\t",
                     row.names = FALSE, quote = FALSE)
}

if (sys.nframe() == 0) {
  main()
}
