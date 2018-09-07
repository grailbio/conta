#!/usr/bin/Rscript

# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' CLI for swap detection
#' 
#' Reads two genotype files and calculates the genotype concordance. Each
#' file should contain a "gt" column. The result is written to the specified
#' output file.
#' 
#' @importFrom conta load_conta_file
#' @importFrom optparse OptionParser
#' @importFrom optparse make_option
#' @importFrom utils write.table
#' @importFrom readr read_tsv
#' @importFrom data.table as.data.table
#' @export

# functions
main <- function() {
  option_list <- list(
    optparse::make_option(c("-a", "--gt1"), type = "character", default = NULL,
                help = "first genotype file, must contain gt column. Required field",
                metavar = "character"),
    optparse::make_option(c("-b", "--gt2"), type = "character", default = NULL,
                help = "second genotype file, must contain gt column. Required field",
                metavar = "character"),
    optparse::make_option(c("-o", "--out"), type = "character", default = NULL,
                help = "out file with swap results, required", metavar = "character"))
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser);
  # Require both gt files and output file
  if (is.null(opt$gt1) | is.null(opt$gt2) | is.null(opt$out)) {
    print_help(opt_parser)
    stop("Both gt input files and output file must be supplied.\n", call. = FALSE)
  }
  # Read gt files to data tables
  dt1 <-data.table::as.data.table(readr::read_tsv(opt$gt1))
  dt2 <-data.table::as.data.table(readr::read_tsv(opt$gt2))
  # Run conta
  df <- conta::genotype_concordance(dt1, dt2)
  # Write results to tsv
  utils::write.table(df, opt$out, sep = "\t", row.names = FALSE, quote = FALSE)
}

if (sys.nframe() == 0) {
  main()
}