#!/usr/bin/Rscript

# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' CLI for swap detection
#' 
#' Reads n genotype files and labels, and calculates the genotype concordance
#' between each sample. Files and labels are paired by index, and the number of
#' files and labels must be equal. Each file should contain a "gt" column. The 
#' result is written as a tsv to the specified output file.
#' 
#' @importFrom optparse OptionParser
#' @importFrom optparse make_option
#' @importFrom utils write.table
#' @importFrom data.table as.data.table
#' @importFrom conta conta_swap

# functions
main <- function() {
  option_list <- list(
    optparse::make_option(c("-f", "--file_paths"), type = "character", default = NULL,
                help = "n comma separated genotype files, n>2. Must contain gt column. Required.",
                metavar = "character"),
    optparse::make_option(c("-l", "--file_labels"), type = "character", default = NULL,
                help = "list of n comma separated labels for genotype files, n>2. Required.",
                metavar = "character"),
    optparse::make_option(c("-o", "--out"), type = "character", default = NULL,
                help = "out file with swap results, required", metavar = "character"))
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser);
  # Require both gt files and output file
  if (is.null(opt$file_paths) | is.null(opt$file_labels) | is.null(opt$out)) {
    stop("Gt input files, file labels, and output file must be supplied.\n", call. = FALSE)
  }
  # Parse comma separated strings into lists of strings
  files <- strsplit(opt$file_paths, ",", fixed=TRUE)
  labels <- strsplit(opt$file_labels, ",", fixed=TRUE)
  # Require that the number of files and labels are equal
  if (length(files[[1]]) != length(labels[[1]])) {
    stop("The number of files and labels must be equal.\n", call = FALSE)
  }
  # Require that there are at least two files
  if (length(files[[1]]) < 2) {
    stop("There must be at least two files.\n", call=FALSE)
  }
  # Run pairwise genotype concordance
  df <- conta::conta_swap(files, labels)
  # Write results to tsv
  utils::write.table(df, opt$out, sep = "\t", row.names = FALSE, quote = FALSE)
}

if (sys.nframe() == 0) {
  main()
}