#!/usr/bin/Rscript

# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

#' Run conta swap detection
#' 
#' Calculates pairwise concordance between n genotype samples.
#' Number of files and number of labels must match and be greater
#' than two. Genotype files must contain gt column. Files and labels
#' will be paired by index.
#' This returns a data frame with two metrics for each pair of samples:
#' concordance and call. Call is true if the concordance is less than
#' or equal to the concordance threshold. Metrics are only reported 
#' once for each pair of samples.
#'
#' @param gt_files list of strings, paths to genotype files
#' @param gt_labels list of strings, labels for genotype files
#' @param concordance_treshold float, threshold for genotype concordance
#' 
#' @return data.frame with all pairwise concordances and swap calls
#' 
#' @importFrom data.table fread
#' @export
conta_swap <- function(files, labels, concordance_threshold) {
  # Create all pairings of samples, without duplicates
  n <- length(files[[1]])
  combination_indices <- combn(n, 2)
  # Compute concordances values and calls for each pair
  num_combn <- ncol(combination_indices)
  concordance_values <- matrix(ncol=4, nrow=num_combn)
  concordance_calls <- matrix(ncol=4, nrow=num_combn)
  for (c in 1:num_combn){
    gt1 <- combination_indices[1, c]
    gt2 <- combination_indices[2, c]
    # Read in the gt files as data tables
    dt1 <- data.table::fread(files[[1]][gt1])
    dt2 <- data.table::fread(files[[1]][gt2])
    # Calculate concordance
    conc <- genotype_concordance(dt1, dt2)
    concordance_values[c, ] <- c(labels[[1]][gt1], labels[[1]][gt2], "concordance", conc)
    # If concordance is less than /equal to the threshold, swap call is true
    call <- conc <= concordance_threshold
    concordance_calls[c, ] <- c(labels[[1]][gt1], labels[[1]][gt2], "call", call)
  }
  # Merge concordance and call metric, return all results
  results <- rbind(concordance_values, concordance_calls)
  df <- data.frame(results)
  colnames(df) <- c("Sample1", "Sample2", "metric_name", "value_name")
  return(df)
}
