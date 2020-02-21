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
#' This returns a list of 2, first element being the results in wide format,
#' and the second element being the results in long format. There are 8 metrics
#' for each pair of samples: 'Sample1_snps', 'Sample2_snps', 'Hom_Ref_to_Hom_Alt',
#' 'Hom_Alt_to_Hom_Ref', 'Het_Change', 'No_Change', 'Concordance', 'Call'.
#' Sample1_snps and Sample2_snps are the total number of SNPs in the
#' respecitive samples gt file. 'Hom_Ref_to_Hom_Alt' and 'Hom_Alt_to_Hom_Ref'
#' are the counts of homozygous gt changes. 'Het_Change' are all genotype changes
#' either from a homozygous to heterogyous genotype or visa versa. 'No_Change'
#' denotes the SNP count where the same genotype is observed between samples.
#' 'Concordance' is the fraction of SNPs that have the same genotype between
#' samples. 'Call' is true if the concordance is less than or equal to the
#' concordance threshold. Metrics are only reported once for each pair of samples.
#'
#' @param gt_files list of strings, paths to genotype files
#' @param gt_labels list of strings, labels for genotype files
#' @param concordance_threshold float, threshold for genotype concordance
#' @param maf_filter float, minimum maf required to include snp in concordance calculations
#' @param depth_fitler float, minimum depth required to include snp in concordance calculations
#'
#' @return list of 2, first element being the results in wide format, and the
#'         second element being the results in long format
#'
#' @importFrom data.table fread
#' @importFrom dplyr mutate
#' @importFrom reshape2 melt
#' @export
conta_swap <- function(files, labels, concordance_threshold, maf_filter, depth_filter) {
  # Create all pairings of samples, without duplicates
  n <- length(files[[1]])
  combination_indices <- combn(n, 2)
  # Compute concordances values and calls for each pair
  num_combn <- ncol(combination_indices)
  concordance_values <- matrix(ncol=4, nrow=num_combn)
  concordance_calls <- matrix(ncol=4, nrow=num_combn)

  results_wide <- NULL
  for (c in 1:num_combn){
    gt1 <- combination_indices[1, c]
    gt2 <- combination_indices[2, c]
    # Read in the gt files as data tables
    dt1 <- data.table::fread(files[[1]][gt1]) %>% conta::filter_gt_file(maf_filter, depth_filter)
    dt2 <- data.table::fread(files[[1]][gt2]) %>% conta::filter_gt_file(maf_filter, depth_filter)
    # Get sample labels and number of snps
    sample_info <- data.frame("Sample1" = labels[[1]][gt1],
                              "Sample2" = labels[[1]][gt2],
                              "Sample1_snps" = nrow(dt1),
                              "Sample2_snps" = nrow(dt2))
    # Calculate concordance metrics
    conc_info <- conta::genotype_concordance(dt1, dt2)
    # If concordance is less than/equal to the threshold, swap call is true
    conc_info <- conc_info %>%
      dplyr::mutate(Call = ifelse((Concordance <= concordance_threshold), "TRUE", "FALSE"))
    conc_results <- cbind(sample_info, conc_info)
    # Generate wide output
    results_wide <- rbind(results_wide, conc_results)
  }
  # Generate long output
  results_long <- reshape2::melt(results_wide,
                                 id.vars = c("Sample1", "Sample2"),
                                 variable.name = "metric_name",
                                 value.name = "metric_value",
                                 value.factor = TRUE)

  return(list(results_wide, results_long))
}
