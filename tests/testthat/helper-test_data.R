# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

# Global variables used by tests (helper-* files are loaded by testthat)

data_dir <- "test_data"
in_dir_source <- paste(data_dir, "source_test", sep = "/")
in_dir_source_single <- paste(data_dir, "source_test_single", sep = "/")
in_dir_source_double <- paste(data_dir, "source_test_double", sep = "/")
out_dir <- "test_output"
out_dir_wgs <- "test_output/wgs"
out_dir_targeted <- "test_output/targeted"
out_dir_sim <- "test_sim"
out_dir_source <- "test_source"
out_source <- paste(out_dir_source, "source.tsv", sep = "/")
out_dir_source_single <- "test_source_single"
out_source_single <- paste(out_dir_source_single, "source.tsv", sep = "/")
out_dir_source_double <- "test_source_double"
out_source_double <- paste(out_dir_source_double, "source.tsv", sep = "/")
out_dir_strand <- "test_strand"

# Remove existing test output directories
unlink(out_dir, recursive = TRUE)
unlink(out_dir_sim, recursive = TRUE)
unlink(out_dir_source, recursive = TRUE)
unlink(out_dir_source_single, recursive = TRUE)
unlink(out_dir_source_double, recursive = TRUE)
unlink(out_dir_strand, recursive = TRUE)

# Create them again
dir.create(out_dir, showWarnings = FALSE)
dir.create(out_dir_wgs, showWarnings = FALSE)
dir.create(out_dir_targeted, showWarnings = FALSE)
dir.create(out_dir_strand, showWarnings = FALSE)

# Input files
tsv_file <- sprintf("%s/test.tsv", data_dir)
pileup_file <- sprintf("%s/test.pileup", data_dir)
male_metrics_file <- sprintf("%s/test.male.bio-metrics.txt", data_dir)
female_metrics_file <- sprintf("%s/test.female.bio-metrics.txt", data_dir)
maf_file <- sprintf("%s/test.maf.tsv", out_dir)
maf_file2 <- sprintf("%s/test.maf2.tsv", out_dir)
maf_file3 <- sprintf("%s/test.maf3.tsv", out_dir)
dbSNP_file <- sprintf("%s/test.dbsnp.vcf", data_dir)
dbSNP_file_art <- sprintf("%s/dbSNP_art_subset.vcf", data_dir)
sim_file <- sprintf("%s/test.regular.maf.tsv", data_dir)
error_tsv <- sprintf("%s/test.error.tsv", data_dir)
supp_tsv <- sprintf("%s/test.supp.maf.tsv", data_dir)
wgs_tsv <- sprintf("%s/test.wgs.maf.tsv", data_dir)
targeted_tsv <- sprintf("%s/test.targeted.maf.tsv", data_dir)
baseline <- sprintf("%s/test.posterior.txt", data_dir)
spam_for_strand <- sprintf("%s/test_strand_for.tsv", data_dir)
spam_for_maf <- sprintf("%s/test_strand_for.maf.tsv", out_dir_strand)
spam_rev_strand <- sprintf("%s/test_strand_rev.tsv", data_dir)
spam_rev_maf <- sprintf("%s/test_strand_rev.maf.tsv", out_dir_strand)
for_dbSNP_file <- sprintf("%s/test_strand_dbsnp_pos.tm.vcf", data_dir)
rev_dbSNP_file <- sprintf("%s/test_strand_dbsnp_neg.tm.vcf", data_dir)

# swap files
swap_test_file_1 <- sprintf("%s/test.sample1_wgs.gt.tsv", data_dir)
swap_test_file_2 <- sprintf("%s/test.sample1_wbc.gt.tsv", data_dir)
swap_test_file_3 <- sprintf("%s/test.sample2_wbc.gt.tsv", data_dir)
