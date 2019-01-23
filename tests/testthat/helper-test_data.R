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
out_dir_context <- "test_output/context"

# Remove existing test output directories
unlink(out_dir, recursive = TRUE)
unlink(out_dir_sim, recursive = TRUE)
unlink(out_dir_source, recursive = TRUE)
unlink(out_dir_source_single, recursive = TRUE)
unlink(out_dir_source_double, recursive = TRUE)
unlink(out_dir_strand, recursive = TRUE)
unlink(out_dir_context, recursive = TRUE)

# Create them again
dir.create(out_dir, showWarnings = FALSE)
dir.create(out_dir_wgs, showWarnings = FALSE)
dir.create(out_dir_targeted, showWarnings = FALSE)
dir.create(out_dir_strand, showWarnings = FALSE)
dir.create(out_dir_context, showWarnings = FALSE)

# Conta input files
dat_tsv <- sprintf("%s/test.dat.tsv", data_dir)
wgs_tsv <- sprintf("%s/test.wgs.maf.tsv", data_dir)
male_metrics_file <- sprintf("%s/test.male.bio-metrics.txt", data_dir)
female_metrics_file <- sprintf("%s/test.female.bio-metrics.txt", data_dir)
targeted_tsv <- sprintf("%s/test.targeted.maf.tsv", data_dir)
# Error input files
error_tsv <- sprintf("%s/test.error.tsv", data_dir)
supp_tsv <- sprintf("%s/test.supp.maf.tsv", data_dir)
baseline_tsv <- sprintf("%s/test.baseline.tsv", data_dir)
# Intersect input files
maf_file <- sprintf("%s/test.maf.tsv", out_dir)
maf_file2 <- sprintf("%s/test.maf2.tsv", out_dir)
maf_file3 <- sprintf("%s/test.maf3.tsv", out_dir)
context_tsv <- sprintf("%s/test.context.tsv", data_dir)
context_out_tsv <- sprintf("%s/test.context.maf.tsv", out_dir)
context_error_tsv <- sprintf("%s/test.error.tsv", out_dir_context)
tsv_file <- sprintf("%s/test.tsv", data_dir)
dbSNP_file <- sprintf("%s/test.dbsnp.vcf", data_dir)
pileup_file <- sprintf("%s/test.pileup", data_dir)
spam_for_maf <- sprintf("%s/test_strand_for.maf.tsv", out_dir_strand)
# Simulation input files
sim_file <- sprintf("%s/test.regular.maf.tsv", data_dir)
# Source input files
dbSNP_file_art <- sprintf("%s/dbSNP_art_subset.vcf", data_dir)
# Strand input files
spam_for_strand <- sprintf("%s/test_strand_for.tsv", data_dir)
spam_rev_strand <- sprintf("%s/test_strand_rev.tsv", data_dir)
spam_rev_maf <- sprintf("%s/test_strand_rev.maf.tsv", out_dir_strand)
rev_dbSNP_file <- sprintf("%s/test_strand_dbsnp_neg.tm.vcf", data_dir)
for_dbSNP_file <- sprintf("%s/test_strand_dbsnp_pos.tm.vcf", data_dir)
# Swap input files
swap_sim_file_1 <- sprintf("%s/test.sim1.gt.tsv", data_dir)
swap_sim_file_2 <- sprintf("%s/test.sim2.gt.tsv", data_dir)
swap_sim_file_3 <- sprintf("%s/test.sim3.gt.tsv", data_dir)
test_genome <- sprintf("%s/test.genome.fa", data_dir)
