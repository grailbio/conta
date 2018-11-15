# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test strand specific")


test_that("conta strand specific run",  {

  expect_true(file.exists(spam_for_strand))
  expect_true(file.exists(spam_rev_strand))

  conta::intersect_snps(spam_for_strand, spam_for_maf, for_dbSNP_file, FALSE)
  conta::intersect_snps(spam_rev_strand, spam_rev_maf, rev_dbSNP_file, FALSE)

  conta_main(tsv_file = spam_for_maf, tsv_rev_file = spam_rev_maf,
             sample = "strand_test", save_dir = out_dir_strand,
             min_depth = 20, cores = 8)
  conta_main(tsv_file = spam_for_maf,
             sample = "strand_test_for", save_dir = out_dir_strand,
             min_depth = 20, cores = 8)
  conta_main(tsv_file = spam_rev_maf,
             sample = "strand_test_rev", save_dir = out_dir_strand,
             min_depth = 20, cores = 8)

  conta_file <- paste(out_dir_strand, "strand_test.conta.tsv", sep = "/")
  conta_file2 <- paste(out_dir_strand, "strand_test_for.conta.tsv", sep = "/")
  conta_file3 <- paste(out_dir_strand, "strand_test_rev.conta.tsv", sep = "/")
  expect_true(file.exists(conta_file))
  expect_true(file.exists(conta_file2))
  expect_true(file.exists(conta_file3))
  result <- read_data_table(conta_file)
  result_for <- read_data_table(conta_file2)
  result_rev <- read_data_table(conta_file3)
  expect_true(result$hom_snps > result_for$hom_snps)
  expect_true(result$hom_snps > result_rev$hom_snps)
  expect_true(result$hom_snps > 0.8 * (result_for$hom_snps + result_rev$hom_snps))
})
