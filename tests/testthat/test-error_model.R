# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test error model")

test_that("error model correctly calculates error rate", {

  expect_true(file.exists(error_tsv))
  fdat <- annotate_and_filter(read_and_prep(error_tsv))
  EE <- calculate_error_model(fdat, out_dir, "errorTest")
  error_out_file <- file.path(out_dir, "errorTest.error.tsv")
  expect_true(file.exists(error_out_file))
  EE <- read_data_table(error_out_file)

  # Correction factor for masked positions (only SNPs are used in this test)
  cof <- 2/3
  avg_error <- mean(c(1 / 200, 1 / 200, 1 / 200, 2 / 300, 1 / 300,
                     1 / 100, 1 / 300) / cof)
  tol <- 1e-4
  expect_equal(EE[ref == "A" & subs == "T", er], 1 / 200 / cof, tolerance = tol)
  expect_equal(EE[ref == "A" & subs == "C", er], 1 / 200 / cof, tolerance = tol)
  expect_equal(EE[ref == "A" & subs == "G", er], avg_error, tolerance = tol)
  expect_equal(EE[ref == "G" & subs == "T", er], 1 / 200 / cof, tolerance = tol)
  expect_equal(EE[ref == "G" & subs == "A", er], avg_error, tolerance = tol)
  expect_equal(EE[ref == "G" & subs == "C", er], avg_error, tolerance = tol)
  expect_equal(EE[ref == "C" & subs == "G", er], 2 / 300 / cof, tolerance = tol)
  expect_equal(EE[ref == "C" & subs == "A", er], 1 / 300 / cof, tolerance = tol)
  expect_equal(EE[ref == "C" & subs == "T", er], 1 / 300 / cof, tolerance = tol)
  expect_equal(EE[ref == "T" & subs == "A", er], 1 / 100 / cof, tolerance = tol)
  expect_equal(EE[ref == "T" & subs == "G", er], avg_error, tolerance = tol)
  expect_equal(EE[ref == "T" & subs == "C", er], avg_error, tolerance = tol)
})

test_that("error model correctly calculates error rate with supp positions", {

  expect_true(file.exists(supp_tsv))
  fdat <- annotate_and_filter(read_and_prep(supp_tsv))
  EE <- calculate_error_model(fdat, out_dir, "suppTest")
  error_out_file <- file.path(out_dir, "suppTest.error.tsv")
  expect_true(file.exists(error_out_file))
  EE <- read_data_table(error_out_file)

  # Correction factor for masked positions not including supplementary positions
  cof <- 2/3 * fdat[, mean(maf > 0)] + fdat[, mean(maf == 0)]
  avg_error <- mean(c(1 / 200, 2 / 200, 3 / 100) / cof)
  tol <- 1e-4
  expect_equal(EE[ref == "A" & subs == "T", er], avg_error, tolerance = tol)
  expect_equal(EE[ref == "A" & subs == "C", er], avg_error, tolerance = tol)
  expect_equal(EE[ref == "A" & subs == "G", er], avg_error, tolerance = tol)
  expect_equal(EE[ref == "G" & subs == "T", er], 1 / 200 / cof, tolerance = tol)
  expect_equal(EE[ref == "G" & subs == "A", er], 2 / 200 / cof, tolerance = tol)
  expect_equal(EE[ref == "G" & subs == "C", er], avg_error, tolerance = tol)
  expect_equal(EE[ref == "C" & subs == "G", er], avg_error, tolerance = tol)
  expect_equal(EE[ref == "C" & subs == "A", er], 3 / 100 / cof, tolerance = tol)
  expect_equal(EE[ref == "C" & subs == "T", er], avg_error, tolerance = tol)
  expect_equal(EE[ref == "T" & subs == "A", er], avg_error, tolerance = tol)
  expect_equal(EE[ref == "T" & subs == "G", er], avg_error, tolerance = tol)
  expect_equal(EE[ref == "T" & subs == "C", er], avg_error, tolerance = tol)
})

test_that("error model correctly calculates tri-base context", {

  conta::intersect_snps(context_tsv, context_out_tsv, dbSNP_file, FALSE,
                        test_genome)
  expect_true(file.exists(context_out_tsv))

  conta::conta_main(tsv_file = context_out_tsv, sample_id = "test",
    save_dir = out_dir_context, context_mode = TRUE, min_maf = 0,
    outlier_frac = 0, min_depth = 1)

  # Read the error out file
  error_out_file <- paste0(out_dir_context, "/test.error.tsv")
  EE <- read_data_table(error_out_file)

  # SNP:0 should not be counted because it is a het. GCC > GAC
  expect_true(EE[ref == "GCC" & subs == "GAC", reads] == 0)
  expect_true(EE[ref == "GCC" & subs == "GAC", denom] == 0)
  expect_true(EE[ref == "GAC" & subs == "GCC", reads] == 0)
  expect_true(EE[ref == "GAC" & subs == "GCC", denom] == 0)

  # SNP-1: CCT > CAT error should not be counted because it is the SNP allele
  # and SNP alles are supposed to be masked from error rate calculation
  expect_true(EE[ref == "CCT" & subs == "CAT", reads] == 0)
  expect_true(EE[ref == "CCT" & subs == "CAT", denom] == 30)

  # SNP-2: GAA > GGA error should be counted
  expect_true(EE[ref == "GAA" & subs == "GGA", reads] == 1)
  expect_true(EE[ref == "GAA" & subs == "GGA", denom] == 31)

  # SNP-3: ACG > ATG error should not be counted.
  # ACG > AAG error should be counted.
  expect_true(EE[ref == "ACG" & subs == "ATG", reads] == 0)
  expect_true(EE[ref == "ACG" & subs == "AAG", reads] == 1)
  expect_true(EE[ref == "ACG" & subs == "ATG", denom] == 7)

  # ATG > ACG error and denom should not be counted.
  expect_true(EE[ref == "ATG" & subs == "ACG", reads] == 0)
  expect_true(EE[ref == "ATG" & subs == "ACG", denom] == 0)

})
