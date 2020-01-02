# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test conta")

test_that("conta test wgs run with biometrics file", {

  expect_true(file.exists(wgs_tsv))
  expect_true(file.exists(male_metrics_file))
  expect_true(file.exists(female_metrics_file))

  sample_id <- "wgs"
  # run conta on dummy sample with female metrics
  conta_main(wgs_tsv, sample_id, out_dir_wgs, metrics_file = female_metrics_file,
             cores = 4)
  conta_out <- file.path(out_dir_wgs, "wgs.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  #expect_true(result[, conta_call])
  expect_equal(result[, y_fraction_covered], 0.0016, tolerance = 1e-4)
  expect_equal(result[, y_count], 21802, tolerance = 0)
  expect_equal(result[, y_norm_count], 0.00004, tolerance = 1e-5)
  expect_equal(result[, x_fraction_covered], 0.74, tolerance = 1e-2)
  expect_equal(result[, x_count], 23445561, tolerance = 0)
  expect_equal(result[, x_norm_count], 0.040, tolerance = 1e-3)

  # run conta on dummy sample with male metrics
  conta_main(wgs_tsv, sample_id, out_dir_wgs, metrics_file = male_metrics_file,
             cores = 4)
  conta_out <- file.path(out_dir_wgs, "wgs.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  #expect_true(result[, conta_call])
  expect_equal(result[, y_fraction_covered], 0.2508, tolerance = 1e-3)
  expect_equal(result[, y_count], 1717521, tolerance = 0)
  expect_equal(result[, y_norm_count], 0.00181, tolerance = 1e-5)
  expect_equal(result[, x_fraction_covered], 0.75, tolerance = 1e-2)
  expect_equal(result[, x_count], 26821414, tolerance = 0)
  expect_equal(result[, x_norm_count], 0.028, tolerance = 1e-3)

  # Check result columns are the same as that for empty results
  expect_equal(colnames(result), colnames(empty_result(sample_id)))
})

test_that("conta test targeted run", {

  expect_true(file.exists(targeted_tsv))

  # run conta on dummy sample
  conta_main(targeted_tsv, "targeted", out_dir_targeted, min_cf = 0.001)
  conta_out <- file.path(out_dir_targeted, "targeted.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  expect_true(result[, conta_call])
  expect_equal(result[, cf], 0.001, tolerance = 2e-4)
})
