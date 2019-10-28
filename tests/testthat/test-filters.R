# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test functions for conta filters")

test_that("test error rate filter", {

  expect_true(file.exists(wgs_tsv))

  conta_main(tsv_file = wgs_tsv,
             sample_id = "error_filter_test", save_dir = out_dir_error_filter,
             min_depth = 15, cores = 8, error_quantile_filter = 1.0,
             outlier_frac = 0)

  gt_loh_file1 <- file.path(out_dir_error_filter, "error_filter_test.gt.loh.tsv")
  expect_true(file.exists(gt_loh_file1))
  gt_loh1 <- read_data_table(gt_loh_file1)

  error_quant_filter <- 0.9
  conta_main(tsv_file = wgs_tsv,
             sample_id = "error_filter_test2", save_dir = out_dir_error_filter,
             min_depth = 15, cores = 8, error_quantile_filter = error_quant_filter,
             outlier_frac = 0)
  gt_loh_file2 <- file.path(out_dir_error_filter, "error_filter_test2.gt.loh.tsv")
  expect_true(file.exists(gt_loh_file2))
  gt_loh2 <- read_data_table(gt_loh_file2)

  expect_true(nrow(gt_loh1) > nrow(gt_loh2))
  expect_true(max(gt_loh1$er) > max(gt_loh2$er))
  expect_true(max(gt_loh2$er) == quantile(gt_loh1$er, error_quant_filter))

})
