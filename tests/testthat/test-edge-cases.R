# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

context("test edge cases")

test_that("conta test run with 3 genotypes", {

  expect_true(file.exists(edge_3_genotypes))

  # run conta on dummy sample
  conta_main(edge_3_genotypes, "edge_3", out_dir_edge, sim_level = 0, cores = 8,
             outlier_frac = 0)
  conta_out <- file.path(out_dir_edge, "edge_3.conta.tsv")
  expect_true(file.exists(conta_out))
  result <- read_data_table(conta_out)
  expect_true(result[, .N] == 1)
  expect_false(result[, conta_call])
})
